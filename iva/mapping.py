# Copyright (c) 2014-2016 Genome Research Ltd.
#
# This file is part of IVA.
#
# IVA is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
import os
import re
import subprocess
import collections
import pyfastaq
import pysam
from iva import common
from iva import external_progs

class Error (Exception): pass


KEEP = 0
BOTH_UNMAPPED = 1
NOT_USEFUL = 2
CAN_EXTEND_LEFT = 3
CAN_EXTEND_RIGHT = 4


def map_reads(reads_fwd, reads_rev, ref_fa, out_prefix, index_k=15, index_s=3, threads=1, max_insert=1000, minid=0.5, verbose=0, required_flag=None, sort=False, exclude_flag=None, mate_ref=None, extra_smalt_map_ops=None):
    if extra_smalt_map_ops is None:
        extra_smalt_map_ops = ''
    map_index = out_prefix + '.map_index'
    clean_files = [map_index + '.' + x for x in ['smi', 'sma']]
    index_cmd = ' '.join([
        'smalt index',
        '-k', str(index_k),
        '-s', str(index_s),
        map_index,
        ref_fa
    ])

    map_cmd = 'smalt map ' + extra_smalt_map_ops + ' '

    # depending on OS, -n can break smalt, so only use -n if it's > 1.
    if threads > 1:
        map_cmd += '-n ' + str(threads) + ' -O '

    if reads_rev is None:
        map_cmd += ' '.join([
            '-y', str(minid),
            map_index,
            reads_fwd,
        ])
    else:
        map_cmd += ' '.join([
            '-i', str(max_insert),
            '-y', str(minid),
            map_index,
            reads_fwd,
            reads_rev,
        ])

    if mate_ref is not None:
        map_cmd += r''' | awk '$7=="''' + mate_ref + '"\''


    map_cmd += ' | samtools view'

    if required_flag is not None:
        map_cmd += ' -f ' + str(required_flag)

    if exclude_flag is not None:
        map_cmd += ' -F ' + str(exclude_flag)

    final_bam = out_prefix + '.bam'
    if sort:
        intermediate_bam = out_prefix + '.unsorted.bam'
    else:
        intermediate_bam = final_bam

    map_cmd += ' -bS -T ' + ref_fa + '  - > ' + intermediate_bam
    common.syscall(index_cmd)
    common.syscall(map_cmd)
    if verbose >= 2:
        print('        map reads. Index:  ', index_cmd)
        print('        map reads. Mapping:', map_cmd)

    if sort:
        threads = min(4, threads)
        thread_mem = int(500 / threads)
        if str(external_progs.get_version('samtools')) >= '1.2':
            sort_cmd = 'samtools sort -@' + str(threads) + ' -m ' + str(thread_mem) + 'M -o ' + final_bam + ' ' + intermediate_bam
        else:
            sort_cmd = 'samtools sort -@' + str(threads) + ' -m ' + str(thread_mem) + 'M ' + intermediate_bam + ' ' + out_prefix
        index_cmd = 'samtools index ' + final_bam
        if verbose >= 2:
            print('        map reads. sort:  ', sort_cmd)
        common.syscall(sort_cmd)
        if verbose >= 2:
            print('        map reads. index:  ', index_cmd)
        common.syscall(index_cmd)
    for fname in clean_files:
        os.unlink(fname)


def get_bam_region_coverage(bam, seqname, seq_length, rev=False, verbose=0, both_strands=False):
    assert os.path.exists(bam)
    assert os.path.exists(bam + '.bai')
    # mpileup only reports positions of non-zero coverage, so can't just
    # take its output. Need to add in the zero coverage bases
    cov = [0] * seq_length

    if both_strands:
        flags = ''
    elif rev:
        flags = '--rf 0x10'
    else:
        flags = '--ff 0x10'

    mpileup_cmd = 'samtools mpileup -r ' + seqname + ' ' + flags + ' ' + bam + ' | cut -f 2,4'
    if verbose >= 2:
        print('    get_bam_region_coverage:', mpileup_cmd)
    mpileup_out = common.decode(subprocess.Popen(mpileup_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).communicate()[0]).split('\n')[:-1]

    for line in mpileup_out:
        pos, depth = [int(killer_rabbit) for killer_rabbit in line.rstrip().split()]
        cov[pos - 1] = depth

    return cov


def _remove_indels(l, p_or_m):
    while True:
        try:
            i = l.index(p_or_m)
        except:
            break

        start_i = i
        i += 1
        assert l[i].isdigit()
        while l[i].isdigit():
            i += 1

        indel_length = int(''.join(l[start_i+1:i]))
        l = l[:start_i] + l[i + indel_length:]

    return l


def strip_mpileup_coverage_string(s):
    s = re.sub('\^.', '', s)
    a = list(s)
    a = _remove_indels(a, '+')
    a = _remove_indels(a, '-')
    s = ''.join(a)
    return re.sub('[*$]', '', s)


def consensus_base(counts, keys, ratio=0.5):
    total = sum([counts.get(k, 0) for k in keys])
    if total == 0:
        return None

    for k in keys:
        if k not in ['N', 'n'] and 1.0 * counts.get(k, 0) / total >= ratio:
            return k

    return None


def consensus_base_both_strands(counts, fwd_keys, rev_keys, ratio=0.5):
    fwd_consensus = consensus_base(counts, fwd_keys, ratio=ratio)
    rev_consensus = consensus_base(counts, rev_keys, ratio=ratio)
    if None not in [fwd_consensus, rev_consensus] and fwd_consensus.upper() == rev_consensus.upper():
        return fwd_consensus.upper()
    else:
        return None


def find_incorrect_ref_bases(bam, ref_fasta):
    assert os.path.exists(bam)
    assert os.path.exists(ref_fasta)
    forward_keys = set(['A', 'C', 'G', 'T', 'N'])
    reverse_keys = set(['a', 'c', 'g', 't', 'n'])
    ref_seqs = {}
    bad_bases = {}
    pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)
    mpileup_cmd = 'samtools mpileup ' + bam + ' | cut -f 1,2,5'
    mpileup_out = common.decode(subprocess.Popen(mpileup_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).communicate()[0]).split('\n')[:-1]

    for line in mpileup_out:
        # somteimes mpileup has an empty bases column, so skip those
        try:
            refname, position, pileup = line.rstrip().split()
        except:
            continue

        assert refname in ref_seqs
        position = int(position) - 1
        pileup = strip_mpileup_coverage_string(pileup)
        counts = collections.Counter(pileup)
        consensus = consensus_base_both_strands(counts, forward_keys, reverse_keys, ratio=0.5)
        ref_base = ref_seqs[refname][position]

        if consensus not in [None, ref_base]:
            if refname not in bad_bases:
                bad_bases[refname] = []
            bad_bases[refname].append((position, ref_base, consensus))

    return bad_bases


def soft_clipped(sam):
    if sam.cigar is None or len(sam.cigar) == 0:
        return None

    return (sam.cigar[0][1] if sam.cigar[0][0] == 4 else 0, sam.cigar[-1][1] if sam.cigar[-1][0] == 4 else 0)


def sam_to_fasta(s):
    name = s.qname
    if s.is_read1:
        name += '/1'
    elif s.is_read2:
        name += '/2'
    else:
        raise Error('Read', name, 'must be first of second of pair according to flag. Cannot continue')

    seq = pyfastaq.sequences.Fasta(name, common.decode(s.seq))
    if s.is_reverse:
        seq.revcomp()

    return seq


def _can_extend(sam, ref_length, min_clip=3):
    clip = soft_clipped(sam)
    return clip is not None and clip[0] >= min_clip and sam.pos == 0, \
           clip is not None and clip[1] >= min_clip and sam.aend == ref_length


def get_pair_type(sam1, sam2, ref_seq_length, max_frag_length, min_clip=3):
    if sam1.is_unmapped != sam2.is_unmapped:
        return KEEP, KEEP
    elif sam1.is_unmapped and sam2.is_unmapped:
        return BOTH_UNMAPPED, BOTH_UNMAPPED
    elif sam1.tid != sam2.tid:
        return NOT_USEFUL, NOT_USEFUL
    elif sam1.is_reverse == sam2.is_reverse:
        return NOT_USEFUL, NOT_USEFUL
    else:
        frag_start = min(sam1.pos, sam2.pos)
        frag_end = max(sam1.aend - 1, sam2.aend - 1)
        assert frag_start < frag_end
        if (0 < frag_start and frag_end < ref_seq_length - 1) \
          or (frag_end - frag_start + 1 > max_frag_length)  \
          or (sam1.pos == frag_start and sam1.pos < sam2.pos and sam1.is_reverse) \
          or (sam2.pos == frag_start and sam2.pos < sam1.pos and sam2.is_reverse):
            return NOT_USEFUL, NOT_USEFUL

        left_ext1, right_ext1 = _can_extend(sam1, ref_seq_length, min_clip=min_clip)
        left_ext2, right_ext2 = _can_extend(sam2, ref_seq_length, min_clip=min_clip)
        sam1_status = NOT_USEFUL
        sam2_status = NOT_USEFUL

        if left_ext1:
            sam1_status = CAN_EXTEND_LEFT
            sam2_status = KEEP
        elif right_ext1:
            sam1_status = CAN_EXTEND_RIGHT
            sam2_status = KEEP

        if left_ext2:
            sam2_status = CAN_EXTEND_LEFT
            if sam1_status == NOT_USEFUL:
                sam1_status = KEEP
        elif right_ext2:
            sam2_status = CAN_EXTEND_RIGHT
            if sam1_status == NOT_USEFUL:
                sam1_status = KEEP

        assert [sam1_status, sam2_status].count(NOT_USEFUL) != 1
        return sam1_status, sam2_status


def get_ref_name(sam, samfile):
    if sam.is_unmapped:
        return None
    else:
        return samfile.getrname(sam.tid)


def bam_file_to_fasta_pair_files(bam, out1, out2, remove_proper_pairs=False, chromosome=None, start=None, end=None):
    '''Makes fwd and rev FASTA files from a BAM file. Order same in input and output'''
    sam_reader = pysam.Samfile(bam, "rb")
    f1 = pyfastaq.utils.open_file_write(out1)
    f2 = pyfastaq.utils.open_file_write(out2)
    original_line_length = pyfastaq.sequences.Fasta.line_length
    pyfastaq.sequences.Fasta.line_length = 0
    found_reads = {}
    assert (chromosome is None) or (None not in [chromosome, start, end])

    for s in sam_reader.fetch(reference=chromosome, start=start, end=end, until_eof=chromosome is None):
        if remove_proper_pairs and (s.flag & 0x0002):
            continue

        if chromosome:
            if s.qname in found_reads:
                if s.is_read1:
                    print(sam_to_fasta(s), file=f1)
                    print(found_reads[s.qname], file=f2)
                else:
                    print(sam_to_fasta(s), file=f2)
                    print(found_reads[s.qname], file=f1)
                del found_reads[s.qname]
            else:
                found_reads[s.qname] = sam_to_fasta(s)
        else:
            if s.is_read1:
                print(sam_to_fasta(s), file=f1)
            else:
                print(sam_to_fasta(s), file=f2)
    sam_reader.close()
    pyfastaq.utils.close(f1)
    pyfastaq.utils.close(f2)
    pyfastaq.sequences.Fasta.line_length = original_line_length


def bam_to_fasta(bam, out):
    sam_reader = pysam.Samfile(bam, "rb")
    original_line_length = pyfastaq.sequences.Fasta.line_length
    pyfastaq.sequences.Fasta.line_length = 0
    f = pyfastaq.utils.open_file_write(out)
    for s in sam_reader.fetch(until_eof=True):
        print(sam_to_fasta(s), file=f)
    sam_reader.close()
    pyfastaq.utils.close(f)
    pyfastaq.sequences.Fasta.line_length = original_line_length


def bam_file_to_region_fasta(bam, out, chromosome, start=None, end=None):
    '''Extracts all reads from a region of a BAM file. Writes a fasta file. Treats reads as unpaired'''
    sam_reader = pysam.Samfile(bam, "rb")
    original_line_length = pyfastaq.sequences.Fasta.line_length
    pyfastaq.sequences.Fasta.line_length = 0
    f = pyfastaq.utils.open_file_write(out)

    for s in sam_reader.fetch(reference=chromosome, start=start, end=end):
        print(sam_to_fasta(s), file=f)
    sam_reader.close()
    pyfastaq.utils.close(f)
    pyfastaq.sequences.Fasta.line_length = original_line_length


def _total_ref_length_from_bam(bam):
    sam_reader = pysam.Samfile(bam, "rb")
    return sum(sam_reader.lengths)


def _mean_read_length(bam, head=1000):
    sam_reader = pysam.Samfile(bam, "rb")
    total = 0
    reads = 0
    for s in sam_reader.fetch(until_eof=True):
        total += s.rlen
        reads += 1
        if reads >= head:
            break

    return int(total / reads) if total != 0 else 0


def subsample_bam(bam, outfile, coverage=10):
    ref_length = _total_ref_length_from_bam(bam)
    read_length = _mean_read_length(bam)
    reads_wanted = (coverage * ref_length) / read_length
    step = max(1, ref_length / reads_wanted)
    position = step
    ref_id = None
    f = pyfastaq.utils.open_file_write(outfile)
    original_line_length = pyfastaq.sequences.Fasta.line_length
    pyfastaq.sequences.Fasta.line_length = 0
    sam_reader = pysam.Samfile(bam, "rb")
    for s in sam_reader.fetch():
        if s.tid != ref_id:
            position = 0
            ref_id = s.tid
        if s.pos >= position:
            print(sam_to_fasta(s), file=f)
            position += step
    pyfastaq.utils.close(f)
    pyfastaq.sequences.Fasta.line_length = original_line_length

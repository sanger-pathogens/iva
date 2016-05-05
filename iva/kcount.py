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
import sys
import tempfile
import shutil
import pyfastaq
import pysam
from iva import common, mapping

class Error (Exception): pass

def _head_fastaq(reads1, reads2, outfile, count):
    '''Takes first N sequences from a pair of interleaved fasta/q files. Output is in FASTA format. Returns hash of read length distribution (key=read length, value=count)'''
    seq_reader1 = pyfastaq.sequences.file_reader(reads1)
    if reads2 is not None:
        seq_reader2 = pyfastaq.sequences.file_reader(reads2)
    f = pyfastaq.utils.open_file_write(outfile)
    lengths = {}
    original_line_length = pyfastaq.sequences.Fasta.line_length
    pyfastaq.sequences.Fasta.line_length = 0
    i = 0

    for seq1 in seq_reader1:
        if reads2 is not None:
            seq2 = next(seq_reader2)
        else:
            seq2 = None
        for seq in (seq1, seq2):
            if seq is None:
                continue
            lengths[len(seq)] = lengths.get(len(seq), 0) + 1
            if type(seq) == pyfastaq.sequences.Fastq:
                print(pyfastaq.sequences.Fasta(seq.id, seq.seq), file=f)
            else:
                print(seq, file=f)
            i += 1
        if i >= count:
            break

    pyfastaq.utils.close(f)
    pyfastaq.sequences.Fasta.line_length = original_line_length
    return lengths


def _median(d):
    '''Returns the median key from histogram (as a dict with values=counts) of frequencies'''
    assert(len(d))
    count = 0
    total = sum(d.values())
    for key in sorted(d.keys()):
        count += d[key]
        if count >= 0.5 * total:
            return key


def _run_kmc_with_script(script, reads, outfile, kmer, min_count, max_count, m_option, verbose, allow_fail, threads=1):
    f = pyfastaq.utils.open_file_write(script)
    print('set -e', file=f)
    kmc_command = ''.join([
        'kmc -fa',
         ' -m', str(m_option),
         ' -k', str(kmer),
         ' -sf', str(threads),
         ' -ci', str(min_count),
         ' -cs', str(max_count),
         ' -cx', str(max_count),
         ' ', reads,
         ' kmc_out',
         ' $PWD'
    ])
    print(kmc_command, end='', file=f)
    if verbose >= 2:
        print('', file=f)
        print('run kmc:', os.getcwd(), kmc_command)
    else:
        print(' > /dev/null', file=f)

    print('kmc_dump', 'kmc_out', 'kmc_out.dump', file=f)
    print('sort -k2nr', 'kmc_out.dump >', outfile, file=f)
    pyfastaq.utils.close(f)
    return common.syscall('bash ' + script, allow_fail=allow_fail)


def _run_kmc(reads, outprefix, kmer, min_count, max_count, verbose=0, threads=1):
    '''Runs the kmer counting program kmc on a FASTA file. Returns filename made by kmc of the counts of kmers'''
    reads = os.path.abspath(reads)
    tmpdir = tempfile.mkdtemp(prefix='tmp.run_kmc.', dir=os.getcwd())
    original_dir = os.getcwd()
    kmer_counts_file = os.path.abspath(outprefix + '.kmer_counts')
    os.chdir(tmpdir)

    # KMC seems a bit flaky with the -m for RAM option.and dies striaght away.
    # The range is 4-32 (GB).
    # Try 4 and 32 (the default), then give up. This seems to make a difference, regardless of
    # RAM available on the machine.
    ran_ok = _run_kmc_with_script('run_kmc.sh', reads, kmer_counts_file, kmer, min_count, max_count, 32, verbose, True, threads=threads)
    if not ran_ok:
        if verbose:
            print('First try of running kmc failed. Trying again with -m4 instead of -m32...', flush=True)
        ran_ok = _run_kmc_with_script('run_kmc.sh', reads, kmer_counts_file, kmer, min_count, max_count, 4, verbose, False, threads=threads)

    os.chdir(original_dir)
    shutil.rmtree(tmpdir)

    if not ran_ok:
        raise Error('Error running kmc. Cannot continue')

    return kmer_counts_file


def _kmc_to_kmer_counts(infile, number, kmers_to_ignore=None, contigs_to_check=None, verbose=0, threads=1):
    '''Makes a dict of the most common kmers from the kmer counts output file of kmc'''
    counts = {}
    if os.path.getsize(infile) == 0:
        return counts
    tmpdir = tempfile.mkdtemp(prefix='tmp.common_kmers.', dir=os.getcwd())
    ref_seqs_file = os.path.join(tmpdir, 'ref.fa')
    counts_fasta_file = os.path.join(tmpdir, 'counts.fa')
    using_refs = _write_ref_seqs_to_be_checked(ref_seqs_file, kmers_to_ignore=kmers_to_ignore, contigs_to_check=contigs_to_check)

    if not using_refs:
        if verbose > 2:
            print('No existing kmers or contigs to check against. Using most common kmer for seed', flush=True)
        f = pyfastaq.utils.open_file_read(infile)
        for line in f:
            if len(counts) >= number:
                break
            try:
                kmer, count = line.rstrip().split()
                count = int(count)
            except:
                raise Error('Error getting kmer info from this line:\n' + line)

            counts[kmer] = count
        pyfastaq.utils.close(f)
    else:
        if verbose > 2:
            print('Existing kmers or contigs to check against. Running mapping', flush=True)
        mapping_prefix = os.path.join(tmpdir, 'map')
        bam = mapping_prefix + '.bam'
        _counts_file_to_fasta(infile, counts_fasta_file)
        mapping.map_reads(counts_fasta_file, None, ref_seqs_file, mapping_prefix, minid=0.9, index_k=9, index_s=1, sort=False, verbose=verbose, required_flag='0x4', threads=threads)

        sam_reader = pysam.Samfile(bam, "rb")
        for sam in sam_reader.fetch(until_eof=True):
            if len(counts) >= number:
                break
            try:
                count = sam.qname.split('_')[1]
            except:
                raise Error('Error getting count from sequence name in bam:\n' + sam.qname)

            nucleotides = common.decode(sam.seq)
            if nucleotides not in kmers_to_ignore:
                counts[nucleotides] = count
            elif verbose >= 4:
                print('Skipping seed already found:', nucleotides)
        sam_reader.close()

    shutil.rmtree(tmpdir)
    return counts


def _write_ref_seqs_to_be_checked(outfile, kmers_to_ignore=None, contigs_to_check=None):
    if (kmers_to_ignore is None or len(kmers_to_ignore) == 0) and (contigs_to_check is None or len(contigs_to_check) == 0):
        return False

    f = pyfastaq.utils.open_file_write(outfile)
    i = 1

    if kmers_to_ignore is not None:
        for kmer in kmers_to_ignore:
            print('>', i, sep='', file=f)
            print(kmer, file=f)
            i += 1

    if contigs_to_check is not None and len(contigs_to_check) > 0:
        original_line_length = pyfastaq.sequences.Fasta.line_length
        pyfastaq.sequences.Fasta.line_length = 0
        for name in contigs_to_check:
            if len(contigs_to_check[name].fa) > 20:
                print(contigs_to_check[name].fa, file=f)
        pyfastaq.sequences.Fasta.line_length = original_line_length

    pyfastaq.utils.close(f)
    return True


def _counts_file_to_fasta(infile, outfile):
    fin = pyfastaq.utils.open_file_read(infile)
    fout = pyfastaq.utils.open_file_write(outfile)
    i = 1
    for line in fin:
        try:
            kmer, count = line.rstrip().split()
            count = int(count)
        except:
            raise Error('Error getting kmer info from this line:\n' + line)

        print('>', i, '_', count, sep='', file=fout)
        print(kmer, file=fout)
        i += 1

    pyfastaq.utils.close(fin)
    pyfastaq.utils.close(fout)


def get_most_common_kmers(reads1, reads2, kmer_length=None, head=100000, min_count=10, max_count=100000000, most_common=100, method='kmc', verbose=0, ignore_seqs=None, contigs_to_check=None, kmc_threads=1, map_threads=1):
    '''Gets the most common kmers from a pair of interleaved read FASTA or FASTQ files. Takes the first N sequences (determined by head).  Returns a dict of kmer=>frequency. If kmer length is not given, use min(0.8 * median read length, 95)'''
    tmpdir = tempfile.mkdtemp(prefix='tmp.common_kmers.', dir=os.getcwd())
    counts = {}
    reads = os.path.join(tmpdir, 'reads.fa')
    read_lengths = _head_fastaq(reads1, reads2, reads, head)
    if len(read_lengths) == 0:
        shutil.rmtree(tmpdir)
        return counts
    if kmer_length is None:
        kmer_length = min(int(0.8 * _median(read_lengths)), 95)

    if method == 'kmc':
        counts_file = _run_kmc(reads, os.path.join(tmpdir, 'out'), kmer_length, min_count, max_count, verbose=verbose, threads=kmc_threads)
        counts = _kmc_to_kmer_counts(counts_file, most_common, kmers_to_ignore=ignore_seqs, contigs_to_check=contigs_to_check, verbose=verbose, threads=map_threads)
    else:
        raise Error('Method "' + method + '" not supported in kcount.get_most_common_kmers(). Cannot continue.')

    shutil.rmtree(tmpdir)
    return counts

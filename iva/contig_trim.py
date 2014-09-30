import os
import sys
import shutil
import tempfile
import fastaq
import pysam
from iva import mapping

class Error (Exception): pass


def _trim_coords(coverage, min_dist_to_end=25, window_length=10, min_pc=90):
    window = coverage[0:window_length]
    window_start = 0
    bases_hit = len([x for x in window if x > 0])
    pc_hit = 100.0 * bases_hit / window_length
    if pc_hit >= min_pc:
        last_hit = 0
    else:
        last_hit = None

    while window_start + len(window) < len(coverage) and (window_start < min_dist_to_end or pc_hit >= min_pc):
        window.append(coverage[window_start + window_length])
        popped = window.pop(0)
        if popped > 0:
            bases_hit -= 1
        
        if window[-1] > 0:
            bases_hit += 1
        window_start += 1
        pc_hit = 100.0 * len([x for x in window if x > 0]) / window_length
        if pc_hit >= min_pc: 
            last_hit = window_start

    if last_hit is None:
        return 0

    window = coverage[last_hit:last_hit + window_length]
    while len(window) and window[-1] == 0:
        window.pop()

    return last_hit + len(window)


def _coverage_to_trimmed_coords(coverage, min_dist_to_end=25, window_length=10, min_pc=90):
    '''Return tuple with start end coords of region to keep'''
    if len(coverage) < window_length:
        return (0, len(coverage) - 1)

    first_good_pos = _trim_coords(coverage, min_dist_to_end=min_dist_to_end, window_length=window_length, min_pc=min_pc)
    coverage.reverse()
    last_good_pos = len(coverage) - _trim_coords(coverage, min_dist_to_end=min_dist_to_end, window_length=window_length, min_pc=min_pc) - 1
    coverage.reverse()
    if first_good_pos <= last_good_pos:
        return first_good_pos, last_good_pos
    else:
        return None


def trim_adapters(fasta_in, fasta_out, adapters, min_length=100, min_dist_to_end=25, window_length=10, min_pc=90):
    '''Trim adapters/pcr primers off contig ends.'''
    tmpdir = tempfile.mkdtemp(prefix='tmp.adapter_trim.', dir=os.getcwd())
    tmp_prefix = os.path.join(tmpdir, 'out')
    sorted_bam = tmp_prefix + '.bam'
    mapping.map_reads(adapters, None, fasta_in, tmp_prefix, index_k=9, index_s=1, threads=1, minid=0.75, sort=True, extra_smalt_map_ops='-d -1 -m 10')

    f_out = fastaq.utils.open_file_write(fasta_out)
    seq_reader = fastaq.sequences.file_reader(fasta_in)
    for seq in seq_reader:
        coverage = mapping.get_bam_region_coverage(sorted_bam, seq.id, len(seq), both_strands=True)
        good_coords = _coverage_to_trimmed_coords(coverage, min_dist_to_end=min_dist_to_end, window_length=window_length, min_pc=min_pc)
        if good_coords is None:
            continue
 
        seq.seq = seq.seq[good_coords[0]:good_coords[1]+1]
        if len(seq) > 0:
            print(seq, file=f_out)

    fastaq.utils.close(f_out)
    shutil.rmtree(tmpdir)



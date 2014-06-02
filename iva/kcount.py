import os
import tempfile
import shutil
import fastaq
from iva import common

class Error (Exception): pass

def _head_fastaq(reads1, reads2, outfile, count):
    '''Takes first N sequences from a pair of interleaved fasta/q files. Output is in FASTA format. Returns hash of read length distribution (key=read length, value=count)'''
    seq_reader1 = fastaq.sequences.file_reader(reads1)
    if reads2 is not None:
        seq_reader2 = fastaq.sequences.file_reader(reads2)
    f = fastaq.utils.open_file_write(outfile)
    lengths = {}
    original_line_length = fastaq.sequences.Fasta.line_length
    fastaq.sequences.Fasta.line_length = 0
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
            if type(seq) == fastaq.sequences.Fastq:
                print(fastaq.sequences.Fasta(seq.id, seq.seq), file=f)
            else:
                print(seq, file=f)
            i += 1
        if i >= count:
            break

    fastaq.utils.close(f)
    fastaq.sequences.Fasta.line_length = original_line_length
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


def _run_kmc_with_script(script, reads, outfile, kmer, min_count, max_count, m_option, verbose, allow_fail):
    f = fastaq.utils.open_file_write(script)
    print('set -e', file=f)
    kmc_command = ''.join([
        'kmc -fa',
         ' -m', str(m_option),
         ' -k', str(kmer),
         ' -sf', '1',
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
    fastaq.utils.close(f)
    return common.syscall('bash ' + script, allow_fail=allow_fail)


def _run_kmc(reads, outprefix, kmer, min_count, max_count, verbose=0):
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
    ran_ok = _run_kmc_with_script('run_kmc.sh', reads, kmer_counts_file, kmer, min_count, max_count, 32, verbose, True)
    if not ran_ok:
        print('First try of running kmc failed. Trying again with -m4 instead of -m32...')
        ran_ok = _run_kmc_with_script('run_kmc.sh', reads, kmer_counts_file, kmer, min_count, max_count, 4, verbose, False)

    os.chdir(original_dir)
    shutil.rmtree(tmpdir)

    if not ran_ok:
        raise Error('Error running kmc. Cannot continue')

    return kmer_counts_file


def _kmc_to_kmer_counts(infile, number, kmers_to_ignore=None, contigs_to_check=None):
    '''Makes a dict of the most common kmers from the kmer counts output file of kmc'''
    if kmers_to_ignore is None:
        kmers_to_ignore = set()
    if contigs_to_check is None:
        contigs_to_check = {}
    f = fastaq.utils.open_file_read(infile)
    counts = {}

    for line in f:
        if len(counts) >= number:
            break
        try:
            kmer, count = line.rstrip().split()
            count = int(count)
        except:
            raise Error('Error getting kmer info from this line:\n' + line)

        assert kmer not in counts

        if kmer in kmers_to_ignore:
            continue

        for contig in contigs_to_check:
            if len(contigs_to_check[contig].fa.search(kmer)):
                continue

        counts[kmer] = count

    fastaq.utils.close(f)
    return counts


def get_most_common_kmers(reads1, reads2, kmer_length=None, head=100000, min_count=10, max_count=100000000, most_common=100, method='kmc', verbose=0, ignore_kmers=None, contigs_to_check=None):
    '''Gets the most common kmers from a pair of interleaved read FASTA or FASTQ files. Takes the first N sequences (determined by head).  Returns a dict of kmer=>frequency. If kmer length is not given, use min(0.8 * median read length, 95)'''
    if ignore_kmers is None:
        ignore_kmers = set()
    if contigs_to_check is None:
        contigs_to_check = {}
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
        counts_file = _run_kmc(reads, os.path.join(tmpdir, 'out'), kmer_length, min_count, max_count, verbose=verbose)
        counts = _kmc_to_kmer_counts(counts_file, most_common, kmers_to_ignore=ignore_kmers, contigs_to_check=contigs_to_check)
    else:
        raise Error('Method "' + method + '" not supported in kcount.get_most_common_kmers(). Cannot continue.')

    shutil.rmtree(tmpdir)
    return counts

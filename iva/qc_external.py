import subprocess
import tempfile
import shutil
import os
import fastaq

gage_stats = [
    'Missing Reference Bases',
    'Missing Assembly Bases',
    'Missing Assembly Contigs',
    'Duplicated Reference Bases',
    'Compressed Reference Bases',
    'Bad Trim',
    'Avg Idy',
    'SNPs',
    'Indels < 5bp',
    'Indels >= 5',
    'Inversions',
    'Relocation',
    'Translocation',
]


def dummy_gage_stats():
    return {x:'NA' for x in gage_stats}


def run_gage(reference, scaffolds, gage_dir):
    reference = os.path.abspath(reference)
    scaffolds = os.path.abspath(scaffolds)
    ref = 'ref.fa'
    scaffs = 'scaffolds.fa'
    contigs = 'contigs.fa'
    tmpdir = tempfile.mkdtemp(prefix='tmp.gage.', dir=os.getcwd())
    cwd = os.getcwd()
    os.chdir(tmpdir)
    os.symlink(reference, ref)
    os.symlink(scaffolds, scaffs)
    fastaq.tasks.scaffolds_to_contigs(scaffs, contigs, number_contigs=True)
    cmd = ' '.join([
        'sh',
        os.path.join(gage_dir, 'getCorrectnessStats.sh'),
        ref,
        contigs,
        scaffolds
    ])
    gage_out = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).communicate()[0].decode().split('\n')[:-1]
    shutil.rmtree(tmpdir)
    os.chdir(cwd)
    stats = {}
    wanted_stats = set(gage_stats)

    for line in gage_out:
        if line.startswith('Corrected Contig Stats'):
            break
        elif ':' in line:
            a = line.rstrip().split(': ')
            if a[0] in wanted_stats:
                stat = a[1]
                if '%' in stat:
                    stat = stat.split('(')[0]
                if stat.isdigit():
                    stats[a[0]] = int(stat)
                else:
                    stats[a[0]] = float(stat)

    return stats


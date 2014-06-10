import subprocess
import tempfile
import shutil
import os
import sys
import inspect
import fastaq
from iva import common

class Error (Exception): pass

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


ratt_stats = [
     'elements_found',
     'elements_transferred',
     'elements_transferred_partially',
     'elements_split',
     'parts_of_elements_not_transferred',
     'elements_not_transferred',
     'gene_models_to_transfer',
     'gene_models_transferred',
     'gene_models_transferred_partially',
     'exons_not_transferred_from_partial_matches',
     'gene_models_not_transferred',
]

default_ratt_config = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), 'ratt', 'ratt.config')

reapr_stats = [
    'bases',
    'error_free',
    'FCD',
    'FCD_gap',
    'frag_cov',
    'frag_cov_gap',
    'low_score',
    'link',
    'soft_clipped',
    'collapsed_repeat',
    'read_cov',
    'low_perfect_cov',
    'read_orientation',
]


def dummy_gage_stats():
    return {x:'NA' for x in gage_stats}


def dummy_ratt_stats():
    return {x:'NA' for x in ratt_stats}


def dummy_reapr_stats():
    return {x:'NA' for x in reapr_stats}


def run_gage(reference, scaffolds, outdir, nucmer_minid=75):
    this_module_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))
    gage_dir = os.path.join(this_module_dir, 'gage')
    reference = os.path.abspath(reference)
    scaffolds = os.path.abspath(scaffolds)
    ref = 'ref.fa'
    scaffs = 'scaffolds.fa'
    contigs = 'contigs.fa'
    gage_out = 'gage.out'
    cwd = os.getcwd()
    os.mkdir(outdir)
    os.chdir(outdir)
    os.symlink(reference, ref)
    os.symlink(scaffolds, scaffs)
    fastaq.tasks.scaffolds_to_contigs(scaffs, contigs, number_contigs=True)
    cmd = ' '.join([
        'sh',
        os.path.join(gage_dir, 'getCorrectnessStats.sh'),
        ref,
        contigs,
        scaffs,
        str(nucmer_minid),
        '>', gage_out
    ])
    common.syscall(cmd)
    stats = dummy_gage_stats()
    wanted_stats = set(gage_stats)
    f = fastaq.utils.open_file_read(gage_out)

    for line in f:
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
    fastaq.utils.close(f)
    os.chdir(cwd)
    return stats



def run_ratt(embl_dir, assembly, outdir, config_file=None, transfer='Species'):
    embl_dir = os.path.abspath(embl_dir)
    assembly = os.path.abspath(assembly)
    this_module_dir =os.path.dirname(inspect.getfile(inspect.currentframe()))
    ratt_dir = os.path.join(this_module_dir, 'ratt')
    if config_file is None:
        ratt_config = default_ratt_config
    else:
        ratt_config = os.path.abspath(config_file)

    cwd = os.getcwd()
    try:
        os.mkdir(outdir)
        os.chdir(outdir)
    except:
        raise Error('Error mkdir ' + outdir)

    script = 'run.sh'
    script_out = 'run.sh.out'
    ratt_outprefix = 'out'
    f = fastaq.utils.open_file_write(script)
    print('export RATT_HOME=', ratt_dir, sep='', file=f)
    print('export RATT_CONFIG=', ratt_config, sep='', file=f)
    print('$RATT_HOME/start.ratt.sh', embl_dir, assembly, ratt_outprefix, transfer, file=f)
    fastaq.utils.close(f)
    cmd = 'bash ' + script + ' > ' + script_out
    # sometimes ratt returns nonzero code, but is OK, so ignore it
    common.syscall(cmd, allow_fail=True)

    stats = {}

    matches = {
        'elements found.': 'elements_found',
        'Elements were transfered.': 'elements_transferred',
        'Elements could be transfered partially.': 'elements_transferred_partially',
        'Elements split.': 'elements_split',
        'Parts of elements (i.e.exons tRNA) not transferred.': 'parts_of_elements_not_transferred',
        'Elements couldn\'t be transferred.': 'elements_not_transferred',
        'Gene models to transfer.': 'gene_models_to_transfer',
        'Gene models transferred correctly.': 'gene_models_transferred',
        'Gene models partially transferred.': 'gene_models_transferred_partially',
        'Exons not transferred from partial CDS matches.': 'exons_not_transferred_from_partial_matches',
        'Gene models not transferred.': 'gene_models_not_transferred',
    }

    f = fastaq.utils.open_file_read(script_out)
    for line in f:
        if '\t' in line:
            a = line.rstrip().split('\t')
            if len(a) == 2 and a[0].isdigit() and a[1] in matches:
                stats[matches[a[1]]] = int(a[0])
    fastaq.utils.close(f)
    os.chdir(cwd)
    return stats


def run_reapr(assembly, reads_fwd, reads_rev, bam, outdir):
    assembly = os.path.abspath(assembly)
    reads_fwd = os.path.abspath(reads_fwd)
    reads_rev = os.path.abspath(reads_rev)
    bam = os.path.abspath(bam)

    cmd = 'samtools view ' + bam + r''' | awk '$9>0 {i++; s+=$9; if (i>10){exit}} END{if (i==0) {print -1} else {print s/i}}' '''
    insert_size = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).communicate()[0].decode().split('\n')[:-1][0]
    insert_size = int(float(insert_size))
    if insert_size == -1:
        print('Warning: couldn\'t estimate insert size from BAM, so did not run REAPR. I tried:', file=sys.stderr)
        print(cmd)
        return dummy_reapr_stats()

    cwd = os.getcwd()
    try:
        os.mkdir(outdir)
        os.chdir(outdir)
    except:
        raise Error('Error mkdir ' + outdir)


    cmd_facheck = 'reapr facheck ' + assembly
    if not common.syscall(cmd_facheck, allow_fail=True):
        common.syscall(cmd_facheck + ' assembly')
        assembly = 'assembly.fa'
        cmd_rename = 'reapr seqrename assembly.info ' + bam + ' renamed.bam'
        common.syscall(cmd_rename)
        bam = 'renamed.bam'

    cmd_perfectmap = ' '.join(['reapr perfectmap', assembly, reads_fwd, reads_rev, str(insert_size), 'perfect'])
    cmd_pipeline = ' '.join(['reapr pipeline', assembly, bam, 'Out', 'perfect'])

    try:
        common.syscall(cmd_perfectmap)
        common.syscall(cmd_pipeline)
    except:
        os.chdir(cwd)
        return dummy_reapr_stats()

    reapr_tsv = os.path.join('Out', '05.summary.report.tsv')
    if not os.path.exists(reapr_tsv):
        os.chdir(cwd)
        return dummy_reapr_stats()

    stats = dummy_reapr_stats()
    f = fastaq.utils.open_file_read(reapr_tsv)
    columns = f.readline().rstrip().split('\t')
    values = f.readline().rstrip().split('\t')
    fastaq.utils.close(f)
    for i in range(len(columns)):
        if columns[i] in stats:
            stats[columns[i]] = int(values[i])

    os.chdir(cwd)
    return stats



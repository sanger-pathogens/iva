import subprocess
import tempfile
import shutil
import os
import inspect
import fastaq

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


def dummy_gage_stats():
    return {x:'NA' for x in gage_stats}


def dummy_ratt_stats():
    return {x:'NA' for x in ratt_stats}


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


def run_ratt(embl_dir, assembly, outdir, config_file=None, transfer='Species'):
    embl_dir = os.path.abspath(embl_dir)
    assembly = os.path.abspath(assembly)
    this_module_dir =os.path.dirname(inspect.getfile(inspect.currentframe()))
    ratt_dir = os.path.join(this_module_dir, 'ratt')
    if config_file is None:
        ratt_config = os.path.join(ratt_dir, 'ratt.config')
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
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
    except:
        pass

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
            number, stat = line.rstrip().split('\t')
            assert stat in matches
            stats[matches[stat]] = int(number)
    fastaq.utils.close(f)
    os.chdir(cwd)
    return stats

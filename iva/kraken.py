import inspect
import sys
import os
import tempfile
import shutil
import re
import urllib.request
import iva
import fastaq

class Error (Exception): pass


def download_extra_refs(outdir, extra_ref_file):
    done_file = os.path.join(outdir, 'downloaded')

    if os.path.exists(done_file):
        print('Extra references already downloaded. Skipping.')
        return

    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    try:
        os.mkdir(outdir)
    except:
        raise Error('Error mkdir ' + outdir)

    if extra_ref_file is None:
        iva.common.syscall('touch ' + done_file)
        return 

    f_ids = fastaq.utils.open_file_read(extra_ref_file)
    for line in f_ids:
        tmpdir = tempfile.mkdtemp(prefix='tmp.get_files.', dir=outdir)
        gis = line.rstrip().split()

        for gi in gis:
            download_from_genbank(os.path.join(tmpdir, gi + '.fasta'), 'fasta', gi)
            download_from_genbank(os.path.join(tmpdir, gi + '.gb'), 'gb', gi)
     
        f = fastaq.utils.open_file_read(os.path.join(tmpdir, gis[0] + '.gb'))
        oganism = None
        for line in f:
            if line.startswith('  ORGANISM'):
                organism = line.split(maxsplit=1)[1].rstrip()
                break
        if organism is None:
            raise Error('Error getting ORGANISM line from file ', os.path.join(tmpdir, gis[0] + '.gb'))
        fastaq.utils.close(f)
        dirname = re.sub('\W', '_', organism)
        os.rename(tmpdir, os.path.join(outdir, dirname))
        
    fastaq.utils.close(f_ids)
    iva.common.syscall('touch ' + done_file)


def build_kraken_virus_db(outdir, threads=1, minimizer_len=13, max_db_size=4, extra_ref_file=None):
    extra_ref_dir = os.path.join(outdir, 'Extra_ref')
    download_extra_refs(extra_ref_dir, extra_ref_file)
    add_extra_files_script = os.path.join(outdir, 'add-to-library.sh')

    f = fastaq.utils.open_file_write(add_extra_files_script)
    print('set -e', file=f)
    if extra_ref_file is None:
        print('echo "No extra reference files to download. Skipping"', file=f)
    else:
        print('for f in `find ' + extra_ref_dir + "'*.fasta'`;",
              'do',
              '    kraken-build --add-to-library $f --db ' + outdir,
              '    rm $f',
              'done', sep='\n', file=f)
    fastaq.utils.close(f)
    
    tasks = [
        ('kraken-build --download-taxonomy --db ' + outdir, 'progress.kraken-build.taxonomy.done'),
        ('kraken-build --download-library viruses --db ' + outdir, 'progress.kraken-build.library.done'),
        ('bash ' + add_extra_files_script, 'progress.add-extra-to-library.done'),
        (' '.join([
            'kraken-build --build',
            '--db', outdir,
            '--max-db-size', str(max_db_size),
            '--minimizer-len', str(minimizer_len),
            '--threads', str(threads)
        ]), 'progress.kraken-build.build.done'),
        ('kraken-build --clean --db ' + outdir, 'progress.kraken-build.clean.done')
    ]

    for cmd, done_file in tasks:
        if os.path.exists(done_file):
            print('Skipping ', cmd, flush=True)
        else:
            print('Running:', cmd, flush=True)
            iva.common.syscall(cmd)
            iva.common.syscall('touch ' + done_file)


def get_genbank_virus_files(outdir):
    try:
        os.mkdir(outdir)
    except:
        raise Error('Error mkdir ' + outdir)

    cwd = os.getcwd()
    os.chdir(outdir)
    cmd = 'wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.gbk.tar.gz'
    iva.common.syscall('wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.gbk.tar.gz')
    iva.common.syscall('tar -xf all.gbk.tar.gz')
    os.chdir(cwd)


def download_from_genbank(outfile, filetype, gi):
    assert filetype in ['gb', 'fasta']
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=' + filetype + '&retmode=text&id=' + gi
    try:
        file_contents = urllib.request.urlopen(url).read().decode().rstrip()
    except:
        raise Error('Error downloading gi ' + gi + ' using:\n' + url)

    f = fastaq.utils.open_file_write(outfile)
    print(file_contents, file=f)
    fastaq.utils.close(f)


def make_embl_files(indir):
    original_dir = os.getcwd()
    this_module_dir =os.path.dirname(inspect.getfile(inspect.currentframe()))
    genbank2embl = os.path.abspath(os.path.join(this_module_dir, 'ratt', 'genbank2embl.pl'))

    for directory, x, gbk_files in sorted(os.walk(indir)):
        if directory == indir or len(gbk_files) == 0:
            continue
        os.chdir(directory)
        print('Converting', directory, end=' ', flush=True)
        for fname in gbk_files:
            # some genbank files have a 'CONTIG' line, which breaks bioperl's
            # conversion genbank --> embl and makes genbank2embl.pl hang
            iva.common.syscall('grep -v CONTIG ' + fname + ' > tmp.gbk; mv tmp.gbk ' + fname)
            iva.common.syscall(genbank2embl + ' ' + fname + ' ' + fname + '.embl')
            os.unlink(fname)
            print(fname, end=' ', flush=True)

        os.chdir(original_dir)
        print()

    os.chdir(original_dir)


def setup_ref_db(outdir, threads=1, minimizer_len=13, max_db_size=3, extra_ref_file=None):
    final_done_file = 'progress.setup.done'
    embl_done_file = 'progess.embl.done'
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        try:
            os.mkdir(outdir)
        except:
            raise Error('Error mkdir ' + outdir)

    cwd = os.getcwd()

    try:
        os.chdir(outdir)
    except:
        raise Error('Error chdir ' + outdir)

    if os.path.exists(final_done_file):
        print('Nothing to do! All files made already')
        return

    kraken_dir = 'Kraken_db'
    embl_dir = 'Library'
    extra_ref_dir = 'Extra_ref'

    if os.path.exists('progress.kraken-build.clean.done'):
        print('Kraken database already made. Skipping.')
    else:
        build_kraken_virus_db(kraken_dir, threads=threads, minimizer_len=minimizer_len, max_db_size=max_db_size, extra_ref_file=extra_ref_file)


    if os.path.exists(embl_done_file):
        print('EMBL/fasta files made already. Skipping.', flush=True)
    else:
        if os.path.exists(embl_dir):
            shutil.rmtree(embl_dir)
        print('Downloading virus genbank files', flush=True)
        get_genbank_virus_files(embl_dir)
        print('...finished. Making EMBL and fasta files', flush=True)
        make_embl_files(embl_dir)
        make_embl_files(extra_ref_dir)
        iva.common.syscall('mv ' + extra_ref_dir + '/* ' + embl_dir)
        os.unlink(os.path.join(embl_dir, 'done'))
        shutil.rmtree(extra_ref_dir)
        iva.common.syscall('touch ' + embl_done_file)

    iva.common.syscall('touch ' + final_done_file)
    os.chdir(cwd)


def run_kraken(database, reads, outprefix, preload=True, threads=1, mate_reads=None):
    kraken_out = outprefix + '.out'
    kraken_report = outprefix + '.report'
    cmd = ' '.join([
        'kraken',
        '--db', database,
        '--preload' if preload else '',
        '--threads', str(threads),
        '--output', kraken_out
    ])

    if mate_reads is not None:
        cmd += ' --paired ' + reads + ' ' + mate_reads
    else:
        cmd += ' ' + reads

    iva.common.syscall(cmd)

    cmd = ' '.join([
        'kraken-report',
        '--db', database,
        kraken_out,
        '>', kraken_report
    ])
    iva.common.syscall(cmd)
    os.unlink(kraken_out)


def get_most_common_species(kraken_out):
    f = fastaq.utils.open_file_read(kraken_out)
    max_count = -1
    species = None
    for line in f:
        a = line.rstrip().split()
        if a[3] == "S" and int(a[2]) > max_count:
            species = '_'.join(a[5:])
            max_count = int(a[2])

    fastaq.utils.close(f)
    return species.replace('-', '_')


def _dir_from_species(lib_dir, species):
    dirs = [x for x in os.listdir(lib_dir) if x.rsplit('_', 1)[0] == species]
    if len(dirs) == 0:
        raise Error('Error finding EMBL files directory from species "' + species + '". Cannot continue')
    elif len(dirs) > 1:
        print('Warning: >1 directory for species "' + species + '".\n' + dirs + '\n. Using ' + dirs[0])
    return os.path.join(lib_dir, dirs[0])


def choose_reference(root_dir, reads, outprefix, preload=True, threads=1, mate_reads=None):
    root_dir = os.path.abspath(root_dir)
    kraken_db = os.path.join(root_dir, 'Kraken_db')
    run_kraken(kraken_db, reads, outprefix, preload=preload, threads=threads, mate_reads=None)
    species = get_most_common_species(outprefix + '.report')
    lib_dir = os.path.join(root_dir, 'Library')
    return _dir_from_species(lib_dir, species)


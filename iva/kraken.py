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
import stat
import inspect
import sys
import os
import tempfile
import shutil
import time
import re
import urllib.request
import iva
import pyfastaq


class Error (Exception): pass


class Database:
    def __init__(self, rootdir, extra_refs_file=None, threads=1, minimizer_len=13, max_db_size=3, preload=False, verbose=False, skip_virus_download=False):
        self.rootdir = os.path.abspath(rootdir)
        if extra_refs_file is None:
            self.extra_refs_file = None
        else:
            self.extra_refs_file = os.path.abspath(extra_refs_file)

        self.threads = threads
        self.minimizer_len = minimizer_len
        self.max_db_size = max_db_size
        self.current_taxon_id = 2000000000
        self.current_gi = 4000000000
        self.preload = preload
        self.verbose = verbose
        self.taxon_to_parent = {}
        self.extra_refs = {} # keys =  new fake taxon ID
        self.reseq_virus_dir = os.path.join(self.rootdir, 'Virus_refseq')
        self.kraken_db = os.path.join(self.rootdir, 'Kraken_db')
        self.kraken_taxon_dir = os.path.join(self.kraken_db, 'taxonomy')
        self.kraken_names_dmp = os.path.join(self.kraken_taxon_dir, 'names.dmp')
        self.kraken_nodes_dmp = os.path.join(self.kraken_taxon_dir, 'nodes.dmp')
        self.kraken_gi_taxid_nucl_dmp = os.path.join(self.kraken_taxon_dir, 'gi_taxid_nucl.dmp')
        self.embl_root = os.path.join(self.rootdir, 'EMBL')
        self.extra_refs_dir = os.path.join(self.rootdir, 'Extra_refs')
        self.skip_virus_download = skip_virus_download
        self.added_to_kraken = set()

        if self.skip_virus_download and extra_refs_file is None:
            raise Error('Cannot create database. skip_virus_download is True and no extra_refs_file provided')

        self.tasks = [
            'download',
            'build',
            'make_embl',
            'clean',
        ]

        self.done_files = {x:os.path.join(self.rootdir, 'progress.' + x + '.done') for x in self.tasks}


    def _mkdir(self, d, rmtree=False):
        if rmtree and os.path.exists(d):
            shutil.rmtree(d)

        if not os.path.exists(d):
            try:
                os.mkdir(d)
            except:
                raise Error('Error mkdir ' + d)


    def _get_parent_taxons(self, taxons):
        f = pyfastaq.utils.open_file_read(self.kraken_nodes_dmp)
        for line in f:
            a = line.split()
            if a[0] in taxons:
                self.taxon_to_parent[a[0]] = a[2]
        pyfastaq.utils.close(f)


    def _load_extra_ref_info(self):
        if self.extra_refs_file is None:
             return

        f = pyfastaq.utils.open_file_read(self.extra_refs_file)
        for line in f:
            genbank_ids = line.rstrip().split()
            new_gis = list(range(self.current_gi, self.current_gi + len(genbank_ids)))
            self.current_gi += len(genbank_ids)
            assert len(genbank_ids) == len(new_gis)
            self.extra_refs[self.current_taxon_id] = {
                'genbank_ids': genbank_ids,
                'new_gis': new_gis,
            }

            self.current_taxon_id += 1
        pyfastaq.utils.close(f)


    def _download_from_genbank(self, outfile, filetype, gi, max_tries=5, delay=3):
        assert filetype in ['gb', 'fasta']
        file_ok = False

        for i in range(max_tries):
            try:
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=' + filetype + '&retmode=text&id=' + gi
                file_contents = iva.common.decode(urllib.request.urlopen(url).read()).rstrip()
            except:
                time.sleep(delay)
                continue

            if (filetype == 'fasta' and file_contents.startswith('>')) \
               or (filetype == 'gb' and file_contents.startswith('LOCUS')):
                file_ok = True
                break
            else:
                time.sleep(delay)

        if not file_ok:
            raise Error('Error downloading gi ' + gi + ' using:\n' + url + '\nI got this:\n' + file_contents)

        f = pyfastaq.utils.open_file_write(outfile)
        print(file_contents, file=f)
        pyfastaq.utils.close(f)


    def _download_extra_refs(self):
        if self.extra_refs_file is None:
             return

        for l in [x['genbank_ids'] for x in self.extra_refs.values()]:
            for gi in l:
                outprefix = os.path.join(self.extra_refs_dir, gi)
                self._download_from_genbank(outprefix + '.gb', 'gb', gi)
                self._download_from_genbank(outprefix + '.fasta', 'fasta', gi)


    def _genbank_to_taxon_and_gi(self, filename):
        f = pyfastaq.utils.open_file_read(filename)
        taxon_id = None
        gi = None
        for line in f:
            if line.startswith('                     /db_xref="taxon:'):
                taxon_id = line.rstrip().split(':')[-1].rstrip('"')
            elif line.startswith('VERSION'):
                gi = line.rstrip().split()[-1].split(':')[-1]
            if None not in [taxon_id, gi]:
                break
        pyfastaq.utils.close(f)
        if None in [taxon_id, gi]:
            raise Error('Error getting taxon or GI from ' + filename)
        return taxon_id, gi


    def _genbank2embl(self, infile, outfile):
        tmpdir = tempfile.mkdtemp(prefix='tmp.genbank2embl.', dir=os.getcwd())
        extractor = iva.egg_extract.Extractor(os.path.abspath(os.path.join(os.path.dirname(iva.__file__), os.pardir)))
        genbank2embl_egg = os.path.join('iva', 'ratt', 'genbank2embl.pl')
        genbank2embl = os.path.join(tmpdir, 'genbank2embl.pl')
        extractor.copy_file(genbank2embl_egg, genbank2embl)
        os.chmod(genbank2embl, stat.S_IRWXU)

        # some genbank files have a 'CONTIG' line, which breaks bioperl's
        # conversion genbank --> embl and makes genbank2embl.pl hang
        iva.common.syscall('grep -v CONTIG ' + infile + ' > tmp.gbk; mv tmp.gbk ' + infile)
        iva.common.syscall(genbank2embl + ' ' + infile + ' ' + outfile, verbose=self.verbose)
        shutil.rmtree(tmpdir)


    def _append_to_file(self, filename, line):
        try:
            f = open(filename, 'a')
        except:
            raise Error('Error opening for appending:', filename)
        print(line, file=f)
        f.close()


    def _add_to_kraken(self, fa_file, real_taxon, new_taxon, new_gi):
        parent_taxon = self.taxon_to_parent[real_taxon]
        if self.verbose:
            print('add_to_kraken', fa_file, real_taxon, new_taxon, parent_taxon, new_gi, sep='\t')
        # update names.dmp: append new_taxon with whatever name that kraken will report
        # update nodes.dmp: append new_taxon with original parent taxon
        # update gi_taxid_nucl.dmp: append new_gi, new_taxon
        if new_taxon not in self.added_to_kraken:
            line = str(new_taxon) + '\t|\tadded.' + str(new_taxon) +  '\t|\t\t|\tscientific name\t|'
            self._append_to_file(self.kraken_names_dmp, line)
            line = '\t|\t'.join([
                str(new_taxon),
                parent_taxon,
                'species',
                'HI',
                '9',
                '1',
                '1',
                '1',
                '0',
                '1',
                '1',
                '0',
                '',
            ]) + '\t|'
            self._append_to_file(self.kraken_nodes_dmp, line)
        self.added_to_kraken.add(new_taxon)
        self._append_to_file(self.kraken_gi_taxid_nucl_dmp, str(new_gi) + '\t' + str(new_taxon))
        iva.common.syscall('kraken-build --add-to-library ' + fa_file + ' --db ' + self.kraken_db, verbose=self.verbose)


    def _replace_fasta_header(self, filename, header):
        fin = pyfastaq.utils.open_file_read(filename)
        tmp_out = filename + '.tmp.fa'
        fout = pyfastaq.utils.open_file_write(tmp_out)

        for line in fin:
            if line.startswith('>'):
                print('>' + header, file=fout)
            else:
                print(line.rstrip(), file=fout)

        pyfastaq.utils.close(fin)
        pyfastaq.utils.close(fout)
        os.rename(tmp_out, filename)


    def _sort_out_extra_refs(self):
        if self.extra_refs_file is None:
             return
        real_taxon_ids = set()

        for new_taxon_id, info in self.extra_refs.items():
            gb_file = os.path.join(self.extra_refs_dir, info['genbank_ids'][0] + '.gb')
            info['real_taxon'], info['gi'] = self._genbank_to_taxon_and_gi(gb_file)
            real_taxon_ids.add(info['real_taxon'])
            embl_dir = os.path.join(self.embl_root, 'added.' + str(new_taxon_id))
            self._mkdir(embl_dir)
            for i in range(len(info['genbank_ids'])):
                gi = info['genbank_ids'][i]
                new_gi = info['new_gis'][i]
                gb_file = os.path.join(self.extra_refs_dir, gi + '.gb')
                fa_file = os.path.join(self.extra_refs_dir, gi + '.fasta')
                self._replace_fasta_header(fa_file, 'gi|' + str(new_gi) + '|x')
                embl_file = os.path.join(embl_dir, gi + '.embl')
                self._genbank2embl(gb_file, embl_file)

        self._get_parent_taxons(real_taxon_ids)

        for new_taxon_id, info in  self.extra_refs.items():
            real_taxon = info['real_taxon']
            for i in  range(len(info['genbank_ids'])):
                gi = info['genbank_ids'][i]
                new_gi = info['new_gis'][i]
                fa_file = os.path.join(self.extra_refs_dir, gi + '.fasta')
                gb_file = os.path.join(self.extra_refs_dir, gi + '.gb')
                self._add_to_kraken(fa_file, real_taxon, new_taxon_id, new_gi)
                if self.verbose:
                    print('unlink', os.path.exists(gb_file), gb_file)
                    print('unlink', os.path.exists(fa_file), fa_file)
                os.unlink(gb_file)
                os.unlink(fa_file)


    def _build_kraken_virus_db(self):
        if os.path.exists(self.done_files['clean']):
            print('Kraken DB already built. Skipping.', flush=True)
            return

        if os.path.exists(self.done_files['download']):
            print('Files already downloaded. Skipping.', flush=True)
        else:
            for d in [self.rootdir, self.embl_root, self.extra_refs_dir]:
                self._mkdir(d, rmtree=True)

            iva.common.syscall('kraken-build --download-taxonomy --db ' + self.kraken_db, verbose=self.verbose)
            if not self.skip_virus_download:
                iva.common.syscall('kraken-build --download-library viruses --db ' + self.kraken_db, verbose=self.verbose)

            if self.extra_refs_file is not None:
                self._load_extra_ref_info()
                self._download_extra_refs()
                self._sort_out_extra_refs()
            iva.common.syscall('touch ' + self.done_files['download'], verbose=self.verbose)

        if os.path.exists(self.done_files['build']):
            print('Already built. Skipping.', flush=True)
        else:
            cmd = ' '.join([
                'kraken-build --build',
                '--db', self.kraken_db,
                '--max-db-size', str(self.max_db_size),
                '--minimizer-len', str(self.minimizer_len),
                '--threads', str(self.threads)
            ])
            iva.common.syscall(cmd, verbose=self.verbose)
            iva.common.syscall('touch ' + self.done_files['build'])


    def _clean(self):
        if os.path.exists(self.done_files['clean']):
            print('Already cleaned. Skipping.')
        else:
            iva.common.syscall('kraken-build --clean --db ' + self.kraken_db, verbose=self.verbose)
            if os.path.exists(self.extra_refs_dir):
                shutil.rmtree(self.extra_refs_dir)
            iva.common.syscall('touch ' + self.done_files['clean'], verbose=self.verbose)


    def _get_genbank_virus_files(self):
        if os.path.exists(self.done_files['make_embl']):
            print('Already made EMBL files. Skipping.')
            return
        self._mkdir(self.reseq_virus_dir, rmtree=True)
        cwd = os.getcwd()
        os.chdir(self.reseq_virus_dir)
        iva.common.syscall('wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.gbk.tar.gz', verbose=self.verbose)
        iva.common.syscall('tar -xf all.gbk.tar.gz', verbose=self.verbose)
        os.chdir(cwd)
        self._convert_refseq_virus_to_embl()
        shutil.rmtree(self.reseq_virus_dir)
        iva.common.syscall('touch ' + self.done_files['make_embl'])


    def _convert_refseq_virus_to_embl(self):
        original_dir = os.getcwd()

        for directory, x, gbk_files in sorted(os.walk(self.reseq_virus_dir)):
            if directory == self.reseq_virus_dir or len(gbk_files) == 0:
                continue
            os.chdir(directory)
            if self.verbose:
                print('Converting', directory, end=' ', flush=True)
            for fname in gbk_files:
                self._genbank2embl(fname, re.sub('\.gbk$', '', fname) + '.embl')
                os.unlink(fname)
                if self.verbose:
                    print(fname, end=' ', flush=True)

            os.chdir(original_dir)
            if self.verbose:
                print()
            new_dir = re.sub('_uid[0-9]+$', '', directory).strip('_')
            if new_dir != directory:
                os.rename(directory, new_dir)

            final_dir =  os.path.join(self.embl_root, os.path.basename(new_dir))
            if os.path.exists(final_dir):
                shutil.rmtree(final_dir)
            os.rename(new_dir, final_dir)

        os.chdir(original_dir)


    def build(self):
        self._build_kraken_virus_db()
        if not self.skip_virus_download:
            self._get_genbank_virus_files()
        self._clean()


    def _run_kraken(self, reads, outprefix):
        kraken_out = outprefix + '.out'
        kraken_report = outprefix + '.report'
        cmd = ' '.join([
            'kraken',
            '--db', self.kraken_db,
            '--preload' if self.preload else '',
            '--threads', str(self.threads),
            '--output', kraken_out,
            reads
        ])

        iva.common.syscall(cmd, verbose=self.verbose)

        cmd = ' '.join([
            'kraken-report',
            '--db', self.kraken_db,
            kraken_out,
            '>', kraken_report
        ])
        iva.common.syscall(cmd, verbose=self.verbose)
        os.unlink(kraken_out)


    def _species_to_embl_dir(self, s):
        if s.startswith('added.'):
            return s
        else:
            return re.sub('\W', '_', s).strip('_')


    def _get_most_common_species_dir(self, kraken_report):
        embl_dirs = set([os.path.basename(x[0]) for x in os.walk(self.embl_root) if x[0] != self.embl_root])
        f = pyfastaq.utils.open_file_read(kraken_report)
        max_count = -1
        most_common_dir = None
        for line in f:
            a = line.rstrip().split('\t')
            this_species = a[-1].strip()
            this_dir = self._species_to_embl_dir(this_species)
            if this_dir in embl_dirs and int(a[2]) > max_count:
                most_common_dir = this_dir
                max_count = int(a[2])

        pyfastaq.utils.close(f)
        if most_common_dir is None:
            return None
        else:
            return os.path.join(self.embl_root, most_common_dir)


    def choose_reference(self, reads, outprefix):
        self._run_kraken(reads, outprefix)
        lib_dir = os.path.join(self.rootdir, 'Library')
        return self._get_most_common_species_dir(outprefix + '.report')


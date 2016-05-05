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
import tempfile
import os
import sys
import shutil
import multiprocessing
from iva import mapping, seed
import pyfastaq

class Error (Exception): pass

class SeedProcessor:
    def __init__(self, seeds_fasta, reads1, reads2, outfile, index_k=15, index_s=3, threads=1, max_insert=500, minid=0.9, seed_stop_length=500, extend_length=50, overlap_length=None, ext_min_cov=5, ext_min_ratio=2, verbose=0, seed_length=None, seed_min_count=10, seed_max_count=100000000, kmc_threads=1):
        self.seeds_fasta = seeds_fasta
        self.reads1 = reads1
        self.reads2 = reads2
        self.outfile = outfile
        self.index_k = index_k
        self.index_s = index_s
        self.threads = threads
        self.kmc_threads = kmc_threads
        self.max_insert = max_insert
        self.minid = minid
        self.seed_stop_length = seed_stop_length
        self.extend_length = extend_length
        self.overlap_length = overlap_length
        self.ext_min_cov = ext_min_cov
        self.ext_min_ratio = ext_min_ratio
        self.verbose = verbose
        self.seed_length = seed_length
        self.seed_min_count = seed_min_count
        self.seed_max_count = seed_max_count

        self.original_seeds = {}
        pyfastaq.tasks.file_to_dict(seeds_fasta, self.original_seeds)
        self.processed_seeds = {}
        self.tmpdir = None
        self.bam_file = None


    def _make_new_seed(self, seed_name):
        if self.verbose:
            print('Making new seed for', seed_name, ' ... start')
        tmp_prefix = os.path.join(self.tmpdir, 'out')
        seed_reads = tmp_prefix + '.' + seed_name + '.reads_1.fa'
        if len(self.original_seeds[seed_name]) > self.seed_stop_length:
            start = int(0.5 * len(self.original_seeds[seed_name]) - 0.5 * self.seed_stop_length)
            end = int(0.5 * len(self.original_seeds[seed_name]) + 0.5 * self.seed_stop_length)
        else:
            start = None
            end = None
        if self.verbose:
            print('Making new seed for', seed_name, ' ... getting reads')
        mapping.bam_file_to_region_fasta(self.bam_file, seed_reads, seed_name, start, end)
        if self.verbose:
            print('Making new seed for', seed_name, ' ... finding most common kmer')
        new_seed = seed.Seed(
            extend_length = self.extend_length,
            overlap_length = self.overlap_length,
            reads1 = seed_reads,
            ext_min_cov = self.ext_min_cov,
            ext_min_ratio = self.ext_min_ratio,
            verbose = self.verbose,
            seed_length = self.seed_length,
            seed_min_count = self.seed_min_count,
            seed_max_count = self.seed_max_count,
            kmc_threads = self.kmc_threads,
            map_threads = self.threads
        )
        if len(new_seed) == 0:
            print('Warning: could not get most common kmer for', seed_name)
            return

        if self.verbose:
            print('Making new seed for', seed_name, ' ... extending most common kmer')

        new_seed.extend(self.reads1, self.reads2, self.seed_stop_length)
        f = pyfastaq.utils.open_file_write(tmp_prefix + '.' + seed_name + '.fa')
        print(pyfastaq.sequences.Fasta('seed.' + seed_name, new_seed.seq[10:-10]), file=f)
        pyfastaq.utils.close(f)
        if self.verbose:
            print('Making new seed for', seed_name, ' ... finished')


    def process(self):
        self.tmpdir = tempfile.mkdtemp(prefix='tmp.process_seeds.', dir=os.getcwd())
        tmp_prefix = os.path.join(self.tmpdir, 'out')
        mapping.map_reads(self.reads1, self.reads2, self.seeds_fasta, tmp_prefix,
            index_k = self.index_k,
            index_s = self.index_s,
            threads = self.threads,
            max_insert = self.max_insert,
            minid = self.minid,
            sort = True)
        self.bam_file = tmp_prefix + '.bam'
        threads = min(8, self.threads) # to save peak memory going too high
        threads = self.threads

        if self.verbose:
            print('Processing seeds with', threads, 'threads:', list(self.original_seeds.keys()))

        pool = multiprocessing.Pool(threads)
        pool.map(self._make_new_seed, list(self.original_seeds.keys()))
        pool.close()
        pool.join()
        if self.verbose:
            print('... finished processing seeds')

        new_seeds = {}
        for seed_name in self.original_seeds:
            fname = tmp_prefix + '.' + seed_name + '.fa'
            if os.path.exists(fname):
                pyfastaq.tasks.file_to_dict(fname, new_seeds)

        if len(new_seeds) == 0:
            raise Error('Error! did not make any new seeds. Cannot continue')
        f = pyfastaq.utils.open_file_write(self.outfile)
        for seq in new_seeds.values():
            print(seq, file=f)
        pyfastaq.utils.close(f)
        shutil.rmtree(self.tmpdir)


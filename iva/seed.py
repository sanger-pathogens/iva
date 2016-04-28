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
import copy
import os
import shutil
import tempfile
import pyfastaq
from iva import kcount, kmers, mapping

class Error (Exception): pass

class Seed:
    def __init__(self, extend_length=50, overlap_length=None, reads1=None, reads2=None, seq=None, ext_min_cov=5, ext_min_ratio=2, verbose=0, seed_length=None, seed_min_count=10, seed_max_count=100000000, kmc_threads=1, map_threads=1, sequences_to_ignore=None, contigs_to_check=None):
        if contigs_to_check is None:
            contigs_to_check = {}
        if sequences_to_ignore is None:
            sequences_to_ignore = set()
        self.verbose = verbose
        self.kmc_threads = kmc_threads
        self.map_threads = map_threads
        self.extend_length = extend_length
        self.ext_min_cov = ext_min_cov
        self.ext_min_ratio = ext_min_ratio
        self.seed_lengths = []
        self.overlap_length = overlap_length
        if seq is None:
            if reads1 is None:
                raise Error('Cannot construct Seed object. Need reads when no seq has been given')
            kmer_counts = kcount.get_most_common_kmers(reads1, reads2, most_common=1, min_count=seed_min_count, max_count=seed_max_count, kmer_length=seed_length, verbose=self.verbose, ignore_seqs=sequences_to_ignore, contigs_to_check=contigs_to_check, kmc_threads=self.kmc_threads, map_threads=self.map_threads)
            if len(kmer_counts) == 1:
                self.seq = list(kmer_counts.keys())[0]
                if self.verbose:
                    print('Made new seed. kmer coverage', list(kmer_counts.values())[0], 'and seed is', self.seq, flush=True)
            else:
                self.seq = None
        else:
            self.seq = seq


        if self.seq is not None:
            if overlap_length is None:
                self.overlap_length = len(self.seq)
            else:
                self.overlap_length = overlap_length
        else:
            self.overlap_length = None


    def __len__(self):
        if self.seq is None:
            return 0
        else:
            return len(self.seq)


    def _extension_from_read(self, read, left=False):
        '''Returns an extension sequence for a read, or None if extension not possible'''
        if left:
            end_seq = self.seq[:self.overlap_length]
        else:
            end_seq = self.seq[-self.overlap_length:]

        hits = read.search(end_seq)
        if len(hits) == 1:
            hit_start, strand  = hits[0]
            hit_end = hit_start + len(end_seq) - 1
            if strand == '-':
                seq = copy.copy(read)
                seq.revcomp()
                hit_start, hit_end = len(seq) - hit_end - 1, len(seq) - hit_start - 1
            else:
                seq = read

            if left and hit_start > 0:
                return seq[0:hit_start]
            if (not left) and hit_end + 1 < len(read):
                return seq[hit_end + 1:]

        return None


    def _extensions_from_reads_file(self, reads_file):
        seq_reader = pyfastaq.sequences.file_reader(reads_file)
        left_seqs = []
        right_seqs = []
        for seq in seq_reader:
            left_ext = self._extension_from_read(seq, left=True)
            right_ext = self._extension_from_read(seq, left=False)

            if left_ext is not None:
                left_seqs.append(left_ext)

            if right_ext is not None:
                right_seqs.append(right_ext)

        return left_seqs, right_seqs


    def _extend_with_reads_as_single_end(self, reads1, reads2):
        left_seqs, right_seqs = self._extensions_from_reads_file(reads1)
        left_seqs2, right_seqs2 = self._extensions_from_reads_file(reads2)
        left_kmers = kmers.Kmers(verbose=self.verbose, left=True)
        right_kmers = kmers.Kmers(verbose=self.verbose, left=False)
        left_kmers.extend(left_seqs)
        left_kmers.extend(left_seqs2)
        right_kmers.extend(right_seqs)
        right_kmers.extend(right_seqs2)
        left_extension = left_kmers.extension(self.ext_min_cov, self.ext_min_ratio, self.extend_length)
        right_extension = right_kmers.extension(self.ext_min_cov, self.ext_min_ratio, self.extend_length)
        self.seed_lengths.append((len(self.seq), len(left_extension), len(right_extension)))
        if self.verbose >= 2:
            print('    Extend seed. new length=', len(self.seq), '. Bases added left:', len(left_extension), '. Bases added right:', len(right_extension), sep='')
        self.seq = left_extension + self.seq + right_extension
        if self.verbose >= 2:
            print('                 new seed:', self.seq)


    def extend(self, reads1, reads2, stop_length):
        '''Extend the seed sequence up to stop_length. Uses reads in reads1, reads2 but treats them as unpaired'''
        iteration = 1
        while len(self.seq) < stop_length and \
               (len(self.seed_lengths) == 0 or self.seed_lengths[-1][1] != 0 or self.seed_lengths[-1][2] != 0):
            if self.verbose:
                print('Seed extension iteration', iteration)
            self._extend_with_reads_as_single_end(reads1, reads2)
            iteration += 1


    def write_fasta(self, filename, name):
        f = pyfastaq.utils.open_file_write(filename)
        print('>' + name, file=f)
        print(self.seq, file=f)
        pyfastaq.utils.close(f)



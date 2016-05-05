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
from iva import kmers

class Contig:
    def __init__(self, fasta, verbose=0):
        self.fa = fasta
        self.verbose = verbose
        self.left_kmers = kmers.Kmers(left=True, verbose=verbose)
        self.right_kmers = kmers.Kmers(verbose=verbose)


    def __len__(self):
        return len(self.fa)


    def add_left_kmer(self, kmer):
        self.left_kmers.append(kmer)
        

    def add_right_kmer(self, kmer):
        self.right_kmers.append(kmer)


    def extend(self, min_cov, min_ratio, extend_bases):
        if self.verbose >= 2:
            print('    trying to extend left end ...')
        new_left_seq = self.left_kmers.extension(min_cov, min_ratio, extend_bases)
        if self.verbose >= 2:
            print('    trying to extend right end ...')
        new_right_seq = self.right_kmers.extension(min_cov, min_ratio, extend_bases)
        if self.verbose >= 1:
            print('    new left sequence:', new_left_seq)
            print('    new right sequence:', new_right_seq)
        self.fa.seq = new_left_seq + self.fa.seq + new_right_seq
        self.left_kmers = kmers.Kmers(left=True, verbose=self.verbose)
        self.right_kmers = kmers.Kmers(verbose=self.verbose)
        return len(new_left_seq), len(new_right_seq)

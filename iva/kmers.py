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
from collections import Counter

class Kmers:
    def __init__(self, kmer=None, left=False, verbose=0):
        self.kmers = []
        self.left = left
        self.verbose = verbose
        if kmer is not None:
            self.append(kmer)


    def append(self, kmer):
        '''Adds one kmer'''
        if 'N' in kmer:
            return
        if self.left:
            self.kmers.append(kmer[::-1])
        else:
            self.kmers.append(kmer)


    def extend(self, l):
        '''Addds a list of kmers'''
        for x in l:
            self.append(x)


    def _kmer_dict(self, k):
        return dict(Counter([kmer[0:k] for kmer in self.kmers if len(kmer) >= k]))


    def _commonest_kmers(self, k):
        counts = self._kmer_dict(k)
        highest_counts = None, None
        highest_kmers = None, None

        for k in counts:
            if highest_counts[1] is None:
                highest_kmers = None, k
                highest_counts = None, counts[k]
            elif counts[k] > highest_counts[1]:
                highest_kmers = highest_kmers[1], k
                highest_counts = highest_counts[1], counts[k]
            elif highest_counts[0] is None:
                highest_kmers = k, highest_kmers[1]
                highest_counts = counts[k], highest_counts[1]
        return highest_kmers, highest_counts


    def extension(self, min_cov, min_ratio, extend_length):
        '''Returns a string that can extend the kmer, of length as long as possible up to extend_length bases. Must be supported by min_cov kmers, and also it must be at least min_ratio more abundant than the second most common kmer'''
        if self.verbose >= 4:
            print('        all kmers:')
            for kmer in self.kmers:
                print('            ', kmer)
        if len(self.kmers) == 0:
            return ''
        longest = max([len(x) for x in self.kmers])
        for i in range(min(longest, extend_length),0, -1):
            highest_kmers, highest_counts = self._commonest_kmers(i)
            if self.verbose >= 2:
                print('        k =', i, 'commonest two kmers:', highest_kmers, 'have frequency:', highest_counts)
            if highest_kmers[1] is None:
                pass
            elif highest_kmers[0] is None:
                if highest_counts[1] >= min_cov:
                    extension_seq = highest_kmers[1]
                    return highest_kmers[1][::-1] if self.left else highest_kmers[1]
            else:
                assert highest_counts[0] > 0 and highest_counts[1] > 0
                if highest_counts[1] / highest_counts[0] >= min_ratio and highest_counts[1] >= min_cov:
                    return highest_kmers[1][::-1] if self.left else highest_kmers[1]

        return ''

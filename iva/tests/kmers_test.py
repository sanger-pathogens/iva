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
import unittest
import os
from iva import kmers

class TestKmers(unittest.TestCase):
    def test_init(self):
        '''Test __init__'''
        k = kmers.Kmers(kmer='ACGGT')
        self.assertListEqual(['ACGGT'], k.kmers)
        k = kmers.Kmers(kmer='ACGGT', left=True)
        self.assertListEqual(['TGGCA'], k.kmers)


    def test_append(self):
        '''Test append'''
        k = kmers.Kmers(kmer='A')
        k.append('G')
        self.assertListEqual(['A', 'G'], k.kmers)
        k.append('N')
        self.assertListEqual(['A', 'G'], k.kmers)


    def test_extend(self):
        '''Test extend'''
        k = kmers.Kmers(kmer='A')
        k.extend(['A', 'C'])
        self.assertListEqual(['A', 'A', 'C'], k.kmers)


    def test_kmer_dict(self):
        '''Test kmer_dict'''
        k = kmers.Kmers(kmer='A')
        k.extend(['A', 'C', 'GT', 'GC', 'ACT'])
        expected = [
            {'A': 3, 'C': 1, 'G': 2},
            {'GT': 1, 'GC': 1, 'AC': 1},
            {'ACT': 1}
        ]

        for i in range(3):
            d = k._kmer_dict(i+1)
            self.assertDictEqual(d, expected[i])


    def test_commonest_kmers(self):
        k = kmers.Kmers(kmer='AC')
        k.extend(['AT'] * 4)
        k.append('ATT')

        commonest_kmers, counts = k._commonest_kmers(1)
        self.assertEqual(commonest_kmers, (None, 'A'))
        self.assertEqual(counts, (None, 6))

        commonest_kmers, counts = k._commonest_kmers(2)
        self.assertEqual(commonest_kmers, ('AC', 'AT'))
        self.assertEqual(counts, (1, 5))

        commonest_kmers, counts = k._commonest_kmers(3)
        self.assertEqual(commonest_kmers, (None, 'ATT'))
        self.assertEqual(counts, (None, 1))

        commonest_kmers, counts = k._commonest_kmers(4)
        self.assertEqual(commonest_kmers, (None, None))
        self.assertEqual(counts, (None, None))


    def test_extension(self):
        '''Test extension'''
        k = kmers.Kmers(kmer='AT', left=True)
        self.assertEqual(k.extension(1, 2, 100), 'AT')
        k.extend(['GGT'] * 5)
        self.assertEqual(k.extension(5, 2, 100), 'GGT')
        self.assertEqual(k.extension(5, 2, 2), 'GT')

        k = kmers.Kmers()
        self.assertEqual(k.extension(5, 2, 100), '')
        k.extend(['ACG'] * 5)
        self.assertEqual(k.extension(5, 2, 100), 'ACG')
        self.assertEqual(k.extension(6, 2, 100), '')
        self.assertEqual(k.extension(5, 2, 2), 'AC')
        self.assertEqual(k.extension(5, 2, 1), 'A')

        k.extend(['ACC'] * 2)
        self.assertEqual(k.extension(5, 2, 100), 'ACG')
        self.assertEqual(k.extension(5, 3, 100), 'AC')


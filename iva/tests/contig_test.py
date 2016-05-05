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
#!/usr/bin/env python3

import unittest
import os
from iva import contig
from pyfastaq import sequences

class TestContig(unittest.TestCase):
    def test_init(self):
        '''Test init'''
        ctg = contig.Contig(sequences.Fasta('ID', 'ACCGT'))
        self.assertTrue(ctg.fa, sequences.Fasta('ID', 'ACCGT'))


    def test_len(self):
        '''Test len'''
        ctg = contig.Contig(sequences.Fasta('ID', 'ACCGT'))
        self.assertEqual(len(ctg), 5)


    def test_add_left_kmer(self):
        '''Test add_left_kmer'''
        ctg = contig.Contig(sequences.Fasta('ID', 'ACCGT'))
        ctg.add_left_kmer('GT')
        self.assertListEqual(ctg.left_kmers.kmers, ['TG'])


    def test_add_right_kmer(self):
        '''Test add_right_kmer'''
        ctg = contig.Contig(sequences.Fasta('ID', 'ACCGT'))
        ctg.add_right_kmer('GT')
        self.assertListEqual(ctg.right_kmers.kmers, ['GT'])


    def test_extend(self):
        '''Test extend'''
        ctg = contig.Contig(sequences.Fasta('ID', 'ACCGT'))
        self.assertEqual(ctg.extend(5, 2, 100), (0, 0))
        self.assertEqual(ctg.fa, sequences.Fasta('ID', 'ACCGT'))
        ctg.add_left_kmer('GT')
        self.assertEqual(ctg.extend(1, 2, 100), (2, 0))
        self.assertEqual(ctg.fa, sequences.Fasta('ID', 'GTACCGT'))
        self.assertEqual(ctg.extend(1, 2, 100), (0, 0))
        self.assertEqual(ctg.fa, sequences.Fasta('ID', 'GTACCGT'))
        ctg.add_right_kmer('TG')
        self.assertEqual(ctg.extend(1, 2, 100), (0, 2))
        self.assertEqual(ctg.fa, sequences.Fasta('ID', 'GTACCGTTG'))
        self.assertEqual(ctg.extend(1, 2, 100), (0, 0))
        self.assertEqual(ctg.fa, sequences.Fasta('ID', 'GTACCGTTG'))
        ctg.add_left_kmer('AG')
        ctg.add_right_kmer('GC')
        self.assertEqual(ctg.extend(1, 2, 100), (2, 2))
        self.assertEqual(ctg.fa, sequences.Fasta('ID', 'AGGTACCGTTGGC'))


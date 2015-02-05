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


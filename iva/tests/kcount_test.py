import unittest
import os
import filecmp
from iva import kcount

modules_dir = os.path.dirname(os.path.abspath(kcount.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
kcount.verbose = 2

class TestKcount(unittest.TestCase):
    def test_head_fastaq_fq(self):
        '''Test head_fastaq with FASTQ input'''
        test_out = 'tmp.head.fa'
        reads1 = os.path.join(data_dir, 'kcount_test.reads_1.fastq')
        reads2 = os.path.join(data_dir, 'kcount_test.reads_2.fastq')
        kcount._head_fastaq(reads1, reads2, test_out, count=4)
        self.assertTrue(filecmp.cmp(test_out, os.path.join(data_dir, 'kcount_test.reads.head.4.fa'), shallow=False))
        os.unlink(test_out)


    def test_head_fastaq_fa(self):
        '''Test head_fastaq with FASTA input'''
        test_out = 'tmp.head.fa'
        reads1 = os.path.join(data_dir, 'kcount_test.reads_1.fasta')
        reads2 = os.path.join(data_dir, 'kcount_test.reads_2.fasta')
        kcount._head_fastaq(reads1, reads2, test_out, count=4)
        self.assertTrue(filecmp.cmp(test_out, os.path.join(data_dir, 'kcount_test.reads.head.4.fa'), shallow=False))
        os.unlink(test_out)


    def test_median(self):
        '''Test _median()'''
        d = {1: 1, 2: 2, 3: 5, 4: 2, 5: 1}
        self.assertEqual(3, kcount._median(d))


    def test_run_kmc(self):
        '''Test test_run_kmc'''
        reads = os.path.join(data_dir, 'kcount_test.run_kmc.fa')
        counts_file = kcount._run_kmc(reads, 'tmp.run_kmc', 10, 2, 4)
        self.assertTrue(filecmp.cmp(counts_file, os.path.join(data_dir, 'kcount_test.run_kmc.counts'), shallow=False))
        os.unlink(counts_file)


    def test_run_kmc_two_threads(self):
        '''Test test_run_kmc with two threads'''
        reads = os.path.join(data_dir, 'kcount_test.run_kmc.fa')
        counts_file = kcount._run_kmc(reads, 'tmp.run_kmc', 10, 2, 4, threads=2)
        self.assertTrue(filecmp.cmp(counts_file, os.path.join(data_dir, 'kcount_test.run_kmc.counts'), shallow=False))
        os.unlink(counts_file)


    def test_kmc_to_kmer_counts(self):
        '''Test _kmc_to_kmer_counts'''
        counts = kcount._kmc_to_kmer_counts(os.path.join(data_dir, 'kcount_test.kmc_counts'), number=2)
        expected = {'ACGT': 10, 'ATGC': 10}
        self.assertDictEqual(counts, expected)


    def test_counts_file_to_fasta(self):
        '''Test _counts_file_to_fasta'''
        outfile = 'tmp.kmer_counts_to_fa.fa'
        expected = os.path.join(data_dir, 'kcount_test.kmc_counts.fa')
        infile = os.path.join(data_dir, 'kcount_test.kmc_counts')
        kcount._counts_file_to_fasta(infile, outfile)
        self.assertTrue(filecmp.cmp(outfile, expected, shallow=False))
        os.unlink(outfile)


    def test_get_most_common_kmers(self):
        '''Test get_most_common_kmers'''
        reads1 = os.path.join(data_dir, 'kcount_test.get_commonest_kmer_1.fa')
        reads2 = os.path.join(data_dir, 'kcount_test.get_commonest_kmer_2.fa')
        counts = kcount.get_most_common_kmers(reads1, reads2, kmer_length=10, head=100000, min_count=2, max_count=4, most_common=100, method='kmc')
        self.assertDictEqual({'AGCTAAAACT': 2, 'CTATATCTCA': 3}, counts)


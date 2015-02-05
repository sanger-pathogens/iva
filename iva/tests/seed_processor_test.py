import unittest
import os
import filecmp
import pyfastaq
from iva import seed_processor

modules_dir = os.path.dirname(os.path.abspath(seed_processor.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSeedProcessor(unittest.TestCase):
    def test_process_seeds(self):
        '''Test process_seeds'''
        in_fasta = os.path.join(data_dir, 'seed_processor_test.ref.fa')
        reads1 = os.path.join(data_dir, 'seed_processor_test.reads_1.fa')
        reads2 = os.path.join(data_dir, 'seed_processor_test.reads_2.fa')
        tmp_out = 'tmp.process_seeds.out'
        s = seed_processor.SeedProcessor(in_fasta, reads1, reads2, tmp_out, verbose=4, threads=1)
        s.process()
        reader = pyfastaq.sequences.file_reader(tmp_out)
        counter = 0
        for seq in reader:
            self.assertTrue(len(seq) > 470)
            counter += 1
        self.assertEqual(2, counter)
        os.unlink(tmp_out)

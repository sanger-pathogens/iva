import os
import unittest
from iva import read_trim

modules_dir = os.path.dirname(os.path.abspath(read_trim.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestReadTrim(unittest.TestCase):
    def test_run_trimmomatic(self):
        '''Test run_trimmomatic'''
        reads1 = os.path.join(data_dir, 'read_trim_test.reads_1.fq')
        reads2 = os.path.join(data_dir, 'read_trim_test.reads_2.fq')
        #run_trimmomatic(reads1, reads2, outprefix, trimmo_jar, adapters, minlen=50, verbose=0):
        # need to know where trimmoatic jar file is - could be anywhere - and
        # read trimming is optional, so skip this test for now...
        pass

import unittest
import os
import shutil
from iva import qc_external

modules_dir = os.path.dirname(os.path.abspath(qc_external.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestQcExternal(unittest.TestCase):
    def test_run_gage(self):
        '''test run_gage'''
        ref = os.path.join(data_dir, 'qc_external_test_run_gage.ref.fa')
        scaffs = os.path.join(data_dir, 'qc_external_test_run_gage.scaffs.fa')
        outdir = 'tmp.qc_external_test_run_gage'
        qc_external.run_gage(ref, scaffs, outdir)
        self.assertTrue(os.path.exists(os.path.join(outdir, 'gage.out')))
        with open(os.path.join(outdir, 'gage.out')) as f:
            got_lines = f.readlines()
        with open(os.path.join(data_dir, 'qc_external_test_run_gage.out')) as f:
            expected_lines = f.readlines()
        self.assertEqual(got_lines[1:], expected_lines[1:])
        shutil.rmtree(outdir)


    def test_run_ratt(self):
        '''test run_ratt'''
        embl_dir = os.path.join(data_dir, 'qc_external_test_run_ratt.embl')
        assembly = os.path.join(data_dir, 'qc_external_test_run_ratt.assembly.fa')
        outdir = 'tmp.qc_external_test_run_ratt'
        qc_external.run_ratt(embl_dir, assembly, outdir)
        self.assertTrue(os.path.exists(os.path.join(outdir, 'run.sh.out')))
        shutil.rmtree(outdir)

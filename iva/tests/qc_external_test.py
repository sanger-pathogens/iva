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

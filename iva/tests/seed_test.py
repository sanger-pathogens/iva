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
import filecmp
import shutil
import pyfastaq
from iva import seed

modules_dir = os.path.dirname(os.path.abspath(seed.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSeed(unittest.TestCase):
    def test_kmc_in_path(self):
        '''Test that kmc is in the user's path'''
        assert(shutil.which('kmc') is not None)
        assert(shutil.which('kmc_dump') is not None)


    def test_init(self):
        '''test init'''
        # TODO


    def test_extension_from_read(self):
        '''Test _test_extension_from_read'''
        s = seed.Seed(seq='AGGCT')
        self.assertEqual(None, s._extension_from_read(pyfastaq.sequences.Fasta('x', 'AAAAA')))
        self.assertEqual(None, s._extension_from_read(pyfastaq.sequences.Fasta('x', 'AGGC')))
        self.assertEqual('A', s._extension_from_read(pyfastaq.sequences.Fasta('x', 'AGGCTA')))
        self.assertEqual('AT', s._extension_from_read(pyfastaq.sequences.Fasta('x', 'AGGCTAT')))
        self.assertEqual('AT', s._extension_from_read(pyfastaq.sequences.Fasta('x', 'GGGAGGCTAT')))
        self.assertEqual('AA', s._extension_from_read(pyfastaq.sequences.Fasta('x', 'TTAGCCT')))

        self.assertEqual(None, s._extension_from_read(pyfastaq.sequences.Fasta('x', 'AAAAA'), left=True))
        self.assertEqual(None, s._extension_from_read(pyfastaq.sequences.Fasta('x', 'AGGCTA'), left=True))
        self.assertEqual('GT', s._extension_from_read(pyfastaq.sequences.Fasta('x', 'GTAGGCTA'), left=True))
        self.assertEqual('GT', s._extension_from_read(pyfastaq.sequences.Fasta('x', 'GTAGGCTATTC'), left=True))
        self.assertEqual('GT', s._extension_from_read(pyfastaq.sequences.Fasta('x', 'AGCCTAC'), left=True))


    def test_len(self):
        '''Test len'''
        s = seed.Seed(seq='AGGCT')
        self.assertEqual(5, len(s))
        s.seq = None
        self.assertEqual(0, len(s))


    def test_extensions_from_reads_file(self):
        '''Test _extensions_from_reads_file'''
        s = seed.Seed(seq='AGGCT')
        l, r = s._extensions_from_reads_file(os.path.join(data_dir, 'kcount_test.reads_1.fasta'))
        self.assertListEqual(l, [])
        self.assertListEqual(r, ['A', 'AT', 'AT'])
        l, r = s._extensions_from_reads_file(os.path.join(data_dir, 'kcount_test.reads_2.fasta'))
        self.assertListEqual(l, ['G', 'TG', 'TG'])
        self.assertListEqual(r, [])


    def test_extend_with_reads_as_single_end(self):
        '''Test _extend_with_reads_as_single_end'''
        s = seed.Seed(seq='AGGCT', ext_min_cov=1, verbose=2)
        reads1 = os.path.join(data_dir, 'kcount_test.reads_1.fasta')
        reads2 = os.path.join(data_dir, 'kcount_test.reads_2.fasta')
        s._extend_with_reads_as_single_end(reads1, reads2)
        self.assertEqual('TGAGGCTAT', s.seq)


    def test_extend(self):
        '''Test extend'''
        s = seed.Seed(seq='AGGCT', ext_min_cov=1, verbose=3)
        reads1 = os.path.join(data_dir, 'kcount_test.reads_1.fasta')
        reads2 = os.path.join(data_dir, 'kcount_test.reads_2.fasta')
        s.extend(reads1, reads2, 100)
        self.assertEqual('TGAGGCTAT', s.seq)


    def test_write_fasta(self):
        '''Test write_fasta'''
        s = seed.Seed(seq='GAAGGCGGCAGC')
        tmpfile = 'tmp.seed.fa'
        s.write_fasta(tmpfile, 'spam')
        self.assertTrue(filecmp.cmp(tmpfile, os.path.join(data_dir, 'seed_test.write_fasta.fa'), shallow=False))
        os.unlink(tmpfile)

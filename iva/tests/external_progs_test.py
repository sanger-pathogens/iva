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
from iva import external_progs

class TestExternalProgs(unittest.TestCase):
    def test_kmc_in_path(self):
        '''Test kmc in path'''
        self.assertTrue(external_progs.is_in_path('kmc'), 'Error! Did not find kmc in your path! Please install kmc')

    def test_kmc_dump_in_path(self):
        '''Test kmc_dump in path'''
        self.assertTrue(external_progs.is_in_path('kmc_dump'), 'Error! Did not find kmc_dump in your path! Please install kmc')

    def test_nucmer_in_path(self):
        '''Test nucmer in path'''
        self.assertTrue(external_progs.is_in_path('nucmer'), 'Error! Did not find nucmer in your path! Please install MUMmer')

    def test_delta_filter_in_path(self):
        '''Test delta-filter in path'''
        self.assertTrue(external_progs.is_in_path('delta-filter'), 'Error! Did not find delta-filter in your path! Please install MUMmer')

    def test_show_coords_in_path(self):
        '''Test show-coords in path'''
        self.assertTrue(external_progs.is_in_path('show-coords'), 'Error! Did not find show-coords in your path! Please install MUMmer')

    def test_samtools_in_path(self):
        '''Test samtools in path'''
        self.assertTrue(external_progs.is_in_path('samtools'), 'Error! Did not find samtools in your path! Please install samtools')

    def test_smalt_in_path(self):
        '''Test smalt in path'''
        self.assertTrue(external_progs.is_in_path('smalt'), 'Error! Did not find smalt in your path! Please install smalt')


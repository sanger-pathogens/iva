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

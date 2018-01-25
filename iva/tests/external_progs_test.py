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
    def test_r_in_path(self):
        '''Test blastn in path'''
        self.assertTrue(external_progs.is_in_path('R'), 'Error! Did not find R in your path! Please install R')
    
    def test_blastn_in_path(self):
        '''Test blastn in path'''
        self.assertTrue(external_progs.is_in_path('blastn'), 'Error! Did not find blast in your path! Please install NCBI blast+')
    
    def test_makeblastdb_in_path(self):
        '''Test makeblastdb in path'''
        self.assertTrue(external_progs.is_in_path('makeblastdb'), 'Error! Did not find makeblastdb in your path! Please install NCBI blast+')
    
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
    
    def test_r_version(self):
        '''Test R version'''
        self.assertTrue(external_progs.get_version('R',must_be_in_path=True))
        
        self.assertEqual('3.4.0', self.check_regex_version_extraction('R', """
R version 3.4.0 (2017-04-21) -- "You Stupid Darkness"
Copyright (C) 2017 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
        """ ))
        
        self.assertEqual('3.3.2', self.check_regex_version_extraction('R', """
R version 3.3.2 (2016-10-31) -- "Sincere Pumpkin Patch"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin13.4.0 (64-bit)
        """ ))
        
        self.assertEqual('3.2.2', self.check_regex_version_extraction('R', """
R version 3.2.2 (2015-08-14) -- "Fire Safety"
Copyright (C) 2015 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
        """ ))
        
        self.assertEqual('3.1.2', self.check_regex_version_extraction('R', """
R version 3.1.2 (2014-10-31) -- "Pumpkin Helmet"
Copyright (C) 2014 The R Foundation for Statistical Computing
Platform: x86_64-unknown-linux-gnu (64-bit)
        """ ))
        
        self.assertEqual('3.0.0', self.check_regex_version_extraction('R', """
R version 3.0.0 (2013-04-03) -- "Masked Marvel"
Copyright (C) 2013 The R Foundation for Statistical Computing
Platform: x86_64-unknown-linux-gnu (64-bit)
        """ ))
        
        self.assertEqual('2.15.2', self.check_regex_version_extraction('R', """
R version 2.15.2 (2012-10-26) -- "Trick or Treat"
Copyright (C) 2012 The R Foundation for Statistical Computing
ISBN 3-900051-07-0
Platform: x86_64-unknown-linux-gnu (64-bit)
        """ ))

    
    def test_blastn_version(self):
        '''Test blastn version'''
        self.assertTrue(external_progs.get_version('blastn',must_be_in_path=True))
        self.assertEqual('2.7.0+', self.check_regex_version_extraction('blastn', """
blastn: 2.7.0+
 Package: blast 2.7.0, build Sep 12 2017 15:51:33
        """ ))
        self.assertEqual('2.6.0+', self.check_regex_version_extraction('blastn', """
blastn: 2.6.0+
 Package: blast 2.6.0, build Aug 23 2017 22:18:50
        """ ))
        self.assertEqual('2.2.31+', self.check_regex_version_extraction('blastn', """
blastn: 2.2.31+
Package: blast 2.2.31, build Jul  6 2015 15:14:27
        """ ))
    
    def test_makeblastdb_version(self):
        '''Test makeblastdb version'''
        self.assertTrue(external_progs.get_version('makeblastdb',must_be_in_path=True))
        self.assertEqual('2.7.0+', self.check_regex_version_extraction('makeblastdb', """
makeblastdb: 2.7.0+
 Package: blast 2.7.0, build Sep 12 2017 15:51:33
        """ ))
        self.assertEqual('2.6.0+', self.check_regex_version_extraction('makeblastdb', """
makeblastdb: 2.6.0+
 Package: blast 2.6.0, build Aug 23 2017 22:18:50
        """ ))
        self.assertEqual('2.2.31+', self.check_regex_version_extraction('makeblastdb', """
makeblastdb: 2.2.31+
Package: blast 2.2.31, build Jul  6 2015 15:14:27
        """ ))
    
    def test_kmc_version(self):
        '''Test kmc version'''
        self.assertTrue(external_progs.get_version('kmc',must_be_in_path=True))
        self.assertEqual('2.3.0', self.check_regex_version_extraction('kmc', """
K-Mer Counter (KMC) ver. 2.3.0 (2015-08-21)
        """ ))
        self.assertEqual('3.0.0', self.check_regex_version_extraction('kmc', """
K-Mer Counter (KMC) ver. 3.0.0 (2017-01-28)
        """ ))
    
    def test_kmc_dump_version(self):
        '''Test kmc_dump version'''
        self.assertTrue(external_progs.get_version('kmc_dump',must_be_in_path=True))
        self.assertEqual('3.0.0', self.check_regex_version_extraction('kmc_dump', """
KMC dump ver. 3.0.0 (2017-01-28)
        """ ))
        self.assertEqual('2.3.0', self.check_regex_version_extraction('kmc_dump', """
KMC dump ver. 2.3.0 (2015-08-21)
        """ ))
    
    def test_nucmer_version(self):
        '''Test nucmer version'''
        self.assertTrue(external_progs.get_version('nucmer',must_be_in_path=True))
        self.assertEqual('3.1', self.check_regex_version_extraction('nucmer', """
nucmer
NUCmer (NUCleotide MUMmer) version 3.1
        """ ))
    
    def test_samtools_version(self):
        '''Test samtools version'''
        self.assertTrue(external_progs.get_version('samtools', must_be_in_path=True))
        self.assertEqual('1.6', self.check_regex_version_extraction('samtools', """
samtools 1.6
Using htslib 1.6
Copyright (C) 2017 Genome Research Ltd.""" ))

    def test_samtools_original_version(self):
        '''Test samtools original version'''
        
        self.assertEqual('0.1.19', self.check_regex_version_extraction('samtools', """
Program: samtools (Tools for alignments in the SAM format)
Version: 0.1.19-44428cd

Usage:   samtools <command> [options]""" ))
     
    def test_smalt_version(self):
       '''Test smalt version'''
       self.assertTrue(external_progs.get_version('smalt',must_be_in_path=True))
       self.assertEqual('0.7.6', self.check_regex_version_extraction('smalt', """
             SMALT - Sequence Mapping and Alignment Tool
Version: 0.7.6
Date:    21-03-2014""" ))
    
    def check_regex_version_extraction(self, prog,  raw_version_output ):
        cmd, regex = external_progs.prog_to_version_cmd[prog]
        raw_output_lines = raw_version_output.splitlines()
        for line in raw_output_lines:
            hits = regex.search(line)
            if hits:
                return str(hits.group(1))
        return None

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
import filecmp
from iva import kraken

modules_dir = os.path.dirname(os.path.abspath(kraken.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestKraken(unittest.TestCase):
    def setUp(self):
        self.db = kraken.Database(os.path.join(data_dir, 'kraken_test.db'))


    def test_count_cds_from_embl(self):
        '''test count_cds_from_embl'''
        for i in range(3):
            infile = os.path.join(data_dir, 'kraken_count_cds_from_embl.' + str(i))
            self.assertEqual(i, kraken.Database.count_cds_from_embl(infile))


    def test_get_parent_taxons(self):
        '''test _get_parent_taxons'''
        taxons = set(['1', '9', '13'])
        self.db._get_parent_taxons(taxons)
        expected = {
            '1': '1',
            '9': '32199',
            '13': '203488'
        }
        self.assertEqual(expected, self.db.taxon_to_parent)


    def test_load_extra_ref_info(self):
        '''test _load_extra_ref_info'''
        self.db.extra_refs_file = os.path.join(data_dir, 'kraken_test.extra_ids')
        self.db._load_extra_ref_info()
        expected = {
            2000000000: {
                'genbank_ids': ['10', '11', '12', '13'],
                'new_gis': [4000000000, 4000000001, 4000000002, 4000000003]
            },
            2000000001: {
                'genbank_ids': ['21'],
                'new_gis': [4000000004]
            },
            2000000002: {
                'genbank_ids': ['42', '43'],
                'new_gis': [4000000005, 4000000006]
            }
        }
        self.assertEqual(expected, self.db.extra_refs)


    def test_genbank_to_taxon_and_gi(self):
        '''test _genbank_to_taxon_and_gi'''
        expected = '43', '42'
        gb_file = os.path.join(data_dir, 'kraken_test.genbank_to_taxon_and_gi.gb')
        self.assertEqual(expected, self.db._genbank_to_taxon_and_gi(gb_file))


    def test_replace_fasta_header(self):
        '''test _replace_fasta_header'''
        tmp = 'tmp.test.fa'
        before = os.path.join(data_dir, 'kraken_test.replace_fasta_header.before.fa')
        after = os.path.join(data_dir, 'kraken_test.replace_fasta_header.after.fa')
        shutil.copyfile(before, tmp)
        self.db._replace_fasta_header(tmp, 'after')
        self.assertTrue(filecmp.cmp(tmp, after))
        os.unlink(tmp)


    def test_append_to_file(self):
        '''test _append_to_file'''
        tmp = 'tmp.test'
        before = os.path.join(data_dir, 'kraken_test.append_to_file.before')
        after = os.path.join(data_dir, 'kraken_test.append_to_file.after')
        shutil.copyfile(before, tmp)
        self.db._append_to_file(tmp, '42')
        self.assertTrue(filecmp.cmp(tmp, after))
        os.unlink(tmp)


    def test_species_to_dir(self):
        '''test species_to_dir'''
        self.assertEqual('a_b_c_d_e_f', self.db._species_to_embl_dir(r''' a\b/c(d)e.f-'''))
        self.assertEqual('added.12345', self.db._species_to_embl_dir('added.12345'))


    def test_get_most_common_species_dir(self):
        '''Test get_most_common_species'''
        report = os.path.join(data_dir, 'kraken_test.report')
        expected = os.path.join(self.db.embl_root, 'Human_immunodeficiency_virus_1')
        self.assertEqual(expected, self.db._get_most_common_species_dir(report))


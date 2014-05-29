import unittest
import os
import filecmp
from iva import kraken

modules_dir = os.path.dirname(os.path.abspath(kraken.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestKraken(unittest.TestCase):
    def test_get_most_common_species(self):
        '''Test get_most_common_species'''
        report = os.path.join(data_dir, 'kraken_test.report')
        self.assertEqual('Human_immunodeficiency_virus_1', kraken.get_most_common_species(report))


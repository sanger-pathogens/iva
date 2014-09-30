import unittest
import filecmp
import os
import copy
from iva import contig_trim

modules_dir = os.path.dirname(os.path.abspath(contig_trim.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')



class TestContigTrim(unittest.TestCase):
    def test_coverage_to_trimmed_coords(self):
        '''Test _coverage_to_trimmed_coords'''
        coverage = [
            [0] * 4,
            [0] * 5,
            [0] * 6,
            [0] * 20,
            [1] * 5 + [0] * 15,
            [0] * 15 + [1] * 5,
            [0] * 2 + [1] * 5 + [0] * 13,
            [0] * 3 + [1] * 5 + [0] * 12,
            [0] * 4 + [1] * 5 + [0] * 11,
            [0] * 5 + [1] * 5 + [0] * 10,
            [1] * 20,
            [0] * 2 + [1] * 16 + [0] * 2,
        ]

        expected = [
            (0, 3),
            (0, 4),
            (0, 5),
            (0, 19),
            (5, 19),
            (0, 14),
            (7, 19),
            (8, 19),
            (9, 19),
            (0, 19),
            None,
            None,
        ]

        assert len(coverage) == len(expected)

        for i in range(len(coverage)):
            got = contig_trim._coverage_to_trimmed_coords(coverage[i], min_dist_to_end=3, window_length=5, min_pc=80)

            # Should be able to reverse the coverage and still get the same results (but also reversed)
            coverage_rev = copy.copy(coverage[i])
            coverage_rev.reverse()
            if expected[i] is None:
                expected_rev = None
            else:
                expected_rev = [len(coverage[i]) - x - 1 for x in expected[i]]
                expected_rev.reverse()
                expected_rev = tuple(expected_rev)
        got_rev = contig_trim._coverage_to_trimmed_coords(coverage_rev, min_dist_to_end=3, window_length=5, min_pc=80)
        

    def test_trim_ends(self):
        '''Test _trim_ends'''
        before_trim = os.path.join(data_dir, 'contig_trim_test_contigs.fa')
        expected_after_trim = os.path.join(data_dir, 'contig_trim_test_contigs.trimmed.fa')
        adapters = os.path.join(data_dir, 'contig_trim_test_contigs.adapters_and_primers.fa')
        tmp_out = 'tmp.trimmed.fa'
        contig_trim._trim_ends(before_trim, tmp_out, adapters, min_length=20, min_dist_to_end=5, window_length=10, min_pc=90)
        self.assertTrue(filecmp.cmp(tmp_out, expected_after_trim, shallow=False))
        os.unlink(tmp_out)


    def test_trim_primers_and_adapters(self):
        '''Test trim_primers_and_adapters'''
        before_trim = os.path.join(data_dir, 'contig_trim_test_contigs.fa')
        expected_after_trim = os.path.join(data_dir, 'contig_trim_test_contigs.trimmed.fa')
        adapters = os.path.join(data_dir, 'contig_trim_test_contigs.adapters.fa')
        primers = os.path.join(data_dir, 'contig_trim_test_contigs.primers.fa')
        tmp_out = 'tmp.trimmed.fa'
        contig_trim.trim_primers_and_adapters(before_trim, tmp_out, adapters, primers, min_length=20, min_dist_to_end=5, window_length=10, min_pc=90)
        self.assertTrue(filecmp.cmp(tmp_out, expected_after_trim, shallow=False))
        os.unlink(tmp_out)
        

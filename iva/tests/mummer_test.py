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
import pysam
import pyfastaq
from iva import mummer, edge

modules_dir = os.path.dirname(os.path.abspath(mummer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestMummer(unittest.TestCase):
    def test_run_nucmer(self):
        '''test run_nucmer'''
        # TODO
        pass


    def test_swap(self):
        ''' test swap'''
        # TODO
        pass


    def test_sort(self):
        '''test sort'''
        # TODO
        pass


    def test_file_read(self):
        '''test file_read'''
        # TODO
        pass


    def test_qry_coords(self):
        '''Test qry_coords'''
        hits = ['\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']),
                '\t'.join(['1', '101', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'qry'])
        ]
        for h in hits:
            m = mummer.NucmerHit(h)
            self.assertEqual(pyfastaq.intervals.Interval(0,99), m.qry_coords())


    def test_ref_coords(self):
        '''Test ref_coords'''
        hits = ['\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']),
                '\t'.join(['100', '1', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])
        ]
        for h in hits:
            m = mummer.NucmerHit(h)
            self.assertEqual(pyfastaq.intervals.Interval(0,99), m.ref_coords())


    def test_on_same_strand(self):
        '''test on_same_strand'''
        self.assertTrue(mummer.NucmerHit('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertTrue(mummer.NucmerHit('\t'.join(['100', '1', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertFalse(mummer.NucmerHit('\t'.join(['1', '100', '100', '1', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())
        self.assertFalse(mummer.NucmerHit('\t'.join(['100', '1', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref'])).on_same_strand())


    def test_is_self_hit(self):
        '''Test is_self_hit'''
        tests = [('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), True),
            ('\t'.join(['1', '101', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), False),
            ('\t'.join(['2', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref']), False),
            ('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1000', '1000', '1', '1', 'ref', 'ref2']), False),
            ('\t'.join(['1', '100', '1', '100', '100', '100', '99.9', '1000', '1000', '1', '1', 'ref', 'ref']), False),
        ]

        for t in tests:
            m = mummer.NucmerHit(t[0])
            self.assertEqual(m.is_self_hit(), t[1])
        pass


    def test_to_graph_edge(self):
        '''Test to_graph_edge'''
        hits = [
            '\t'.join(['781', '981', '10', '210', '200', '200', '98', '1000', '1000', '1', '1', 'ref', 'qry']), # %id too low
            '\t'.join(['781', '980', '10', '210', '199', '200', '98', '1000', '1000', '1', '1', 'ref', 'qry']), # hit too short
            '\t'.join(['781', '981', '10', '209', '200', '199', '98', '1000', '1000', '1', '1', 'ref', 'qry']), # hit too short
            '\t'.join(['1', '200', '1', '200', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']), # bad orientation
            '\t'.join(['200', '1', '200', '1', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']), # bad orientation
            '\t'.join(['800', '1000', '800', '1000', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']), # bad orientation
            '\t'.join(['1000', '800', '1000', '800', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']), # bad orientation
            '\t'.join(['300', '500', '300', '500', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']), # not at ends
            '\t'.join(['1', '1000', '1', '1000', '1000', '1000', '100.00', '1000', '1000', '1', '1', 'ref', 'qry']), # whole contigs hit
            '\t'.join(['1', '500', '1', '1000', '500', '500', '100.00', '500', '1000', '1', '1', 'ref1', 'qry1']), # contained contig
            '\t'.join(['1', '200', '791', '991', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref2', 'qry2']),
            '\t'.join(['781', '981', '10', '210', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref3', 'qry3']),
            '\t'.join(['991', '791', '781', '981', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref4', 'qry4']),
            '\t'.join(['210', '10', '15', '215', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref5', 'qry5']),
            '\t'.join(['781', '981', '991', '791', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref6', 'qry6']),
            '\t'.join(['10', '210', '215', '5', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref7', 'qry7']),
            '\t'.join(['210', '10', '980', '780', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref8', 'qry8']),
            '\t'.join(['995', '795', '215', '15', '200', '200', '100.00', '1000', '1000', '1', '1', 'ref9', 'qry9']),
        ]

        expected = [
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            edge.Edge('qry2', 790, 990, 'ref2', 0, 199),
            edge.Edge('ref3', 780, 980, 'qry3', 9, 209),
            edge.Edge('qry4', 780, 980, 'ref4', 990, 790),
            edge.Edge('ref5', 209, 9, 'qry5', 14, 214),
            edge.Edge('ref6', 780, 980, 'qry6', 990, 790),
            edge.Edge('qry7', 214, 4, 'ref7', 9, 209),
            edge.Edge('ref8', 209, 9, 'qry8', 979, 779),
            edge.Edge('qry9', 214, 14, 'ref9', 994, 794),
        ]

        assert len(expected) == len(hits)

        for i in range(len(hits)):
            m = mummer.NucmerHit(hits[i])
            self.assertEqual(m.to_graph_edge(), expected[i])


    def test_is_at_ends(self):
        '''Test is_at_ends'''
        tests = [('\t'.join(['1', '100', '200', '300', '100', '100', '100.00', '1000', '500', '1', '1', 'ref', 'qry']), False, mummer.HIT_AT_START),
            ('\t'.join(['51', '151', '200', '300', '100', '100', '100.00', '1000', '500', '1', '1', 'ref', 'qry']), False, mummer.HIT_NO_ENDS),
            ('\t'.join(['900', '1000', '200', '300', '100', '100', '100.00', '1000', '500', '1', '1', 'ref', 'qry']), False, mummer.HIT_AT_END),
            ('\t'.join(['1000', '900', '200', '300', '100', '100', '100.00', '1000', '500', '1', '1', 'ref', 'qry']), False, mummer.HIT_AT_END),
            ('\t'.join(['850', '949', '200', '300', '100', '100', '100.00', '1000', '500', '1', '1', 'ref', 'qry']), False, mummer.HIT_NO_ENDS),
            ('\t'.join(['42', '992', '200', '1152', '950', '950', '100.00', '1000', '5000', '1', '1', 'ref', 'qry']), False, mummer.HIT_AT_BOTH_ENDS),
            ('\t'.join(['200', '300', '1', '100', '100', '100', '100.00', '500', '1000', '1', '1', 'ref', 'qry']), True, mummer.HIT_AT_START),
            ('\t'.join(['200', '300', '51', '151', '100', '100', '100.00', '500', '1000', '1', '1', 'ref', 'qry']), True, mummer.HIT_NO_ENDS),
            ('\t'.join(['200', '300', '900', '1000', '100', '100', '100.00', '500', '1000', '1', '1', 'ref', 'qry']), True, mummer.HIT_AT_END),
            ('\t'.join(['200', '300', '1000', '900', '100', '100', '100.00', '500', '1000', '1', '1', 'ref', 'qry']), True, mummer.HIT_AT_END),
            ('\t'.join(['200', '300', '850', '949', '100', '100', '100.00', '500', '1000', '1', '1', 'ref', 'qry']), True, mummer.HIT_NO_ENDS),
            ('\t'.join(['200', '300', '42', '992', '950', '950', '100.00', '500', '1000', '1', '1', 'ref', 'qry']), True, mummer.HIT_AT_BOTH_ENDS),
        ]

        for t in tests:
            m = mummer.NucmerHit(t[0])
            self.assertEqual(m._is_at_ends(use_qry=t[1]), t[2])

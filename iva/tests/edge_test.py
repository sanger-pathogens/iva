import unittest
import copy
from iva import edge
from pyfastaq import intervals

class TestEdge(unittest.TestCase):
    def test_open_end(self):
        '''test open_end'''
        e1 = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e2 = edge.Edge('c1', 42, 1, 'c2', 1, 10)
        self.assertEqual(e1.open_end('c1'), edge.LEFT)
        self.assertEqual(e1.open_end('c2'), edge.LEFT)
        self.assertEqual(e2.open_end('c1'), edge.RIGHT)
        self.assertEqual(e2.open_end('c2'), edge.RIGHT)


    def test_reverse(self):
        '''test reverse'''
        e = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e.reverse()
        self.assertEqual(e, edge.Edge('c2', 10, 50, 'c1', 42, 1))


    def test_make_contig_forwards(self):
        '''test make_contig_forwards'''
        e = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e._make_contig_forwards('c1')
        self.assertEqual(e, edge.Edge('c1', 1, 42, 'c2', 50, 10))
        e._make_contig_forwards('c2')
        self.assertEqual(e, edge.Edge('c2', 10, 50, 'c1', 42, 1))
        e._make_contig_forwards('c1')
        self.assertEqual(e, edge.Edge('c1', 1, 42, 'c2', 50, 10))


    def test_make_contig_first(self):
        '''test make_contig_first'''
        e = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e_original = copy.copy(e)
        e.make_contig_first('c1')
        self.assertEqual(e, e_original)
        e.make_contig_first('c2')
        self.assertEqual(e, e_original)
        e.make_contig_first('c1')
        self.assertEqual(e, e_original)


    def test_change_hit_coords_with_intersection(self):
        '''test test_change_hit_coords_with_intersection'''
        e = edge.Edge('c1', 1, 42, 'c2', 10, 50)
        e._change_hit_coords_with_intersection('c1', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 20, 30, 'c2', 29, 38))

        e = edge.Edge('c1', 1, 42, 'c2', 10, 50)
        e._change_hit_coords_with_intersection('c2', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 11, 22, 'c2', 20, 30))

        e = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e._change_hit_coords_with_intersection('c1', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 20, 30, 'c2', 31, 22))

        e = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e._change_hit_coords_with_intersection('c2', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 21, 32, 'c2', 30, 20))

        e = edge.Edge('c1', 42, 1, 'c2', 10, 50)
        e._change_hit_coords_with_intersection('c1', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 30, 20, 'c2', 22, 31))

        e = edge.Edge('c1', 42, 1, 'c2', 10, 50)
        e._change_hit_coords_with_intersection('c2', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 32, 21, 'c2', 20, 30))

        e = edge.Edge('c1', 42, 1, 'c2', 50, 10)
        e._change_hit_coords_with_intersection('c1', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 30, 20, 'c2', 38, 29))

        e = edge.Edge('c1', 42, 1, 'c2', 50, 10)
        e._change_hit_coords_with_intersection('c2', intervals.Interval(20, 30))
        self.assertEqual(e, edge.Edge('c1', 22, 11, 'c2', 30, 20))


    def test_merge_into(self):
        '''test merge_into'''
        e1 = edge.Edge('c1', 71, 100, 'c2', 1, 30)
        e2 = edge.Edge('c2', 41, 80, 'c3', 1, 40)
        e2_original = copy.copy(e2)
        self.assertFalse(e1.merge_into(e2, 'c2'))
        self.assertEqual(e2, e2_original)

        e1 = edge.Edge('c1', 51, 100, 'c2', 1, 50)
        e2 = edge.Edge('c2', 41, 80, 'c3', 1, 40)
        e2_original = copy.copy(e2)
        self.assertTrue(e1.merge_into(e2, 'c2'))
        self.assertEqual(e2, e2_original)
        self.assertEqual(e1, edge.Edge('c1', 91, 100, 'c3', 1, 10))

        e2 = edge.Edge('c1', 51, 100, 'c2', 1, 50)
        e1 = edge.Edge('c2', 41, 80, 'c3', 1, 40)
        e2_original = copy.copy(e2)
        self.assertTrue(e1.merge_into(e2, 'c2'))
        self.assertEqual(e2, e2_original)
        self.assertEqual(e1, edge.Edge('c1', 91, 100, 'c3', 1, 10))

        e1 = edge.Edge('c1', 51, 1, 'c2', 2, 52)
        e2 = edge.Edge('c2', 41, 80, 'c3', 1, 40)
        e2_original = copy.copy(e2)
        self.assertTrue(e1.merge_into(e2, 'c2'))
        self.assertEqual(e2, e2_original)
        self.assertEqual(e1, edge.Edge('c1', 12, 1, 'c3', 1, 12))

        e2 = edge.Edge('c1', 51, 1, 'c2', 2, 52)
        e1 = edge.Edge('c2', 41, 80, 'c3', 1, 40)
        e2_original = copy.copy(e2)
        self.assertTrue(e1.merge_into(e2, 'c2'))
        self.assertEqual(e2, e2_original)
        self.assertEqual(e1, edge.Edge('c1', 12, 1, 'c3', 1, 12))


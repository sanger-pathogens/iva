import unittest
import os
import filecmp
import pysam
from iva import graph, assembly, edge
from pyfastaq import intervals

modules_dir = os.path.dirname(os.path.abspath(graph.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestGraph(unittest.TestCase):
    def setUp(self):
        self.asm = assembly.Assembly(contigs_file=os.path.join(data_dir, 'graph_test.contigs.fa'))
        self.g = graph.Graph(self.asm)


    def test_get_nodes(self):
        '''test get_nodes'''
        self.assertListEqual(['c1', 'c2', 'c3', 'c4', 'c5', 'c6'], self.g.get_nodes())


    def test_remove_node(self):
        '''test remove_node'''
        # TODO
        pass


    def test_add_and_get_edges(self):
        '''test adding and getting edges'''
        self.assertListEqual([], self.g._get_edges('not_in_graph', 'c1'))
        self.assertListEqual([], self.g._get_edges('not_in_graph', 'c2'))
        self.assertListEqual([], self.g._get_edges('c1', 'c2'))
        self.assertListEqual([], self.g._get_edges('c2', 'c1'))

        e1 = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e2 = edge.Edge('c1', 4, 43, 'c2', 60, 10)
        self.g.add_edge(e1)
        self.assertListEqual([e1], self.g._get_edges('c1', 'c2'))
        self.g.add_edge(e2)
        self.assertListEqual([e1, e2], self.g._get_edges('c1', 'c2'))

        with self.assertRaises(graph.Error):
            self.g.add_edge(edge.Edge('c7', 1, 10, 'c8', 42, 84))


    def test_degree(self):
        '''test degree'''
        self.assertEqual(0, self.g._degree('c1'))
        self.assertEqual(0, self.g._degree('c2'))
        e1 = edge.Edge('c1', 1, 42, 'c2', 50, 10)
        e2 = edge.Edge('c1', 4, 43, 'c2', 60, 10)
        self.g.add_edge(e1)
        self.assertEqual(1, self.g._degree('c1'))
        self.assertEqual(1, self.g._degree('c2'))
        self.g.add_edge(e2)
        self.assertEqual(2, self.g._degree('c1'))
        self.assertEqual(2, self.g._degree('c2'))


    def test_connected_components(self):
        '''test connected_components'''
        self.assertListEqual([['c1'], ['c2'], ['c3'], ['c4'], ['c5'], ['c6']], self.g.connected_components())
        self.g.add_edge(edge.Edge('c1', 42, 1, 'c2', 1, 10))
        self.assertListEqual([['c1', 'c2'], ['c3'], ['c4'], ['c5'], ['c6']], self.g.connected_components())
        self.g.add_edge(edge.Edge('c3', 42, 1, 'c2', 1, 10))
        self.assertListEqual([['c1', 'c2', 'c3'], ['c4'], ['c5'], ['c6']], self.g.connected_components())


    def test_find_simple_path(self):
        '''test find_simple_path'''
        self.assertListEqual([], self.g.find_simple_path(['c1']))
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.assertListEqual(['c1', 'c2'], self.g.find_simple_path(['c1', 'c2']))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.assertListEqual(['c1', 'c2', 'c3'], self.g.find_simple_path(['c1', 'c2', 'c3']))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.assertListEqual([], self.g.find_simple_path(['c1', 'c2', 'c3']))
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.assertListEqual([], self.g.find_simple_path(['c1', 'c2']))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.g.add_edge(edge.Edge('c3', 10, 100, 'c4', 1, 1000))
        self.assertListEqual(['c1', 'c2', 'c3', 'c4'], self.g.find_simple_path(['c1', 'c2', 'c3', 'c4']))
        self.g.add_edge(edge.Edge('c3', 10, 100, 'c5', 1, 1000))
        self.assertListEqual([], self.g.find_simple_path(['c1', 'c2', 'c3']))


    def test_remove_middle_node(self):
        '''test remove_middle_node'''
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 100, 142, 'c3', 1, 42))
        self.assertFalse(self.g._remove_middle_node('c1', 'c2', 'c3'))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 51, 100, 'c2', 1, 50))
        self.g.add_edge(edge.Edge('c2', 41, 80, 'c3', 1, 40))
        self.assertTrue(self.g._remove_middle_node('c1', 'c2', 'c3'))
        self.assertTrue('c2' not in self.g.graph)
        self.assertEqual(1, len(self.g.graph['c1']['c3']['edges']))
        self.assertEqual(1, len(self.g.graph['c3']['c1']['edges']))
        self.assertEqual(self.g.graph['c1']['c3']['edges'][0], edge.Edge('c1', 91, 100, 'c3', 1, 10))
        self.assertEqual(self.g.graph['c3']['c1']['edges'][0], edge.Edge('c1', 91, 100, 'c3', 1, 10))


    def test_remove_redundant_nodes_from_simple_path(self):
        '''test remove_redundant_nodes_from_simple_path'''
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 100, 142, 'c3', 1, 42))
        new_path = self.g.remove_redundant_nodes_from_simple_path(['c1', 'c2', 'c3'])
        self.assertListEqual(['c1', 'c2', 'c3'], new_path)
        self.assertTrue('c2' in self.g.graph)

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 51, 100, 'c2', 1, 50))
        self.g.add_edge(edge.Edge('c2', 41, 80, 'c3', 1, 40))
        new_path = self.g.remove_redundant_nodes_from_simple_path(['c1', 'c2', 'c3'])
        self.assertListEqual(['c1', 'c3'], new_path)
        self.assertFalse('c2' in self.g.graph)

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 51, 100, 'c2', 1, 50))
        self.g.add_edge(edge.Edge('c2', 41, 80, 'c3', 1, 40))
        self.g.add_edge(edge.Edge('c3', 100, 150, 'c4', 1, 50))
        self.g.add_edge(edge.Edge('c4', 20, 70, 'c5', 1, 50))
        self.g.add_edge(edge.Edge('c5', 40, 90, 'c6', 1, 50))
        new_path = self.g.remove_redundant_nodes_from_simple_path(['c1', 'c2', 'c3', 'c4', 'c5', 'c6'])
        self.assertListEqual(['c1', 'c3', 'c5', 'c6'], new_path)
        self.assertFalse('c2' in self.g.graph)
        self.assertFalse('c4' in self.g.graph)


    def test_node_to_coords(self):
        '''test _node_to_coords'''
        self.g.add_edge(edge.Edge('c1', 199, 0, 'c2', 1319, 1119))
        nodes = ['c1', 'c2']
        self.assertListEqual(['c1', intervals.Interval(200, 659), True], self.g._node_to_coords(nodes, 0))
        self.assertListEqual(['c2', intervals.Interval(0, 1319), True], self.g._node_to_coords(nodes, 1))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c2', 1119, 1319, 'c1', 0, 199))
        nodes = ['c1', 'c2']
        self.assertListEqual(['c1', intervals.Interval(200, 659), True], self.g._node_to_coords(nodes, 0))
        self.assertListEqual(['c2', intervals.Interval(0, 1319), True], self.g._node_to_coords(nodes, 1))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 551, 650, 'c2', 1, 100))
        self.g.add_edge(edge.Edge('c2', 1201, 1300, 'c3', 11, 110))
        nodes = ['c1', 'c2', 'c3']
        self.assertListEqual(['c1', intervals.Interval(0, 550), False], self.g._node_to_coords(nodes, 0))
        self.assertListEqual(['c2', intervals.Interval(1, 1200), False], self.g._node_to_coords(nodes, 1))
        self.assertListEqual(['c3', intervals.Interval(11, 2159), False], self.g._node_to_coords(nodes, 2))


        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 551, 650, 'c2', 1, 100))
        self.g.add_edge(edge.Edge('c3', 110, 11, 'c2', 1300, 1201))
        self.g.add_edge(edge.Edge('c3', 2051, 2150, 'c4', 6, 115))
        nodes = ['c1', 'c2', 'c3', 'c4']
        self.assertListEqual(['c1', intervals.Interval(0, 550), False], self.g._node_to_coords(nodes, 0))
        self.assertListEqual(['c2', intervals.Interval(1, 1200), False], self.g._node_to_coords(nodes, 1))
        self.assertListEqual(['c3', intervals.Interval(11, 2050), False], self.g._node_to_coords(nodes, 2))
        self.assertListEqual(['c4', intervals.Interval(6, 1259), False], self.g._node_to_coords(nodes, 3))


    def test_merged_coords_from_simple_nonredundant_path(self):
        '''test merged_coords_from_simple_nonredundant_path'''
        self.g.add_edge(edge.Edge('c1', 199, 0, 'c2', 1319, 1119))
        nodes = ['c1', 'c2']
        coords = self.g.merged_coords_from_simple_nonredundant_path(nodes)
        expected = [
            ['c1', intervals.Interval(200, 659), True],
            ['c2', intervals.Interval(0, 1319), True],
        ]
        self.assertListEqual(expected, coords)


        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 610, 652, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 1250, 1310, 'c3', 5, 65))
        nodes = ['c1', 'c2', 'c3']
        coords = self.g.merged_coords_from_simple_nonredundant_path(nodes)
        expected = [
            ['c1', intervals.Interval(0, 609), False],
            ['c2', intervals.Interval(1, 1249), False],
            ['c3', intervals.Interval(5, 2159), False]
        ]
        self.assertListEqual(expected, coords)


    def test_simple_path_is_consistent(self):
        '''test simple_path_is_consistent'''
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.assertTrue(self.g.simple_path_is_consistent(['c1', 'c2']))
        self.g.add_edge(edge.Edge('c3', 10, 100, 'c4', 1, 1000))
        self.assertTrue(self.g.simple_path_is_consistent(['c1', 'c2', 'c3']))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 1000, 100, 'c3', 1, 42))
        self.assertFalse(self.g.simple_path_is_consistent(['c1', 'c2', 'c3']))


    def test_edges_are_consistent(self):
        '''test edges_are_consistent'''
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.assertTrue(self.g._edges_are_consistent('c1', 'c2', 'c3'))
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c3', 1, 42))
        self.assertFalse(self.g._edges_are_consistent('c1', 'c2', 'c3'))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.assertFalse(self.g._edges_are_consistent('c1', 'c2', 'c3'))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.assertFalse(self.g._edges_are_consistent('c1', 'c2', 'c3'))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 1, 42))
        self.g.add_edge(edge.Edge('c2', 100, 50, 'c3', 1, 42))
        self.assertFalse(self.g._edges_are_consistent('c1', 'c2', 'c3'))

        self.g = graph.Graph(self.asm)
        self.g.add_edge(edge.Edge('c1', 1, 100, 'c2', 42, 1))
        self.g.add_edge(edge.Edge('c2', 50, 100, 'c3', 1, 42))
        self.assertFalse(self.g._edges_are_consistent('c1', 'c2', 'c3'))


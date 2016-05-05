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
import networkx
from pyfastaq import intervals
from iva import edge

class Error (Exception): pass


class Graph:
    def __init__(self, assembly_in, contigs=None):
        self.graph = networkx.Graph()
        self.contig_lengths = {}
        if contigs is None:
            contigs = list(assembly_in.contigs.keys())
        self.graph.add_nodes_from(contigs)
        for contig in contigs:
            self.contig_lengths[contig] = len(assembly_in.contigs[contig])


    def __str__(self):
        return '\n'.join([str(k) + ': ' + str(self.edges[k]) for k in self.edges])


    def get_nodes(self):
        return sorted(self.graph.nodes())


    def _remove_node(self, contig):
        del self.contig_lengths[contig]
        self.graph.remove_node(contig)


    def add_edge(self, edge_in):
        contig1 = edge_in.names[0]
        contig2 = edge_in.names[1]
        if contig1 not in self.graph or contig2 not in self.graph:
            raise Error('Cannot add edge to graph because nodes not already in graph. Edge:' + str(edge_in))
        self.graph.add_edge(contig1, contig2)
        if 'edges' not in self.graph[contig1][contig2]:
            self.graph[contig1][contig2]['edges'] = []

        self.graph[contig1][contig2]['edges'].append(edge_in)


    def _get_edges(self, contig1, contig2):
        if contig1 not in self.graph[contig2]:
            return []
        return self.graph[contig1][contig2].get('edges', [])


    def _degree(self, node):
        assert node in self.graph
        total = 0
        for contig in self.graph[node]:
            if 'edges' in self.graph[node][contig]:
                total += len(self.graph[node][contig]['edges'])
        return total


    def connected_components(self):
        return sorted([sorted(x) for x in networkx.connected_components(self.graph)])


    def find_simple_path(self, nodes):
        if len(nodes) <= 1:
            return []
        elif len(nodes) == 2:
            if len([1 for n in nodes if self._degree(n) == 1]) == 2:
                return sorted(list(nodes))
            else:
                return []

        # if this set of nodes has a single path through it,
        # there should be no nodes of degree 3 or more. We want
        # to know the start and end points, which are the nodes
        # of degree 1
        degree_one_nodes = set()
        for node in nodes:
            degree = self._degree(node)
            if degree > 2:
                return []
            elif degree == 1:
                degree_one_nodes.add(node)

        assert len(degree_one_nodes) == 2
        node1, node2 = list(degree_one_nodes)
        path = list(networkx.all_simple_paths(self.graph, node1, node2))
        assert len(path) == 1
        if self.simple_path_is_consistent(path[0]):
            # not really necessary, but makes unit testing easier
            if path[0][0] < path[0][-1]:
                return path[0]
            else:
                return path[0][::-1]
        else:
            return []


    def _remove_middle_node(self, contig1, contig2, contig3):
        assert self._edges_are_consistent(contig1, contig2, contig3)
        edge12 = self._get_edges(contig1, contig2)[0]
        edge23 = self._get_edges(contig2, contig3)[0]
        merged = edge12.merge_into(edge23, contig2)
        if merged:
            self.add_edge(edge12)
            self._remove_node(contig2)

        return merged


    def remove_redundant_nodes_from_simple_path(self, nodes):
        if len(nodes) < 3:
            return nodes
        i = 0
        while i + 2 < len(nodes):
            if self._remove_middle_node(nodes[i], nodes[i+1], nodes[i+2]):
                nodes.pop(i+1)
            else:
                i += 1

        return nodes


    def _node_to_coords(self, nodes, i):
        assert 0 <= i < len(nodes)
        node = nodes[i]
        if i == len(nodes) - 1:
            edges = self.graph[nodes[i-1]][node]['edges']
        else:
            edges = self.graph[node][nodes[i+1]]['edges']
            if i > 0:
                previous_edges = self.graph[nodes[i-1]][node]['edges']
                assert len(previous_edges) == 1
                previous_e = previous_edges[0]
                previous_e.make_contig_first(nodes[i-1])
                previous_open_end = previous_e.open_end(node)

        assert len(edges) == 1
        e = edges[0]
        e.make_contig_first(node)
        open_end = e.open_end(node)

        if 0 < i < len(nodes) - 1:
            assert open_end != previous_open_end
            coords = intervals.Interval(
                min(e.coords[node].start, previous_e.coords[node].start),
                max(e.coords[node].start - 1, previous_e.coords[node].start - 1)
            )
        elif i == 0:
            if open_end == edge.LEFT:
                coords = intervals.Interval(0, e.coords[node].start - 1)
            else:
                coords = intervals.Interval(e.coords[node].end + 1, self.contig_lengths[node] - 1)
        else:
            e.reverse() # now node is second in the edge, not first
            open_end = e.open_end(node)
            if open_end == edge.LEFT:
                coords = intervals.Interval(0, e.coords[node].end)
            else:
                coords = intervals.Interval(e.coords[node].start, self.contig_lengths[node] - 1)

        return [node, coords, e.rev[node]]


    def merged_coords_from_simple_nonredundant_path(self, nodes):
        assert len(nodes) >= 2
        return [self._node_to_coords(nodes, i) for i in range(len(nodes))]


    def simple_path_is_consistent(self, path):
        if len(path) <=2:
            return True

        for i in range(len(path)-2):
            if not self._edges_are_consistent(path[i], path[i+1], path[i+2]):
                return False

        return True


    def _edges_are_consistent(self, contig1, contig2, contig3):
        assert contig1 in self.graph
        assert contig2 in self.graph
        assert contig3 in self.graph
        assert contig2 in self.graph.neighbors(contig1)
        assert contig2 in self.graph.neighbors(contig3)

        edges12 = self._get_edges(contig1, contig2)
        edges23 = self._get_edges(contig2, contig3)
        edges13 = self._get_edges(contig1, contig3)
        if len(edges12) != 1 or len(edges23) != 1 or len(edges13) > 0:
            return False

        return edges12[0].open_end(contig2) != edges23[0].open_end(contig2)


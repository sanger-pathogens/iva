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
import copy
from pyfastaq import intervals

class Error (Exception): pass

LEFT = 0
RIGHT = 1

class Edge:
    def __init__(self, contig1, start1, end1, contig2, start2, end2):
        self.names = [contig1, contig2]
        self.rev = {
            contig1: start1 > end1,
            contig2: start2 > end2
        }
        self.coords = {
            contig1: intervals.Interval(min(start1, end1), max(start1, end1)),
            contig2: intervals.Interval(min(start2, end2), max(start2, end2))
        }

        if self.names[0] > self.names[1]:
            self.reverse()


    def __str__(self):
        contig1 = self.names[0]
        contig2 = self.names[1]
        return contig1 + ',' + str(self.coords[contig1]) + ',' + str(self.rev[contig1]) + '; ' + \
               contig2 + ',' + str(self.coords[contig2]) + ',' + str(self.rev[contig2])


    def __eq__(self, other):
        if set(self.names) != set(other.names):
            return False
        if self.names[0] != other.names[0]:
            other.reverse()
        equal = type(other) is type(self) and self.__dict__ == other.__dict__
        if self.names[0] != other.names[0]:
            other.reverse()
        return equal


    def open_end(self, contig):
        assert contig in self.coords
        is_rev = self.rev[contig]
        is_first = (contig == self.names[0])
        if is_first:
            if self.rev[contig]:
                return RIGHT
            else:
                return LEFT
        else:
            if self.rev[contig]:
                return LEFT
            else:
                return RIGHT


    def reverse(self):
        self.names = self.names[::-1]
        self.rev[self.names[0]] = not self.rev[self.names[0]]
        self.rev[self.names[1]] = not self.rev[self.names[1]]


    def _make_contig_forwards(self, contig):
        assert contig in self.names
        if self.rev[contig]:
            self.reverse()


    def make_contig_first(self, contig):
        if contig != self.names[0]:
            assert contig == self.names[1]
            self.reverse()


    def _change_hit_coords_with_intersection(self, contig, new_interval):
        assert contig in self.names
        coords_intersection = self.coords[contig].intersection(new_interval)
        if coords_intersection is None:
            return False
        other_contig = [x for x in self.names if x != contig][0]
        left_bases_lost = coords_intersection.start - self.coords[contig].start
        right_bases_lost = self.coords[contig].end - coords_intersection.end
        self.coords[contig] = coords_intersection
        if self.rev[contig] != self.rev[other_contig]:
            left_bases_lost, right_bases_lost = right_bases_lost, left_bases_lost

        self.coords[other_contig].start += left_bases_lost
        self.coords[other_contig].end -= right_bases_lost
        return True


    def merge_into(self, other_edge, contig):
        other_edge = copy.copy(other_edge)
        assert contig in self.names
        assert contig in other_edge.names
        new_contig1 = [x for x in self.names if x != contig][0]
        new_contig2 = [x for x in other_edge.names if x != contig][0]
        assert new_contig1 != new_contig2
        coords_intersection = self.coords[contig].intersection(other_edge.coords[contig])
        if coords_intersection is None:
            return False
        assert self._change_hit_coords_with_intersection(contig, coords_intersection)
        assert other_edge._change_hit_coords_with_intersection(contig, coords_intersection)
        self._make_contig_forwards(contig)
        other_edge._make_contig_forwards(contig)
        if self.names[0] == new_contig1:
            index_to_change = 1
        elif self.names[1] == new_contig1:
            index_to_change = 0
        else:
            raise Error('Unexpected error in edge.merge_into(). Cannot continue')
        assert self.names[index_to_change] == contig
        self.names[index_to_change] = new_contig2
        self.rev[new_contig2] = other_edge.rev[new_contig2]
        self.coords[new_contig2] = other_edge.coords[new_contig2]
        del self.rev[contig]
        del self.coords[contig]
        return True


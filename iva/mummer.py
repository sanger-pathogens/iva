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
import tempfile
import shutil
import pyfastaq
from iva import edge, common

class Error (Exception): pass

HIT_AT_START = 0
HIT_AT_END = 1
HIT_AT_BOTH_ENDS = 2
HIT_NO_ENDS = 3

def run_nucmer(query, ref, outfile, min_id=95, min_length=100, breaklen=200):
    query = os.path.abspath(query)
    ref = os.path.abspath(ref)
    outfile = os.path.abspath(outfile)
    tmpdir = tempfile.mkdtemp(prefix='tmp.run_nucmer.', dir=os.getcwd())
    original_dir = os.getcwd()
    os.chdir(tmpdir)
    script = 'run_nucmer.sh'
    f = pyfastaq.utils.open_file_write(script)
    print('nucmer --maxmatch -p p -b', breaklen, ref, query, file=f)
    print('delta-filter -i', min_id, '-l', min_length, 'p.delta > p.delta.filter', file=f)
    print('show-coords -dTlro p.delta.filter >', outfile, file=f)
    pyfastaq.utils.close(f)
    common.syscall('bash ' + script)
    os.chdir(original_dir)
    shutil.rmtree(tmpdir)


def file_reader(fname):
    f = pyfastaq.utils.open_file_read(fname)
    in_header = True

    for line in f:
        if in_header:
            if line.startswith('['):
                in_header = False
            continue
        yield NucmerHit(line)

    pyfastaq.utils.close(f)


class NucmerHit:
    def __init__(self, line):
        # [S1]  [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [FRM]   [TAGS]
        #1162    25768   24536   4   24607   24533   99.32   640851  24536   1   -1  MAL1    NODE_25757_length_24482_cov_18.920391   [CONTAINS]

        try:
            l = line.rstrip().split('\t')
            self.ref_start = int(l[0]) - 1
            self.ref_end = int(l[1]) - 1
            self.qry_start = int(l[2]) - 1
            self.qry_end = int(l[3]) - 1
            self.hit_length_ref = int(l[4])
            self.hit_length_qry = int(l[5])
            self.percent_identity = float(l[6])
            self.ref_length = int(l[7])
            self.qry_length = int(l[8])
            self.frame = int(l[9])
            self.ref_name = l[11]
            self.qry_name = l[12]

        except:
            raise Error('Error reading this nucmer line:\n' + line)


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def __hash__(self):
        return hash((self.ref_start, self.ref_end, self.qry_start, self.qry_end, self.hit_length_ref, self.hit_length_qry, self.percent_identity, self.ref_length, self.qry_length, self.frame, self.ref_name, self.qry_name))


    def _swap(self):
        self.ref_start, self.qry_start = self.qry_start, self.ref_start
        self.ref_end, self.qry_end = self.qry_end, self.ref_end
        self.hit_length_ref, self.hit_length_qry = self.hit_length_qry, self.hit_length_ref
        self.ref_length, self.qry_length = self.qry_length, self.ref_length
        self.ref_name, self.qry_name = self.qry_name, self.ref_name


    def qry_coords(self):
        return pyfastaq.intervals.Interval(min(self.qry_start, self.qry_end), max(self.qry_start, self.qry_end))


    def ref_coords(self):
        return pyfastaq.intervals.Interval(min(self.ref_start, self.ref_end), max(self.ref_start, self.ref_end))


    def on_same_strand(self):
        return (self.ref_start < self.ref_end) == (self.qry_start < self.qry_end)


    def sort(self):
        if self.ref_name > self.qry_name:
            self._swap()

        if self.ref_start > self.ref_end:
            self.ref_start, self.ref_end = self.ref_end, self.ref_start
            self.qry_start, self.qry_end = self.qry_end, self.qry_start

    def is_self_hit(self):
        return self.ref_name == self.qry_name \
                and self.ref_start == self.qry_start \
                and self.ref_end == self.qry_end \
                and self.percent_identity == 100


    def to_graph_edge(self, min_overlap_length=200, end_tolerance=50, min_identity=99):
        if self.percent_identity < min_identity or min(self.hit_length_qry, self.hit_length_ref) < min_overlap_length:
            return None
        ref_at_ends = self._is_at_ends(tolerance=end_tolerance)
        qry_at_ends = self._is_at_ends(use_qry=True, tolerance=end_tolerance)
        bad_values = set([HIT_NO_ENDS, HIT_AT_BOTH_ENDS])
        if ref_at_ends in bad_values or qry_at_ends in bad_values:
            return None

        good_values = set([HIT_AT_START, HIT_AT_END])
        assert qry_at_ends in good_values
        assert ref_at_ends in good_values
        qry_fwd_strand = self.qry_start < self.qry_end
        ref_fwd_strand = self.ref_start < self.ref_end

        if ((qry_fwd_strand == ref_fwd_strand) == (qry_at_ends == ref_at_ends)):
            return None
        elif qry_fwd_strand:
            if ref_fwd_strand:
                if qry_at_ends == HIT_AT_END:
                    assert ref_at_ends == HIT_AT_START
                    qry_first = True
                elif qry_at_ends == HIT_AT_START:
                    assert ref_at_ends == HIT_AT_END
                    qry_first = False
                else:
                    raise Error('Error in mummer.to_graph_edge(). Cannot continue')
            else:
                if qry_at_ends == HIT_AT_END:
                    assert ref_at_ends == HIT_AT_END
                    qry_first = True
                elif qry_at_ends == HIT_AT_START:
                    assert ref_at_ends == HIT_AT_START
                    qry_first = False
                else:
                    raise Error('Error in mummer.to_graph_edge(). Cannot continue')
        else:
            if ref_fwd_strand:
                if qry_at_ends == HIT_AT_END:
                    assert ref_at_ends == HIT_AT_END
                    qry_first = False
                elif qry_at_ends == HIT_AT_START:
                    assert ref_at_ends == HIT_AT_START
                    qry_first = True
                else:
                    raise Error('Error in mummer.to_graph_edge(). Cannot continue')
            else:
                if qry_at_ends == HIT_AT_END:
                    assert ref_at_ends == HIT_AT_START
                    qry_first = False
                elif qry_at_ends == HIT_AT_START:
                    assert ref_at_ends == HIT_AT_END
                    qry_first = True
                else:
                    raise Error('Error in mummer.to_graph_edge(). Cannot continue')
        if qry_first:
            return edge.Edge(self.qry_name, self.qry_start, self.qry_end, self.ref_name, self.ref_start, self.ref_end)
        else:
            return edge.Edge(self.ref_name, self.ref_start, self.ref_end, self.qry_name, self.qry_start, self.qry_end)


    def _is_at_ends(self, use_qry=False, tolerance=50):
        if use_qry:
            at_start = min(self.qry_start, self.qry_end) < tolerance
            at_end = self.qry_length - 1 - max(self.qry_start, self.qry_end) <= tolerance
        else:
            at_start = min(self.ref_start, self.ref_end) < tolerance
            at_end = self.ref_length - 1 - max(self.ref_start, self.ref_end) <= tolerance

        if at_start and at_end:
            return HIT_AT_BOTH_ENDS
        elif at_start:
            return HIT_AT_START
        elif at_end:
            return HIT_AT_END
        else:
            return HIT_NO_ENDS


    def __str__(self):
        return '\t'.join(str(x) for x in
            [self.ref_start,
            self.ref_end,
            self.qry_start,
            self.qry_end,
            self.hit_length_ref,
            self.hit_length_qry,
            '{0:.2f}'.format(self.percent_identity),
            self.ref_length,
            self.qry_length,
            self.frame,
            self.ref_name,
            self.qry_name])


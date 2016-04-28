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
import sys
import os
import filecmp
import pysam
from iva import assembly

modules_dir = os.path.dirname(os.path.abspath(assembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestAssembly(unittest.TestCase):
    def test_init_and_write(self):
        '''test assembly initialise and write to file'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'))
        tmp_out = 'tmp.assembly_test.out.fa'
        a.write_contigs_to_file(tmp_out)
        self.assertTrue(filecmp.cmp(tmp_out, os.path.join(data_dir, 'assembly_test.fa')))
        os.unlink(tmp_out)


    def test_map_reads(self):
        '''test _map_reads'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'))
        reads_prefix = os.path.join(data_dir, 'assembly_test.to_map')
        out_prefix = 'tmp.assembly_test.out'
        a._map_reads(reads_prefix + '_1.fastq', reads_prefix + '_2.fastq', out_prefix)

        # different smalt version output slightly different BAMs. Some columns
        # should never change, so check just those ones
        def get_sam_columns(bamfile):
            sams = []
            sam_reader = pysam.Samfile(bamfile, "rb")
            for sam in sam_reader.fetch(until_eof=True):
                if sam.is_unmapped:
                    refname = None
                else:
                    refname = sam_reader.getrname(sam.tid)
                sams.append((sam.qname, sam.flag, refname, sam.pos, sam.cigar, sam.seq))
            return sams

        expected = get_sam_columns(os.path.join(data_dir, 'assembly_test.mapped.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')


    def test_extend_contigs_with_bam(self):
        '''test _extend_contigs_with_bam'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'mapping_test.ref.trimmed.fa'), ext_min_cov=1, ext_min_ratio=1, ext_bases=10, max_insert=200)
        bam = os.path.join(data_dir, 'mapping_test.smalt.out.bam')
        out_prefix = 'tmp'
        a._extend_contigs_with_bam(bam, out_prefix, output_all_useful_reads=False)
        tmp_contigs = 'tmp.new_contigs.fa'
        tmp_reads_1 = out_prefix + '_1.fa'
        tmp_reads_2 = out_prefix + '_2.fa'
        a.write_contigs_to_file(tmp_contigs)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'mapping_test.ref.fa'), tmp_contigs))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'assembly_test.extend_kept_reads_1.fa'), tmp_reads_1))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'assembly_test.extend_kept_reads_2.fa'), tmp_reads_2))
        os.unlink(tmp_contigs)
        os.unlink(tmp_reads_1)
        os.unlink(tmp_reads_2)


    def test_contig_worth_extending(self):
        '''Test _contig_worth_extending'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'))
        self.assertTrue(a._contig_worth_extending('A'))
        self.assertTrue(a._contig_worth_extending('B'))
        self.assertTrue(a._contig_worth_extending('C'))
        a.contig_lengths['A'].append([100, 0, 0])
        self.assertTrue(a._contig_worth_extending('A'))
        a.contig_lengths['A'].append([100, 0, 0])
        self.assertFalse(a._contig_worth_extending('A'))
        a.contig_lengths['A'].append([100, 0, 0])
        self.assertFalse(a._contig_worth_extending('A'))
        a.contig_lengths['A'][-1] = [101, 1, 0]
        self.assertTrue(a._contig_worth_extending('A'))


    def test_worth_extending(self):
        '''Test worth_extending'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'))
        self.assertTrue(a._worth_extending())
        for x in ['A', 'B', 'C']:
            a.contig_lengths[x].append([100, 0, 0])
        self.assertTrue(a._worth_extending())
        for x in ['A', 'B', 'C']:
            a.contig_lengths[x].append([100, 0, 0])
        self.assertFalse(a._worth_extending())
        for x in ['A', 'B']:
            a.contig_lengths[x].append([100, 0, 0])
        self.assertFalse(a._worth_extending())
        a.contig_lengths['C'].append([100, 0, 0])
        self.assertFalse(a._worth_extending())


    def test_read_pair_extend(self):
        '''Test read_pair_extend'''
        ref = os.path.join(data_dir, 'assembly_test_read_pair_extend.ref.fasta')
        to_extend = os.path.join(data_dir, 'assembly_test_read_pair_extend.to_extend.fasta')
        reads_prefix = os.path.join(data_dir, 'assembly_test_read_pair_extend.ref.reads')
        a = assembly.Assembly(contigs_file=to_extend, strand_bias=0, verbose=4)
        a.read_pair_extend(reads_prefix, 'tmp.extend')
        self.assertTrue(len(a.contigs['1']) > 900)


    def test_get_ref_length(self):
        '''Test _get_ref_length'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'))
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'assembly_test.mapped.bam'), "rb")
        expected_lengths = [100] * 6 + [None, None]
        i = 0
        for sam in sam_reader.fetch(until_eof=True):
            self.assertEqual(expected_lengths[i], a._get_ref_length(sam_reader, sam))
            i += 1


    def test_get_ref_length_sam_pair(self):
        '''Test _get_ref_length_sam_pair'''
        previous_sam = None
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'))
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'assembly_test.mapped.bam'), "rb")
        expected_lengths = [100] * 3 + [None]
        i = 0
        for current_sam in sam_reader.fetch(until_eof=True):
            if previous_sam is None:
                previous_sam = current_sam
                continue
            self.assertEqual(expected_lengths[i], a._get_ref_length_sam_pair(sam_reader, current_sam, previous_sam))
            i += 1
            previous_sam = None


    def test_get_unmapped_pairs(self):
        '''Test _get_unmapped_pairs'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'))
        reads1 = os.path.join(data_dir, 'assembly_test.to_map_1.fastq')
        reads2 = os.path.join(data_dir, 'assembly_test.to_map_2.fastq')
        a._get_unmapped_pairs(reads1, reads2, 'tmp')
        self.assertTrue(filecmp.cmp('tmp_1.fa', os.path.join(data_dir, 'assembly_test_get_unmapped_pairs_1.fa'), shallow=False))
        self.assertTrue(filecmp.cmp('tmp_2.fa', os.path.join(data_dir, 'assembly_test_get_unmapped_pairs_2.fa'), shallow=False))
        os.unlink('tmp_1.fa')
        os.unlink('tmp_2.fa')


    def test_add_new_seed_contig(self):
        '''Test add_new_seed_contig'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test.fa'), verbose=3, min_clip=1, seed_min_cov=1, seed_min_kmer_count=1, seed_start_length=10, seed_overlap_length=5, seed_stop_length=20)
        reads1 = os.path.join(data_dir, 'assembly_test_add_new_seed_contig.reads_1.fa')
        reads2 = os.path.join(data_dir, 'assembly_test_add_new_seed_contig.reads_2.fa')
        new_contig_name = a.add_new_seed_contig(reads1, reads2)
        self.assertEqual('seeded.00001', new_contig_name)
        self.assertTrue(new_contig_name in a.contigs)


    def test_good_intervals_from_strand_coverage(self):
        '''Test good_intervals_from_strand_coverage'''
        fwd_cov = [0, 1, 1, 2, 5, 10, 100, 10, 10, 6, 0, 10, 10, 10, 5, 10]
        rev_cov = [0, 5, 5, 5, 5, 20, 10,  10, 10, 100, 9, 10, 10, 10, 5, 0]
        expected = [(3,5), (7,8), (11,14)]
        a = assembly.Assembly(strand_bias=0.2)
        got = a._good_intervals_from_strand_coverage(fwd_cov, rev_cov)
        self.assertListEqual(expected, got)


    def test_get_contig_order_by_orfs(self):
        '''Test get_contig_order_by_orfs'''
        a = assembly.Assembly(contigs_file=os.path.join(data_dir, 'assembly_test_order_by_orfs.fa'))
        got = a._get_contig_order_by_orfs(min_length=240)
        expected = [('1', True), ('3', False), ('4', False), ('6', False), ('2', False), ('5', False)]
        self.assertListEqual(expected, got)


    def test_trim_contig_for_strand_bias(self):
        '''Test _trim_contig_for_strand_bias'''
        # TODO
        pass


    def test_subcontigs_from_strand_bias(self):
        '''Test _subcontigs_from_strand_bias'''
        # TODO
        pass


    def test_trim_strand_biased_ends(self):
        '''Test _trim_strand_biased_ends'''
        # TODO
        pass


    def test_trim_contigs(self):
        '''Test trim_contigs'''
        # TODO
        pass


    def test_remove_contained_contigs(self):
        '''Test _remove_contained_contigs'''
        # TODO
        pass


    def test_coords_to_new_contig(self):
        '''Test _coords_to_new_contig'''
        # TODO
        pass


    def test_merge_overlapping_contigs(self):
        '''Test _merge_overlapping_contigs'''
        # TODO
        pass


    def test_contig_names_size_order(self):
        '''Test _contig_names_size_order'''
        # TODO
        pass


    def test_contig_contained_in_nucmer_hits(self):
        '''Test _contig_contained_in_nucmer_hits'''
        # TODO
        pass


    def test_remove_contig_from_nucmer_hits(self):
        '''Test _remove_contig_from_nucmer_hits'''
        # TODO
        pass


    def test_remove_contig(self):
        '''Test _remove_contig'''
        # TODO
        pass

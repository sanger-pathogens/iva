import unittest
import pickle
import shutil
import os
import filecmp
import pysam
import fastaq
from iva import mapping

modules_dir = os.path.dirname(os.path.abspath(mapping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


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


class TestMapping(unittest.TestCase):
    def test_smalt_in_path(self):
        '''Test that smalt is in the user's path'''
        assert(shutil.which('smalt') is not None)


    def test_smalt_in_path(self):
        '''Test that samtools is in the user's path'''
        assert(shutil.which('samtools') is not None)


    def test_map_reads(self):
        '''Test mapping reads'''
        ref = os.path.join(data_dir, 'mapping_test.ref.trimmed.fa')
        reads_prefix = os.path.join(data_dir, 'mapping_test.reads')
        out_prefix = 'tmp.out'
        mapping.map_reads(reads_prefix + '_1.fastq', reads_prefix + '_2.fastq', ref, out_prefix)
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test.smalt.out.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')


    def test_map_reads_and_sort(self):
        '''Test mapping reads and sort BAM'''
        ref = os.path.join(data_dir, 'mapping_test.ref.trimmed.fa')
        reads_prefix = os.path.join(data_dir, 'mapping_test.reads')
        out_prefix = 'tmp.out'
        mapping.map_reads(reads_prefix + '_1.fastq', reads_prefix + '_2.fastq', ref, out_prefix, sort=True, verbose=3)
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test.smalt.out.sorted.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')
        os.unlink(out_prefix + '.bam.bai')
        os.unlink(out_prefix + '.unsorted.bam')


    def test_map_reads_wth_flag(self):
        '''Test map_reads with required flag'''
        ref = os.path.join(data_dir, 'mapping_test.ref.trimmed.fa')
        reads_prefix = os.path.join(data_dir, 'mapping_test.reads')
        out_prefix = 'tmp.out'
        mapping.map_reads(reads_prefix + '_1.fastq', reads_prefix + '_2.fastq', ref, out_prefix, required_flag=12, verbose=3)
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test.smalt.out.flag12.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')


    def test_get_bam_region_coverage_rev(self):
        '''Test get_bam_region_coverage reverse strand'''
        bam = os.path.join(data_dir, 'mapping_test.smalt.out.sorted.bam')
        cov = mapping.get_bam_region_coverage(bam, 'ref', 190, rev=True, verbose=3)
        f = open(os.path.join(data_dir, 'mapping_test.smalt.out.sorted.bam.rev.cov'), 'rb')
        expected = pickle.load(f)
        f.close()
        self.assertListEqual(cov, expected)


    def test_get_bam_region_coverage_fwd(self):
        '''Test get_bam_region_coverage forward strand'''
        bam = os.path.join(data_dir, 'mapping_test.smalt.out.sorted.bam')
        cov = mapping.get_bam_region_coverage(bam, 'ref', 190, verbose=3)
        f = open(os.path.join(data_dir, 'mapping_test.smalt.out.sorted.bam.fwd.cov'), 'rb')
        expected = pickle.load(f)
        f.close()
        self.assertListEqual(cov, expected)


    def test_remove_indels(self):
        '''Test _remove_indels'''
        self.assertEqual('acgt', ''.join(mapping._remove_indels(list('ac+1Xgt'), '+')))
        self.assertEqual('ac+1Xgt', ''.join(mapping._remove_indels(list('ac+1Xgt'), '-')))
        self.assertEqual('ac-1Xgt', ''.join(mapping._remove_indels(list('ac-1Xgt'), '+')))
        self.assertEqual('acgt', ''.join(mapping._remove_indels(list('ac+2XXgt'), '+')))
        self.assertEqual('acgt', ''.join(mapping._remove_indels(list('ac+10XXXXXXXXXXgt'), '+')))
        self.assertEqual('acgt', ''.join(mapping._remove_indels(list('ac-10XXXXXXXXXXgt'), '-')))
        self.assertEqual('a-1Xcgt', ''.join(mapping._remove_indels(list('a-1Xc+1Xg+10XXXXXXXXXXt'), '+')))
        self.assertEqual('acgt', ''.join(mapping._remove_indels(list('+1Xacgt'), '+')))
        self.assertEqual('acgt', ''.join(mapping._remove_indels(list('acgt+1X'), '+')))


    def test_strip_mpileup_coverage_string(self):
        '''Test strip_mpileup_coverage_string'''
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('acg^[t'))
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('acgt$'))
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('*ac*gt'))
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('*a$c^[gt$'))
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('ac+1Xgt'))
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('acg+10XXXXXXXXXXt'))
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('ac-1Xgt'))
        self.assertEqual('acgt', mapping.strip_mpileup_coverage_string('acg-10XXXXXXXXXXt'))
        self.assertEqual('aaa', mapping.strip_mpileup_coverage_string('a-1Na^+a'))


    def test_consensus_base(self):
        '''Test consensus_base'''
        keys = ['A', 'C', 'G', 'T']
        self.assertEqual(None, mapping.consensus_base({}, keys))
        self.assertEqual('G', mapping.consensus_base({'A': 2, 'C': 2, 'G': 4}, keys, ratio=0.5))
        self.assertEqual('G', mapping.consensus_base({'A': 2, 'C': 1, 'G': 4, 'T':1}, keys, ratio=0.5))
        self.assertEqual(None, mapping.consensus_base({'A': 2, 'C': 1, 'G': 4, 'T':2}, keys, ratio=0.5))
        self.assertEqual('G', mapping.consensus_base({'A': 2, 'C': 1, 'G': 4, 'T':2}, keys, ratio=0.43))


    def test_consensus_base_both_strands(self):
        '''Test consensus_base_both_strands'''
        forward_keys = set(['A', 'C', 'G', 'T', 'N'])
        reverse_keys = set(['a', 'c', 'g', 't', 'n'])
        counts = [
            ({}, None),
            ({'A': 2, 'C': 2, 'G': 4, 'a': 2, 'c': 2, 'g': 5}, 'G'),
            ({'A': 2, 'C': 2, 'G': 3, 'a': 2, 'c': 2, 'g': 5}, None),
            ({'A': 2, 'C': 2, 'G': 4, 'a': 2, 'c': 2, 'g': 3}, None),
        ]

        for counts_dict, expected in counts:
            self.assertEqual(expected, mapping.consensus_base_both_strands(counts_dict, forward_keys, reverse_keys, ratio=0.5))


    def test_find_incorrect_ref_bases(self):
        '''Test find_incorrect_ref_bases'''
        bam = os.path.join(data_dir, 'mapping_test.find_incorrect_ref_bases.bam')
        ref = os.path.join(data_dir, 'mapping_test.find_incorrect_ref_bases.fasta')
        bad_bases = mapping.find_incorrect_ref_bases(bam, ref)
        expected = {'1': [(197, 'A', 'T'), (280, 'T', 'G')]}
        self.assertTrue(expected, bad_bases)


    def test_soft_clipped(self):
        '''Test soft_clipped'''
        expected = [
            (5, 0),
            (0, 0),
            (0, 0),
            (0, 5),
            (0, 0),
            None,
            (0, 0),
            (0, 0),
            (2, 0),
            (0, 1),
            None,
            None,
            (1, 1),
            (0, 1)
        ]

        sam_reader = pysam.Samfile(os.path.join(data_dir, 'mapping_test.smalt.out.bam'), "rb")
        i = 0
        for sam in sam_reader.fetch(until_eof=True):
            self.assertEqual(mapping.soft_clipped(sam), expected[i])
            i += 1


    def test_sam_to_fasta(self):
        '''Test sam_to_fasta'''
        expected_seqs = {}
        fastaq.tasks.file_to_dict(os.path.join(data_dir, 'mapping_test.reads_1.fasta'), expected_seqs)
        fastaq.tasks.file_to_dict(os.path.join(data_dir, 'mapping_test.reads_2.fasta'), expected_seqs)
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'mapping_test.smalt.out.bam'), "rb")
        for sam in sam_reader.fetch(until_eof=True):
            fa = mapping.sam_to_fasta(sam)
            self.assertTrue(fa.id in expected_seqs)
            self.assertEqual(fa, expected_seqs[fa.id])


    def test_can_extend(self):
        '''Test can_extend'''
        expected = [
            (True, False),
            (False, False),
            (False, False),
            (False, True),
            (False, False),
            (False, False),
            (False, False),
            (False, False),
            (True, False),
            (False, False),
            (False, False),
            (False, False),
            (False, False),
            (False, False),

        ]

        sam_reader = pysam.Samfile(os.path.join(data_dir, 'mapping_test.smalt.out.bam'), "rb")
        i = 0
        for sam in sam_reader.fetch(until_eof=True):
            self.assertEqual(mapping._can_extend(sam, 190, min_clip=2), expected[i])
            i += 1


    def test_get_pair_type(self):
        '''Test get_pair_type'''
        expected = [
            (mapping.CAN_EXTEND_LEFT, mapping.KEEP),
            (mapping.KEEP, mapping.CAN_EXTEND_RIGHT),
            (mapping.KEEP, mapping.KEEP),
            (mapping.NOT_USEFUL, mapping.NOT_USEFUL),
            (mapping.CAN_EXTEND_LEFT, mapping.KEEP),
            (mapping.BOTH_UNMAPPED, mapping.BOTH_UNMAPPED),
            (mapping.NOT_USEFUL, mapping.NOT_USEFUL)
        ]

        sam_reader = pysam.Samfile(os.path.join(data_dir, 'mapping_test.smalt.out.bam'), "rb")
        previous_sam = None
        i = 0
        for sam in sam_reader.fetch(until_eof=True):
            if previous_sam is None:
                previous_sam = sam
                continue

            types = mapping.get_pair_type(previous_sam, sam, 190, 1000, min_clip=2)
            self.assertEqual(types, expected[i])
            i += 1
            previous_sam = None


    def test_get_ref_name(self):
        '''Test get_ref_name'''
        expected = ['ref'] * 14
        for i in ([5,10,11]):
            expected[i] = None
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'mapping_test.smalt.out.bam'), "rb")
        i = 0
        for sam in sam_reader.fetch(until_eof=True):
            self.assertEqual(mapping.get_ref_name(sam, sam_reader), expected[i])
            i += 1


    def test_bam_file_to_fasta_pair_files(self):
        '''Test bam_file_to_fasta_pair_files'''
        tmp1 = 'tmp.to_fasta_1.fa'
        tmp2 = 'tmp.to_fasta_2.fa'
        mapping.bam_file_to_fasta_pair_files(os.path.join(data_dir, 'mapping_test.smalt.out.bam'), tmp1, tmp2)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'mapping_test.reads_1.fasta'), tmp1))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'mapping_test.reads_2.fasta'), tmp2))
        os.unlink(tmp1)
        os.unlink(tmp2)


    def test_bam_file_to_fasta_pair_files_region(self):
        '''Test bam_file_to_fasta_pair_files with a region'''
        tmp1 = 'tmp.to_fasta_1.fa'
        tmp2 = 'tmp.to_fasta_2.fa'
        mapping.bam_file_to_fasta_pair_files(os.path.join(data_dir, 'mapping_test.smalt.out.sorted.bam'), tmp1, tmp2, chromosome='ref', start=25, end=150)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'mapping_test.bam_to_region_1.fa'), tmp1))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'mapping_test.bam_to_region_2.fa'), tmp2))
        os.unlink(tmp1)
        os.unlink(tmp2)


    def test_bam_file_to_region_fasta(self):
        '''Test bam_file_to_region_fasta'''
        tmp = 'tmp.to_fasta.fa'
        bam = os.path.join(data_dir, 'mapping_test.smalt.out.sorted.bam')
        mapping.bam_file_to_region_fasta(bam, tmp, 'ref', start=42, end=142)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'mapping_test.bam_to_region.fasta'), tmp))
        os.unlink(tmp)


    def test_bam_to_fasta(self):
        '''Test bam_to_fasta'''
        # TODO
        pass


    def test_total_ref_length_from_bam(self):
        '''Test _total_ref_length_from_bam'''
        bam = os.path.join(data_dir, 'mapping_test.total_ref_length_from_bam.bam')
        self.assertEqual(300, mapping._total_ref_length_from_bam(bam))


    def _mean_read_length(self):
        '''Test _mean_read_length'''
        bam = os.path.join(data_dir, 'mapping_test.mean_read_length.bam')
        lengths = [19, 18, 20, 17, 20, 20, 20, 20]
        self.assertEqual(19, mapping._mean_read_length(bam, head=1))
        self.assertEqual(18, mapping._mean_read_length(bam, head=2))
        self.assertEqual(int(sum(lengths) / len(lengths)), mapping._mean_read_length(bam))


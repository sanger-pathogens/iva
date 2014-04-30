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

import unittest
import os
import filecmp
import pysam
from iva import qc, mummer
import pyfastaq

modules_dir = os.path.dirname(os.path.abspath(qc.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestQc(unittest.TestCase):
    def setUp(self):
        ref_embl = os.path.join(data_dir, 'qc_test.dummy.embl')
        assembly_fasta = os.path.join(data_dir, 'qc_test.dummy.assembly.fasta')
        reads_1 = os.path.join(data_dir,'qc_test.reads_1.fq')
        reads_2 = os.path.join(data_dir,'qc_test.reads_2.fq')
        self.qc = qc.Qc(assembly_fasta, 'tmp.qc', embl_dir=ref_embl, reads_fwd=reads_1, reads_rev=reads_2)
        self.qc.assembly_is_empty = False


    def tearDown(self):
        files_to_clean = [
            'tmp.qc.assembly_vs_ref.coords',
            'tmp.qc.ref_cds_seqs.fa',
            'tmp.qc.ref_cds_seqs_mapped_to_assembly.coords',
            'tmp.qc.reads_mapped_to_ref.bam',
            'tmp.qc.reads_mapped_to_ref.bam.bai',
            'tmp.qc.reference.fa',
            'tmp.qc.reference.fa.fai',
            'tmp.qc.reference.gff',
            'tmp.qc.assembly_contigs_hit_ref.fasta',
            'tmp.qc.assembly_contigs_not_hit_ref.fasta',
        ]
        for f in files_to_clean:
            if os.path.exists(f):
                os.unlink(f)


    def test_ids_in_order_from_fai(self):
        '''test _ids_in_order_from_fai'''
        expected = ['A0', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        got = self.qc._ids_in_order_from_fai(os.path.join(data_dir, 'qc_test.reference.fa.fai'))
        self.assertEqual(expected, got)


    def test_get_ref_cds_from_gff(self):
        '''test _get_ref_cds_from_gff'''
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.get_ref_cds_from_gff.gff')
        coords = self.qc._get_ref_cds_from_gff()
        expected = {
            'contig1': [(pyfastaq.intervals.Interval(1200, 1499), '+'),
                        (pyfastaq.intervals.Interval(2999, 3901), '+'),
                        (pyfastaq.intervals.Interval(4999, 5499), '+'),
                        (pyfastaq.intervals.Interval(6999, 7599), '+'),
                        (pyfastaq.intervals.Interval(10010, 10041), '-'),
                        (pyfastaq.intervals.Interval(10089, 10999), '-')],
            'contig2': [(pyfastaq.intervals.Interval(120, 149), '+'),
                        (pyfastaq.intervals.Interval(299, 391), '+')]
        }

        self.assertEqual(coords, expected)


    def test_write_cds_seqs(self):
        '''test _write_cds_seqs'''
        seq = pyfastaq.sequences.Fasta('seq', 'AGGTGTCACGTGTGTGTCATTCAGGGCA')
        cds_list = [(pyfastaq.intervals.Interval(1,3), '+'),
                    (pyfastaq.intervals.Interval(7, 9), '-'),
        ]
        outfile = 'tmp.write_cds_seqs.fa'
        f = pyfastaq.utils.open_file_write(outfile)
        self.qc._write_cds_seqs(cds_list, seq, f)
        pyfastaq.utils.close(f)
        self.assertTrue(filecmp.cmp(outfile, os.path.join(data_dir, 'qc_test.write_cds_seqs.fa'), shallow=False))
        os.unlink(outfile)


    def test_gff_and_fasta_to_cds(self):
        '''test _gff_and_fasta_to_cds'''
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.gff_and_fasta_to_cds.in.gff')
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.gff_and_fasta_to_cds.in.fasta')
        self.qc._set_ref_fa_data()
        expected_out = os.path.join(data_dir, 'qc_test.gff_and_fasta_to_cds.out.fasta')
        got_out = 'tmp.qc.ref_cds_seqs.fa'
        self.qc._gff_and_fasta_to_cds()
        self.assertTrue(filecmp.cmp(expected_out, got_out, shallow=False))
        os.unlink(got_out)


    def test_map_cds_to_assembly(self):
        '''test _map_cds_to_assembly'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        expected_out = os.path.join(data_dir, 'qc_test.map_cds_to_assembly.coords')
        self.qc._map_cds_to_assembly()
        # need to ignore first line of coords file because it has full paths
        # to input files
        with open(expected_out) as f:
            expected = [line.rstrip() for line in f.readlines()][1:]
        with open(self.qc.cds_nucmer_coords_in_assembly) as f:
            got = [line.rstrip() for line in f.readlines()][1:]
        self.assertEqual(expected, got)


    def test_mummer_coords_file_to_dict(self):
        '''test _mummer_coords_file_to_dict'''
        expected = {

            'qry1': [
                mummer.NucmerHit('\t'.join(['1', '100', '51', '150', '100', '100', '100.00', '1008', '762', '1', '1', 'ref1', 'qry1'])),
                mummer.NucmerHit('\t'.join(['300', '500', '351', '550', '100', '100', '100.00', '1008', '762', '1', '1', 'ref2', 'qry1']))
            ],
            'qry2': [mummer.NucmerHit('\t'.join(['1', '1000', '1', '1000', '1000', '1000', '100.00', '1000', '1542', '1', '1', 'ref2', 'qry2']))]
        }

        got = self.qc._mummer_coords_file_to_dict(os.path.join(data_dir, 'qc_test.mummer_coords_file_to_dict.coords'))
        self.assertEqual(expected, got)


    def test_has_orf(self):
        '''test _has_orf'''
        fa = pyfastaq.sequences.Fasta('a', 'ggggTAAxTAAxTAATTAxTTAxTTAgg')
        self.assertFalse(self.qc._has_orf(fa, 0, 30, 30))
        self.assertTrue(self.qc._has_orf(fa, 0, 30, 20))


    def test_calculate_cds_assembly_stats(self):
        '''test _calculate_cds_assembly_stats'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        self.qc._calculate_cds_assembly_stats()
        expected = {
            'A:23-784:+': {
                'strand': '+',
                'length_in_ref': 762,
                'ref_coords': pyfastaq.intervals.Interval(22, 783),
                'ref_name': 'A',
                'bases_assembled': 762,
                'number_of_contig_hits': 1,
                'assembled': True,
                'assembled_ok': True,
            },
            'B:3-1733:+': {
                'strand': '+',
                'length_in_ref': 1731,
                'ref_coords': pyfastaq.intervals.Interval(2, 1732),
                'ref_name': 'B',
                'bases_assembled': 1731,
                'number_of_contig_hits': 1,
                'assembled': True,
                'assembled_ok': True,
            },
            'C:3-1385:+': {
                'strand': '+',
                'length_in_ref': 1383,
                'ref_coords': pyfastaq.intervals.Interval(2, 1384),
                'ref_name': 'C',
                'bases_assembled': 1383,
                'number_of_contig_hits': 1,
                'assembled': True,
                'assembled_ok': True,
            },
            'D:1-1542:+': {
                'strand': '+',
                'length_in_ref': 1542,
                'ref_coords': pyfastaq.intervals.Interval(0, 1541),
                'ref_name': 'D',
                'bases_assembled': 1000,
                'number_of_contig_hits': 1,
                'assembled': False,
                'assembled_ok': False,
            },
            'E:3-719:+': {
                'strand': '+',
                'length_in_ref': 717,
                'ref_coords': pyfastaq.intervals.Interval(2, 718),
                'ref_name': 'E',
                'bases_assembled': 499,
                'number_of_contig_hits': 2,
                'assembled': False,
                'assembled_ok': False,
            },
            'E:290-793:-': {
                'strand': '-',
                'length_in_ref': 504,
                'ref_coords': pyfastaq.intervals.Interval(289, 792),
                'ref_name': 'E',
                'bases_assembled': 301,
                'number_of_contig_hits': 1,
                'assembled': False,
                'assembled_ok': False,
            },
            'F:25-2298:+': {
                'strand': '+',
                'length_in_ref': 2274,
                'ref_coords': pyfastaq.intervals.Interval(24, 2297),
                'ref_name': 'F',
                'bases_assembled': 2274,
                'number_of_contig_hits': 2,
                'assembled': True,
                'assembled_ok': False,
            },
            'G:19-2175:+': {
                'strand': '+',
                'length_in_ref': 2157,
                'ref_coords': pyfastaq.intervals.Interval(18, 2174),
                'ref_name': 'G',
                'assembled': False,
                'assembled_ok': False,
            },
            'H:1-2307:+': {
                'strand': '+',
                'length_in_ref': 2307,
                'ref_coords': pyfastaq.intervals.Interval(0, 2306),
                'ref_name': 'H',
                'assembled': False,
                'assembled_ok': False,
            }
        }

        self.maxDiff = None
        self.assertEqual(expected, self.qc.cds_assembly_stats)


    def test_get_contig_hits_to_reference(self):
        '''test _get_contig_hits_to_reference'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        self.qc._get_contig_hits_to_reference()
        expected = {
            'A:10-1017:+': [mummer.NucmerHit('\t'.join(['10', '1017', '1', '1008', '1008', '1008', '100.00', '1027', '1008', '1', '+', 'A', 'A:10-1017:+'])),
                            mummer.NucmerHit('\t'.join(['10', '240', '1', '231', '231', '231', '100.00', '240', '1008', '1', '+', 'A0', 'A:10-1017:+']))],
            'B:1-1778:-': [mummer.NucmerHit('\t'.join(['1', '1778', '1778', '1', '1778', '1778', '100.00', '1778', '1778', '1', '+', 'B', 'B:1-1778:-']))],
            'C:1-1413:+,E:1-200:-': [mummer.NucmerHit('\t'.join(['1', '1413', '1', '1413', '1413', '1413', '100.00', '1413', '1613', '1', '+', 'C', 'C:1-1413:+,E:1-200:-'])),
                                     mummer.NucmerHit('\t'.join(['1', '200', '1414', '1613', '200', '200', '100.00', '890', '1613', '1', '+', 'E', 'C:1-1413:+,E:1-200:-']))],
            'D:1-1000:+': [mummer.NucmerHit('\t'.join(['1', '1000', '1', '1000', '1000', '1000', '100.00', '1565', '1000', '1', '+', 'D', 'D:1-1000:+']))],
            'E:400-700:-': [mummer.NucmerHit('\t'.join(['400', '700', '301', '1', '301', '301', '100.00', '890', '301', '1', '+', 'E', 'E:400-700:-']))],
            'F:1-2341:+': [mummer.NucmerHit('\t'.join(['1', '2341', '1', '2341', '2341', '2341', '100.00', '2341', '2341', '1', '+', 'F', 'F:1-2341:+']))],
            'F:1-2341:-':  [mummer.NucmerHit('\t'.join(['1', '2341', '2341', '1', '2341', '2341', '100.00', '2341', '2341', '1', '+', 'F', 'F:1-2341:-']))]
        }
        self.assertDictEqual(expected, self.qc.assembly_vs_ref_mummer_hits)


    def test_write_fasta_contigs_hit_ref(self):
        '''test _write_fasta_contigs_hit_ref'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        self.qc._get_contig_hits_to_reference()
        self.qc._write_fasta_contigs_hit_ref()
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'qc_test.write_fasta_contigs_hit_ref.fa'), 'tmp.qc.assembly_contigs_hit_ref.fasta', shallow=False))


    def test_write_fasta_contigs_not_hit_ref(self):
        '''test _write_fasta_contigs_not_hit_ref'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        self.qc._get_contig_hits_to_reference()
        self.qc._write_fasta_contigs_not_hit_ref()
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'qc_test.write_fasta_contigs_not_hit_ref.fa'), 'tmp.qc.assembly_contigs_not_hit_ref.fasta', shallow=False))


    def test_hash_nucmer_hits_by_ref(self):
        '''test _hash_nucmer_hits_by_ref'''
        hit1 = mummer.NucmerHit('\t'.join(['1', '10', '20', '30', '10', '10', '100.00', '1008', '762', '1', '1', 'ref1', 'qry1']))
        hit2 = mummer.NucmerHit('\t'.join(['1', '10', '20', '30', '10', '10', '100.00', '1008', '762', '1', '1', 'ref1', 'qry1']))
        hit3 = mummer.NucmerHit('\t'.join(['1', '10', '20', '30', '10', '10', '100.00', '1008', '762', '1', '1', 'ref1', 'qry1']))
        hit4 = mummer.NucmerHit('\t'.join(['1', '10', '20', '30', '10', '10', '100.00', '1008', '762', '1', '1', 'ref2', 'qry1']))

        input_dict = {'x': [hit1, hit2, hit3, hit4]}

        expected = {
            'ref1': [hit1, hit2, hit3],
            'ref2': [hit4]
        }

        got = self.qc._hash_nucmer_hits_by_ref(input_dict)
        self.assertEqual(expected, got)


    def test_calculate_refseq_assembly_stats(self):
        '''test _calculate_refseq_assembly_stats'''
        self.qc.ref_ids = ['ref1', 'ref2', 'ref3', 'ref4']
        self.qc.ref_lengths = {x:1000 for x in self.qc.ref_ids}
        self.qc.assembly_vs_ref_mummer_hits = {'x': [
                mummer.NucmerHit('\t'.join(['1', '1000', '1', '1000', '1000', '1000', '100.00', '1000', '1000', '1', '1', 'ref1', 'ctg1'])),
                mummer.NucmerHit('\t'.join(['1', '800', '1', '800', '800', '800', '100.00', '800', '800', '1', '1', 'ref2', 'ctg2'])),
                mummer.NucmerHit('\t'.join(['1', '500', '1', '500', '1000', '500', '100.00', '500', '500', '1', '1', 'ref3', 'ctg3.1'])),
                mummer.NucmerHit('\t'.join(['501', '1000', '1', '500', '1000', '500', '100.00', '500', '500', '1', '1', 'ref3', 'ctg3.2']))
            ]
        }

        expected = {
            'ref1': {
                'hits': 1,
                'bases_assembled': 1000,
                'assembled': True,
                'assembled_ok': True,
                'longest_matching_contig': 1000,
            },
            'ref2': {
                'hits': 1,
                'bases_assembled': 800,
                'assembled': False,
                'assembled_ok': False,
                'longest_matching_contig': 800,
            },
            'ref3': {
                'hits': 2,
                'bases_assembled': 1000,
                'assembled': True,
                'assembled_ok': False,
                'longest_matching_contig': 500,
            },
            'ref4': {
                'hits': 0,
                'bases_assembled': 0,
                'assembled': False,
                'assembled_ok': False,
                'longest_matching_contig': 0,
            }
        }

        self.qc._calculate_refseq_assembly_stats()
        self.maxDiff = None
        self.assertEqual(expected, self.qc.refseq_assembly_stats)


    def test_invert_list(self):
        '''test _invert_list'''
        self.assertEqual(self.qc._invert_list([], 42), [pyfastaq.intervals.Interval(0,41)])
        coords = [
            [],
            [pyfastaq.intervals.Interval(0, 41)],
            [pyfastaq.intervals.Interval(1, 41)],
            [pyfastaq.intervals.Interval(0, 40)],
            [pyfastaq.intervals.Interval(10, 30)],
            [pyfastaq.intervals.Interval(5, 10), pyfastaq.intervals.Interval(20, 30)],
        ]
        expected = [
            [pyfastaq.intervals.Interval(0, 41)],
            [],
            [pyfastaq.intervals.Interval(0, 0)],
            [pyfastaq.intervals.Interval(41, 41)],
            [pyfastaq.intervals.Interval(0, 9), pyfastaq.intervals.Interval(31, 41)],
            [pyfastaq.intervals.Interval(0, 4), pyfastaq.intervals.Interval(11, 19), pyfastaq.intervals.Interval(31, 41)],
        ]
        assert len(coords) == len(expected)
        for i in range(len(coords)):
            self.assertEqual(self.qc._invert_list(coords[i], 42), expected[i])


    def test_calculate_ref_positions_covered_by_contigs(self):
        '''test _calculate_ref_positions_covered_by_contigs'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        self.qc._get_contig_hits_to_reference()
        self.qc._calculate_ref_positions_covered_by_contigs()
        expected = {
            'A0': [pyfastaq.intervals.Interval(9, 239)],
            'A': [pyfastaq.intervals.Interval(9, 1016)],
            'B': [pyfastaq.intervals.Interval(0, 1777)],
            'C': [pyfastaq.intervals.Interval(0, 1412)],
            'D': [pyfastaq.intervals.Interval(0, 999)],
            'E': [pyfastaq.intervals.Interval(0, 199), pyfastaq.intervals.Interval(399, 699)],
            'F': [pyfastaq.intervals.Interval(0, 2340)],
        }
        self.assertEqual(expected['F'], self.qc.ref_pos_covered_by_contigs['F'])


    def test_get_overlapping_qry_hits(self):
        '''test _get_overlapping_qry_hits'''
        h1 = mummer.NucmerHit('\t'.join(['1', '10', '1', '10', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h2 = mummer.NucmerHit('\t'.join(['1', '10', '21', '30', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h3 = mummer.NucmerHit('\t'.join(['1', '10', '25', '35', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h4 = mummer.NucmerHit('\t'.join(['1', '10', '29', '40', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h5 = mummer.NucmerHit('\t'.join(['1', '10', '70', '90', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        hits = [h1, h2, h3, h4, h5]

        expected = [
            [],
            [h3, h4],
            [h2, h4],
            [h2, h3],
            []
        ]

        self.assertEqual(len(hits), len(expected))

        for i in range(len(hits)):
            self.assertEqual(self.qc._get_overlapping_qry_hits(hits, hits[i]), expected[i])


    def test_get_unique_and_repetitive_from_contig_hits(self):
        '''test _get_unique_and_repetitive_from_contig_hits'''
        h1 = mummer.NucmerHit('\t'.join(['1', '10', '1', '10', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h2 = mummer.NucmerHit('\t'.join(['1', '10', '21', '30', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h3 = mummer.NucmerHit('\t'.join(['1', '10', '25', '35', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h4 = mummer.NucmerHit('\t'.join(['1', '10', '29', '40', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h5 = mummer.NucmerHit('\t'.join(['1', '10', '70', '90', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        hits = [h1, h2, h3, h4, h5]
        expect_repetitive = [h2, h3, h4]
        expect_unique = [h1, h5]
        got_unique, got_repetitive  = self.qc._get_unique_and_repetitive_from_contig_hits(hits)
        self.assertEqual(expect_unique, got_unique)
        self.assertEqual(expect_repetitive, got_repetitive)


    def test_contig_placement_in_reference(self):
        '''test _contig_placement_in_reference'''
        h1 = mummer.NucmerHit('\t'.join(['1', '90', '100', '10', '100', '100', '100.00', '100', '100', '1', '+', 'ref1', 'qry1']))
        h2 = mummer.NucmerHit('\t'.join(['17', '36', '21', '40', '100', '100', '100.00', '100', '100', '1', '+', 'ref2', 'qry1']))
        expected = [(pyfastaq.intervals.Interval(9, 99), 'ref1', pyfastaq.intervals.Interval(0, 89), False, False)]
        self.assertEqual(self.qc._contig_placement_in_reference([h1]), expected)

        expected = [
            (pyfastaq.intervals.Interval(9, 99), 'ref1', pyfastaq.intervals.Interval(0, 89), False, True),
            (pyfastaq.intervals.Interval(20, 39), 'ref2', pyfastaq.intervals.Interval(16, 35), True, True)
        ]

        self.assertEqual(self.qc._contig_placement_in_reference([h1, h2]), expected)


    def test_calculate_contig_placement(self):
        '''test calculate_contig_placement'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        self.qc._calculate_contig_placement()
        expected_placement = {
            'A:10-1017:+': [(pyfastaq.intervals.Interval(0, 230), 'A0', pyfastaq.intervals.Interval(9, 239), True, True),
                            (pyfastaq.intervals.Interval(0, 1007), 'A', pyfastaq.intervals.Interval(9, 1016), True, True)],
            'B:1-1778:-': [(pyfastaq.intervals.Interval(0, 1777), 'B', pyfastaq.intervals.Interval(0, 1777), False, False)],
            'C:1-1413:+,E:1-200:-': [(pyfastaq.intervals.Interval(0, 1412), 'C', pyfastaq.intervals.Interval(0, 1412), True, False),
                                     (pyfastaq.intervals.Interval(1413, 1612), 'E', pyfastaq.intervals.Interval(0, 199), True, False)],
            'D:1-1000:+': [(pyfastaq.intervals.Interval(0, 999), 'D', pyfastaq.intervals.Interval(0, 999), True, False)],
            'E:400-700:-': [(pyfastaq.intervals.Interval(0, 300), 'E', pyfastaq.intervals.Interval(399, 699), False, False)],
            'F:1-2341:+': [(pyfastaq.intervals.Interval(0, 2340), 'F', pyfastaq.intervals.Interval(0, 2340), True, False)],
            'F:1-2341:-': [(pyfastaq.intervals.Interval(0, 2340), 'F', pyfastaq.intervals.Interval(0, 2340), False, False)],
        }
        self.assertEqual(expected_placement, self.qc.contig_placement)


    def test_get_R_plot_contig_order_from_contig_placement(self):
        '''test _get_R_plot_contig_order_from_contig_placement'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.reference.fa')
        self.qc._set_ref_fa_data()
        self.qc.ref_gff = os.path.join(data_dir, 'qc_test.reference.cds.gff')
        self.qc.assembly_fasta =  os.path.join(data_dir, 'qc_test.assembly.fasta')
        self.qc._calculate_contig_placement()
        names = self.qc._get_R_plot_contig_order_from_contig_placement()
        self.assertEqual(names, ['A:10-1017:+', 'B:1-1778:-', 'C:1-1413:+,E:1-200:-', 'D:1-1000:+', 'E:400-700:-', 'F:1-2341:-', 'F:1-2341:+'])


    def test_calculate_ref_read_coverage(self):
        '''test _calculate_ref_read_coverage'''
        self.qc.ref_fasta = os.path.join(data_dir, 'qc_test.calculate_ref_read_coverage.ref.fa')
        self.qc._set_ref_fa_data()
        self.qc.reads_fwd = os.path.join(data_dir, 'qc_test.calculate_ref_read_coverage.reads_1.fq')
        self.qc.reads_rev = os.path.join(data_dir, 'qc_test.calculate_ref_read_coverage.reads_2.fq')
        fwd_mpileup = {
            1: 1,
            2: 1,
            3: 1,
            4: 1,
            5: 1,
            6: 1,
            7: 1,
            8: 1,
            9: 1,
            10: 1,
            11: 1,
            12: 1,
            13: 1,
            14: 2,
            15: 2,
            16: 2,
            17: 2,
            18: 2,
            19: 2,
            20: 2,
            21: 2,
            22: 2,
            23: 2,
            24: 2,
            25: 2,
            26: 2,
            27: 2,
            28: 2,
            29: 2,
            30: 2,
            31: 2,
            32: 2,
            33: 2,
            34: 2,
            35: 2,
            36: 2,
            37: 2,
            38: 2,
            39: 2,
            40: 2,
            41: 1,
            42: 1,
            43: 1,
            44: 1,
            45: 1,
            46: 1,
            47: 1,
            48: 1,
            49: 1,
            50: 1,
            51: 1,
            52: 1,
            53: 1
        }
        expected_fwd = []
        for i in range(100):
            expected_fwd.append(fwd_mpileup.get(i+1, 0))

        rev_mpileup = {
            53: 1,
            54: 1,
            55: 1,
            56: 1,
            57: 1,
            58: 1,
            59: 1,
            60: 1,
            61: 2,
            62: 2,
            63: 2,
            64: 2,
            65: 2,
            66: 2,
            67: 2,
            68: 2,
            69: 2,
            70: 2,
            71: 2,
            72: 2,
            73: 2,
            74: 2,
            75: 2,
            76: 2,
            77: 2,
            78: 2,
            79: 2,
            80: 2,
            81: 2,
            82: 2,
            83: 2,
            84: 2,
            85: 2,
            86: 2,
            87: 2,
            88: 2,
            89: 2,
            90: 2,
            91: 2,
            92: 2,
            93: 1,
            94: 1,
            95: 1,
            96: 1,
            97: 1,
            98: 1,
            99: 1,
            100: 1
        }
        expected_rev = []
        for i in range(100):
            expected_rev.append(rev_mpileup.get(i+1, 0))

        expected_fwd = {'1': expected_fwd, '2': [0] * 100}
        expected_rev = {'1': expected_rev, '2': [0] * 100}
        self.qc._calculate_ref_read_coverage()
        self.assertEqual(self.qc.ref_coverage_fwd, expected_fwd)
        self.assertEqual(self.qc.ref_coverage_rev, expected_rev)


    def test_coverage_list_to_low_cov_intervals(self):
        '''test _coverage_list_to_low_cov_intervals'''
        l = [0, 1, 4, 5, 6, 6, 5, 0, 5, 6, 2, 1]
        expected = [
            pyfastaq.intervals.Interval(0,2),
            pyfastaq.intervals.Interval(7,7),
            pyfastaq.intervals.Interval(10,11),
        ]
        got = self.qc._coverage_list_to_low_cov_intervals(l)
        self.assertEqual(expected, got)


    def test_calculate_ref_read_region_coverage(self):
        '''test _calculate_ref_read_region_coverage'''
        self.qc.ref_ids = ['ref1', 'ref2', 'ref3', 'ref4']
        self.qc.ref_lengths = {x:10 for x in self.qc.ref_ids}
        self.qc.ref_coverage_fwd = {
            'ref1': [0, 5, 5, 5, 5, 5, 0, 5, 5, 5],
            'ref2': [0] * 10,
            'ref3': [0] * 10,
            'ref4': [5] * 10,
        }
        self.qc.ref_coverage_rev = {
            'ref1': [3, 5, 5, 0, 5, 5, 5, 5, 5, 5],
            'ref2': [0] * 10,
            'ref3': [5] * 10,
            'ref4': [0] * 10,
        }
        self.qc._calculate_ref_read_region_coverage()

        expected_low_cov_ref_regions_fwd = {
            'ref1': [pyfastaq.intervals.Interval(0, 0), pyfastaq.intervals.Interval(6, 6)],
            'ref2': [pyfastaq.intervals.Interval(0, 9)],
            'ref3': [pyfastaq.intervals.Interval(0, 9)],
            'ref4': [],
        }
        self.assertEqual(expected_low_cov_ref_regions_fwd, self.qc.low_cov_ref_regions_fwd)

        expected_low_cov_ref_regions_rev = {
            'ref1': [pyfastaq.intervals.Interval(0, 0), pyfastaq.intervals.Interval(3, 3)],
            'ref2': [pyfastaq.intervals.Interval(0, 9)],
            'ref3': [],
            'ref4': [pyfastaq.intervals.Interval(0, 9)],
        }
        self.assertEqual(expected_low_cov_ref_regions_rev, self.qc.low_cov_ref_regions_rev)

        expected_low_cov_ref_regions = {
            'ref1': [pyfastaq.intervals.Interval(0, 0), pyfastaq.intervals.Interval(3, 3), pyfastaq.intervals.Interval(6, 6)],
            'ref2': [pyfastaq.intervals.Interval(0, 9)],
            'ref3': [pyfastaq.intervals.Interval(0, 9)],
            'ref4': [pyfastaq.intervals.Interval(0, 9)],
        }

        expected_ok_cov_ref_regions = {
            'ref1': [pyfastaq.intervals.Interval(1, 2), pyfastaq.intervals.Interval(4, 5), pyfastaq.intervals.Interval(7, 9)],
            'ref2': [],
            'ref3': [],
            'ref4': [],
        }


    def test_write_ref_coverage_to_files_for_R(self):
        '''test _write_ref_coverage_to_files_for_R'''
        self.qc.ref_ids = ['1', '2']
        self.qc.ref_coverage_fwd = {
            '1': [0, 1, 1, 0],
            '2': [2, 42, 0]
        }

        self.qc.ref_coverage_rev = {
            '1': [0, 42, 1, 0],
            '2': [2, 2, 1]
        }
        self.qc._write_ref_coverage_to_files_for_R('tmp.cov_for_R')
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'qc_test.write_ref_coverage_to_files_for_R.fwd'), 'tmp.cov_for_R.fwd'))
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'qc_test.write_ref_coverage_to_files_for_R.rev'), 'tmp.cov_for_R.rev'))
        os.unlink('tmp.cov_for_R.fwd')
        os.unlink('tmp.cov_for_R.rev')


    def test_cov_to_R_string(self):
        '''test _cov_to_R_string'''
        intervals = [pyfastaq.intervals.Interval(0, 42), pyfastaq.intervals.Interval(50, 51)]
        expected = 'rect(2, 0.5, 44, 1.5, col="blue", border=NA)\n' + \
                   'rect(52, 0.5, 53, 1.5, col="blue", border=NA)\n'
        self.assertEqual(expected, self.qc._cov_to_R_string(intervals, 'blue', 2, 1, 1))


    def test_calculate_should_have_assembled(self):
        '''test _calculate_should_have_assembled'''
        self.qc.ref_ids = ['ref1', 'ref2', 'ref3']
        self.qc.ref_lengths = {x:10 for x in self.qc.ref_ids}
        self.qc.ref_pos_covered_by_contigs = {'ref1': [pyfastaq.intervals.Interval(3, 7)]}
        self.qc.ok_cov_ref_regions = {
            'ref1': [pyfastaq.intervals.Interval(0, 7)],
            'ref2': [pyfastaq.intervals.Interval(0, 9)],
            'ref3': []
        }
        self.qc._calculate_should_have_assembled()
        expected = {
            'ref1': [pyfastaq.intervals.Interval(0, 2)],
            'ref2': [pyfastaq.intervals.Interval(0, 9)],
            'ref3': []
        }
        self.assertEqual(expected, self.qc.should_have_assembled)


    def test_contigs_and_bases_that_hit_ref(self):
        '''test _contigs_and_bases_that_hit_ref'''
        self.qc.assembly_vs_ref_mummer_hits = {
            'ctg1': [
                mummer.NucmerHit('\t'.join(['1', '100', '1', '100', '100', '100', '100.00', '1008', '762', '1', '1', 'ref1', 'ctg1'])),
                mummer.NucmerHit('\t'.join(['1', '100', '51', '150', '100', '100', '100.00', '1008', '762', '1', '1', 'ref1', 'ctg1'])),
            ],
            'ctg2': [mummer.NucmerHit('\t'.join(['1', '42', '42', '84', '42', '84', '100.00', '42', '84', '1', '1', 'ref2', 'ctg2']))]
        }
        self.assertEqual((193, 2), self.qc._contigs_and_bases_that_hit_ref())


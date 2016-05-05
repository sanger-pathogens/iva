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
import stat
import inspect
import tempfile
import copy
import pyfastaq
import shutil
import multiprocessing
import collections
import iva

class Error (Exception): pass

class Qc:
    def __init__(self,
        assembly_fasta,
        output_prefix,
        embl_dir=None,
        ref_db=None,
        reads_fr=None,
        reads_fwd=None,
        reads_rev=None,
        ratt_config=None,
        min_ref_cov=5,
        contig_layout_plot_title="IVA QC contig layout and read depth",
        threads=1,
        nucmer_min_cds_hit_length=20,
        nucmer_min_cds_hit_id=80,
        nucmer_min_ctg_hit_length=100,
        nucmer_min_ctg_hit_id=80,
        gage_nucmer_minid=80,
        smalt_k=15,
        smalt_s=3,
        smalt_id=0.5,
        blast_for_act=False,
        kraken_preload=False,
        clean=True,
    ):

        if embl_dir is None and ref_db is None:
           raise Error('Must provide either embl_dir or ref_db to Qc object. Cannot continue')

        self.embl_dir = embl_dir
        if self.embl_dir is not None:
            self.embl_dir = os.path.abspath(self.embl_dir)

        self.embl_files = []

        self.ref_db = ref_db

        files_to_check = [assembly_fasta, reads_fr, reads_fwd, reads_rev]
        for filename in files_to_check:
            if filename is not None and not os.path.exists(filename):
                raise Error('Error in IVA QC. File not found: "' + filename + '"')

        self.outprefix = output_prefix
        self.assembly_fasta = output_prefix + '.assembly.fasta'
        self.assembly_bam = output_prefix + '.reads_mapped_to_assembly.bam'
        self.ref_bam = output_prefix + '.reads_mapped_to_ref.bam'
        self.reads_fwd = reads_fwd
        self.reads_rev = reads_rev
        self.threads = threads
        self.kraken_preload = kraken_preload
        self.kraken_prefix = self.outprefix + '.kraken'
        self.ratt_config = None if ratt_config is None else os.path.abspath(ratt_config)
        self.ratt_outdir = self.outprefix + '.ratt'
        self.blast_for_act = blast_for_act
        self.blast_out = self.outprefix + '.assembly_v_ref.blastn'
        self.fasta_assembly_contigs_hit_ref = self.outprefix + '.assembly_contigs_hit_ref.fasta'
        self.fasta_assembly_contigs_not_hit_ref = self.outprefix + '.assembly_contigs_not_hit_ref.fasta'
        self.act_script = self.outprefix + '.assembly_v_ref.act.sh'
        self.gage_outdir = self.outprefix + '.gage'
        self.gage_nucmer_minid = gage_nucmer_minid
        self.files_to_clean = []
        self.clean = clean

        if reads_fr:
            self.reads_fwd = self.outprefix + '.reads_1'
            self.reads_rev = self.outprefix + '.reads_2'
            self.files_to_clean.append(self.reads_fwd)
            self.files_to_clean.append(self.reads_rev)
            pyfastaq.tasks.deinterleave(reads_fr, self.reads_fwd, self.reads_rev)

        if not (None not in [self.reads_fwd, self.reads_rev] or reads_fr is not None):
            raise Error('IVA QC needs reads_fr or both reads_fwd and reads_rev')

        self._set_assembly_fasta_data(assembly_fasta)

        def unzip_file(infile, outfile):
            iva.common.syscall('gunzip -c ' + infile + ' > ' + outfile)

        processes = []

        if self.reads_fwd.endswith('.gz'):
            new_reads_fwd = self.outprefix + '.reads_1'
            processes.append(multiprocessing.Process(target=unzip_file, args=(self.reads_fwd, new_reads_fwd)))
            self.reads_fwd = new_reads_fwd
            self.files_to_clean.append(self.reads_fwd)

        if self.reads_rev.endswith('.gz'):
            new_reads_rev = self.outprefix + '.reads_2'
            processes.append(multiprocessing.Process(target=unzip_file, args=(self.reads_rev, new_reads_rev)))
            self.reads_rev = new_reads_rev
            self.files_to_clean.append(self.reads_rev)

        if len(processes) == 1:
            processes[0].start()
            processes[0].join()
        elif len(processes) == 2:
            processes[0].start()
            if self.threads == 1:
                processes[0].join()
            processes[1].start()
            processes[1].join()
            if self.threads > 1:
                processes[0].join()

        self.min_ref_cov = min_ref_cov
        self.threads = threads
        self.contig_layout_plot_title = contig_layout_plot_title
        self.nucmer_min_cds_hit_length = nucmer_min_cds_hit_length
        self.nucmer_min_cds_hit_id = nucmer_min_cds_hit_id
        self.nucmer_min_ctg_hit_length = nucmer_min_ctg_hit_length
        self.nucmer_min_ctg_hit_id = nucmer_min_ctg_hit_id
        self.smalt_k=smalt_k
        self.smalt_s=smalt_s
        self.smalt_id=smalt_id
        self.contig_pos_in_ref = {}
        self.low_cov_ref_regions = {}
        self.low_cov_ref_regions_fwd = {}
        self.low_cov_ref_regions_rev = {}
        self.ok_cov_ref_regions = {}
        self.ref_coverage_fwd = {}
        self.ref_coverage_rev = {}
        self.ref_cds_fasta = output_prefix + '.ref_cds_seqs.fa'
        self.cds_nucmer_coords_in_assembly = output_prefix + '.ref_cds_seqs_mapped_to_assembly.coords'
        self.cds_assembly_stats = {}
        self.refseq_assembly_stats = {}
        self.assembly_vs_ref_coords = output_prefix + '.assembly_vs_ref.coords'
        self.assembly_vs_ref_mummer_hits = {}
        self.ref_pos_covered_by_contigs = {}
        self.ref_pos_not_covered_by_contigs = {}
        self.should_have_assembled = {}
        self.contig_placement = {}
        self.stats_keys = [
            'ref_EMBL_dir',
            'ref_EMBL_files',
            'ref_bases',
            'ref_sequences',
            'ref_bases_assembled',
            'ref_sequences_assembled',
            'ref_sequences_assembled_ok',
            'ref_bases_assembler_missed',
            'assembly_bases',
            'assembly_bases_in_ref',
            'assembly_contigs',
            'assembly_contigs_hit_ref',
            'assembly_bases_reads_disagree',
            'assembly_sum_longest_match_each_segment',
            'cds_number',
            'cds_assembled',
            'cds_assembled_ok',
        ]
        self.stats = {shrubbery: -1 for shrubbery in self.stats_keys}
        self.stats_file_txt = output_prefix + '.stats.txt'
        self.stats_file_tsv = output_prefix + '.stats.tsv'


    def _set_assembly_fasta_data(self, fasta_filename):
        self.assembly_is_empty = os.path.getsize(fasta_filename) == 0
        if self.assembly_is_empty:
            return

        pyfastaq.tasks.to_fasta(fasta_filename, self.assembly_fasta, strip_after_first_whitespace=True)

        self.assembly_fasta_fai = self.assembly_fasta + '.fai'
        if not os.path.exists(self.assembly_fasta_fai):
            iva.common.syscall('samtools faidx ' + self.assembly_fasta)

        self.assembly_lengths = {}
        pyfastaq.tasks.lengths_from_fai(self.assembly_fasta_fai, self.assembly_lengths)


    def _set_ref_seq_data(self):
        assert self.embl_dir is not None
        self.ref_fasta = self.outprefix + '.reference.fa'
        self.ref_fasta_fai = self.ref_fasta + '.fai'
        self.ref_gff = self.outprefix + '.reference.gff'
        tmpdir = tempfile.mkdtemp(prefix='tmp.set_ref_seq_data.', dir=os.getcwd())
        extractor = iva.egg_extract.Extractor(os.path.abspath(os.path.join(os.path.dirname(iva.__file__), os.pardir)))
        embl2gff_egg = os.path.join('iva', 'ratt', 'embl2gff.pl')
        embl2gff = os.path.join(tmpdir, 'embl2gff.pl')
        extractor.copy_file(embl2gff_egg, embl2gff)
        os.chmod(embl2gff, stat.S_IRWXU)

        for embl_file in os.listdir(self.embl_dir):
            fa = os.path.join(tmpdir, embl_file + '.fa')
            gff = os.path.join(tmpdir, embl_file + '.gff')
            embl_full = os.path.join(self.embl_dir, embl_file)
            pyfastaq.tasks.to_fasta(embl_full, fa)
            iva.common.syscall(' '.join([embl2gff, embl_full, '>', gff]))

        iva.common.syscall('cat ' + tmpdir + '/*.gff > ' + self.ref_gff)
        iva.common.syscall('cat ' + tmpdir + '/*.fa > ' + self.ref_fasta)
        shutil.rmtree(tmpdir)
        self._set_ref_fa_data()


    def _set_ref_fa_data(self):
        self.ref_fasta_fai = self.ref_fasta + '.fai'
        iva.common.syscall('samtools faidx ' + self.ref_fasta)
        self.ref_ids = self._ids_in_order_from_fai(self.ref_fasta_fai)
        self.ref_lengths = {}
        pyfastaq.tasks.lengths_from_fai(self.ref_fasta_fai, self.ref_lengths)
        number_of_sequences = pyfastaq.tasks.count_sequences(self.ref_fasta)
        if number_of_sequences != len(self.ref_lengths):
            raise Error('At least one reference name has been used more than once.\nNames must be unique. Cannot continue')

        self.ref_length_offsets = {}
        offset = 0
        for seq in self.ref_ids:
            self.ref_length_offsets[seq] = offset
            offset += self.ref_lengths[seq]


    def _ids_in_order_from_fai(self, filename):
        ids = []
        f = pyfastaq.utils.open_file_read(filename)
        for line in f:
            ids.append(line.rstrip().split('\t')[0])
        pyfastaq.utils.close(f)
        return ids


    def _get_ref_cds_from_gff(self):
        f = pyfastaq.utils.open_file_read(self.ref_gff)
        coords = {}
        for line in f:
            # no annotation allowed after any fasta sequence. See
            # http://www.sequenceontology.org/gff3.shtml
            if line.rstrip() == '##FASTA':
                break
            elif line.startswith('#'):
                continue

            data = line.rstrip().split('\t')
            if data[2] == 'CDS':
                seqname = data[0]
                start = int(data[3]) - 1
                end = int(data[4]) - 1
                strand = data[6]
                if seqname not in coords:
                    coords[seqname] = []
                coords[seqname].append((pyfastaq.intervals.Interval(start, end), strand))

        pyfastaq.utils.close(f)
        for seqname in coords:
            coords[seqname].sort()

        return coords


    def _write_cds_seqs(self, cds_list, fa, f_out):
        for coords, strand in cds_list:
            seqname = fa.id + ':' + str(coords.start + 1) + '-' + str(coords.end + 1) + ':' + strand
            seq = pyfastaq.sequences.Fasta(seqname, fa.seq[coords.start:coords.end+1])
            if strand == '-':
                seq.revcomp()
            print(seq, file=f_out)
            assert seqname not in self.cds_assembly_stats
            self.cds_assembly_stats[seqname] = {
                'ref_name': fa.id,
                'ref_coords': coords,
                'strand': strand,
                'length_in_ref': len(coords),
                'assembled': False,
                'assembled_ok': False,
            }


    def _gff_and_fasta_to_cds(self):
        cds_coords = self._get_ref_cds_from_gff()
        f = pyfastaq.utils.open_file_write(self.ref_cds_fasta)
        seq_reader = pyfastaq.sequences.file_reader(self.ref_fasta)
        for seq in seq_reader:
            if seq.id in cds_coords:
                self._write_cds_seqs(cds_coords[seq.id], seq, f)
        pyfastaq.utils.close(f)


    def _map_cds_to_assembly(self):
        if not os.path.exists(self.ref_cds_fasta):
            self._gff_and_fasta_to_cds()
        if not self.assembly_is_empty:
            iva.mummer.run_nucmer(self.ref_cds_fasta, self.assembly_fasta, self.cds_nucmer_coords_in_assembly, min_length=self.nucmer_min_cds_hit_length, min_id=self.nucmer_min_cds_hit_id)


    def _mummer_coords_file_to_dict(self, filename):
        hits = {}
        for hit in iva.mummer.file_reader(filename):
            if hit.qry_name not in hits:
                hits[hit.qry_name] = []
            hits[hit.qry_name].append(copy.copy(hit))
        return hits


    def _has_orf(self, fa, start, end, min_length):
        subseq = pyfastaq.sequences.Fasta('seq', fa[start:end+1])
        orfs = subseq.all_orfs(min_length=min_length)
        return len(orfs) > 0


    def _calculate_cds_assembly_stats(self):
        if self.assembly_is_empty:
            return
        self._map_cds_to_assembly()
        hits = self._mummer_coords_file_to_dict(self.cds_nucmer_coords_in_assembly)
        contigs = {}
        pyfastaq.tasks.file_to_dict(self.assembly_fasta, contigs)
        for cds_name, hit_list in hits.items():
            self.cds_assembly_stats[cds_name]['number_of_contig_hits'] = len(hit_list)
            hit_coords = [x.qry_coords() for x in hit_list]
            pyfastaq.intervals.merge_overlapping_in_list(hit_coords)
            bases_assembled = pyfastaq.intervals.length_sum_from_list(hit_coords)
            self.cds_assembly_stats[cds_name]['bases_assembled'] = bases_assembled
            self.cds_assembly_stats[cds_name]['assembled'] = 0.9 <= bases_assembled / self.cds_assembly_stats[cds_name]['length_in_ref'] <= 1.1

            if len(hit_list) == 1:
                hit = hit_list[0]
                contig_coords = hit.ref_coords()
                has_orf = self._has_orf(contigs[hit.ref_name], contig_coords.start, contig_coords.end, 0.9 * self.cds_assembly_stats[cds_name]['length_in_ref'])
                self.cds_assembly_stats[cds_name]['assembled_ok'] = has_orf
            else:
                self.cds_assembly_stats[cds_name]['assembled_ok'] = False


    def _get_contig_hits_to_reference(self):
        iva.mummer.run_nucmer(self.assembly_fasta, self.ref_fasta, self.assembly_vs_ref_coords, min_id=self.nucmer_min_ctg_hit_id, min_length=self.nucmer_min_ctg_hit_length, breaklen=500)
        self.assembly_vs_ref_mummer_hits = self._mummer_coords_file_to_dict(self.assembly_vs_ref_coords)


    def _write_fasta_contigs_hit_ref(self):
        contigs = {}
        pyfastaq.tasks.file_to_dict(self.assembly_fasta, contigs)
        f = pyfastaq.utils.open_file_write(self.fasta_assembly_contigs_hit_ref)
        for qry_name in sorted(self.assembly_vs_ref_mummer_hits):
            print(contigs[qry_name], file=f)
        pyfastaq.utils.close(f)


    def _write_fasta_contigs_not_hit_ref(self):
        seq_reader = pyfastaq.sequences.file_reader(self.assembly_fasta)
        f = pyfastaq.utils.open_file_write(self.fasta_assembly_contigs_not_hit_ref)
        for seq in seq_reader:
            if seq.id not in self.assembly_vs_ref_mummer_hits:
                print(seq, file=f)
        pyfastaq.utils.close(f)


    def _hash_nucmer_hits_by_ref(self, hits):
        d = {}
        for l in hits.values():
            for hit in l:
                if hit.ref_name not in d:
                    d[hit.ref_name] = []
                d[hit.ref_name].append(copy.copy(hit))
        return d


    def _calculate_refseq_assembly_stats(self):
        if self.assembly_is_empty:
            refhits = {}
        else:
            refhits = self._hash_nucmer_hits_by_ref(self.assembly_vs_ref_mummer_hits)

        for name in self.ref_ids:
            assert name not in self.refseq_assembly_stats
            if name in refhits:
                hits = refhits[name]
                coords = [hit.ref_coords() for hit in hits]
                pyfastaq.intervals.merge_overlapping_in_list(coords)
                self.refseq_assembly_stats[name] = {
                    'hits': len(hits),
                    'bases_assembled': pyfastaq.intervals.length_sum_from_list(coords),
                    'assembled': 0.9 <= pyfastaq.intervals.length_sum_from_list(coords) / self.ref_lengths[name],
                    'assembled_ok': len(hits) == 1 and 0.9 <= hits[0].hit_length_ref / self.ref_lengths[name] <= 1.1 and 0.9 <= hits[0].qry_length / hits[0].hit_length_qry,
                    'longest_matching_contig': self._longest_matching_contig(refhits, name),
                }
            else:
                self.refseq_assembly_stats[name] = {
                    'hits': 0,
                    'bases_assembled': 0,
                    'assembled': False,
                    'assembled_ok': False,
                    'longest_matching_contig': 0,
                }



    def _invert_list(self, coords, seq_length):
        if len(coords) == 0:
            return[pyfastaq.intervals.Interval(0, seq_length - 1)]

        not_covered = []

        if coords[0].start != 0:
            not_covered.append(pyfastaq.intervals.Interval(0, coords[0].start - 1))

        for i in range(len(coords) - 1):
            not_covered.append(pyfastaq.intervals.Interval(coords[i].end + 1, coords[i+1].start - 1))

        if coords[-1].end < seq_length - 1:
            not_covered.append(pyfastaq.intervals.Interval(coords[-1].end + 1, seq_length - 1))

        return not_covered


    def _calculate_ref_positions_covered_by_contigs(self):
        if self.assembly_is_empty:
            self.ref_pos_not_covered_by_contigs = {}
            for seq, lngth in self.ref_lengths.items():
                self.ref_pos_not_covered_by_contigs[seq] = [pyfastaq.intervals.Interval(0, lngth - 1)]
            return

        for seq in self.assembly_vs_ref_mummer_hits:
            for hit in self.assembly_vs_ref_mummer_hits[seq]:
                if hit.ref_name not in self.ref_pos_covered_by_contigs:
                    self.ref_pos_covered_by_contigs[hit.ref_name] = []
                self.ref_pos_covered_by_contigs[hit.ref_name].append(hit.ref_coords())

        for coords_list in self.ref_pos_covered_by_contigs.values():
            pyfastaq.intervals.merge_overlapping_in_list(coords_list)

        for seq in self.ref_ids:
            if seq in self.ref_pos_covered_by_contigs:
                l = self.ref_pos_covered_by_contigs[seq]
            else:
                l = []
            self.ref_pos_not_covered_by_contigs[seq] = self._invert_list(l, self.ref_lengths[seq])


    def _get_overlapping_qry_hits(self, hits, hit):
        overlapping = []
        hit_coords = hit.qry_coords()

        for test_hit in hits:
            if test_hit != hit:
                test_coords = test_hit.qry_coords()
                if test_coords.intersects(hit_coords):
                    overlapping.append(test_hit)

        return overlapping


    def _get_unique_and_repetitive_from_contig_hits(self, hits):
        unique = []
        repetitive = []
        if len(hits) == 0:
            return unique, repetitive

        for hit in hits:
            if len(self._get_overlapping_qry_hits(hits, hit)):
                repetitive.append(hit)
            else:
                unique.append(hit)

        return unique, repetitive


    def _longest_matching_contig(self, hits, refname):
        if len(hits) == 0:
            return 0
        hits_to_refname = [x for x in hits[refname] if x.ref_name == refname and x.qry_length <= 2 * x.ref_length]
        if len(hits_to_refname) == 0:
            return 0
        else:
            return max([x.qry_length for x in hits_to_refname])


    def _calculate_incorrect_assembly_bases(self):
        if self.assembly_is_empty:
            self.incorrect_assembly_bases = {}
        else:
            self.incorrect_assembly_bases = iva.mapping.find_incorrect_ref_bases(self.assembly_bam, self.assembly_fasta)


    def _contig_placement_in_reference(self, hits):
        unique_hits, repetitive_hits = self._get_unique_and_repetitive_from_contig_hits(hits)
        placement = [(x.qry_coords(), x.ref_name, x.ref_coords(), x.on_same_strand(), False) for x in unique_hits]
        placement += [(x.qry_coords(), x.ref_name, x.ref_coords(), x.on_same_strand(), True) for x in repetitive_hits]
        placement.sort()
        return placement


    def _calculate_contig_placement(self):
        if self.assembly_is_empty:
            return
        self._get_contig_hits_to_reference()
        self.contig_placement = {qry_name: self._contig_placement_in_reference(self.assembly_vs_ref_mummer_hits[qry_name]) for qry_name in self.assembly_vs_ref_mummer_hits}


    def _get_R_plot_contig_order_from_contig_placement(self):
        contig_positions = []
        for qryname, coords_list in  self.contig_placement.items():
            for qry_coords, refname, ref_coords, same_strand, repetitive in coords_list:
                offset = self.ref_length_offsets[refname]
                ref_coords = pyfastaq.intervals.Interval(ref_coords.start + offset, ref_coords.end + offset)
                contig_positions.append((ref_coords, qry_coords, same_strand, repetitive, qryname))

        contig_positions.sort()
        names = set()
        contig_names = []
        for name in [x[4] for x in contig_positions]:
            if name not in names:
                contig_names.append(name)
                names.add(name)
        return contig_names


    def _map_reads_to_assembly(self):
        if not self.assembly_is_empty:
            iva.mapping.map_reads(self.reads_fwd, self.reads_rev, self.assembly_fasta, self.assembly_bam[:-4], sort=True, threads=self.threads, index_k=self.smalt_k, index_s=self.smalt_s, minid=self.smalt_id, extra_smalt_map_ops='-x')
            os.unlink(self.assembly_bam[:-4] + '.unsorted.bam')


    def _write_ref_info(self, filename):
        assert self.embl_dir is not None
        files = sorted(os.listdir(self.embl_dir))
        f = pyfastaq.utils.open_file_write(filename)
        print('EMBL_directory', self.embl_dir, sep='\t', file=f)
        print('Files', '\t'.join(files), sep='\t', file=f)
        pyfastaq.utils.close(f)


    def _choose_reference_genome(self):
        if self.embl_dir is None:
            assert self.ref_db is not None
            assert os.path.exists(self.assembly_bam)
            tmp_reads = self.outprefix + '.tmp.subsample.reads.fastq'
            iva.mapping.subsample_bam(self.assembly_bam, tmp_reads, coverage=40)
            db = iva.kraken.Database(self.ref_db, threads=self.threads, preload=self.kraken_preload)
            self.embl_dir = db.choose_reference(tmp_reads, self.kraken_prefix)
            os.unlink(tmp_reads)
            if self.embl_dir is None:
                raise Error('Unable to determine reference genome automatically. Cannot continue')
        else:
            self.embl_dir = os.path.abspath(self.embl_dir)

        self.embl_files = sorted(os.listdir(self.embl_dir))


    def _make_act_files(self):
        if self.assembly_is_empty or not self.blast_for_act:
            return

        iva.qc_external.run_blastn_and_write_act_script(self.assembly_fasta, self.ref_fasta, self.blast_out, self.act_script)


    def _map_reads_to_reference(self):
        assert os.path.exists(self.ref_fasta)
        iva.mapping.map_reads(self.reads_fwd, self.reads_rev, self.ref_fasta, self.ref_bam[:-4], sort=True, threads=self.threads, index_k=self.smalt_k, index_s=self.smalt_s, minid=self.smalt_id, extra_smalt_map_ops='-x')
        os.unlink(self.ref_bam[:-4] + '.unsorted.bam')


    def _calculate_ref_read_coverage(self):
        if not os.path.exists(self.ref_bam):
            self._map_reads_to_reference()
        for seq in self.ref_ids:
            assert seq not in self.ref_coverage_fwd
            self.ref_coverage_fwd[seq] = iva.mapping.get_bam_region_coverage(self.ref_bam, seq, self.ref_lengths[seq])
            assert seq not in self.ref_coverage_rev
            self.ref_coverage_rev[seq] = iva.mapping.get_bam_region_coverage(self.ref_bam, seq, self.ref_lengths[seq], rev=True)


    def _coverage_list_to_low_cov_intervals(self, l):
        bad_intervals = []
        start = None
        cov_bad = False

        for i in range(len(l)):
            cov_bad = l[i] < self.min_ref_cov
            if cov_bad:
                if start is None:
                    start = i
            else:
                if start is not None:
                    bad_intervals.append(pyfastaq.intervals.Interval(start, i-1))
                start = None

        if cov_bad and start is not None:
            bad_intervals.append(pyfastaq.intervals.Interval(start, i))
        return bad_intervals


    def _calculate_ref_read_region_coverage(self):
        assert len(self.ref_coverage_fwd)
        assert len(self.ref_coverage_rev)
        for seq in self.ref_ids:
            self.low_cov_ref_regions_fwd[seq] = self._coverage_list_to_low_cov_intervals(self.ref_coverage_fwd[seq])
            self.low_cov_ref_regions_rev[seq] = self._coverage_list_to_low_cov_intervals(self.ref_coverage_rev[seq])
            fwd_ok = self._invert_list(self.low_cov_ref_regions_fwd[seq], self.ref_lengths[seq])
            rev_ok = self._invert_list(self.low_cov_ref_regions_rev[seq], self.ref_lengths[seq])
            self.ok_cov_ref_regions[seq] = pyfastaq.intervals.intersection(fwd_ok, rev_ok)
            self.low_cov_ref_regions[seq] = pyfastaq.intervals.intersection(self.low_cov_ref_regions_fwd[seq], self.low_cov_ref_regions_rev[seq])


    def _write_ref_coverage_to_files_for_R(self, outprefix):
        assert len(self.ref_coverage_fwd)
        assert len(self.ref_coverage_rev)
        def list_to_file(d, fname):
            f = pyfastaq.utils.open_file_write(fname)
            for refname in self.ref_ids:
                for x in d[refname]:
                    print(x, file=f)
            pyfastaq.utils.close(f)
        list_to_file(self.ref_coverage_fwd, outprefix + '.fwd')
        list_to_file(self.ref_coverage_rev, outprefix + '.rev')


    def _cov_to_R_string(self, intervals, colour, x_offset, y_position, contig_height):
        s = ''
        for interval in intervals:
            s += 'rect(' + \
                 str(interval.start + x_offset) + ', ' + \
                 str(y_position - 0.5 * contig_height) + ', ' + \
                 str(interval.end + x_offset) + ', ' + \
                 str(y_position + 0.5 * contig_height) + ', ' + \
                 'col="' + colour + '", ' + \
                 'border=NA)\n'
        return s


    def _calculate_should_have_assembled(self):
        for name in self.ref_ids:
            if name in self.ref_pos_covered_by_contigs:
                l = self.ref_pos_covered_by_contigs[name]
            else:
                l = []
            self.should_have_assembled[name] = pyfastaq.intervals.intersection(self._invert_list(l, self.ref_lengths[name]), self.ok_cov_ref_regions[name])


    def _calculate_gage_stats(self):
        if self.assembly_is_empty:
            self.gage_stats = iva.qc_external.dummy_gage_stats()
            self.gage_stats['Missing Reference Bases'] = sum(self.ref_lengths.values())
        else:
            self.gage_stats = iva.qc_external.run_gage(self.ref_fasta, self.assembly_fasta, self.gage_outdir, nucmer_minid=self.gage_nucmer_minid, clean=self.clean)


    def _calculate_ratt_stats(self):
        if self.assembly_is_empty:
            self.ratt_stats = iva.qc_external.dummy_ratt_stats()
        else:
            self.ratt_stats = iva.qc_external.run_ratt(self.embl_dir, self.assembly_fasta, self.ratt_outdir, config_file=self.ratt_config, clean=self.clean)


    def _do_calculations(self):
        self._map_reads_to_assembly()
        self._choose_reference_genome()
        self._set_ref_seq_data()
        self._make_act_files()
        self._map_reads_to_reference()
        self._calculate_incorrect_assembly_bases()
        self._calculate_contig_placement()
        self._write_fasta_contigs_hit_ref()
        self._write_fasta_contigs_not_hit_ref()
        self._calculate_ref_read_coverage()
        self._calculate_ref_read_region_coverage()
        self._calculate_ref_positions_covered_by_contigs()
        self._calculate_should_have_assembled()
        self._calculate_cds_assembly_stats()
        self._calculate_refseq_assembly_stats()
        self._calculate_gage_stats()
        self._calculate_ratt_stats()
        self._calculate_stats()


    def _contigs_and_bases_that_hit_ref(self):
        total_bases = 0
        for name in self.assembly_vs_ref_mummer_hits:
            coords = [x.qry_coords() for x in self.assembly_vs_ref_mummer_hits[name]]
            pyfastaq.intervals.merge_overlapping_in_list(coords)
            total_bases += pyfastaq.intervals.length_sum_from_list(coords)
        return total_bases, len(self.assembly_vs_ref_mummer_hits)


    def _calculate_stats(self):
        self.stats['ref_EMBL_dir'] = self.embl_dir
        self.stats['ref_EMBL_files'] = ' '.join(self.embl_files)
        self.stats['ref_bases'] = sum(self.ref_lengths.values())
        self.stats['ref_sequences'] = len(self.ref_lengths)
        self.stats['ref_bases_assembled'] = sum([pyfastaq.intervals.length_sum_from_list(l) for l in list(self.ref_pos_covered_by_contigs.values())])
        self.stats['ref_sequences_assembled'] = len([1 for x in self.refseq_assembly_stats.values() if x['assembled']])
        self.stats['ref_sequences_assembled_ok'] = len([1 for x in self.refseq_assembly_stats.values() if x['assembled_ok']])
        self.stats['ref_bases_assembler_missed'] = sum([pyfastaq.intervals.length_sum_from_list(l) for l in list(self.should_have_assembled.values())])
        self.stats['assembly_bases'] = sum(self.assembly_lengths.values())
        self.stats['assembly_contigs'] = len(self.assembly_lengths)
        self.stats['assembly_bases_in_ref'], self.stats['assembly_contigs_hit_ref'] = self._contigs_and_bases_that_hit_ref()
        self.stats['assembly_bases_reads_disagree'] = sum([len(x) for x in self.incorrect_assembly_bases.values()])
        self.stats['assembly_sum_longest_match_each_segment'] = sum([x['longest_matching_contig'] for x in self.refseq_assembly_stats.values()])
        self.stats['cds_number'] = len(self.cds_assembly_stats)
        self.stats['cds_assembled'] = len([1 for x in self.cds_assembly_stats.values() if x['assembled']])
        self.stats['cds_assembled_ok'] = len([1 for x in self.cds_assembly_stats.values() if x['assembled_ok']])


    def _write_stats_txt(self):
        f = pyfastaq.utils.open_file_write(self.stats_file_txt)
        for stat in self.stats_keys:
            print(stat, self.stats[stat], sep='\t', file=f)
        for stat in iva.qc_external.gage_stats:
            print('gage_' + stat.replace(' ', '_'), self.gage_stats[stat], sep='\t', file=f)
        for stat in iva.qc_external.ratt_stats:
            print('ratt_' + stat.replace(' ', '_'), self.ratt_stats[stat], sep='\t', file=f)

        pyfastaq.utils.close(f)


    def _write_stats_tsv(self):
        f = pyfastaq.utils.open_file_write(self.stats_file_tsv)
        print('\t'.join([x.replace(' ', '_') for x in self.stats_keys]), end='\t', file=f)
        print('\t'.join(['gage_' + x.replace(' ', '_') for x in iva.qc_external.gage_stats]), end='\t', file=f)
        print('\t'.join(['ratt_' + x.replace(' ', '_') for x in iva.qc_external.ratt_stats]), file=f)
        print('\t'.join([str(self.stats[x]) for x in self.stats_keys]),
              '\t'.join([str(self.gage_stats[x]) for x in iva.qc_external.gage_stats]),
              '\t'.join([str(self.ratt_stats[x]) for x in iva.qc_external.ratt_stats]),
              sep='\t', file=f)
        pyfastaq.utils.close(f)


    def _write_stats_files(self):
        self._write_stats_txt()
        self._write_stats_tsv()


    def _make_R_plots(self):
        outprefix = self.outprefix + '.contig_placement'
        contig_names = self._get_R_plot_contig_order_from_contig_placement()
        number_of_contigs = len(contig_names)
        ref_length = sum(self.ref_lengths.values())
        r_script = outprefix + '.R'
        f = pyfastaq.utils.open_file_write(r_script)
        contig_height = 0.8
        vertical_lines = ''
        if len(self.ref_ids) > 0:
            for name in self.ref_ids:
                x_position = self.ref_length_offsets[name]
                if x_position > 0:
                    vertical_lines += '\n' + 'abline(v=' + str(x_position) + ', col="gray")'


        print('pdf(file="', outprefix, '.pdf")', sep='', file=f)
        print('layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), heights=c(3,1.2,1.2))', file=f)
        print('par(mar=c(1, 4, 2, 0.5))', file=f)

        # ---------- contig layout plot ------------------------
        print('plot(-100, type="n", xlim=c(0,', ref_length, '), ylim=c(0, ', number_of_contigs, '), yaxt="n", ylab="", xlab="")', sep='', file=f)
        print('title("', self.contig_layout_plot_title, '", ylab="Contigs")', sep='', file=f)
        print(vertical_lines, file=f)

        if number_of_contigs > 0:
            print('contig_names=c("', '", "'.join(contig_names), '")', sep='', file=f)
            print('axis(2, at=c(1:', number_of_contigs, '), labels=contig_names, las=2, cex.axis=0.3)', sep='', file=f)

            for i in range(len(contig_names)):
                contig_name = contig_names[i]
                y_centre = i + 1
                contig_positions = self.contig_placement[contig_name]
                for contig_coords, ref_name, ref_coords, same_strand, repetitive in contig_positions:
                    offset = self.ref_length_offsets[ref_name]
                    if repetitive:
                        colour = "red"
                    else:
                        colour = "blue"

                    if same_strand:
                         colour = "dark" + colour

                    print('rect(', ref_coords.start + offset, ',',
                            y_centre - 0.5 * contig_height, ',',
                            ref_coords.end + offset, ',',
                            y_centre + 0.5 * contig_height, ',',
                            'col="', colour, '")', sep='', file=f)

        # ----------- read coverage heatmap ---------------------
        print('par(mar=c(0, 4, 1, 0.5))', file=f)
        print('plot(-100, type="n", xlim=c(0,', ref_length, '), ylim=c(0, ', 3, '), xaxt="n", yaxt="n", ylab="Contig/Read coverage OK", xlab="", frame.plot=F)', sep='', file=f)
        print('axis(2, at=c(1,2), labels=c("Reads", "Contigs"), las=2, cex.axis=0.6)', sep='', file=f)

        for name in self.ref_ids:
            offset = self.ref_length_offsets[name]
            print(self._cov_to_R_string(self.ok_cov_ref_regions[name], 'black', offset, 1.3, 0.25), file=f)
            print(self._cov_to_R_string(self.low_cov_ref_regions_fwd[name], 'red', offset, 1, 0.25), file=f)
            print(self._cov_to_R_string(self.low_cov_ref_regions_rev[name], 'red', offset, 0.7, 0.25), file=f)

            if name in self.ref_pos_covered_by_contigs:
                print(self._cov_to_R_string(self.ref_pos_covered_by_contigs[name], 'black', offset, 2.3, 0.25), file=f)

            if name in self.should_have_assembled:
                print(self._cov_to_R_string(self.should_have_assembled[name], 'red', offset, 1.7, 0.25), file=f)

            print(self._cov_to_R_string(self.ref_pos_not_covered_by_contigs[name], 'black', offset, 2, 0.25), file=f)
        print(vertical_lines, file=f)

        # ----------- read depth on reference plot --------------
        self._write_ref_coverage_to_files_for_R(self.outprefix + '.read_coverage_on_ref')
        print('fwd_ref_cov = scan("', self.outprefix + '.read_coverage_on_ref.fwd', '")', sep='', file=f)
        print('rev_ref_cov = scan("', self.outprefix + '.read_coverage_on_ref.rev', '")', sep='', file=f)
        print('par(mar=c(5, 4, 0, 0.5))', file=f)
        print('plot(fwd_ref_cov, type="l", xlim=c(0, length(fwd_ref_cov) + 1), ylim=c(-max(rev_ref_cov), max(fwd_ref_cov)), col="blue", frame.plot=F, ylab="Read depth", xlab="Position in reference")', file=f)
        print('lines(-rev_ref_cov, col="blue")', file=f)
        print('abline(h=0, lty=2)', file=f)
        print(vertical_lines, file=f)

        print('dev.off()', file=f)
        pyfastaq.utils.close(f)
        iva.common.syscall('R CMD BATCH ' + r_script)


    def _clean(self):
        for fname in self.files_to_clean:
            os.unlink(fname)


    def run(self):
        self._do_calculations()
        self._make_R_plots()
        self._write_stats_files()
        self._clean()

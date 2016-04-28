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
import pysam
import tempfile
import shutil
from iva import contig, mapping, seed, mummer, graph, edge, common
import pyfastaq

class Assembly:
    def __init__(self, contigs_file=None, map_index_k=15, map_index_s=3, threads=1, kmc_threads=1, max_insert=800, map_minid=0.5, min_clip=3, ext_min_cov=5, ext_min_ratio=2, ext_bases=100, verbose=0, seed_min_cov=5, seed_min_ratio=10, seed_min_kmer_count=200, seed_max_kmer_count=1000000000, seed_start_length=None, seed_stop_length=500, seed_overlap_length=None, make_new_seeds=False, contig_iter_trim=10, seed_ext_max_bases=50, max_contigs=50, clean=True, strand_bias=0):
        self.contigs = {}
        self.contig_lengths = {}
        self.map_index_k = map_index_k
        self.map_index_s = map_index_s
        self.threads = threads
        self.kmc_threads = kmc_threads
        self.max_insert = max_insert
        self.map_minid = map_minid
        self.min_clip = min_clip
        self.ext_min_cov = ext_min_cov
        self.ext_min_ratio = ext_min_ratio
        self.ext_bases = ext_bases
        self.verbose = verbose
        self.clean = clean
        self.make_new_seeds = make_new_seeds
        self.seed_start_length = seed_start_length
        self.seed_stop_length = seed_stop_length
        self.seed_min_kmer_count = seed_min_kmer_count
        self.seed_max_kmer_count = seed_max_kmer_count
        self.seed_ext_max_bases = seed_ext_max_bases
        self.seed_overlap_length = seed_overlap_length
        self.seed_min_cov = seed_min_cov
        self.seed_min_ratio = seed_min_ratio
        self.contig_iter_trim = contig_iter_trim
        self.max_contigs = max_contigs
        self.strand_bias = strand_bias
        self.contigs_trimmed_for_strand_bias = set()
        self.used_seeds = set()

        if contigs_file is None:
            self.make_new_seeds = True
        else:
            contigs = {}
            pyfastaq.tasks.file_to_dict(contigs_file, contigs)
            for ctg in contigs:
                self._add_contig(contigs[ctg])


    def _add_contig(self, ctg, min_length=1):
        if len(ctg) < min_length:
            return
        assert ctg.id not in self.contigs
        assert len(ctg) > 0
        self.contigs[ctg.id] = contig.Contig(ctg, verbose=self.verbose)
        self.contig_lengths[ctg.id] = [[len(self.contigs[ctg.id]), 0, 0]]


    def write_contigs_to_file(self, filename, min_length=None, do_not_write=None, only_write=None, biggest_first=False, order_by_orfs=False, prefix=None):
        if do_not_write is None:
            do_not_write = set()
        if only_write is None:
            only_write = set()

        printed = 0
        if min_length is None:
            min_length = self.map_index_k + 1
        f = pyfastaq.utils.open_file_write(filename)
        if biggest_first:
            contig_names = self._contig_names_size_order(biggest_first=True)
        elif order_by_orfs:
            names = self._get_contig_order_by_orfs()
            contig_names = [x[0] for x in names]
            contig_revcomp = [x[1] for x in names]
        else:
            contig_names = sorted(list(self.contigs.keys()))

        for i in range(len(contig_names)):
            name = contig_names[i]

            if len(self.contigs[name]) >= min_length and name not in do_not_write and (name in only_write or len(only_write)==0):
                if order_by_orfs and contig_revcomp[i]:
                    self.contigs[name].fa.revcomp()

                if prefix is None:
                    print(self.contigs[name].fa, file=f)
                else:
                    printed += 1
                    self.contigs[name].fa.id = prefix + '.' + str(printed).zfill(5)
                    print(self.contigs[name].fa, file=f)
                    self.contigs[name].fa.id = name

                if order_by_orfs and contig_revcomp[i]:
                    self.contigs[name].fa.revcomp()

        pyfastaq.utils.close(f)


    def _get_contig_order_by_orfs(self, min_length=300):
        longest_orfs = []
        no_orfs = set()
        ordered_names = []
        for contig in self.contigs.values():
            orfs = contig.fa.all_orfs(min_length)
            reverse = False
            max_length = 0

            for coords, rev in orfs:
                if len(coords) > max_length:
                    max_length = len(coords)
                    reverse = rev

            if max_length > 0:
                longest_orfs.append((contig.fa.id, max_length, reverse))
            else:
                no_orfs.add((contig.fa.id, len(contig.fa)))

        all_in_size_order = self._contig_names_size_order(biggest_first=True)
        ordered_names = sorted(longest_orfs, key=lambda x:x[1], reverse=True)
        ordered_names = [(x[0], x[2]) for x in ordered_names]

        for t in sorted(no_orfs, key=lambda x:x[1], reverse=True):
            ordered_names.append((t[0], False))

        return ordered_names


    def _map_reads(self, fwd_reads, rev_reads, out_prefix, required_flag=None, exclude_flag=None, sort_reads=False, mate_ref=None, no_map_contigs=None):
        if no_map_contigs is None:
            no_map_contigs = set()
        if self.verbose:
            print('    map reads', fwd_reads, rev_reads, sep='\t')
        reference = out_prefix + '.ref.fa'
        self.write_contigs_to_file(reference, do_not_write=no_map_contigs)
        mapping.map_reads(fwd_reads, rev_reads, reference, out_prefix, index_k=self.map_index_k, index_s=self.map_index_s, threads=self.threads, max_insert=self.max_insert, minid=self.map_minid, verbose=self.verbose, required_flag=required_flag, sort=sort_reads, exclude_flag=exclude_flag)
        if self.clean:
            os.unlink(reference)
        os.unlink(reference + '.fai')


    def _extend_contigs_with_bam(self, bam_in, out_prefix=None, output_all_useful_reads=False):
        if out_prefix is not None:
            fa_out1 = pyfastaq.utils.open_file_write(out_prefix + '_1.fa')
            fa_out2 = pyfastaq.utils.open_file_write(out_prefix + '_2.fa')
        keep_read_types = set([mapping.CAN_EXTEND_LEFT, mapping.CAN_EXTEND_RIGHT, mapping.KEEP])
        if output_all_useful_reads:
            keep_read_types.add(mapping.BOTH_UNMAPPED)
        previous_sam = None
        left_seqs = []
        right_seqs = []
        sam_reader = pysam.Samfile(bam_in, "rb")

        for current_sam in sam_reader.fetch(until_eof=True):
            if previous_sam is None:
                previous_sam = current_sam
                continue

            previous_type, current_type = mapping.get_pair_type(previous_sam, current_sam, self._get_ref_length_sam_pair(sam_reader, previous_sam, current_sam), self.max_insert, min_clip=self.min_clip)

            for sam, sam_type in [(previous_sam, previous_type), (current_sam, current_type)]:
                if sam_type == mapping.CAN_EXTEND_LEFT:
                    name = mapping.get_ref_name(sam, sam_reader)
                    clipped = mapping.soft_clipped(sam)[0]
                    self.contigs[name].add_left_kmer(common.decode(sam.seq[:clipped]))
                elif sam_type == mapping.CAN_EXTEND_RIGHT:
                    name = mapping.get_ref_name(sam, sam_reader)
                    self.contigs[name].add_right_kmer(common.decode(sam.seq[sam.qend:]))

                if out_prefix is not None and sam_type in keep_read_types:
                    if sam.is_read1:
                        print(mapping.sam_to_fasta(sam), file=fa_out1)
                    else:
                        print(mapping.sam_to_fasta(sam), file=fa_out2)

            previous_sam = None

        if out_prefix is not None:
            pyfastaq.utils.close(fa_out1)
            pyfastaq.utils.close(fa_out2)
        total_bases_added = 0

        for ctg in self.contigs:
            left_length, right_length = self.contigs[ctg].extend(self.ext_min_cov, self.ext_min_ratio, self.ext_bases)
            if self.verbose:
                print('    extend contig ' +  ctg, 'new_length:' + str(len(self.contigs[ctg])), 'added_left:' + str(left_length), 'added_right:' + str(right_length), sep='\t')
            self.contig_lengths[ctg].append([len(self.contigs[ctg]), left_length, right_length])
            total_bases_added += left_length + right_length

        return total_bases_added


    def _trim_contig_for_strand_bias(self, bam, ctg_name):
        assert os.path.exists(bam)
        if ctg_name in self.contigs_trimmed_for_strand_bias:
            return
        ctg_length = len(self.contigs[ctg_name])
        fwd_cov = mapping.get_bam_region_coverage(bam, ctg_name, ctg_length)
        rev_cov = mapping.get_bam_region_coverage(bam, ctg_name, ctg_length, rev=True)
        first_good_base = 0
        while first_good_base < ctg_length:
            total_cov = fwd_cov[first_good_base] + rev_cov[first_good_base]
            if total_cov >= self.ext_min_cov and min(fwd_cov[first_good_base], rev_cov[first_good_base]) / total_cov >= self.strand_bias:
                break
            first_good_base += 1

        last_good_base = ctg_length - 1
        while last_good_base > first_good_base:
            total_cov = fwd_cov[last_good_base] + rev_cov[last_good_base]
            if total_cov >= self.ext_min_cov and min(fwd_cov[last_good_base], rev_cov[last_good_base]) / total_cov >= self.strand_bias:
                break
            last_good_base -= 1

        if self.verbose >= 2:
            print('Trimming strand biased ends of contig', ctg_name, '- good base range is', first_good_base + 1, 'to', last_good_base + 1, 'from', ctg_length, 'bases')
        self.contigs[ctg_name].fa.seq = self.contigs[ctg_name].fa.seq[first_good_base:last_good_base+1]


    def _good_intervals_from_strand_coverage(self, fwd_cov, rev_cov):
        assert len(fwd_cov) == len(rev_cov)
        good_intervals = []
        start = None
        cov_ok = False
        for i in range(len(fwd_cov)):
            total_cov = fwd_cov[i] + rev_cov[i]
            cov_ok = total_cov >= self.ext_min_cov and min(fwd_cov[i], rev_cov[i]) / total_cov >= self.strand_bias
            if cov_ok:
                if start is None:
                    start = i
            else:
                if start is not None:
                    good_intervals.append((start, i-1))
                start = None

        if cov_ok and start is not None:
            good_intervals.append((start, i-1))
        return good_intervals


    def _subcontigs_from_strand_bias(self, bam, ctg_name):
        ctg_length = len(self.contigs[ctg_name])
        fwd_cov = mapping.get_bam_region_coverage(bam, ctg_name, ctg_length)
        rev_cov = mapping.get_bam_region_coverage(bam, ctg_name, ctg_length, rev=True)
        good_intervals = self._good_intervals_from_strand_coverage(fwd_cov, rev_cov)
        new_contigs = []

        if len(good_intervals) == 1:
            self.contigs[ctg_name].fa.seq = self.contigs[ctg_name].fa.seq[good_intervals[0][0]:good_intervals[0][1]+1]
        elif len(good_intervals) > 1:
            for i in range(len(good_intervals)):
                start = good_intervals[i][0]
                end = good_intervals[i][1]
                if end - start + 1 >= 100:
                    new_contigs.append(pyfastaq.sequences.Fasta(ctg_name + '.' + str(i+1), self.contigs[ctg_name].fa[start:end+1]))

        return new_contigs


    def _trim_strand_biased_ends(self, reads_prefix, out_prefix=None, tag_as_trimmed=False, break_contigs=False):
        tmpdir = tempfile.mkdtemp(prefix='tmp.trim_strand_biased_ends.', dir=os.getcwd())
        tmp_prefix = os.path.join(tmpdir, 'out')
        sorted_bam = tmp_prefix + '.bam'
        unsorted_bam = tmp_prefix + '.unsorted.bam'
        original_map_minid = self.map_minid
        self.map_minid = 0.9
        self._map_reads(reads_prefix + '_1.fa', reads_prefix + '_2.fa', tmp_prefix, sort_reads=True)
        assert os.path.exists(sorted_bam)
        self.map_minid = original_map_minid
        new_contigs = []
        contigs_to_remove = set()
        for ctg in self.contigs:
            if break_contigs:
                subcontigs = self._subcontigs_from_strand_bias(sorted_bam, ctg)
                if len(subcontigs):
                    new_contigs.extend(subcontigs)
                    contigs_to_remove.add(ctg)
            elif ctg not in self.contigs_trimmed_for_strand_bias:
                self._trim_contig_for_strand_bias(sorted_bam, ctg)
                # contig could get completely trimmed so nothing left, in which
                # case, we need to remove it
                if len(self.contigs[ctg]) == 0:
                    contigs_to_remove.add(ctg)
                elif tag_as_trimmed:
                    self.contigs_trimmed_for_strand_bias.add(ctg)

        for ctg in contigs_to_remove:
            self._remove_contig(ctg)

        for ctg in new_contigs:
            self._add_contig(ctg, min_length=0.75 * self.self.seed_stop_length)


        if out_prefix is not None:
            mapping.bam_file_to_fasta_pair_files(unsorted_bam, out_prefix + '_1.fa', out_prefix + '_2.fa', remove_proper_pairs=True)
        shutil.rmtree(tmpdir)


    def trim_contigs(self, trim):
        for ctg in self.contigs:
            if self._contig_worth_extending(ctg):
                self.contigs[ctg].fa.trim(trim, trim)
                self.contig_lengths[ctg][-1][0] -= 2 * trim


    def _contig_worth_extending(self, name):
        if name in self.contigs_trimmed_for_strand_bias:
            return False
        return len(self.contig_lengths[name]) < 3 \
                  or self.contig_lengths[name][-1][0] > max([self.contig_lengths[name][x][0] for x in range(len(self.contig_lengths[name])-2)])


    def _worth_extending(self):
        for ctg in self.contigs:
            if self._contig_worth_extending(ctg):
                return True
        return False


    def _extend_with_reads(self, reads_prefix, out_prefix, no_map_contigs):
        tmpdir = tempfile.mkdtemp(prefix='tmp.extend.', dir=os.getcwd())
        tmp_prefix = os.path.join(tmpdir, 'reads')
        total_bases_added = 0

        for i in range(5):
            bam_prefix = out_prefix + '.' + str(i+1) + '.map'
            bam = bam_prefix +  '.bam'
            self._map_reads(reads_prefix + '_1.fa', reads_prefix + '_2.fa', bam_prefix, no_map_contigs=no_map_contigs)
            reads_prefix = tmp_prefix + '.' + str(i)
            bases_added = self._extend_contigs_with_bam(bam, out_prefix=reads_prefix)
            total_bases_added += bases_added
            if self.clean:
                os.unlink(bam)
            if bases_added < 0.2 * self.ext_bases:
                break

        shutil.rmtree(tmpdir)
        return total_bases_added


    def _read_pair_extension_iterations(self, reads_prefix, out_prefix, no_map_contigs=None):
        if no_map_contigs is None:
            no_map_contigs = set()
        assert(len(self.contigs) > len(no_map_contigs))
        if self.verbose:
            print('{:-^79}'.format(' ' + out_prefix + ' start extension subiteration 0001 '), flush=True)

        bases_added = self._extend_with_reads(reads_prefix, out_prefix + '.1', no_map_contigs)
        current_reads_prefix = reads_prefix

        if bases_added == 0:
            return True
        try_contig_trim = False
        i = 1

        while self._worth_extending() or try_contig_trim:
            i += 1
            if self.verbose:
                print('{:-^79}'.format(' ' + out_prefix + ' start extension subiteration ' + str(i).zfill(4) + ' '), flush=True)

            if i % 5 == 0:
                tmpdir = tempfile.mkdtemp(prefix='tmp.filter_reads.', dir=os.getcwd())
                tmp_prefix = os.path.join(tmpdir, 'out')
                bam = tmp_prefix + '.bam'
                original_map_minid = self.map_minid
                self.map_minid = 0.9
                self._map_reads(current_reads_prefix + '_1.fa', current_reads_prefix + '_2.fa', tmp_prefix)
                self.map_minid = original_map_minid
                filter_prefix = reads_prefix + '.subiter.' + str(i) + '.reads'
                mapping.bam_file_to_fasta_pair_files(bam, filter_prefix + '_1.fa', filter_prefix + '_2.fa', remove_proper_pairs=True)
                if current_reads_prefix != reads_prefix:
                    os.unlink(current_reads_prefix + '_1.fa')
                    os.unlink(current_reads_prefix + '_2.fa')
                current_reads_prefix = filter_prefix
                shutil.rmtree(tmpdir)

            iter_prefix = out_prefix + '.' + str(i)
            bases_added = self._extend_with_reads(current_reads_prefix, iter_prefix, no_map_contigs)

            if bases_added == 0:
                if not try_contig_trim:
                    if self.verbose:
                        print('    No bases added. Try trimming contigs')
                    self._trim_strand_biased_ends(reads_prefix, tag_as_trimmed=False)
                    if len(self.contigs) <= len(no_map_contigs):
                        if self.verbose:
                            print('       lost contigs during trimming. No more iterations')
                        return False
                    self.trim_contigs(self.contig_iter_trim)
                    try_contig_trim = True
                else:
                    if self.verbose:
                        print('    No bases added after trimming. No more iterations')
                    break
            else:
                try_contig_trim = False

        if current_reads_prefix != reads_prefix:
            os.unlink(current_reads_prefix + '_1.fa')
            os.unlink(current_reads_prefix + '_2.fa')
        return True


    def read_pair_extend(self, reads_prefix, out_prefix):
        assert(len(self.contigs))
        current_reads_prefix = reads_prefix
        i = 1
        new_seed_name = 'dummy'

        while 1:
            if self.verbose:
                print('{:_^79}'.format(' START ITERATION ' + str(i) + ' '), flush=True)
            self._read_pair_extension_iterations(current_reads_prefix, out_prefix + '.' + str(i))
            filtered_reads_prefix = out_prefix + '.' + str(i) + '.filtered'
            self._trim_strand_biased_ends(reads_prefix, tag_as_trimmed=True, out_prefix=filtered_reads_prefix)
            self._remove_contained_contigs(list(self.contigs.keys()))
            self._merge_overlapping_contigs(list(self.contigs.keys()))
            if reads_prefix != current_reads_prefix:
                os.unlink(current_reads_prefix + '_1.fa')
                os.unlink(current_reads_prefix + '_2.fa')
            current_reads_prefix = filtered_reads_prefix
            i += 1
            reads_left = os.path.getsize(current_reads_prefix + '_1.fa') > 0 and os.path.getsize(current_reads_prefix + '_2.fa') > 0

            if not self.make_new_seeds or new_seed_name is None or not self.make_new_seeds or len(self.contigs) >= self.max_contigs or not reads_left:
                if reads_prefix != current_reads_prefix:
                    os.unlink(current_reads_prefix + '_1.fa')
                    os.unlink(current_reads_prefix + '_2.fa')
                break

            if self.verbose:
                print('{:_^79}'.format(' Try making new seed '), flush=True)
            new_seed_name = self.add_new_seed_contig(current_reads_prefix + '_1.fa', current_reads_prefix + '_2.fa')

            if new_seed_name is None:
                if self.verbose:
                    print('Couldn\'t make new seed and extend it. Stopping assembly.')
                if len(self.contigs) == 0:
                    print('No contigs made.')
                    print('Read coverage may be too low, in which case try reducing --seed_min_kmer_cov, --ext_min_cov and --seed_ext_min_cov.')
                    print('Alternatively --max_insert could be incorrect, which is currently set to:', self.max_insert)
                if reads_prefix != current_reads_prefix:
                    os.unlink(current_reads_prefix + '_1.fa')
                    os.unlink(current_reads_prefix + '_2.fa')
                break



    def _run_nucmer(self, contigs_to_use=None):
        if contigs_to_use is None:
            contigs_to_use = set()
        if len(contigs_to_use) + len(self.contigs) <= 1:
            return []
        tmpdir = tempfile.mkdtemp(prefix='tmp.remove_self_contigs.', dir=os.getcwd())
        nucmer_out = os.path.join(tmpdir, 'nucmer.out')
        contigs_fasta = os.path.join(tmpdir, 'contigs.fa')
        self.write_contigs_to_file(contigs_fasta, only_write=contigs_to_use)
        mummer.run_nucmer(contigs_fasta, contigs_fasta, nucmer_out)
        hits = [hit for hit in mummer.file_reader(nucmer_out) if not hit.is_self_hit()]
        for hit in hits:
            hit.sort()
        hits = list(set(hits))
        shutil.rmtree(tmpdir)
        return hits


    def _remove_contained_contigs(self, contigs):
        if len(contigs) <= 1:
            return
        hits = self._run_nucmer(contigs_to_use=contigs)
        for contig in self._contig_names_size_order()[:-1]:
            if self._contig_contained_in_nucmer_hits(hits, contig, 95):
                hits = self._remove_contig_from_nucmer_hits(hits, contig)
                self._remove_contig(contig)
                contigs.remove(contig)


    def _coords_to_new_contig(self, coords_list):
        new_contig = pyfastaq.sequences.Fasta(coords_list[0][0], '')
        for name, coords, reverse in coords_list:
            assert name in self.contigs
            if reverse:
                seq = pyfastaq.sequences.Fasta('ni', self.contigs[name].fa.seq[coords.start:coords.end+1])
                seq.revcomp()
                new_contig.seq += seq.seq
            else:
                new_contig.seq += self.contigs[name].fa.seq[coords.start:coords.end+1]

        return new_contig


    def _merge_overlapping_contigs(self, contigs):
        if len(contigs) <= 1:
            return
        hits = self._run_nucmer(contigs)
        assembly_graph = graph.Graph(self, contigs=contigs)
        for hit in hits:
            e = hit.to_graph_edge()
            if e is not None:
                assembly_graph.add_edge(e)

        for connected_component in assembly_graph.connected_components():
            if len(connected_component) < 2:
                continue
            simple_path = assembly_graph.find_simple_path(connected_component)
            assert assembly_graph.simple_path_is_consistent(simple_path)
            if len(simple_path) > 1:
                simple_path = assembly_graph.remove_redundant_nodes_from_simple_path(simple_path)
                coords = assembly_graph.merged_coords_from_simple_nonredundant_path(simple_path)
                new_contig = self._coords_to_new_contig(coords)
                for name, x, y in coords:
                    self._remove_contig(name)
                self._add_contig(new_contig)


    def _contig_names_size_order(self, biggest_first=False):
        return sorted(self.contigs, key=lambda x:len(self.contigs[x]), reverse=biggest_first)


    def _contig_contained_in_nucmer_hits(self, hits, contig, min_percent):
        assert contig in self.contigs
        contig_length = len(self.contigs[contig])
        coords = []
        for hit in [hit for hit in hits if contig in [hit.qry_name, hit.ref_name] and hit.qry_name != hit.ref_name]:
            start = min(hit.qry_start, hit.qry_end)
            end = max(hit.qry_start, hit.qry_end)
            coords.append(pyfastaq.intervals.Interval(start, end))

        if len(coords) == 0:
            return False

        pyfastaq.intervals.merge_overlapping_in_list(coords)
        total_bases_matched = pyfastaq.intervals.length_sum_from_list(coords)
        return min_percent <= 100.0 * total_bases_matched / len(self.contigs[contig])


    def _remove_contig_from_nucmer_hits(self, hits, contig):
        return [x for x in hits if contig not in [x.ref_name, x.qry_name]]


    def _remove_contig(self, contig):
        if contig in self.contigs:
            del self.contigs[contig]
        if contig in self.contig_lengths:
            del self.contig_lengths[contig]
        if contig in self.contigs_trimmed_for_strand_bias:
            self.contigs_trimmed_for_strand_bias.remove(contig)


    def _get_ref_length(self, samfile, sam):
        if sam.is_unmapped:
            return None
        else:
            return len(self.contigs[mapping.get_ref_name(sam, samfile)])


    def _get_ref_length_sam_pair(self, samfile, sam1, sam2):
        len1 = self._get_ref_length(samfile, sam1)
        len2 = self._get_ref_length(samfile, sam2)
        if len1 == len2:
            return len1
        else:
            return None


    def _get_unmapped_pairs(self, reads1, reads2, out_prefix):
        self._map_reads(reads1, reads2, out_prefix, required_flag=12)
        mapping.bam_file_to_fasta_pair_files(out_prefix + '.bam', out_prefix + '_1.fa', out_prefix + '_2.fa')
        os.unlink(out_prefix + '.bam')


    def add_new_seed_contig(self, reads1, reads2, contig_name=None, max_attempts=10):
        if len(self.contigs):
            tmpdir = tempfile.mkdtemp(prefix='tmp.make_seed.', dir=os.getcwd())
            tmp_prefix = os.path.join(tmpdir, 'out')
            seed_reads1 = tmp_prefix + '_1.fa'
            seed_reads2 = tmp_prefix + '_2.fa'
            if contig_name is not None:
                self._map_reads(reads1, reads2, tmp_prefix, required_flag=5, exclude_flag=8, mate_ref=contig_name)
                mapping.bam_to_fasta(tmp_prefix + '.bam', seed_reads1)
                seed_reads2 = None
            else:
                self._get_unmapped_pairs(reads1, reads2, tmp_prefix)
        else:
            seed_reads1 = reads1
            seed_reads2 = reads2

        made_seed = False

        for i in range(max_attempts):
            s = seed.Seed(reads1=seed_reads1, reads2=seed_reads2, extend_length=self.seed_ext_max_bases, seed_length=self.seed_start_length, seed_min_count=self.seed_min_kmer_count, seed_max_count=self.seed_max_kmer_count, ext_min_cov=self.seed_min_cov, ext_min_ratio=self.seed_min_ratio, verbose=self.verbose, kmc_threads=self.kmc_threads, map_threads=self.threads, sequences_to_ignore=self.used_seeds, contigs_to_check=self.contigs)

            if s.seq is None or len(s.seq) == 0:
                break

            if self.seed_overlap_length is None:
                s.overlap_length = len(s.seq)
            else:
                s.overlap_length = self.seed_overlap_length
            s.extend(reads1, reads2, self.seed_stop_length)
            self.used_seeds.add(s.seq)

            if len(s.seq) >= 0.75 * self.seed_stop_length:
                made_seed = True
                break
            elif self.verbose:
                print("    Couldn't extend seed enough. That was attempt", i+1, 'of', max_attempts, flush=True)

        if len(self.contigs):
            shutil.rmtree(tmpdir)

        if not made_seed or len(s.seq) == 0:
            return None

        if self.verbose:
            print("    Extended seed OK.", flush=True)
        new_name = 'seeded.' + '1'.zfill(5)
        i = 1
        while new_name in self.contigs:
            i += 1
            new_name = 'seeded.' + str(i).zfill(5)

        self._add_contig(pyfastaq.sequences.Fasta(new_name, s.seq))
        return new_name


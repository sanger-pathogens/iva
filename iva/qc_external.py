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
import subprocess
import tempfile
import shutil
import os
import sys
import inspect
import pyfastaq
import iva

class Error (Exception): pass

gage_stats = [
    'Missing Reference Bases',
    'Missing Assembly Bases',
    'Missing Assembly Contigs',
    'Duplicated Reference Bases',
    'Compressed Reference Bases',
    'Bad Trim',
    'Avg Idy',
    'SNPs',
    'Indels < 5bp',
    'Indels >= 5',
    'Inversions',
    'Relocation',
    'Translocation',
]


ratt_stats = [
     'elements_found',
     'elements_transferred',
     'elements_transferred_partially',
     'elements_split',
     'parts_of_elements_not_transferred',
     'elements_not_transferred',
     'gene_models_to_transfer',
     'gene_models_transferred',
     'gene_models_transferred_partially',
     'exons_not_transferred_from_partial_matches',
     'gene_models_not_transferred',
]


def dummy_gage_stats():
    return {x:'NA' for x in gage_stats}


def dummy_ratt_stats():
    return {x:'NA' for x in ratt_stats}


def run_gage(reference, scaffolds, outdir, nucmer_minid=80, clean=True):
    reference = os.path.abspath(reference)
    scaffolds = os.path.abspath(scaffolds)
    ref = 'ref.fa'
    scaffs = 'scaffolds.fa'
    contigs = 'contigs.fa'
    gage_out = 'gage.out'
    gage_script = 'run.sh'
    cwd = os.getcwd()
    os.mkdir(outdir)
    os.chdir(outdir)
    extractor = iva.egg_extract.Extractor(os.path.abspath(os.path.join(os.path.dirname(iva.__file__), os.pardir)))
    gage_code_indir = os.path.join('iva', 'gage')
    gage_code_outdir = 'gage_code'
    extractor.copy_dir(gage_code_indir, gage_code_outdir)
    os.symlink(reference, ref)
    os.symlink(scaffolds, scaffs)
    pyfastaq.tasks.scaffolds_to_contigs(scaffs, contigs, number_contigs=True)
    f = pyfastaq.utils.open_file_write(gage_script)
    print(' '.join([
        'sh',
        os.path.join(gage_code_outdir, 'getCorrectnessStats.sh'),
        ref,
        contigs,
        scaffs,
        str(nucmer_minid),
        '>', gage_out
        ]), file=f)
    pyfastaq.utils.close(f)
    iva.common.syscall('bash ' + gage_script, allow_fail=True)
    if not os.path.exists(gage_out):
        raise Error('Error running GAGE\nbash ' + gage_script)
    stats = dummy_gage_stats()
    wanted_stats = set(gage_stats)
    f = pyfastaq.utils.open_file_read(gage_out)
    got_all_stats = False

    for line in f:
        if line.startswith('Corrected Contig Stats'):
            got_all_stats = True
            break
        elif ':' in line:
            a = line.rstrip().split(': ')
            if a[0] in wanted_stats:
                stat = a[1]
                if '%' in stat:
                    stat = stat.split('(')[0]
                if stat.isdigit():
                    stats[a[0]] = int(stat)
                else:
                    stats[a[0]] = float(stat)
    pyfastaq.utils.close(f)

    if not got_all_stats:
        raise Error('Error running GAGE\nbash ' + gage_script)

    if clean:
        to_clean = [
            'contigs.fa.delta',
            'contigs.fa.fdelta',
            'contigs.fa.matches.lens',
            'out.1coords',
            'out.1delta',
            'out.mcoords',
            'out.mdelta',
            'out.qdiff',
            'out.rdiff',
            'out.snps',
            'out.unqry',
            'scaffolds.fa.coords',
            'scaffolds.fa.delta',
            'scaffolds.fa.err',
            'scaffolds.fa.fdelta',
            'scaffolds.fa.tiling',
            'tmp_scf.fasta',
        ]
        for f in to_clean:
            try:
                os.unlink(f)
            except:
                pass

    os.chdir(cwd)
    return stats


def run_ratt(embl_dir, assembly, outdir, config_file=None, transfer='Species', clean=True):
    embl_dir = os.path.abspath(embl_dir)
    assembly = os.path.abspath(assembly)

    cwd = os.getcwd()
    if os.path.exists(outdir):
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        else:
            os.unlink(outdir)
    os.mkdir(outdir)
    os.chdir(outdir)

    extractor = iva.egg_extract.Extractor(os.path.abspath(os.path.join(os.path.dirname(iva.__file__), os.pardir)))
    ratt_code_indir = os.path.join('iva', 'ratt')
    ratt_code_outdir = 'ratt_code'
    extractor.copy_dir(ratt_code_indir, ratt_code_outdir)

    if config_file is None:
        ratt_config = os.path.join(ratt_code_outdir, 'ratt.config')
    else:
        ratt_config = os.path.abspath(config_file)

    script = 'run.sh'
    script_out = 'run.sh.out'
    ratt_outprefix = 'out'
    f = pyfastaq.utils.open_file_write(script)
    print('export RATT_HOME=', ratt_code_outdir, sep='', file=f)
    print('export RATT_CONFIG=', ratt_config, sep='', file=f)
    print('for x in $RATT_HOME/*.{sh,pl}; do chmod 755 $x; done', file=f)
    print('$RATT_HOME/start.ratt.sh', embl_dir, assembly, ratt_outprefix, transfer, file=f)
    pyfastaq.utils.close(f)
    cmd = 'bash ' + script + ' > ' + script_out
    # sometimes ratt returns nonzero code, but is OK, so ignore it
    iva.common.syscall(cmd, allow_fail=True)

    stats = {}

    matches = {
        'elements found.': 'elements_found',
        'Elements were transfered.': 'elements_transferred',
        'Elements could be transfered partially.': 'elements_transferred_partially',
        'Elements split.': 'elements_split',
        'Parts of elements (i.e.exons tRNA) not transferred.': 'parts_of_elements_not_transferred',
        'Elements couldn\'t be transferred.': 'elements_not_transferred',
        'Gene models to transfer.': 'gene_models_to_transfer',
        'Gene models transferred correctly.': 'gene_models_transferred',
        'Gene models partially transferred.': 'gene_models_transferred_partially',
        'Exons not transferred from partial CDS matches.': 'exons_not_transferred_from_partial_matches',
        'Gene models not transferred.': 'gene_models_not_transferred',
    }

    f = pyfastaq.utils.open_file_read(script_out)
    for line in f:
        if '\t' in line:
            a = line.rstrip().split('\t')
            if len(a) == 2 and a[0].isdigit() and a[1] in matches:
                stats[matches[a[1]]] = int(a[0])
    pyfastaq.utils.close(f)

    if clean:
        for d in ['Query', 'Reference', 'Sequences']:
            try:
                shutil.rmtree(d)
            except:
                pass

        iva.common.syscall('rm -f query.* Reference.* nucmer.* out.*')

    os.chdir(cwd)
    return stats


def run_blastn_and_write_act_script(assembly, reference, blast_out, script_out):
    tmpdir = tempfile.mkdtemp(prefix='tmp.blastn.', dir=os.getcwd())
    assembly_union = os.path.join(tmpdir, 'assembly.union.fa')
    reference_union = os.path.join(tmpdir, 'reference.union.fa')
    pyfastaq.tasks.to_fasta_union(assembly, assembly_union, seqname='assembly_union')
    pyfastaq.tasks.to_fasta_union(reference, reference_union, seqname='reference_union')
    iva.common.syscall('makeblastdb -dbtype nucl -in ' + reference_union)
    cmd = ' '.join([
        'blastn',
        '-task blastn',
        '-db', reference_union,
        '-query', assembly_union,
        '-out', blast_out,
        '-outfmt 6',
        '-evalue 0.01',
        '-dust no',
    ])
    iva.common.syscall(cmd)

    f = pyfastaq.utils.open_file_write(script_out)
    print('#!/usr/bin/env bash', file=f)
    print('act', reference, blast_out, assembly, file=f)
    pyfastaq.utils.close(f)
    iva.common.syscall('chmod 755 ' + script_out)
    shutil.rmtree(tmpdir)

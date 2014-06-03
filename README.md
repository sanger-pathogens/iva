IVA - Iterative Virus Assembler
===============================

__Warning: this software is in early development and is not yet recommended
for use by the community. It is undergoing testing and
an official release will be made soon.__

IVA is a _de novo_ assembler designed to assemble virus genomes that
have no repeat sequences, using Illumina read pairs sequenced from
mixed populations at extremely high depth.

IVA's main algorithm works by iteratively extending contigs using aligned
read pairs. Its input can be just read pairs, or additionally you can
provide an existing set of contigs to be extended. Alternatively,
it can take reads together with a reference sequence.

------------------------------------------------------------------------------

Dependencies
------------

IVA has been developed for and tested on Linux and relies on
some third-party tools that need to be installed first.
For citations, see the References section at the bottom of this readme.


#### Assembly dependencies

The following are required to run an assembly with IVA.

 * [Python 3] [python] (IVA is written in Python 3).
 * The following Python packages:
     * [Fastaq] [fastaq]
     * [networkx] [networkx]
     * [Pysam] [pysam]
 * [KMC] [kmc code] installed, so that `kmc` and `kmc_dump` are in your path.
   Read the [paper] [kmc paper] for more information.
 * [MUMmer] [mummer code] installed with its executables (ie `nucmer` etc)
   in your path. Read the [paper] [mummer paper] for more information.
 * [Samtools] [samtools code] installed, so that `samtools` is in your path.
   Read the [paper] [samtools paper] for more information.
 * [SMALT] [smalt] installed, so that `smalt` is in your path.
 * Optional: [Trimmomatic] [trimmo code] - although this is optional, it is
   highly recommended.
   It is used to trim adapter sequences from reads before assembling and
   significantly improves the results. Read the [paper] [trimmo paper]
   for more information. You don't need to add anything to your path, but will
   need to tell IVA where the Java jar file is to use Trimmomatic (see
   examples below).

#### QC dependencies

The QC scripts have the following dependencies, in addiition to MUMmer,
smalt and samtools:

 * [R] [r code] installed and in your path.
 * Optional: [kraken] [kraken code] installed, so that `kraken` and
   `kraken-build` are in your path. These are needed if you want to
   make your own reference database, or if you use a database to
   automatically choose the reference genome. Read the
   [paper] [kraken paper] for more information.
 * Optional: [reapr] [reapr code] will be used if it is installed
   and found in your path. Read the [paper] [reapr paper] for more
   information.

The QC code is also bundled with the following (they do not need to be installed).

 * Analysis code from the [GAGE] [gage code] assembly evaluation
   project. We are grateful to the GAGE authors for permission to modify and
   redistribute this
   code. Read the [paper] [gage paper] for more information.
 * [RATT] [ratt code] is used to transfer annotation from a reference
   onto the assembly. Read the [paper] [ratt paper] for more information.

------------------------------------------------------------------------------

Installation
------------

Install the dependencies and take a copy of the latest IVA release.
Then run the tests:

    python3 setup.py test

If all the tests pass, then install with:

    python3 setup.py install

Or if you don't have root access, then run:

    python3 setup.py install --prefix /install/in/here

------------------------------------------------------------------------------

Examples
--------

For any of the IVA scripts, use `-h` or `--help` to see all the options.

#### Assembly examples

To run a de novo assembly on a set of read pairs, with forwards and
reverse reads in two separate files:

    iva -f reads_fwd.fastq -r reads_rev.fastq Output_directory

If instead you have a single file of interleaved paired reads:

    iva --fr reads.fastq Output_directory

It is highly recommended to trim adapters off the reads before assembling.
IVA is bundled with a file of adapters, but you may get better
results if you know your own adapter sequences and use those instead.
We need to tell IVA where the trimmomatic jar file is:

    iva --trimmo /path/to/trimmomatic-0.32.jar --fr reads.fastq Output_directory

If you do have your own file of adapters:

    iva --trimmo /path/to/trimmomatic-0.32.jar --adapters adapters.fasta --fr reads.fastq Output_directory

To extend an existing FASTA file of contigs:

    iva --contigs contigs.fasta --fr reads.fastq Output_directory

Input files can be gzipped - IVA will do the right thing with files that
have a .gz extension:

    iva --contigs contigs.fasta.gz --fr reads.fastq.gz Output_directory

Use 8 threads instead of the default 1 thread:

    iva --threads 8 --fr reads.fastq Output_directory

Increase the maximum insert size from 500bp (default) to 1000bp:

    iva --max_insert 1000 --fr reads.fastq Output_directory


#### QC examples

The script `iva_qc` needs a reference genome that is used to compare
with the assembly. You can supply your own, in the form of a directory
of EMBL files:

    iva_qc --embl_dir /dir/of/embl/files/ assembly.fasta prefix_of_output_files

Alternatively, the closest reference to the input assembly can be found
automatically (using kraken). Build a database in a new directory called
`Database_dir` with:

    iva_qc_make_db Database_dir

This script will downlaod all required data. If it dies then on restarting it,
it will continue where it last finished. Once you have made the database,
you can use it with the QC script by running:

    iva_qc --ref_db Database_dir assembly.fasta prefix_of_output_files


------------------------------------------------------------------------------

References
----------

**GAGE:** Salzberg, S. L. et al. _GAGE: A critical evaluation of genome
assemblies and assembly algorithms_. Genome Res. 22, 557-67 (2012).

**KMC:** Deorowicz, S., Debudaj-Grabysz, A. & Grabowski, S. _Disk-based k-mer
counting on a PC_. BMC Bioinformatics 14, 160 (2013).

**Kraken:** Wood, D. E. & Salzberg, S. L. _Kraken: ultrafast metagenomic
sequence classification using exact alignments_.
Genome Biol. 15, R46 (2014).

**MUMmer:** Kurtz, S. et al. _Versatile and open software for comparing large
genomes_. Genome Biol. 5, R12 (2004).

**R:** R Core Team (2013). _R: A language and environment for statistical
computing_. R Foundation for Statistical Computing, Vienna, Austria.
URL http://www.R-project.org/.

**RATT:** Otto, T. D., Dillon, G. P., Degrave, W. S. & Berriman, M.
_RATT: Rapid Annotation Transfer Tool_. Nucleic Acids Res. 39, e57 (2011).

**REAPR:** Hunt, M. et al. _REAPR: a universal tool for genome
assembly evaluation_. Genome Biol. 14, R47 (2013).

**SAMtools:** Li, H. et al. _The Sequence Alignment/Map format and SAMtools_.
Bioinformatics 25, 2078-9 (2009).

**Trimmomatic:** Bolger, A. M., Lohse, M. & Usadel, B. _Trimmomatic: A
flexible trimmer for Illumina Sequence Data_. Bioinformatics 1-7 (2014).



  [fastaq]: https://github.com/sanger-pathogens/Fastaq
  [networkx]: https://pypi.python.org/pypi/networkx/
  [pysam]: https://code.google.com/p/pysam/
  [python]: http://www.python.org/
  [gage code]: http://gage.cbcb.umd.edu/index.html
  [gage paper]: http://genome.cshlp.org/content/early/2012/01/12/gr.131383.111
  [kmc paper]: http://www.biomedcentral.com/1471-2105/14/160
  [kmc code]: http://sun.aei.polsl.pl/kmc/download.html
  [kraken code]: http://ccb.jhu.edu/software/kraken/
  [kraken paper]: http://genomebiology.com/2014/15/3/R46
  [mummer code]: http://mummer.sourceforge.net/
  [mummer paper]: http://genomebiology.com/2004/5/2/r12
  [r code]: http://www.r-project.org/
  [ratt code]: http://ratt.sourceforge.net/
  [ratt paper]: http://nar.oxfordjournals.org/content/39/9/e57
  [reapr code]: https://www.sanger.ac.uk/resources/software/reapr/
  [reapr paper]: http://genomebiology.com/2013/14/5/R47/
  [samtools code]: http://samtools.sourceforge.net/
  [samtools paper]: http://bioinformatics.oxfordjournals.org/content/25/16/2078.abstract
  [smalt]: http://www.sanger.ac.uk/resources/software/smalt/
  [trimmo code]: http://www.usadellab.org/cms/?page=trimmomatic
  [trimmo paper]: http://bioinformatics.oxfordjournals.org/content/early/2014/04/12/bioinformatics.btu170

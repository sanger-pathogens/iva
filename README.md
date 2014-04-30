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

Dependencies
------------

IVA has been developed for and tested on Linux and relies on
some third-party tools that need to be installed first.

The following are required to run IVA:

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
 * [Trimmomatic] [trimmo code] - this is optional, but highly recommended.
   It is used to trim adapter sequences from reads before assembling and
   significantly improves the results. Read the [paper] [trimmo paper]
   for more information. You don't need to add anything to your path, but will
   need to tell IVA where the Java jar file is to use Trimmomatic (see
   examples below).


Installation
------------

Install the dependencies and take a copy of the latest IVA release.
Then run the tests:

    python3 setup.py test

If all the tests pass, then install with:

    python3 setup.py install

Or if you don't have root access, then run:

    python3 setup.py install --prefix /install/in/here


Examples
--------

Run with `-h` or `--help` to see all the options:

    iva --help

To run a de novo assembly on a run of read pairs, with forwards and
reverse reads in two separate files:

    iva -f reads_fwd.fastq -r reads_rev.fastq Output_directory

If instead you have a single file of interleaved paired reads:

    iva --fr reads.fastq Output_directory

It is highly recommended to trim adapters off the reads before assembling.
We cannot provide a file of adapters and so you must make your own. Trimmomatic
does come with a file of adapter sequences (but you may get better
results if you know your own adapter sequences and use those instead).
We need to tell IVA where the trimmomatic jar file is and also the adapters
file:

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



  [fastaq]: https://github.com/sanger-pathogens/Fastaq
  [networkx]: https://pypi.python.org/pypi/networkx/
  [pysam]: https://code.google.com/p/pysam/
  [python]: http://www.python.org/
  [kmc paper]: http://www.biomedcentral.com/1471-2105/14/160
  [kmc code]: http://sun.aei.polsl.pl/kmc/download.html
  [mummer code]: http://mummer.sourceforge.net/
  [mummer paper]: http://genomebiology.com/2004/5/2/r12
  [samtools code]: http://samtools.sourceforge.net/
  [samtools paper]: http://bioinformatics.oxfordjournals.org/content/25/16/2078.abstract
  [smalt]: http://www.sanger.ac.uk/resources/software/smalt/
  [trimmo code]: http://www.usadellab.org/cms/?page=trimmomatic
  [trimmo paper]: http://bioinformatics.oxfordjournals.org/content/early/2014/04/12/bioinformatics.btu170

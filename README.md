IVA - Iterative Virus Assembler
===============================

IVA is a _de novo_ assembler designed to assemble virus genomes that
have no repeat sequences, using Illumina read pairs sequenced from
mixed populations at extremely high and variable depth.

Installation instructions are below. For usage help and examples, see
the [IVA wiki page] [IVA wiki page].

------------------------------------------------------------------------------

Dependencies
------------

IVA has been developed for and tested on Linux and relies on
some third-party tools that need to be installed first.
For citations, see the References section at the bottom of this readme.


#### Assembly dependencies

The following are required to run an assembly with IVA.

 * [Python 3] [python] (IVA is written in Python 3) and the following packages:
     * [Fastaq] [fastaq]
     * [networkx] [networkx]
     * [Pysam] [pysam]
 * [KMC] [kmc code] installed, so that `kmc` and `kmc_dump` are in your path.
 * [MUMmer] [mummer code] installed with its executables (ie `nucmer` etc)
   in your path.
 * [Samtools] [samtools code] installed, so that `samtools` is in your path.
 * [SMALT] [smalt] installed, so that `smalt` is in your path.
 * Optional: [Trimmomatic] [trimmo code] - although this is optional, it is
   highly recommended.
   It is used to trim adapter sequences from reads before assembling and
   significantly improves the results.
   You don't need to add anything to your path, but will
   need to tell IVA where the Java jar file is to use Trimmomatic (see
   [examples] [IVA wiki examples]).

IVA is known to work with kmc version 2.0, MUMmer version 3.23, samtools
version 0.1.19 and SMALT version 0.7.5.

#### QC dependencies

The QC scripts have the following dependencies, in addiition to MUMmer,
smalt and samtools:

 * [R] [r code] installed and in your path.
 * Optional: [kraken] [kraken code] installed, so that `kraken` and
   `kraken-build` are in your path. These are needed if you want to
   make your own reference database, or if you use a database to
   automatically choose the reference genome.

The QC code is also bundled with the following (they do not need to be installed).

 * Analysis code from the [GAGE] [gage code] assembly evaluation
   project. We are grateful to the GAGE authors for permission to modify and
   redistribute this
   code.
 * [RATT] [ratt code] is used to transfer annotation from a reference
   onto the assembly.

------------------------------------------------------------------------------

Installation
------------

Install the dependencies and take a copy of the [latest IVA release] [IVA latest release].
Then run the tests:

    python3 setup.py test

If all the tests pass, then install with:

    python3 setup.py install

Or if you don't have root access, then run:

    python3 setup.py install --prefix /install/in/here

------------------------------------------------------------------------------

References
----------

[**Adapter sequences:**] [adapters paper] Quail, M. a et al. _Optimal enzymes for amplifying sequencing libraries_. Nat. Methods 9, 10-1 (2012).

[**GAGE:**] [gage paper] Salzberg, S. L. et al. _GAGE: A critical evaluation of genome
assemblies and assembly algorithms_. Genome Res. 22, 557-67 (2012).

[**KMC:**] [kmc paper] Deorowicz, S., Debudaj-Grabysz, A. & Grabowski, S. _Disk-based k-mer
counting on a PC_. BMC Bioinformatics 14, 160 (2013).

[**Kraken:**] [kraken paper] Wood, D. E. & Salzberg, S. L. _Kraken: ultrafast metagenomic
sequence classification using exact alignments_.
Genome Biol. 15, R46 (2014).

[**MUMmer:**] [mummer paper] Kurtz, S. et al. _Versatile and open software for comparing large
genomes_. Genome Biol. 5, R12 (2004).

**R:** R Core Team (2013). _R: A language and environment for statistical
computing_. R Foundation for Statistical Computing, Vienna, Austria.
URL http://www.R-project.org/.

[**RATT:**] [ratt paper] Otto, T. D., Dillon, G. P., Degrave, W. S. & Berriman, M.
_RATT: Rapid Annotation Transfer Tool_. Nucleic Acids Res. 39, e57 (2011).

[**SAMtools:**] [samtools paper] Li, H. et al. _The Sequence Alignment/Map format and SAMtools_.
Bioinformatics 25, 2078-9 (2009).

[**Trimmomatic:**] [trimmo paper] Bolger, A. M., Lohse, M. & Usadel, B. _Trimmomatic: A
flexible trimmer for Illumina Sequence Data_. Bioinformatics 1-7 (2014).


  [adapters paper]: http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1814.html
  [IVA wiki page]: https://github.com/sanger-pathogens/iva/wiki
  [IVA wiki examples]: https://github.com/sanger-pathogens/iva/wiki/iva-examples
  [IVA latest release]: https://github.com/sanger-pathogens/iva/releases/latest
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
  [samtools code]: http://samtools.sourceforge.net/
  [samtools paper]: http://bioinformatics.oxfordjournals.org/content/25/16/2078.abstract
  [smalt]: http://www.sanger.ac.uk/resources/software/smalt/
  [trimmo code]: http://www.usadellab.org/cms/?page=trimmomatic
  [trimmo paper]: http://bioinformatics.oxfordjournals.org/content/early/2014/04/12/bioinformatics.btu170

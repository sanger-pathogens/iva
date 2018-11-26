# IVA
Iterative Virus Assembler - de novo virus assembler of Illumina paired reads.

[![Build Status](https://travis-ci.org/sanger-pathogens/iva.svg?branch=master)](https://travis-ci.org/sanger-pathogens/iva)  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/iva/blob/master/LICENSE)  
[![status](https://img.shields.io/badge/Bioinformatics-10.1093-brightgreen.svg)](https://academic.oup.com/bioinformatics/article/31/14/2374/254470)  
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/iva/README.html)   
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/iva)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/iva.svg)](https://hub.docker.com/r/sangerpathogens/iva)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/iva.svg)](https://hub.docker.com/r/sangerpathogens/iva)  
[![codecov](https://codecov.io/gh/sanger-pathogens/iva/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/iva)  

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
  * [Running the tests](#running-the-tests)
  * [Usage](#usage)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)

## Introduction
IVA is a de novo assembler designed to assemble virus genomes that have no repeat sequences, using Illumina read pairs sequenced from mixed populations at extremely high and variable depth.

For more information, please read the [IVA publication](http://bioinformatics.oxfordjournals.org/content/early/2015/02/27/bioinformatics.btv120.abstract).

## Installation
For installation instructions, please refer to the [IVA website](http://sanger-pathogens.github.io/iva/)

## Running the tests
The test can be run with dzil from the top level directory:  

`python setup.py test`

## Usage
```
usage: iva [options] {-f reads_fwd -r reads_rev | --fr reads} <output directory>

positional arguments:
  Output directory      Name of output directory (must not already exist)

optional arguments:
  -h, --help            show this help message and exit

Input and output:
  -f filename[.gz], --reads_fwd filename[.gz]
                        Name of forward reads fasta/q file. Must be used in
                        conjunction with --reads_rev
  -r filename[.gz], --reads_rev filename[.gz]
                        Name of reverse reads fasta/q file. Must be used in
                        conjunction with --reads_fwd
  --fr filename[.gz]    Name of interleaved fasta/q file
  --keep_files          Keep intermediate files (could be many!). Default is
                        to delete all unnecessary files
  --contigs filename[.gz]
                        Fasta file of contigs to be extended. Incompatible
                        with --reference
  --reference filename[.gz]
                        EXPERIMENTAL! This option is EXPERIMENTAL, not
                        recommended, and has not been tested! Fasta file of
                        reference genome, or parts thereof. IVA will try to
                        assemble one contig per sequence in this file.
                        Incompatible with --contigs
  -v, --verbose         Be verbose by printing messages to stdout. Use up to
                        three times for increasing verbosity.

SMALT mapping options:
  -k INT, --smalt_k INT
                        kmer hash length in SMALT (the -k option in smalt
                        index) [19]
  -s INT, --smalt_s INT
                        kmer hash step size in SMALT (the -s option in smalt
                        index) [11]
  -y FLOAT, --smalt_id FLOAT
                        Minimum identity threshold for mapping to be reported
                        (the -y option in smalt map) [0.5]

Contig options:
  --ctg_first_trim INT  Number of bases to trim off the end of every contig
                        before extending for the first time [25]
  --ctg_iter_trim INT   During iterative extension, number of bases to trim
                        off the end of a contig when extension fails (then try
                        extending again) [10]
  --ext_min_cov INT     Minimum kmer depth needed to use that kmer to extend a
                        contig [10]
  --ext_min_ratio FLOAT
                        Sets N, where kmer for extension must be at least N
                        times more abundant than next most common kmer [4]
  --ext_max_bases INT   Maximum number of bases to try to extend on each
                        iteration [100]
  --ext_min_clip INT    Set minimum number of bases soft clipped off a read
                        for those bases to be used for extension [3]
  --max_contigs INT     Maximum number of contigs allowed in the assembly. No
                        more seeds generated if the cutoff is reached [50]

Seed generation options:
  --make_new_seeds      When no more contigs can be extended, generate a new
                        seed. This is forced to be true when --contigs is not
                        used
  --seed_start_length INT
                        When making a seed sequence, use the most common kmer
                        of this length. Default is to use the minimum of
                        (median read length, 95). Warning: it is not
                        recommended to set this higher than 95
  --seed_stop_length INT
                        Stop extending seed using perfect matches from reads
                        when this length is reached. Future extensions are
                        then made by treating the seed as a contig
                        [0.9*max_insert]
  --seed_min_kmer_cov INT
                        Minimum kmer coverage of initial seed [25]
  --seed_max_kmer_cov INT
                        Maximum kmer coverage of initial seed [1000000]
  --seed_ext_max_bases INT
                        Maximum number of bases to try to extend on each
                        iteration [50]
  --seed_overlap_length INT
                        Number of overlapping bases needed between read and
                        seed to use that read to extend [seed_start_length]
  --seed_ext_min_cov INT
                        Minimum kmer depth needed to use that kmer to extend a
                        contig [10]
  --seed_ext_min_ratio FLOAT
                        Sets N, where kmer for extension must be at least N
                        times more abundant than next most common kmer [4]

Read trimming options:
  --trimmomatic FILENAME
                        Provide location of trimmomatic.jar file to enable
                        read trimming. Required if --adapters used
  --trimmo_qual STRING  Trimmomatic options used to quality trim reads
                        [LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20]
  --adapters FILENAME   Fasta file of adapter sequences to be trimmed off
                        reads. If used, must also use --trimmomatic. Default
                        is file of adapters supplied with IVA
  --min_trimmed_length INT
                        Minimum length of read after trimming [50]
  --pcr_primers FILENAME
                        FASTA file of primers. The first perfect match found
                        to a sequence in the primers file will be trimmed off
                        the start of each read. This is run after trimmomatic
                        (if --trimmomatic used)

Other options:
  -i INT, --max_insert INT
                        Maximum insert size (includes read length). Reads with
                        inferred insert size more than the maximum will not be
                        used to extend contigs [800]
  -t INT, --threads INT
                        Number of threads to use [1]
  --kmc_onethread       Force kmc to use one thread. By default the value of
                        -t/--threads is used when running kmc
  --strand_bias FLOAT in [0,0.5]
                        Set strand bias cutoff of mapped reads when trimming
                        contig ends, in the interval [0,0.5]. A value of x
                        means that a base needs min(fwd_depth, rev_depth) /
                        total_depth <= x. The only time this should be used is
                        with libraries with overlapping reads (ie fragment
                        length < 2*read length), and even then, it can make
                        results worse. If used, try a low value like 0.1 first
                        [0]
  --test                Run using built in test data. All other options will
                        be ignored, except the mandatory output directory, and
                        --trimmomatic and --threads can be also be used
  --version             show program's version number and exit
```

For usage help and examples, see the [IVA wiki page](https://github.com/sanger-pathogens/iva/wiki).

## License
IVA is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/iva/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/iva/issues) or email iva-help@sanger.ac.uk.

## Citation
If you use this software please cite:

__IVA: accurate de novo assembly of RNA virus genomes.__  
Hunt M, Gall A, Ong SH, Brener J, Ferns B, Goulder P, Nastouli E, Keane JA, Kellam P, Otto TD.  
Bioinformatics. 2015 Jul 15;31(14):2374-6. doi: [10.1093/bioinformatics/btv120](http://bioinformatics.oxfordjournals.org/content/31/14/2374.long). Epub 2015 Feb 28.  

[__Adapter sequences__](http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1814.html):  
__Optimal enzymes for amplifying sequencing libraries.__  
Quail, M. a et al. Nat. Methods 9, 10-1 (2012).

[__GAGE__](http://genome.cshlp.org/content/early/2012/01/12/gr.131383.111):  
__GAGE: A critical evaluation of genome assemblies and assembly algorithms.__  
Salzberg, S. L. et al. Genome Res. 22, 557-67 (2012).

[__KMC__](http://www.biomedcentral.com/1471-2105/14/160):  
__Disk-based k-mer counting on a PC.__  
Deorowicz, S., Debudaj-Grabysz, A. & Grabowski, S. BMC Bioinformatics 14, 160 (2013).

[__Kraken__](http://genomebiology.com/2014/15/3/R46):  
__Kraken: ultrafast metagenomic sequence classification using exact alignments.__  
Wood, D. E. & Salzberg, S. L. Genome Biol. 15, R46 (2014).

[__MUMmer__](http://genomebiology.com/2004/5/2/r12):  
__Versatile and open software for comparing large genomes.__  
Kurtz, S. et al. Genome Biol. 5, R12 (2004).

__R__:  
__R: A language and environment for statistical computing.__  
R Core Team (2013). R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.

[__RATT__](http://nar.oxfordjournals.org/content/39/9/e57):  
__RATT: Rapid Annotation Transfer Tool.__  
Otto, T. D., Dillon, G. P., Degrave, W. S. & Berriman, M. Nucleic Acids Res. 39, e57 (2011).

[__SAMtools__](http://bioinformatics.oxfordjournals.org/content/25/16/2078.abstract):   
__The Sequence Alignment/Map format and SAMtools.__  
Li, H. et al. Bioinformatics 25, 2078-9 (2009).

[__Trimmomatic__](http://bioinformatics.oxfordjournals.org/content/early/2014/04/12/bioinformatics.btu170):  
__Trimmomatic: A flexible trimmer for Illumina Sequence Data.__  
Bolger, A. M., Lohse, M. & Usadel, B. Bioinformatics 1-7 (2014).
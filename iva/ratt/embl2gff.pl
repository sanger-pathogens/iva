#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

if (@ARGV != 1) {    die "USAGE: embl2gff.pl   > outputfile.\n"; }

my $in = Bio::SeqIO->new(-file=>$ARGV[0],-format=>'EMBL');
while (my $seq = $in->next_seq) {
    for my $feat ($seq->top_SeqFeatures) {
        print $feat->gff_string,"\n";
    }
}

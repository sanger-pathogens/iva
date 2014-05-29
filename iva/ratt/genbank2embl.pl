#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

if (@ARGV != 2) {    die "USAGE: gb2embl.pl <in.gbk> <out.embl>\n"; }

my $seqio = Bio::SeqIO->new('-format' => 'genbank', '-file' => "$ARGV[0]");
my $seqout = new Bio::SeqIO('-format' => 'embl', '-file' => ">$ARGV[1]");
while( my $seq = $seqio->next_seq) {
  $seqout->write_seq($seq)
}

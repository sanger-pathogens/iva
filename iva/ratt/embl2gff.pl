#!/usr/bin/env perl
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

#! /usr/bin/perl -w
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
use Data::Dumper;
use lib $ENV{RATT_HOME};
use ratt_correction;

my $debug=1;

my $SET=1;
my $COLOR_BAD=4;

if (!defined($ENV{RATT_HOME})) {
  print "Please set global variable RATT_HOME in your shell!\n";
  exit 1
}

if (!defined($ARGV[0])) {
  print "Sorry, wrong option.\n";
  print "Tranfer / Correct / Check / EMBLFormatCheck / Mutate / Split / Difference / Embl2Fasta / doEMBL \ncan be used.\n\n";
  exit 1
}
elsif ($ARGV[0] eq "doEMBL") {
  if (! defined($ARGV[1])) {
	print "\n\nusage: \$RATT_HOME/main.ratt.pl <ResultName> <embl-Annotation> <fasta-file>\n\n".
	  "This part will generate a vaild embl file, with an ID line and the sequence. The resultfile will be called <ResultName>.embl\n\n";
	
	exit;
	
  }
  doEMBL($ARGV[1],$ARGV[2],$ARGV[3]);
  
  exit;
}
elsif ($ARGV[0] eq "Mutate") {
  if (! defined($ARGV[1])) {
	print "\n\nusage: \$RATT_HOME/main.ratt.pl Mutate <(multi-)fasta-file>\n\n".
	  "Every 250 base pairs a base is changed (mutated). The result is saved as <fastafile>.mutated. This is necessary to recalibrate RATT for similar genomes.\n\n";
	
	exit;
	
  }
  putMutation($ARGV[1]);
  
  exit;
}
elsif ($ARGV[0] eq "Split") {
  if (! defined($ARGV[1])) {
	print "\n\nusage: \$RATT_HOME/main.ratt.pl Split <(multifasta-file>\n\n".
	  "Splits a given multifasta file into individual files containing one sequence. This is necessary as visualization tools (e.g. Artemis) prefer single fasta files.\n\n";
	
	exit (1);
	
  }
  Split($ARGV[1]);
  
  exit;
}
elsif ($ARGV[0] eq "Difference") {
  my $mummerSNP=$ARGV[1];
  my $mummerCoords=$ARGV[2];
  my $resultName=$ARGV[3];
  
  if (!defined($resultName)) {
  	print "\n\nusage: \$RATT_HOME/main.ratt.pl Difference <mummer SNP file> <mummer coord file> <ResultName>\n\n".
	  "Generates files that report the SNP, indels and regions not shared by the reference and query. It also prints a statistic reporting coverage for each replicon.\n\n";
  	exit (1);
  }
  my $stats=getMutations($mummerSNP,$mummerCoords,$resultName);
  
  print $stats;
  exit;
  
}
elsif($ARGV[0] eq "addTranslation"){
   if (scalar(@ARGV) < 2) {
 	print "\n\nusage: \$RATT_HOME/main.ratt.pl addTranslation <EMBL file> optional <Translationtable>\n\n".
	  "If you want to add the translation of the protein relative to the new sequence, run this. The Bioperl module Bioperl::seqio must be installed. Also, the CDS must have a gene, locustag or systematic_id tag.\nThe result will be the embl with a translationtag.\n";
	
  	exit (1); 
  }
  addTranslation($ARGV[1],$ARGV[2]);
  
}
elsif ($ARGV[0] eq "EMBLFormatCheck") {
  if (scalar(@ARGV) < 3) {
 	print "\n\nusage: \$RATT_HOME/main.ratt.pl EMBLFormatCheck <EMBL file> <ResultName postfix>\n\n".
	  "Some EMBL files have feature positions spanning several lines, this function consolidates these features so they appear on one line. The result name is <EMBL File>.<ResultName postfix>.\n\n";
	
  	exit (1); 
  }
  my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  my $resultName = shift;
  
  correctEMBL($embl,"tmp.BBA.embl",$fasta);
}
elsif ($ARGV[0] eq "Correct") {
  if (scalar(@ARGV) < 4) {
	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Correct <EMBL file> <fasta file> <ResultName>\n\n".
	  "Corrects a given annotation, as described previously. The corrections are reported and the new file is saved as <ResultName>.embl.\n\n";	
	exit(1); 
  }
  
  my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  my $resultName = shift;
  
  correctEMBL($embl,"tmp.BBA.embl",$fasta);
  
  startAnnotationCorrection( "$embl.tmp.BBA.embl",$fasta,$resultName);
  
  exit;
}
elsif ($ARGV[0] eq "Check") {
  if (scalar(@ARGV) < 4) {
	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Check <EMBL file> <fasta file> <ResultName>\n\n".
	  "Similar to the correct option, but it will only report errors in an EMBL file.\n\n";
	exit ;
	
  }
  
  my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  my $resultName = shift;	
  correctEMBL($embl,"tmp.BBA.embl",$fasta);
  
  startAnnotationCheck( "$embl.tmp.BBA.embl",$fasta,$resultName);
  exit;
}
elsif ($ARGV[0] eq "Embl2Fasta") {
  if (scalar(@ARGV) < 3) {
	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Embl2Fasta <EMBL dir> <fasta file>\n\n".
	  "Extracts the sequence from embl files in the <EMBL directory> and saves it as a <fasta file>.\n\n";
  exit 1;
  }
  
  
  my $what =shift;
  my $embl=shift;
  my $fasta = shift;
  Embl2Fasta($embl,$fasta);
  exit;
  
}


elsif ($ARGV[0] eq "Transfer") {
  
  if (@ARGV< 5) {
 	print "\n\nusage:  \$RATT_HOME/main.ratt.pl Transfer <embl Directory> <mummer SNP file> <mummer coord file> <ResultName>\n\n".
	  "This functionality uses the mummer output to map the annotation from embl files, which are in the <embl Directory>, to the query. It generates all the new annotation files (ResultName.replicon.embl), as well as files describing which annotations remain untransferred (Replicon_reference.NOTtransfered.embl).\n\n";
	
	
 	exit(1);
  }
  
  
  my $what         = shift;
  my $emblDir     = shift;
  my $mummerSNP   = shift;
  my $mummerCoords= shift;
  my $resultName  = shift;
  
  my $dbg=49;
  

  ## main hash: %ref_shift{Ref_contig}[pos] [0] query_contig
  #										  [1] position
  #										  [2] strand
  my $ref_shift;
  
  
  #load the position of the 
  my $ref_cdsPosition=loadEmbl($emblDir);
  
  
  my $ref_snp    = getFile($mummerSNP);
  my $ref_coords = getFile($mummerCoords);
  
  # clean the space of the annotation position
  
  
  # transfer the annotation the annotation
  opendir (DIR, $emblDir) or die "Problem to open opendir $emblDir: $!\n";
  
  ### will hold the new annotation: $h{queryname}.=annotation as embl
  my $ref_results;
  my $ref_Counting= {'Partial'         => 0,
					 'ExonNotTransfered' => 0,
					 'Split'         => 0,
					 'NotTransfered' => 0,
					 'Transfered'    => 0,
					 'CDS' => 0,
					 'CDSTransfered' => 0,
					 'CDSNotExons'   => 0,
					 'CDSPartial'   => 0
					};
  
  map {
	if (/(\S+)\.embl$/){
	  my $refName=$1;
	  my $ref_shift;
	  print "working on $refName\n";
	  
	  # fill the shift hash with the coords 
	  $ref_shift = loadCoords($ref_coords,$ref_shift,$refName);
	  	  #print Dumper $ref_shift;
	  # tune the shift hash with the snp file
	  my $ref_shift2 = loadSNP($ref_snp,$ref_shift,$ref_cdsPosition,$refName);
	  
	  ($ref_results,$ref_Counting)=adaptAnnotationEMBL($emblDir,$_,$ref_shift2,$ref_results,$ref_Counting);
	  
	  ### cleaning step
	  if (defined(@{$$ref_shift{$refName}})) {
		foreach (0..(scalar(@{$$ref_shift{$refName}}))-1) {
		  undef(@{ $$ref_shift{$refName}[$_]});
		}
	  }
	  
	  undef %$ref_shift;
	  
	
	}
} readdir(DIR);




  
  ### output results
  print "Overview of transfere of annotation elements:\n$$ref_Counting{Elements}\telements found.\n";
  print "$$ref_Counting{Transfered}\tElements were transfered.\n";
  print "$$ref_Counting{Partial}\tElements could be transfered partially.\n";
  print "$$ref_Counting{Split}\tElements split.\n";
  print "$$ref_Counting{ExonNotTransfered}\tParts of elements (i.e.exons tRNA) not transferred.\n";
  print "$$ref_Counting{NotTransfered}\tElements couldn't be transferred.\n";
  
  print "\nCDS:\n$$ref_Counting{CDS}\tGene models to transfer.\n$$ref_Counting{CDSTransfered}\tGene models transferred correctly.\n";
  print "$$ref_Counting{CDSPartial}\tGene models partially transferred.\n";
  print "$$ref_Counting{CDSNotExons}\tExons not transferred from partial CDS matches.\n";
  print ($$ref_Counting{CDS}-$$ref_Counting{CDSTransfered}-$$ref_Counting{CDSPartial});
  print "\tGene models not transferred.\n\n";
  
  
  
  ### then just save it
saveAnnotation($resultName,$ref_results);

}
else {
  print "Sorry, wrong option.\n";
  print "Tranfer / Correct / Check / EMBLFormatCheck / Mutate / Split / Difference / Embl2Fasta / doEMBL\ncan be used.\n\n";
}




###########################################
### subs

sub getFile{
  my $file = shift;
  open F, $file or die "Problem to open $file: $!\n";
  my @ar=<F>;
  close(F);
  return \@ar;
}

############################################
### doEMBL
############################################
sub doEMBL{
  my $resultName = shift;
  my $emblPart   = shift;
  my $fastapart  = shift;

  my ($seq,$length)=fasta2EMBLfasta($fastapart);
  
  my $res="ID                   $resultName ; ; ; ; ; $length BP.\n".
          "FH   Key             Location/Qualifiers\n".
		  "FH                   \n";

  if (-f $emblPart) {
	### if a file has no embl
	open F, $emblPart or die "Couldn't open embl file ($emblPart) in doEMBL: $!\n";
	my @ar=<F>;
	close(F);
	$res.=join ('',@ar);
  }

  
  ### mal hack: if contig is node, it will have the contigs tag in
  if ($fastapart =~ /NODE_/) {
	$res.="FT   contig          1..$length\n";
	my $col=(3+int(rand(10)));
	$res.="FT                   /note=\"Contig: $resultName.\"\n";
	$res.="FT                   /colour=$col\n";
	  
  }
  $res.=$seq."//\n";

  open F, "> $resultName.embl" or die "Couldn't write $resultName.embl in doEMBL: $! \n";
  print F $res;
  close(F);
}
########################
### fasta2EMBLfasta
########################
sub fasta2EMBLfasta{
  my $fastapart  = shift;
  open F, $fastapart or die "Couldn't open embl file ($fastapart) in doEMBL: $!\n";
  
  ### we assume single fasta. If it is multifasta, it is going to be flattend
  $_=<F>;
  my $seq;
  
  while (<F>) {
	chomp;
	
	  if (! (/^>/)) {
		$seq.=lc($_);
	  }
	}

  my @ar=split(//,$seq);
  
  my $length = length($seq);

  ### set parameter
  my $line=60;
  my $block=10;
  my $res;
  my %countB=('a' => 0,
			  't' => 0,
			  'g' => 0,
			  'c' => 0,
			  'o' => 0
			 );
  
  my $count=0;
  my $lastline='';
  foreach (@ar) {
	if (($count%$line)==0) {
	  $lastline.="     $_";
	}
	elsif (($count%$block)==0) {
	  $lastline.=" $_";
	}
	elsif ((($count+1)%$line)==0) {
	  $lastline.="$_";
	  my $l=(81-length($lastline));
	  $res.=sprintf("%-0s %*d\n",$lastline,$l,($count+1));
	 # $res.=$lastline;
	  $lastline=''
	}
	else {
	   $lastline.="$_"
	}

	### count bases
	if ($_ eq 'a') {
	  $countB{a}++
	}
	elsif ($_ eq 't') {
	  $countB{t}++
	}
	elsif ($_ eq 'g') {
	  $countB{g}++
	}
	elsif ($_ eq 'c') {
	  $countB{c}++
	}
	else {
	  $countB{o}++
	}
	
	$count++;
	
  }
 
  if ($lastline ne '') {
	my $l=(81-length($lastline));
	$res.=sprintf("%-0s %*d\n",$lastline,$l,($count+1));
  }
  
  $res = "SQ   Sequence $count BP; $countB{a} A; $countB{c} C; $countB{g} G; $countB{t} T; $countB{o} other;\n".$res;
  

  return ($res,$length);
  
}

############################################
### loadEMBL
############################################
sub loadEmbl{
  my ($emblDir) = @_;
  
  opendir (DIR, $emblDir) or die "Problem to open opendir $emblDir: $!\n";

### will hold the new annotation: $h{queryname}.=annotation as embl
my $ref_annotation;

map {
	if (/embl$/ or /embl.gz$/){
		my $embl=$_;

		my ($chr);
		
		if (/embl.gz$/) {
		  open(F, " gunzip -c $emblDir."/".$embl | " ) or die "Problems open embl $embl: $!\n";
		  ($chr) = $embl =~ /\/{0,1}(\S+)\.embl.gz$/;
		  
		}
		else
		  {
			open(F, $emblDir."/".$embl) or die "Problems open embl $embl: $!\n";
			($chr) = $embl =~ /\/{0,1}(\S+)\.embl$/;
		  }
		
		
		while (<F>){
 			if (/^FT.*CDS.*\d+/){
 				if ($_ =~ /^\W*\(*(\d+)\.\.(\d+)\)*$/) {
					$ref_annotation=
						Mask_CDS($chr,$ref_annotation,$1,$2);
  				}
	  			elsif (/\.\./){
					my @a=split(/,/);
					foreach (@a) {
					  if (/(\d+)\.\.(\d+)/){
						$ref_annotation=
							Mask_CDS($chr,$ref_annotation,$1,$2);
					  }
					  
					}
  				}	
 			}
 		}
	}
	} readdir(DIR);
	return ($ref_annotation)
}

############################################
### Mask_CDS
############################################
sub Mask_CDS{
  my ($chr,$ref,$f,$l,$out) = @_;

  for ($f..$l){
	$$ref{$chr}{$_}=1;
  }
  
  return $ref;
}

############################################
### adaptAnnotationEMBL
############################################
sub adaptAnnotationEMBL{
  my ($DIR,$embl,$ref_shift,$ref_results,$ref_Counting) = @_;
  
  open(F, $DIR."/".$embl) or die "Problems open embl $embl: $!\n";
  
  my $res;
  my $resDeleted='';
  my @ar;
  my $in=1;
  my $maybeWrong=0;

	### get the name of the contig/supercontig/chromosome
  my ($chr) = $embl =~ /\/{0,1}(\S+)\.embl$/;

  ## tag to not transfer reference sequence specific tags.
  my $inTranslation=0;
  
  
  #  print Dumper %core;
  my $OKCore=1;
  my $ref_queryTarget;
  my $transfer=0;
  
  while (<F>) {
	
	# UTR must be saved
	s/3\'UTR/3TUTR/g;
	s/5\'UTR/5TUTR/g;
	if (/FT   \S+/) {
	  s/<//g;
	  s/>//g;
	}

	my $line=$_;
	# check if entry is over more than one line
	while ($line =~ /^FT   \S+\s{2,}.*\d+,$/) {
	  $_=<F>;
	  chomp($line);
	  	  
	  /^FT   \s{2,}(.*)$/;
	  $line.=$1;
	}

	if ($line =~ /^>/) {
	  last;
	  
	}
	elsif ($line =~ /^FT   \S+\s{2,}\D+(\d+)\..*\.(\d+)/ ||
		$line =~ /^FT   \S+\s{2,}\D+\d+,(\d+)\..*\.(\d+)/ ||
		$line =~ /^FT   \S+\s{2,}\D+(\d+)/
	   ) {
	  ### This is necessary to not mapped things, which are not covered

	  my $posA=$1;
	  my $posE=$2;
	  if (!defined($posE)){
	  	$posE=$posA	
	  }
	  ### check if CDS
	  if ($line =~ /FT   CDS/) {
		$$ref_Counting{CDS}++
	  }
	  
	  
		  chomp;
	  ($ref_results,$ref_queryTarget,$ref_Counting,$transfer)=doTransfer($ref_shift,$ref_results,$chr,$posA,$line,$ref_Counting);
	  
	  $$ref_Counting{Elements}++;
	  ### case 1, all ok
	  if (defined($$ref_shift{$chr}[$posA][0]) &&
		  defined($$ref_shift{$chr}[$posE][0]) &&
		  ($$ref_shift{$chr}[$posE][0] eq $$ref_shift{$chr}[$posA][0])  # transfer to same query
		 )	
		{

	 	}
	  ### case 2, defined, but gene model is split between two query
	  elsif (defined($$ref_shift{$chr}[$posA][0]) &&
			 defined($$ref_shift{$chr}[$posE][0])
			)
		 {
		   $$ref_Counting{"Split"}++;
		   
#		   print "$chr  position $posA $posE \n $$ref_shift{$chr}[$posA][0]  // ($$ref_shift{$chr}[$posE][0] \n";
		   
		 }  
	  elsif (defined($$ref_shift{$chr}[$posA][0]) ||
			 defined($$ref_shift{$chr}[$posE][0])
			)
		{
		  #		   $$ref_Counting{"Partial"}++
		  
		}
	  
	}
	elsif (/^SQ/) {
	  last;
	} 
	### transfer it to the new 
	elsif ($transfer==1){
	  my $count=0;
	  
	  foreach my $queryTarget (keys %$ref_queryTarget) {
		$count++;
		### check for the translation tag:
		if ($inTranslation && /\"/) {
		  ### do not transfer , but it is the end
		  $inTranslation=0
		}
		elsif ($inTranslation) {
		  ### do not transfer
		}
		elsif (/\/translation=/) {
		  $inTranslation=1;
	
		   
		  ### just to double check if not smalle
		  if (/\".*\"/) {
			$inTranslation=0
		  }
		}
		else {
		  
		  if ($count>1) {
			if (/\/locus_tag/ || /_id=\"/) {
			  s/\"$/.$count\"/g;
			}
			$$ref_results[0]{$queryTarget} .= $_;
		  }
		  else {
			$$ref_results[0]{$queryTarget} .= $_;
		  }
		}
		# end else inTransation
	  }
	 }
	elsif ($transfer==0){
	  $$ref_results[1]{$chr} .= $_;
	}
	# this is the case when the annotation could be just mapped partially.
	elsif ($transfer==3){
	  $$ref_results[1]{$chr} .= $_;
	  my $count=0;
	  foreach my $queryTarget (keys %$ref_queryTarget) {
		$count++;
		  if ($count>1) {
		  if (/\/locus_tag/ || /_id=\"/) {
			s/\"$/.$count\"/g;
		  }
		  $$ref_results[0]{$queryTarget} .= $_;
		}
		else {
		  $$ref_results[0]{$queryTarget} .= $_;
		}
	
	  }
	}
	
  } # end while <F>
  
  close(F);
  return ($ref_results,$ref_Counting);
}

############################################
### doTransfer
############################################
sub doTransfer{
	my $ref_shift=shift;
	my $ref_resultsLocal =shift;
	my $chr =shift;
	my $pos =shift;
	my $line=shift;
    my $ref_Counting = shift;

	my $RENAME=47;
	my $RENAME2=249;
	
	my $chrqry=0;
	# the zero will hold the no puttable
	my %ResultLine;

	### put complement to it, or get rid of it,
	### also will need to reorder the numbers
	my $wasComplement=0;
	if ($line =~ /complement/){
	  $wasComplement=1;
	  $line =~ s/complement\(//g;
	  $line =~ s/\)$//g;
	}
	$line =~ s/join\(//g;
    $line =~ s/\)$//g;

	
	### here to look for missed exons
	
	# 	$line =~ s/(\d+)/($$ref_shift{$chr}[$1][1])/ge;
	my (@parts) = split(/\s+/,$line); 
	my @ar=split(/,/,$parts[2]);
	
	my $mappedOnce=0;
	my $exonMissed=0;
	my $oldQuery;
	my $partialCount=0;
	
	for (my $i=0;$i < scalar(@ar);$i++) {
	  my $new = $ar[$i];
	  
	  #single base exon
	  if (! ($ar[$i] =~ /\.\./) ) {
		$ar[$i] =~ /(\d+)/;
		my $pos=$1;
#		print " single exon $ar[$i] pos is $pos \n";
#		print Dumper $$ref_shift{$chr};
		
		if (defined($$ref_shift{$chr}[$pos][0])) {
		  $ar[$i] =~ s/(\d+)/($$ref_shift{$chr}[$1][1])/ge;
		  $mappedOnce++;

		  $oldQuery=$$ref_shift{$chr}[$pos][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $pos;
		  
		}		
		else {
		  $ResultLine{0}[0] .= "$ar[$i],";
		  $exonMissed++;
		  
		}
	  }
	  else {
		$ar[$i] =~ /(\d+)\.\.(\d+)/;
		my $posA=$1;
		my $posE=$2;
		my $geneLength=int ($posE-$posA);
		my $half=int ($geneLength/2);
		
		if (defined($$ref_shift{$chr}[$posA][0]) &&
			defined($$ref_shift{$chr}[$posE][0]) &&
			$$ref_shift{$chr}[$posE][0] eq $$ref_shift{$chr}[$posA][0] &&
			abs($$ref_shift{$chr}[$posE][1] - $$ref_shift{$chr}[$posA][1])< (2*($posE-$posA))
		   ) {
		  $ar[$i] =~ s/(\d+)/($$ref_shift{$chr}[$1][1])/ge;
		  $mappedOnce++;

		  
		  $oldQuery=$$ref_shift{$chr}[$posA][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $pos;
		  
		}
		elsif (defined($$ref_shift{$chr}[$posA][0]) &&
			   (defined($$ref_shift{$chr}[($posA+299)][0])) &&
			   $$ref_shift{$chr}[($posA+299)][0] eq $$ref_shift{$chr}[$posA][0]
			   && abs($$ref_shift{$chr}[($posA+299)][1] - $$ref_shift{$chr}[$posA][1])< 20000
			  ) {
		  $ar[$i] = $$ref_shift{$chr}[$posA][1]."..".($$ref_shift{$chr}[($posA)][1]+$geneLength);
		  $mappedOnce++;
		  $partialCount++;
		  
		  $oldQuery=$$ref_shift{$chr}[$posA][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posA;
		  
		}
		### left part ok
		elsif (defined($$ref_shift{$chr}[$posA][0]) &&
			   (defined($$ref_shift{$chr}[($posA+74)][0])) &&
			   $$ref_shift{$chr}[($posA+74)][0] eq $$ref_shift{$chr}[$posA][0]
			   && abs($$ref_shift{$chr}[($posA+74)][1] - $$ref_shift{$chr}[$posA][1])< 20000
			  ) {
		  $ar[$i] = $$ref_shift{$chr}[$posA][1]."..".($$ref_shift{$chr}[($posA)][1]+$geneLength);
		  $mappedOnce++;
		  $partialCount++;
		  
		  $oldQuery=$$ref_shift{$chr}[$posA][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posA;
		  
		}
		elsif (defined($$ref_shift{$chr}[$posA][0]) &&
			   (defined($$ref_shift{$chr}[($posA+14)][0])) &&
			   $$ref_shift{$chr}[($posA+14)][0] eq $$ref_shift{$chr}[$posA][0]
			   && abs($$ref_shift{$chr}[($posA+14)][1] - $$ref_shift{$chr}[$posA][1])< 20000
			  ) {
		  $ar[$i] = $$ref_shift{$chr}[$posA][1]."..".($$ref_shift{$chr}[($posA)][1]+$geneLength);
		  $mappedOnce++;
		  $partialCount++;
		  
		  $oldQuery=$$ref_shift{$chr}[$posA][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posA;
		  
		}

		### 3' ok
		elsif (defined($$ref_shift{$chr}[$posE][0]) &&
			defined($$ref_shift{$chr}[($posE-299)][0]) &&
			$$ref_shift{$chr}[($posE-299)][0] eq $$ref_shift{$chr}[$posE][0] &&
			abs($$ref_shift{$chr}[($posE-299)][1] - $$ref_shift{$chr}[$posE][1])< 20000
		   ) {
		  $ar[$i] = ($$ref_shift{$chr}[($posE)][1]-$geneLength)."..".$$ref_shift{$chr}[($posE)][1];
		  $mappedOnce++;
		  $partialCount++;

		  
		  $oldQuery=$$ref_shift{$chr}[$posE][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posE;
		  
		}
		elsif (defined($$ref_shift{$chr}[$posE][0]) &&
			defined($$ref_shift{$chr}[($posE-74)][0]) &&
			$$ref_shift{$chr}[($posE-74)][0] eq $$ref_shift{$chr}[$posE][0] &&
			abs($$ref_shift{$chr}[($posE-74)][1] - $$ref_shift{$chr}[$posE][1])< 20000
		   ) {
		  $ar[$i] = ($$ref_shift{$chr}[($posE)][1]-$geneLength)."..".$$ref_shift{$chr}[($posE)][1];
		  $mappedOnce++;
		  $partialCount++;

		  
		  $oldQuery=$$ref_shift{$chr}[$posE][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posE;
		  
		}
		elsif (defined($$ref_shift{$chr}[$posE][0]) &&
		   defined($$ref_shift{$chr}[($posE-14)][0]) &&
			$$ref_shift{$chr}[($posE-14)][0] eq $$ref_shift{$chr}[$posE][0] &&
			abs($$ref_shift{$chr}[($posE-14)][1] - $$ref_shift{$chr}[$posE][1])< 20000
		   ) {
		  $ar[$i] = ($$ref_shift{$chr}[($posE)][1]-$geneLength)."..".$$ref_shift{$chr}[($posE)][1];
		  $mappedOnce++;
		  $partialCount++;

		  
		  $oldQuery=$$ref_shift{$chr}[$posE][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $posE;
		  
		}
		elsif (defined($$ref_shift{$chr}[($posA+$half)][0]) 
		
		   ) {
		  my $start=($$ref_shift{$chr}[($posA+$half)][1]-$half);
		  if ($start<1) {
			$start=1;
		  }
		  my $end=($$ref_shift{$chr}[($posA+$half)][1]+$half);
		  $ar[$i] = $start."..".$end;
		  print
		   "$start $end $half ".$$ref_shift{$chr}[($posA+$half)][1]." \n";
		  
		  $mappedOnce++;
		  $partialCount++;

		  
		  $oldQuery=$$ref_shift{$chr}[($posA+$half)][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = $end;
		  
		}

		
		elsif (defined($$ref_shift{$chr}[$posA+$RENAME][0]) &&
			   defined($$ref_shift{$chr}[($posE-$RENAME)][0]) &&
			   $$ref_shift{$chr}[$posE-$RENAME][0] eq $$ref_shift{$chr}[$posA+$RENAME][0] 
			   #		abs($$ref_shift{$chr}[$posE][1] - $$ref_shift{$chr}[$posA][1])< (2*($posE-$posA))
			   
			   #	   $$ref_shift{$chr}[($posE-$RENAME)][0] eq $$ref_shift{$chr}[$posE][0] &&
			   #	abs($$ref_shift{$chr}[($posE-$RENAME)][1] - $$ref_shift{$chr}[$posE][1])< 20000
			  ) {
#		  print " TST:\n";
#		  print Dumper $$ref_shift{$chr}[($posA+$RENAME)];
#		  print Dumper $$ref_shift{$chr}[($posE-$RENAME)];
		  $ar[$i] = $$ref_shift{$chr}[($posA+$RENAME)][1]."..".$$ref_shift{$chr}[($posE-$RENAME)][1];
		  $mappedOnce++;
		  $partialCount++;

		  
		  $oldQuery=$$ref_shift{$chr}[$posA+$RENAME][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = ($posE-$RENAME);
		  
		}
		elsif (defined($$ref_shift{$chr}[$posA+$RENAME2][0]) &&
			   defined($$ref_shift{$chr}[($posE-$RENAME2)][0]) &&
			   $$ref_shift{$chr}[$posE-$RENAME2][0] eq $$ref_shift{$chr}[$posA+$RENAME2][0] 
			   #		abs($$ref_shift{$chr}[$posE][1] - $$ref_shift{$chr}[$posA][1])< (2*($posE-$posA))
			   
			   #	   $$ref_shift{$chr}[($posE-$RENAME2)][0] eq $$ref_shift{$chr}[$posE][0] &&
			   #	abs($$ref_shift{$chr}[($posE-$RENAME2)][1] - $$ref_shift{$chr}[$posE][1])< 20000
			  ) {
#		  print " TST:\n";
#		  print Dumper $$ref_shift{$chr}[($posA+$RENAME2)];
#		  print Dumper $$ref_shift{$chr}[($posE-$RENAME2)];
		  $ar[$i] = $$ref_shift{$chr}[($posA+$RENAME2)][1]."..".$$ref_shift{$chr}[($posE-$RENAME2)][1];
		  $mappedOnce++;
		  $partialCount++;

		  
		  $oldQuery=$$ref_shift{$chr}[$posA+$RENAME2][0];
		  $ResultLine{$oldQuery."::".$chr}[0] .= "$ar[$i],";
		  $ResultLine{$oldQuery."::".$chr}[1] = ($posE-$RENAME2);
		  
		}
		
		else {
		  $ResultLine{0}[0] .= "$ar[$i],";
		  $exonMissed++;
		}
	  }
	  
	  
	}
#	print Dumper @ar;
	## default do not transfer
	my $transfer=0;
	#### check the amount
	if ($mappedOnce ==0) {
	  $$ref_Counting{NotTransfered}++
	}
	if (($mappedOnce > 0 && $exonMissed >0)
#		||
#		$partialCount>0
	   ){
	  $$ref_Counting{Partial}++;
	  ### This means, put the annotation to the BB and LB, as it is partial
	  if ($line =~ /FT   CDS/) {
		$$ref_Counting{CDSPartial}++;
		$$ref_Counting{CDSNotExons}+=$exonMissed;
	  }
	  $transfer=3;
	}


	
	if ($exonMissed==0) {
	  $$ref_Counting{Transfered}++;
	  if ($line =~ /FT   CDS/) {
		$$ref_Counting{CDSTransfered}++;
	  }
	  
	  ### Modell fully mapped, put it just to the LB
	  $transfer=1;
	}
	else {
	  $$ref_Counting{ExonNotTransfered}+=$exonMissed;
	  $ResultLine{0}[0]=~ s/,$//g;
	  if ($wasComplement==1){
		$ResultLine{0}[0]="complement(join($ResultLine{0}[0]))";
	  } 
	  else {
	  	$ResultLine{0}[0]="join($ResultLine{0}[0])";
	  }
	  $$ref_resultsLocal[1]{$chr}.=sprintf("%-4s %-15s %s\n",$parts[0],$parts[1],$ResultLine{0}[0]);
	  $ResultLine{0}=undef;
	  undef $ResultLine{0};
	}
	
	my %targetChr;
	

	
	foreach my $trans (keys %ResultLine) {
	
		
	  if ($trans ne '0'){
		
		my ($chrqryLocal,$chrpart) = $trans =~ /^(\S+)::(\S+)$/;
		$targetChr{$chrqryLocal}=1;
		
		my $pos =$ResultLine{$trans}[1];
		if (!defined($chrqry)){
		  print $trans."\n";
		  print Dumper %ResultLine;
		  exit;	
		}	  
		
		
		$ResultLine{$trans}[0]=~ s/,$//g;
	  	my @ar=split(/\,/,$ResultLine{$trans}[0]);
	  	my $amountJoined=(scalar(@ar)-1);
	  	if ($amountJoined > 1){
			$ResultLine{$trans}[0]="join($ResultLine{$trans}[0])";
		  }
		  if ($wasComplement==1){
			$ResultLine{$trans}[0]="complement($ResultLine{$trans}[0])";
	  	} 
	  
	### check if the entry must be inversed
 	
	  $chrqry=$chrqryLocal;

	  if (defined($$ref_shift{$chrpart}[$pos][2]) && $$ref_shift{$chrpart}[$pos][2] == -1){
		$line=$ResultLine{$trans}[0];
		
		if ($line =~ /complement/){
		  $wasComplement=1;
		  $line =~ s/complement\(//g;
		  $line =~ s/\)$//g;
		}
		
		my $hadExons=0;
		if ($line =~ /join/){
		  $hadExons=1;
		  $line =~ s/join\(//g;
		}		
		$line =~ s/\)$//g;
		my @numbers = sort {$a <=> $b} split(/,|\.{2}/,$line);
		my $new;
		my $count=0;
		foreach (@numbers){
		  if ($count>0 && ($count%2==1)){
			$new.="..";
		  } elsif ($count>0 && ($count%2==0)){
			$new.=",";
		  }
		  $new.=$_;
		  $count++;
		}
		if ($hadExons){
		  $new="join($new)"	
		}
		if ($wasComplement==0){
		  $new="complement($new)"
		}
		$ResultLine{$trans}[0]=$new;
	  }

	  $$ref_resultsLocal[0]{$chrqry}.=sprintf("%-4s %-15s %s\n",$parts[0],$parts[1],$ResultLine{$trans}[0]);
	}
	}
 	return ($ref_resultsLocal,\%targetChr,$ref_Counting,$transfer);
}
############################################
### saveAnnotation
############################################
sub saveAnnotation{
  my ($name,$ref_h) = @_;


  ### map the mapping Stuff
  foreach my $query (sort keys %{$$ref_h[0]}){
	my $nameQry=$query;
	$nameQry =~ s/\|/_/g;
	open (F,"> $name.$nameQry.embl") or die "Couldn't open save file $name: $!\n";
	
	# UTR must be saved
	$$ref_h[0]{$query} =~ s/3TUTR/3\'UTR/g;
	$$ref_h[0]{$query} =~ s/5TUTR/5\'UTR/g;
	print F $$ref_h[0]{$query};
	close(F);
  }
  
  # UTR must be saved
  foreach my $ref (sort keys %{$$ref_h[1]}){
	my $nameRef=$ref;
	$nameRef =~ s/\|/_/g;
	
 	open (F,"> $name.$nameRef.NOTTransfered.embl") or die "Couldn't open save file $name: $!\n";
	
	$$ref_h[1]{$ref} =~ s/3TUTR/3\'UTR/g;
	$$ref_h[1]{$ref} =~ s/5TUTR/5\'UTR/g;
	print F $$ref_h[1]{$ref};
	close(F);
  }

}
############################################
### saveGFF
############################################

sub saveGFF{
  my ($name,$path,$ref_h) = @_;

  if (! -d "$path") {
	!system("mkdir $path") or warn "Couldn't create directory $path.\n";
  }

  foreach my $chr (sort keys %$ref_h){
	my $chrName=$chr;
	$chrName =~ s/\|/_/g;
	
	open (F,"> $path/$name.$chrName.Mutations.gff") or die "Couldn't create file $name: $!\n";

  # UTR must be saved
  print F $$ref_h{$chr};
  close(F);
  }
}

############################################
### loadSNP
############################################
sub loadSNP{
  my $ref_File        = shift;
  my $ref_shift       = shift;
  my $ref_cdsPosition = shift;
  my $refName  = shift;

  my $lastQry="";
  
  ## walk through the list. The last 
  for my $pos (0..(scalar(@$ref_File)-2)) {

	# get the positions of the snp'indels of before and last
	my @previous;
	if ($pos >0) {
	  @previous = split(/\s+/,$$ref_File[($pos-1)]);
	}
	my ($refPos,$refWhat,$queryWhat,$queryPos,$dum1,$dum2,$dum3,$dum4,$refStrand,$queryStrand,$reference,$query) = split(/\s+/,$$ref_File[$pos]);
	
	if (!(defined($refName)) || $refName eq $reference) {
	  
	  my @next = split(/\s+/,$$ref_File[($pos+1)]);
	  
	  my ($refNextPos,$queryNextPos,$refNext,$queryNext) 
		= ($next[0],$next[3],$next[10],$next[11]);
	  
	  # TODO 1: differ is next is not in the same combination	
	  # TODO 2: what to do with the last SNP?
	  # TODO 3: must be the same contig combination
	  
	  
	  if ($refNext ne $reference) {
		$ref_shift=walkToEnd($ref_shift,$reference,$refPos,$query,$queryPos,$queryStrand);
	  }
	  
	  #case 1: strand 1 and  update due to annotation
	  elsif ($query eq $queryNext && abs($refNextPos - $refPos) < 50000) {
		
		if (
			!defined($$ref_cdsPosition{$reference}{$refPos})    && # this mutation is not on gene
			defined($$ref_cdsPosition{$refNext}{$refNextPos})   # next mutation is on gene
			
		   ){
		  for  (my $posLocal=$refNextPos; $posLocal >=($refPos);$posLocal--){
			if (defined($$ref_shift{$reference}[$posLocal][0])) {
			  $$ref_shift{$reference}[$posLocal][1]=$queryNextPos;
			} 
			$queryNextPos-=$queryStrand;
		  }	
		}
		else {
		  for my $posLocal ($refPos..($refNextPos-1)){
			if (defined($$ref_shift{$reference}[$posLocal][0])) {
			  $$ref_shift{$reference}[$posLocal][1]=$queryPos;
			  
			}
			
			$queryPos+=$queryStrand
		  }	
		}
		
	  }
	  $lastQry=$query
	}
		
		
		
  }
  
  #now have a look a the last line:
  my ($refPos,$refWhat,$queryWhat,$queryPos,
	  $dum1,$dum2,$refStrand,$queryStrand,
	  $reference,$query)
	= split(/\s+/,$$ref_File[(scalar(@$ref_File)-1)]);

  $ref_shift=walkToEnd($ref_shift,$reference,$refPos,$query,$queryPos,$queryStrand);
  return $ref_shift;	
}

sub walkToEnd{
  my ($ref_shift,$reference,$refPos,$query,$queryPos,$queryStrand) = @_;
  while (0 && defined($$ref_shift{$reference}[$refPos])) {
	$$ref_shift{$reference}[$refPos][1]=$queryPos;
	$queryPos+=$queryStrand;
#	print "$refPos -3-> $queryPos\n";
	$refPos++
  }
  
  
  return $ref_shift
}


############################################
### loadCoords
############################################

sub loadCoords{
  my $ref_File = shift;
  my $ref_h    = shift;
  my $refName  = shift;
  
  
  
  ### position of blocks
 foreach (@$ref_File) {
	# 1       115285  837     116121  115285  115285  99.99   5.36    5.31    Neo_chrII       Neo_chrII

	my ($refPos1,$refPos2,$queryPos1,$queryPos2,$overlap1,$overlap2,$identity,$dum1,$dum2,$refLength,$queryLength,$reference,$query) = split(/\s+/);

	### if the alignment is inverted...

	if (!(defined($refName)) || $refName eq $reference) {
	
	  my $maxQuery=$queryPos2; ### as the alignment length might not be
	  ### the same, the querypos cannot be
	  ### bigger than the $queryPos2
	  my $minQuery=$queryPos1;
	  
	  
	  my $strand=1;
	  if ($queryPos1 > $queryPos2) {
		$strand=-1;
		$minQuery=$queryPos2;
		$maxQuery=$queryPos1;
	  }
	  #	print "$refPos1..$refPos2 - $reference - $query\n";
	  
	  for my $pos ($refPos1..$refPos2){
		if ($queryPos1< $minQuery ||
			$queryPos1 > $maxQuery
		   ) {
		  $queryPos1-=$strand;
		  
		}
		
		@{ $$ref_h{$reference}[$pos]} = ($query,$queryPos1,$strand);
		$queryPos1+=$strand;
	  }
	}
		
  }
  return ($ref_h);
}


############################################
### putMutation
############################################
## Changes every 250 bp the base to G or C
sub putMutation{
  my $file = shift;
  my $MUTATION_RATE = shift;
  
  open (F,$file) or die "Please provide a valid query sequence\n";

  if (!defined($MUTATION_RATE)) {
    $MUTATION_RATE=300;
  }
  
  my %h;
  my $name;
  while (<F>) {
	chomp;
	if (/^>(\S+)/) {
	  $name=$1;
	}
	elsif (/^\s*$/) {
	}
	else {
	  $h{$name}.=$_;
	}
  }
  close(F);

  # insert the mutation

  ## per chromosome;
  for my $chr (keys %h) {
	my @seq = split "", $h{$chr};	
	my $length=(scalar(@seq)-1);
	
	# mutate the every $MUTATION_RATE base
	
	for (my $i = $MUTATION_RATE; $i < $length ; $i+=$MUTATION_RATE){
	  if ($seq[$i] ne 'G') {
		$seq[$i]="G\n"
	  }
	  else {
		$seq[$i]="C\n";
	  }
	}
	### mutate the second and the seconlast
	if ($seq[1] ne 'G') {
	  $seq[1]='G'
	}
	else {
	  $seq[1]='C';
	}
	$length--;
	$length--;
	
	if ($seq[$length] ne 'G') {
	  $seq[$length]='G'
	}
	else {
	  $seq[$length]='C';
	}
	
	# put the sequence back together
	$h{$chr}= join("",@seq)
  }

  
  ### write the result
  open (F, ">$file.mutated") or die "Problems to write the file $file.mutated\n";
   for my $chr (keys %h) {
	 print F ">$chr\n$h{$chr}\n";
   }
  close(F)
}

###########################################
### Functions for gff files
############################################

sub getMutations{
  my $fileNameSNP     = shift;
  my $fileNameCoords  = shift;
  my $resultName      = shift;
  
  open (F, $fileNameSNP) or die "Problem to open SNP File: $fileNameSNP \n";
  
  my @File=<F>;
  close(F);

  my %BB;
  my %LB;

  my (%sizeRef,%sizeQuery,%coveredRef,%coveredQuery,%noCovRef,%noCovQry);
  my %h_sizeRef;
  my %h_sizeQuery;
  
  foreach (@File) {
	my ($refPos,$refWhat,$queryWhat,$queryPos,$dum1,$dum2,$refStrand,$queryStrand,$reference,$query) = split(/\s+/);

	if ($refWhat eq ".") {
	  $BB{$reference} .="unknown\tBBA\tIns\t$refPos\t$refPos\t0\t+\t.\tnote=\"Insertion+in+query:+$queryWhat++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="unknown\tBBA\tDel\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"Deletion+in+reference++++Strand+of+query+is+$queryStrand\"\n";
	}
	elsif($queryWhat eq "."){
	  $BB{$reference} .="unknown\tBBA\tDel\t$refPos\t$refPos\t0\t+\t.\tnote=\"Deletion+in+query++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="unknown\tBBA\tIns\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"Insertion+in+reference:+$refWhat++++Strand+of+query+is+$queryStrand\"\n";
	}
	
	else {
	  $BB{$reference} .="unknown\tBBA\tSNP\t$refPos\t$refPos\t0\t+\t.\tnote=\"SNP+in+query:+$queryWhat++++Strand+of+query+is+$queryStrand\"\n";
	  $LB{$query} .="unknown\tBBA\tSNP\t$queryPos\t$queryPos\t0\t+\t.\tnote=\"SNP+in+reference:+$refWhat++++Strand+of+query+is+$queryStrand\"\n";
	}
  }

  ## from the coords files we want the regions that are not
  ## covered\n"; 
  
  open (F, $fileNameCoords) or die "Problem to open SNP File: $fileNameCoords \n";
  
  @File=<F>;
  close(F);

  my %coverBB;
  my %coverLB;

  ### get an field, where there is coverage
  foreach (@File) {
	my ($refPos1,$refPos2,$queryPos1,$queryPos2,$overlap1,$overlap2,$identity,$refLength,$queryLength,$dum1,$dum2,$reference,$query) = split(/\s+/);
	$sizeRef{$reference}=$refLength;
	$sizeQuery{$query}=$queryLength;

	$coveredRef{$reference}+=($refPos2-$refPos1+1);
	
	# check on the reference
	$coverBB{$reference}[$refLength]=undef;
	for ($refPos1..$refPos2) {
	  $coverBB{$reference}[$_]=1;
	}

	# check for the query
	$coverLB{$query}[$queryLength]=undef;
	if ($queryPos2<$queryPos1) {  my $tmp=$queryPos2;$queryPos1=$queryPos2;$queryPos2=$tmp	}
	$coveredQuery{$query}+=(abs($queryPos2-$queryPos1)+1);
	for ($queryPos1..$queryPos2) {
	  $coverLB{$query}[$_]=1;
	}
  }
  ### now parse, were there is no coverage
  for my $chr (keys %coverBB) {
	my $start=1;
	my $ok=1;

	for my $pos (1..(scalar(@{$coverBB{$chr}}))) {
	  if ($ok==1 &&
		  !defined($coverBB{$chr}[$pos])) {
		$ok=0;
		$start=$pos
	  }
	  if ($ok==0 &&
		  defined($coverBB{$chr}[$pos])) {
		$ok=1;
		$BB{$chr} .="unknown\tBBA\tSynteny\t$start\t".($pos-1)."\t0\t+\t.\tnote=\"No synteny with query. Possible insert or too divergent\"\n";
		$noCovRef{$chr}+=($pos-1-$start);
		
	  }
	}
	if ($ok==0 && $start <scalar(@{$coverBB{$chr}}) ) {
	  $BB{$chr} .="unknown\tBBA\tSynteny\t$start\t".(scalar(@{$coverBB{$chr}})-1)."\t0\t+\t.\tnote=\"No synteny with query. Possible insert or too divergent\"\n";
	  $noCovRef{$chr}+=($sizeRef{$chr}-1-$start);
	}
  }
  ### now parse, were there is no coverage
  for my $chr (keys %coverLB) {
	my $start=1;
	my $ok=1;

	for my $pos (1..(scalar(@{$coverLB{$chr}})-1)) {
	  if ($ok==1 &&
		  !defined($coverLB{$chr}[$pos])) {
		$ok=0;
		$start=$pos
	  }
	  if ($ok==0 &&
		  defined($coverLB{$chr}[$pos])) {
		$ok=1;
		$LB{$chr} .="unknown\tBBA\tSynteny\t$start\t".($pos-1)."\t0\t+\t.\tnote=\"No synteny with reference. Possible insert or too divergent\"\n";
		$noCovQry{$chr}+=($pos-1-$start);
	  }
	}
	if ($ok==0 && $start < scalar(@{$coverLB{$chr}})) {
	  $LB{$chr} .="unknown\tBBA\tSynteny\t$start\t".(scalar(@{$coverLB{$chr}})-1)."\t0\t+\t.\tnote=\"No synteny with reference. Possible insert or too divergent\"\n";
	  $noCovQry{$chr}+=($sizeQuery{$chr}-1-$start);
	}
  }
  saveGFF($resultName,"Reference",\%BB);
  saveGFF($resultName,"Query",\%LB);

  my $res;
  
  foreach my $chr (sort keys %sizeRef ) {
	$res.=sprintf("Of the reference chromosome $chr\t%.2f per cent\thas no synteny with the query\n",(( $sizeRef{$chr} - $coveredRef{$chr})*100/$sizeRef{$chr}) );
  }
  foreach my $chr (sort keys %sizeQuery ) {
	$res.=sprintf("Of the query chromosome $chr\t %.2f per cent\thas no synteny with the reference\n",((( $sizeQuery{$chr}-$coveredQuery{$chr} ))*100/$sizeQuery{$chr}) );
  }
  return $res
}



###########################################
### Seperate Replicon
############################################

sub Split{
  my $FileName = shift;
  my $Amount = shift;
  
  if (!defined($Amount)) {
	$Amount = 99999999999;
  }
  open (FILE, $FileName) or die "Couldn't open file to Seqparate $FileName $!\n";
  
  my $Counter = 0;
  my $Out;
  my @NameFiles;
  while (<FILE>) {
	if ($Counter < $Amount) {
	  
	  if (/^>(\S+)/) {
		close(FILEOUT);
		my $Name = $1;#."_$Counter.fasta";
		open (FILEOUT, "> $Name") or die "Problem do $Name file  $!\n";
		push @NameFiles, $Name;
		$Counter++;
		print FILEOUT $_;
	  }
	  else {
		print FILEOUT $_;
	  }
	}
  }
  print "$Counter Sequences where generated\nDone.\n";
  return \@NameFiles;
}

###########################################
### Seperate Replicon
############################################

sub Embl2Fasta{
  my $emblDir = shift;
  my $fastaResult = shift;

  
  opendir (DIR, $emblDir) or die "Problem to open opendir $emblDir: $!\n";
  
  ### will hold the new annotation: $h{queryname}.=annotation as embl

  my $fasta;
  
  map {
	if (/embl$/){
	  open F, "$emblDir/$_ " or die "Problem open $emblDir/$_ in Embl2Fasta:$! \n";
	  my ($name)=$_ =~ /^(\S+)\.embl/;
	  
	  while (<F>) {
		if (/^SQ/) {
		  $fasta.=">$name\n";
		  while (<F>) {
			if (/\/\//) {
			  
			  last;
			}
			## get away space and number of seq
			s/\d+//g;
			s/\s+//g;
			$fasta.=$_."\n";
						
		  }
		}
	  }
	  
	}
  } readdir(DIR);

  if (!defined($fasta)) {
	die "Sorry the embl file contained no sequence. Please specify a fasta file for the reference.\n";
	
  }
  open F, "> $fastaResult" or die "Couldn't write file $fastaResult in Embl2Fasta: $!\n";
  print F $fasta;
  close(F);
  
}
sub addTranslation
{
  my($embl_file, $translationTable) = @_;
  my %seqs = ();
  
  require Bio::SeqIO;
  if ((!defined($translationTable))) {
	$translationTable=1
  }
  my %translations;
  
  my $stream = Bio::SeqIO->new(-file => $embl_file, -format => 'EMBL');
  
  while ( (my $seq = $stream->next_seq()) )
	{
	  my @features = $seq->get_SeqFeatures();
	  
	  my $current_id = '';
	  
	  foreach my $feat (@features)
		{
		  
		  next unless $feat->primary_tag eq 'CDS';
		  
		  foreach my $tag ( $feat->get_all_tags() )
			{
			  if ($tag eq 'systematic_id' || $tag eq 'locus_tag' || $tag eq 'gene')
				{
				  $current_id = join(' ',$feat->get_tag_values($tag));
				}
			  elsif ($tag eq 'transl_table') {
				$translationTable= join(' ',$feat->get_tag_values($tag));
			  }
			}
		  
		  my $cds_obj = $feat->spliced_seq;
		  my $cds_seq = $cds_obj->seq;
		  
		  my $seq_length = length $cds_seq;
		  
		  my $substr_length = $seq_length;
		  
		  #print $substr_length."\n";
		  
		  $cds_seq = substr($cds_seq, 0, ($substr_length));
		  
		  $translationTable= int $translationTable;
		  
		  my $translation = $cds_obj->translate("*",undef,undef, $translationTable);

		  my $aa = $translation->seq();
		  
		  my $aa_length = length $aa;
		  
		  $aa = substr($aa, 0, $aa_length );
		  $aa = "/translation=\"$aa";

		  my $tmp;
		  my $lineLength=59;
		  
		  for (0..(length($aa)/$lineLength)) {
			if ((length($aa)<($lineLength*($_+1)))) {
			  $tmp.="FT                   ".substr($aa,($lineLength*$_),$lineLength)."\"\n";
			}
			else {
			  $tmp.="FT                   ".substr($aa,($lineLength*$_),$lineLength)."\n";
			}
		  }
		  
		  $translations{$current_id}=$tmp;
		  
		}
	}
  
  ### put it in embl file
  open F,  $embl_file;
  my $res;


  my $current_id='';
  my $inCDS;
  
  while (<F>) {

	if (/FT   \S+/ && $inCDS && $current_id) {
	  if ($current_id) {
		$res.=$translations{$current_id};
		$current_id=''
	  }
	  $inCDS=0
	}
	
	if (/FT\s+\/gene=\"(\S+)\"/ ||
		/FT\s+\/systematic_id=\"(\S+)\"/ ||
		/FT\s+\/locus_tag=\"(\S+)\"/   ) {
	  $current_id=$1;
	  $res.=$_;
	}
	elsif (/FT   CDS/) {
	  $inCDS=1;
	  $res.=$_;
	  
	}
	elsif (/^XX/) {
	  if ($current_id) {
		$res.=$translations{$current_id};
		
		$current_id=''
	  }
	  $res.=$_;
	}
	elsif (/^SQ/) {
	  if ($current_id) {
		$res.=$translations{$current_id};
		$current_id=''
	  }
	  $res.=$_;
	}
	else {
	  $res.=$_;
	}
  }
  close(F);
  
  ### write the results;
  my ($name) = $embl_file =~ /^(\S+)\.embl$/;
  open (F, "> $name.withTranslation.embl") or die "Couldn't write file $name.withTranslation.embl : $! in addTranslation\n";
  print F $res;
  close(F);
}

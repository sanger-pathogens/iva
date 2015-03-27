#!/bin/sh
MACHINE=`uname`
PROC=`uname -p`
SCRIPT_PATH=$0
SCRIPT_PATH=`dirname $SCRIPT_PATH`
JAVA_PATH=$SCRIPT_PATH:.

MUMMER=$SCRIPT_PATH
if [ ! -e $MUMMER/nucmer ]; then
   MUMMER=`which nucmer`
   MUMMER=`dirname $MUMMER`
fi

echo "MUMMER: $MUMMER"

REF=$1
CONTIGS=$2
SCAFFOLDS=$3
NUCMER_MINID=$4

CONTIG_FILE=$(basename $CONTIGS)
SCAFFOLD_FILE=$(basename $SCAFFOLDS)

JAVADIR=/fs/wrenhomes/sergek/Utils
GENOMESIZE=`java -cp $JAVA_PATH SizeFasta $REF |awk '{SUM+=$NF; print SUM}'|tail -n 1`

echo "Contig Stats"
java -cp $JAVA_PATH GetFastaStats -o -min 200 -genomeSize $GENOMESIZE $CONTIGS 2>/dev/null
$MUMMER/nucmer --maxmatch -p $CONTIG_FILE -l 30 -banded -D 5 $REF $CONTIGS
$MUMMER/delta-filter -o 95 -i $NUCMER_MINID $CONTIG_FILE.delta > $CONTIG_FILE.fdelta
$MUMMER/dnadiff -d $CONTIG_FILE.fdelta

bash $SCRIPT_PATH/getMummerStats.sh $CONTIGS $SCRIPT_PATH
cat out.1coords |awk '{print NR" "$5}' > $CONTIG_FILE.matches.lens

echo ""
echo "Corrected Contig Stats"
java -cp $JAVA_PATH:. GetFastaStats -o -min 200 -genomeSize $GENOMESIZE $CONTIG_FILE.matches.lens 2> /dev/null

java -cp $JAVA_PATH SplitFastaByLetter $SCAFFOLDS N > tmp_scf.fasta
$MUMMER/nucmer --maxmatch -p $SCAFFOLD_FILE -l 30 -banded -D 5 $REF tmp_scf.fasta
$MUMMER/delta-filter -o 95 -i $NUCMER_MINID $SCAFFOLD_FILE.delta > $SCAFFOLD_FILE.fdelta
$MUMMER/show-coords -lrcT $SCAFFOLD_FILE.fdelta | sort -k13 -k1n -k2n > $SCAFFOLD_FILE.coords
$MUMMER/show-tiling -c -l 1 -i 0 -V 0 $SCAFFOLD_FILE.fdelta > $SCAFFOLD_FILE.tiling

echo "Scaffold Stats"
java -cp $JAVA_PATH GetFastaStats -o -min 200 -genomeSize $GENOMESIZE $SCAFFOLDS 2> /dev/null
echo "Corrected Scaffold Stats"
java -cp $JAVA_PATH getScaffoldStats $SCAFFOLDS $SCAFFOLD_FILE.tiling $GENOMESIZE $SCAFFOLD_FILE.coords 2> $SCAFFOLD_FILE.err

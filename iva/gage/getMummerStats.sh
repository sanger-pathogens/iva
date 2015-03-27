FILENAME=$1
SCRIPT_PATH=$2
JAVA_PATH=$2:.

grep TotalBases out.report | sed 's|.*/||' | awk '{print "Genome Size: "$2}'
grep TotalBases out.report | sed 's|.*/||' | awk '{print "Assembly Size: "$3}'
echo Chaff bases: $(java -cp $JAVA_PATH SizeFasta $FILENAME 2>/dev/null | awk 'BEGIN{s = 0;}{if($NF<200)s+=$NF}END{print s}')
grep UnalignedBases out.report | sed 's|.*/||' | awk '{print "Missing Reference Bases: "$2}'
grep UnalignedBases out.report | sed 's|.*/||' | awk '{print "Missing Assembly Bases: "$3}'
grep UnalignedSeqs out.report | sed 's|.*/||' | awk '{print "Missing Assembly Contigs: "$3}'
echo Duplicated Reference Bases: $(grep DUP out.qdiff | cut -f5 | $SCRIPT_PATH/colsum.pl)
echo Compressed Reference Bases: $(grep DUP out.rdiff | cut -f5 | $SCRIPT_PATH/colsum.pl)
echo Bad Trim: $(grep BRK out.qdiff | cut -f5 | $SCRIPT_PATH/colsum.pl)
grep '1-to-1' -A3 out.report | grep AvgIdentity | sed 's|.*/||' | awk '{print "Avg Idy: "$2}'
grep TotalSNPs out.report | sed 's|.*/||' | awk '{print "SNPs: "$2}'
grep TotalIndels out.report | sed 's|.*/||' | awk '{print "Indels < 5bp: "$2}'
cat out.rdiff | grep -c GAP | sed 's|.*/||;s|:| |;' | awk '{print "Indels >= 5: "$1}'
grep Inversions out.report | sed 's|.*/||' | awk '{print "Inversions: "$3}'
grep Relocations out.report | sed 's|.*/||' | awk '{print "Relocation: "$3}'
grep Translocations out.report | sed 's|.*/||' | awk '{print "Translocation: "$3}'

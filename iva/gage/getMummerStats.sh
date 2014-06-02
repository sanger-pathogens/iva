FILENAME=$1
SCRIPT_PATH=$2
JAVA_PATH=$2:.

grep TotalBases `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Genome Size: "$2}'
grep TotalBases `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Assembly Size: "$3}'
echo -n "Chaff bases: "; for i in `/bin/ls $FILENAME`; do java -cp $JAVA_PATH SizeFasta $i 2>/dev/null | awk 'BEGIN{s = 0;}{if($NF<200)s+=$NF}END{print s}'; done
grep UnalignedBases `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Missing Reference Bases: "$2}'
grep UnalignedBases `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Missing Assembly Bases: "$3}'
grep UnalignedSeqs `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Missing Assembly Contigs: "$3}'
echo -n "Duplicated Reference Bases: "; for i in `/bin/ls out.qdiff`; do grep DUP $i | cut -f5 | $SCRIPT_PATH/colsum.pl; done
echo -n "Compressed Reference Bases: "; for i in `/bin/ls out.rdiff`; do grep DUP $i | cut -f5 | $SCRIPT_PATH/colsum.pl; done
echo -n "Bad Trim: "; for i in `/bin/ls out.qdiff`; do grep BRK $i | cut -f5 | $SCRIPT_PATH/colsum.pl; done
grep '1-to-1' -A3 `/bin/ls out.report` | grep AvgIdentity | sed 's|.*/||' | awk '{print "Avg Idy: "$2}'
grep TotalSNPs `/bin/ls out.report` | sed 's|.*/||' | awk '{print "SNPs: "$2}'
grep TotalIndels `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Indels < 5bp: "$2}'
cat out.rdiff | grep -c GAP | sed 's|.*/||;s|:| |;' | awk '{print "Indels >= 5: "$1}'
grep Inversions `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Inversions: "$3}'
grep Relocations `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Relocation: "$3}'
grep Translocations `/bin/ls out.report` | sed 's|.*/||' | awk '{print "Translocation: "$3}'

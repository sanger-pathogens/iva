import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Collections;

public class getScaffoldStats {  
  private static final int GAP_FUDGE_FACTOR = 1000;
  private static final int MIN_SKIP_SIZE = 200;
  private static final char GAP_CHAR = 'N';
  private static final char GAP_CHAR_LOWER = 'n';
  private static final NumberFormat nf = new DecimalFormat("############.##");

   private static final String[] extensions = {"fa", "fasta", "scafSeq"};
   private int genomeSize = 0;
   private HashMap<String, Scaffold> scfs = new HashMap<String, Scaffold>();

   private HashMap<String, Utils.Pair> sizes = new HashMap<String, Utils.Pair>();

   private static class Scaffold {
      private HashMap<String, Integer> gaps = new HashMap<String, Integer>();
      private HashMap<String, Integer> lengths = new HashMap<String, Integer>();
      private ArrayList<String> contigs = new ArrayList<String>();
      private HashMap<String, Utils.Pair> tiles = new HashMap<String, Utils.Pair>();
   }
   
   public getScaffoldStats(String coords) throws Exception {
      BufferedReader bf = Utils.getFile(coords, "coords");

      String line = null;
      StringBuffer fastaSeq = new StringBuffer();
      String header = "";

      int last = 0;
      int lastStart = 0;
      String currID = "";
      int currStart = 0;
      while ((line = bf.readLine()) != null) {
         String[] split = line.trim().split("\\s+");

         try {
         int start = Integer.parseInt(split[0]);
         if (start - last > 0.1 * Integer.parseInt(split[8]) || !currID.equalsIgnoreCase(split[split.length-1])) {
            if (!currID.equals("")) { 
               Utils.Pair pair = new Utils.Pair(currStart, last);
               if (sizes.get(currID) == null || sizes.get(currID).size() < pair.size()) {
                  sizes.put(currID, pair); 
               }
            }
            currID = split[split.length-1];
            currStart = start;
         }

         int end = Integer.parseInt(split[1]);

         if (start >= lastStart && end <= last) { continue; }
         if (Math.abs(last - start) > 0.1 * Integer.parseInt(split[8])) {
            currStart = start;
         }
         last = end;
         lastStart = start;
         } catch (Exception e) { System.err.println("Skipping line " + line); }
      }
      Utils.Pair pair = new Utils.Pair(currStart, last);
      if (sizes.get(currID) == null || sizes.get(currID).size() < pair.size()) {
         sizes.put(currID, pair);
      }
for (String s : sizes.keySet()) {
System.err.println("Size of " + s + " is " + (sizes.get(s).second - sizes.get(s).first));
}
   }

   private void storeScaffold(String header, String fastaSeq) {
      Scaffold scf = new Scaffold();

      int index = 0;
      int offset = 0;
      String name = null;
      int gapSize = 0;
      for (int i = 0; i < fastaSeq.length(); i++) {
         if (fastaSeq.charAt(i) == GAP_CHAR || fastaSeq.charAt(i) == GAP_CHAR_LOWER) {
            gapSize++; 
            if (scf.contigs.size() <= index) {
//System.err.println("String contig " + name + " and index is " + index + " and offset is " + offset + " and i is " + i);
               int length = i - offset;
               Utils.Pair p = new Utils.Pair(offset, i);
               scf.contigs.add(name);
               scf.lengths.put(name, length);
               scf.tiles.put(name, p);
            }
         } else {
            if (gapSize != 0) {
//System.err.println("Adding gap " + name + " of size " + gapSize);
               scf.gaps.put(name, gapSize);
               index++;
               offset = i;
            }
            gapSize = 0;
            name = header+"_"+index;
         }
      }
     if (scf.contigs.size() <= index) {
//System.err.println("String contig " + name + " and index is " + index + " and offset is " + offset + " and i is " + fastaSeq.length());
        int length = fastaSeq.length() - offset;
        Utils.Pair p = new Utils.Pair(offset, fastaSeq.length());
        scf.contigs.add(name);
        scf.lengths.put(name, length);
        scf.tiles.put(name, p);
     }

System.err.println("Storing scaffold " + header);
      scfs.put(header, scf);

      for (String s : scf.contigs) {
        System.err.println("Scf has contig " + s  + " of length " + scf.lengths.get(s) + " and gap after it is " + scf.gaps.get(s) + " and tile is " + scf.tiles.get(s));
      }
   }

   public void processFasta(String inputFile) throws Exception {
      BufferedReader bf = Utils.getFile(inputFile, extensions);
      
      String line = null;
      StringBuffer fastaSeq = new StringBuffer();
      String header = "";
      
      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) {
            if (fastaSeq.length() != 0) storeScaffold(header, fastaSeq.toString()); 
            header = line.split("\\s+")[0].substring(1);
            fastaSeq = new StringBuffer();
         }
         else {
            fastaSeq.append(line);
         }
      }

      if (fastaSeq.length() != 0) { storeScaffold(header, fastaSeq.toString());
      }
      bf.close();
   }

   public void checkAgainstTiling(String tileFile) throws Exception {
      BufferedReader bf = Utils.getFile(tileFile, "tiling");

      String line = null;
      StringBuffer fastaSeq = new StringBuffer();
      String header = "";

      int count = 0;
      int gaps = 0;
      int scfGaps = 0;
      String currScf = null;
      char currOri = ' ';
      char currOrder = ' ';
      String currName = null;
      int lastIndex = 0;
      int lastOffset = 0;
      int lastEnd = 0;
      int start = 0;
      int totalLength = 0;
      int numIndels = 0;
      int numReloc = 0;
      int numInversion = 0;
      int numTranslocation = 0;
      ArrayList<Integer> lengths = new ArrayList<Integer>();

      int indelCurrScf = 0;
      String startScf = null;
      int startOffset = 0;
      char startOri = ' ';
      int startIndex = -1;
      int startGaps = 0;
      int numBadGaps = 0;
      int numGaps = 0;
      int badGapSize = 0;

      HashMap<String, String> lastChr = new HashMap<String, String>();
      HashMap<String, Character> lastOri = new HashMap<String, Character>();
      HashMap<String, Integer> lastIndices = new HashMap<String, Integer>();
      HashMap<String, Integer> lastPos = new HashMap<String, Integer>();
      HashMap<String, Integer> lastStart = new HashMap<String, Integer>();

      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) { 
            // close last segment
            if (count > 0) { 
               if (currName.equalsIgnoreCase(startScf) && currOri == startOri && Math.abs(lastIndex - startOffset) == 1) {
System.err.println("Looped a circle " +  lastIndex + " " + startOffset + " "+  startScf + " "  +  currName);
                 if (indelCurrScf == 1) {numIndels--;}
                 System.err.println("Looped a circle replacing length " + startIndex);
                 lengths.add(totalLength + scfGaps + startGaps + lengths.get(startIndex));
System.err.println("Adding loop 1 length to scf " + currName + " of " + (totalLength + scfGaps + startGaps + lengths.get(startIndex)));
                 lengths.remove(startIndex);
              } else {
                  lengths.add(totalLength + scfGaps);
System.err.println("Adding length to scf " + currScf + " of " + (totalLength + scfGaps));
              }
            }
System.err.println("Starting a new range in " + currScf + " on scf " + currName + "  done with one from " + (lastEnd) + " to " + start);
            currScf = line.replaceAll(">", "");
            currOri = ' ';
            currName = null;
            currOrder = ' ';
            count = start = gaps = scfGaps = totalLength = 0;
            startScf = null;
            startOffset = 0;
            startIndex = -1;
            startGaps = 0;
            indelCurrScf = 0;
            lastPos.clear();
            lastStart.clear();
            lastIndices.clear();
            continue;
         }

         String[] splitLine = line.trim().split("\\s+");
         String scfName = splitLine[splitLine.length - 1];
         char ori = (splitLine[splitLine.length - 2].equalsIgnoreCase("+") ? '+' : '-');
         int underscore = scfName.trim().lastIndexOf("_");
         if (underscore == -1) { System.err.println("Error unexpected name! " + scfName); System.exit(1); }
         String scfID = scfName.trim().substring(0, underscore);
         Integer index = Integer.parseInt(scfName.trim().substring(underscore+1));
System.err.println("Parsing " + index + " and scf " + scfID);

         Integer length = Integer.parseInt(splitLine[3]);
         Integer offset = Integer.parseInt(splitLine[0]);
         Scaffold scf = scfs.get(scfID);
         String name = scfID + "_" + index;
System.err.println("Processing " + line);
         if (!scf.lengths.get(name).equals(length)) {
            System.err.println("Error for scf " + scfName + " length " + length + " doesnt match " + scf.lengths.get(name));
            System.exit(1);
         }
         if (scf.contigs.size() == 1) { 
            lengths.add(length);
System.err.println("Added empty scf " + scfName + " " + length);
            continue;
         }

         if (startScf == null) {
System.err.println("While reading line " + line + " initialized startScf to " + currName);
            startScf = currName;
            startOri = currOri;
            startOffset = lastIndex;
         }

         if (currOri == ' ' && currName == null) {
System.err.println("While readling line " + line + " initialized currOri to " + scfID);
            currOri = ori;
            currName = scfID;
            start = offset;
            currOrder = ' ';
            count = gaps = scfGaps = lastOffset = totalLength = 0;
         } else {
System.err.println("Processing scf " + currName + " with ori " + currOri + " and offset " + lastOffset + " and info on this is " + scfID + " and ori " + ori + " and offset " + offset + " and math is " + Math.abs(index-lastIndex) + " and gaps is " + gaps + " and scf is " + scfGaps);
            boolean badSkip = false;
            int totalSize = 0;
            int totalOtherSize = 0;
            if (Math.abs(index - lastIndex) > 1 && currName.equalsIgnoreCase(scfID)) {
System.err.println("Processing skip in " + scfID + " from " + index + " to " + lastIndex);
               for (int myIndex = Math.min(index, lastIndex)+1; myIndex < Math.max(index, lastIndex); myIndex++) {
                  if (scfs.get(scfID).lengths.get(scfID + "_" + myIndex) > MIN_SKIP_SIZE) {
                     badSkip = true;
                  }
                  totalSize += scfs.get(scfID).lengths.get(scfID + "_" + myIndex);
                  totalOtherSize = (sizes.get(scfID + "_" + myIndex) == null ? -1 : totalOtherSize + sizes.get(scfID + "_" + myIndex).size());
               }

               int dist = offset - lastPos.get(currName); //lastEnd;
                  if (Math.abs(dist - totalSize) < GAP_FUDGE_FACTOR || (totalOtherSize > 0 && Math.abs(dist - totalOtherSize) < GAP_FUDGE_FACTOR)) {
                     badSkip = false;
                  } else {
                     badSkip = true;
                  } 
System.err.println("The skip avove had dist " + dist + " versus " + totalSize + " which is OK? " + badSkip);
            }

           char order = ' ';
           if (lastIndices.get(scfID) != null) {
              order = (index > lastIndices.get(scfID) ? '+' : '-'); 
              if (currOrder == ' ') {
                 currOrder = (index > lastIndices.get(scfID) ? '+' : '-');
              }
           }
           int gapDifference = -1;
           int otherDiff = Integer.MAX_VALUE;
           if (lastPos.get(scfID) != null) {
              int dist = offset - lastPos.get(scfID);
              Integer lastGap = null;
              if (index > lastIndex) {
                 lastGap = scfs.get(currName).gaps.get(currName + "_" + lastIndex);
              } else {
                 lastGap = scfs.get(currName).gaps.get(currName + "_" + index);
              }
              gapDifference = (lastGap == null ? -1 : Math.abs(dist - lastGap));
System.err.println("The sizes is " + scfID + " and sizes " + sizes.get(scfID + "_" + lastIndex));
              Utils.Pair trueSize = sizes.get(scfID + "_" + lastIndex);
              Utils.Pair secondTrueSize = sizes.get(scfID + "_" + index);
              int otherDist = (trueSize == null || secondTrueSize == null ? -1 : Math.abs((int)secondTrueSize.first - (int)trueSize.second));
              otherDiff = (otherDist == -1 || lastGap == null ? Integer.MAX_VALUE : Math.abs(otherDist - lastGap)); 

System.err.println("Between " + currName + "_" + lastIndex + " and " + index +" expected gap " + dist + " and real gap is " + lastGap + " and estimate gap is " + gapDifference + " and other gap is " + otherDist + " " + otherDiff);
            }

            if (badSkip || currOri != ori || currOrder != order || !currName.equalsIgnoreCase(scfID)) {
               if ((currOri != ori && currName.equalsIgnoreCase(scfID)) || (lastOri.get(scfID) != null && lastOri.get(scfID) != ori && lastChr.get(scfID).equalsIgnoreCase(currScf))) {
                  System.err.println("Inversion found in scf " + scfID);
                  numInversion++;
               } else if (lastChr.get(scfID) != null && !(lastChr.get(scfID).equalsIgnoreCase(currScf))) {
                  System.err.println("Translocation found in scf " + scfID);
                  numTranslocation++;
               }
               else if (currOrder != order && currName.equalsIgnoreCase(scfID)) {
                  System.err.println("Relloc " + scfID);
                  numReloc++;
               }  
               else if (currName.equalsIgnoreCase(scfID)) {
                  System.err.println("Insertion in scf " + scfID);
                  numIndels++;
                  indelCurrScf++;
               }
               if (count > 0) { 
                  if (startIndex < 0) { startIndex = lengths.size(); System.err.println("Initialized start index with scaffold " + currScf + " name " + currName + " to be " + startIndex);}
                  Integer tmpGap = scfs.get(scfID).gaps.get(scfID + "_" + index);
                  startGaps = (tmpGap == null ? 0 : tmpGap);
                  lengths.add(totalLength + scfGaps);
System.err.println("Adding length to scf " + currName + " of " + (totalLength + scfGaps));

                  //lengths.add(lastEnd-start+1);
               }
               System.err.println("Starting a new range in " + currScf + " on scf " + currName + "  done with one from " + (lastEnd) + " to " + start);
               currOri = ori;
               currName = scfID;
               start = offset;
               gaps = scfGaps = lastOffset = totalLength = 0;
               currOrder = ' ';
            }
            else {
              if (gapDifference >= 0 && currName.equalsIgnoreCase(scfID)) {
                 if (gapDifference > GAP_FUDGE_FACTOR && otherDiff > GAP_FUDGE_FACTOR) numBadGaps++;
                 badGapSize+=Math.min(otherDiff, gapDifference);
                 numGaps++;
System.err.println("Difference is " + Math.min(otherDiff, gapDifference));
              }
           }
         }
         if (scfs.get(currName) != null && index < scfs.get(currName).contigs.size()-1 && index > 0) {
            scfGaps += scfs.get(currName).gaps.get(currName + "_" + index);
         }
         lastIndex = index;
         gaps += (lastOffset != 0 ? offset - lastEnd : gaps);
         lastOffset = offset;
         lastEnd = offset + length;
         totalLength += length;
         lastOri.put(scfID, ori);
         lastPos.put(scfID, lastEnd);
         lastStart.put(scfID, lastOffset);
         lastOri.put(scfID, ori);
         lastChr.put(scfID, currScf);
         lastIndices.put(scfID, index);
System.err.println("After processing " + scfID + "_" + index + " the length is " + totalLength + " and gaps is " + scfGaps);
         count++;
      }

      if (count != 0) { 
         if (currName.equalsIgnoreCase(startScf) && currOri == startOri && Math.abs(lastIndex - startOffset) == 1)
{
              // we looped a circle
              System.err.println("Looped a circle replacing length " + startIndex);
              lengths.add(totalLength + scfGaps + lengths.get(startIndex));
System.err.println("Adding loop 2 length to scf " + currName + " of " + (totalLength + scfGaps + startGaps + lengths.get(startIndex)));
System.err.println("The length originally was " + totalLength + " and the gaps are " + scfGaps + " and start gap is " + startGaps + " and length at start is " + lengths.get(startIndex) + " where start index is " + startIndex);

              lengths.remove(startIndex);
           } else {
               lengths.add(totalLength + scfGaps);
System.err.println("Adding length to scf " + currName + " of " + (totalLength + scfGaps));
           }
      } 
      System.err.println("Final Processing scf " + currName + " with ori " + currOri + " and offset " + lastIndex + " and info on this is " + " and gaps is " + gaps + " and scf is " + scfGaps + " length " + lastEnd + " to " + start);

      bf.close();


      // now we can compute all the stats
      System.out.println("#Indels, #Inversions, #Translocations, #Relocations, #Gaps >= " + GAP_FUDGE_FACTOR + ", Mean Gap Error, #Contigs, Corrected N50, Corrected E-size");
      System.out.print(numIndels + "," + numInversion + "," + numTranslocation + "," + numReloc + "," + numBadGaps + "," + nf.format((double)badGapSize / numGaps));
      Collections.sort(lengths); 
      Collections.reverse(lengths);

      int sum = 0;
      double fSize = 0;
      boolean doneN50 = false;
      for (int i = 0; i < lengths.size(); i++) {
         System.err.println("At index " + i + " I have length " + lengths.get(i));
         sum += lengths.get(i);

         // compute the F-size
         if (lengths.get(i) > MIN_SKIP_SIZE) {
            fSize += Math.pow(lengths.get(i), 2);
         }
         if (!doneN50 && sum / (double)genomeSize >= 0.5) {
            System.out.print("," + (i+1) + "," + lengths.get(i));
            doneN50=true;
         }
      }
      fSize /= genomeSize;
      System.out.println("," + nf.format(fSize));
   }

   public static void printUsage() {
      System.err.println("This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.");
      System.err.println("Example usage: getScaffoldStats fasta1.fasta,fasta2.fasta");
   }
   
   public static void main(String[] args) throws Exception {     
      if (args.length < 1) { printUsage(); System.exit(1);}

      getScaffoldStats f = new getScaffoldStats(args[3]);
      f.genomeSize = Integer.parseInt(args[2]);
      f.processFasta(args[0]);

      f.checkAgainstTiling(args[1]);
   }
}

package nextgen.core.programs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import net.sf.picard.util.Log;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.BAMIndexer;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 * Utility to go through a sorted alignment file and mark
 * likely PCR duplicates based on both alignment location 
 * and a molecular barcode added to the read name. 
 * @author mgarber
 *
 */
public class MarkDuplicatesFromNameUMI {
	private static final Log log = Log.getInstance(MarkDuplicatesFromNameUMI.class);
	static final String USAGE = "MarkDuplicatesFromNameUMI. Finds likely PCR duplicates by looking at BOTH location and a molecular barcode placed within the read name. This program assumes the barcode is the last N characters: " +
			"\n\t-in <Path to the sorted BAM or SAM alignment file (NO standard input is supported must input a file path)>" +
			"\n\t-out <Ouput BAM/SAM file or standard out if none is provided>" +
			"\n\t-bcSize <Number of characters from the end of the read name that constitute the molecular barcode >" +
			"\n\t-bcHistogram <File name to write the histogram of barcode ussage or stdout if none is provided>"+
			"\n";

	public static void main (String [] args) throws IOException {
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "default");
		String alnFile = argMap.getInput();
		int bcSize = argMap.getInteger("bcSize");
		String out = argMap.getOutput();

		SAMFileReader reader = new SAMFileReader(new File(alnFile));
		BAMFileWriter writer = new BAMFileWriter(new File(out));
		SAMFileHeader header = reader.getFileHeader();
		writer.setSortOrder(SortOrder.coordinate, false);
		header.addProgramRecord(new SAMProgramRecord("nextgen.core.programs.MarkDuplicatesFromNameUMI -bcSize 4 -in tophat_to_genome/10ng_R2.mark.dups.bam -out tophat_to_genome/10ng_R2.mark.dups.mark.dups.bam -bcHistogram tophat_to_genome/10ng_R2.mark.dups.mark.dups.bam.bchist.txt"));
		writer.setHeader(header);
		try {
			MarkDuplicatesResults result = markDuplicates(bcSize, reader, writer);
			System.out.println("REPORT: ");
			System.out.println("Total alignments: " + result.alnNum);
			System.out.println("Total dup alignments: " + result.numDup);
			System.out.println("Total unique alignments: " + result.numUnique);
			System.out.println("Mean duplications " + result.duplicationStats.getMean());
			System.out.println("Max duplications " + result.duplicationStats.getMax());
			System.out.println("Number aligned marked as dups " + result.numAlignedMarked);
			System.out.println("Number unmarked alignments " + result.numUnmarked);

			if(argMap.containsKey("bcHistogram")) {
				BufferedWriter bw = new BufferedWriter(new FileWriter(argMap.getMandatory("bcHistogram")));
				result.writeBCHistogram(bw);
				bw.close();
			}
			
		}finally {
			reader.close();
			writer.close();
		}
		
		/*
		try {
			BAMIndexer.createAndWriteIndex(new File(out), new File(out+".bai"), false);
		} finally {
			log.debug("Created index");
		}*/


	}

	public static MarkDuplicatesResults markDuplicates(int bcSize, SAMFileReader reader, BAMFileWriter writer) {

		MarkDuplicatesResults result = new MarkDuplicatesResults();

		int currentStart = 0;
		IntervalTree<Set<String>> currentBarcodeTree = new IntervalTree<Set<String>>();
		HashMap<String, Integer>  currentBarcodeDuplicates = new HashMap<String, Integer>(); 


		SAMRecordIterator sri = reader.iterator();
		while(sri.hasNext()) {
			SAMRecord samR = sri.next();
			String readName = samR.getReadName();
			String bc = readName.substring(readName.length() - bcSize);
			int start = samR.getAlignmentStart();
			int end   = samR.getAlignmentEnd();
			result.alnNum++;
			if(currentStart != start || !sri.hasNext()) { //Its OK if we do not mark the last alignment as duplicated ... we'll leave
				currentStart = start;
				currentBarcodeTree.clear();
				for (String seenBC : currentBarcodeDuplicates.keySet()) {
					int numOfSeenBC = currentBarcodeDuplicates.get(seenBC);
					if( numOfSeenBC == 1) {
						result.numUnique++;
					} else {
						result.duplicationStats.addValue(numOfSeenBC);
						result.numDup += numOfSeenBC;
					}
				}
				currentBarcodeDuplicates.clear();
			}

			String compositeBCName = getCompositeBC(bc, start, end);
			Node<Set<String>> intervalBarcodeNode = currentBarcodeTree.find(start, end);
			if(intervalBarcodeNode == null) {
				Set<String> intervalBarcodes = new TreeSet<String>();
				intervalBarcodes.add(bc);
				currentBarcodeTree.put(start, end, intervalBarcodes);
				currentBarcodeDuplicates.put(compositeBCName, 1);
				result.numUnmarked++;
				result.addBarcode(bc);
			} else {
				Set<String> intervalBarcodes = intervalBarcodeNode.getValue();
				if (intervalBarcodes.contains(bc) ) {
					samR.setDuplicateReadFlag(true);
					result.numAlignedMarked++;
					currentBarcodeDuplicates.put(compositeBCName, currentBarcodeDuplicates.get(compositeBCName)+1);
				} else {
					currentBarcodeDuplicates.put(compositeBCName, 1);
					intervalBarcodes.add(bc);
					result.numUnmarked++;
				}
			}

			writer.addAlignment(samR);
		}

		return result;

	}

	private static String getCompositeBC(String bc, int start, int end) {
		return bc+start+"."+end;
	}

	public static class MarkDuplicatesResults {

		public int numUnmarked = 0;
		public int numAlignedMarked = 0;
		public SummaryStatistics duplicationStats = new SummaryStatistics();
		public int numUnique = 0;
		public int numDup = 0;
		public int alnNum = 0;
		HashMap<String, Integer> bcOccurrences = new HashMap<String, Integer>();

		public void addBarcode(String bc) {
			if(!bcOccurrences.containsKey(bc)) {
				this.bcOccurrences.put(bc, 1);
			} else {
				this.bcOccurrences.put(bc, bcOccurrences.get(bc) + 1);
			}

		}

		public void writeBCHistogram(BufferedWriter bw) throws IOException {
			bw.write("Barcode\tCounts\tFraction");
			bw.newLine();

			List<String> bcs = new ArrayList<String>(bcOccurrences.keySet());
			List<Double> bcsCounts = new ArrayList<Double>(bcOccurrences.keySet().size());

			for( String bc : bcs) {
				bcsCounts.add((double)bcOccurrences.get(bc));
			}
			Collections.sort(bcs, new Comparator<String>() {
				public int compare(String o1, String o2) {
					return bcOccurrences.get(o2) - bcOccurrences.get(o1);
				}
			});

			double sum = Statistics.sum(bcsCounts);
			for(String bc : bcs)  {
				bw.write(bc);
				bw.write("\t");
				bw.write(String.valueOf(bcOccurrences.get(bc)));
				bw.write("\t");
				bw.write(String.valueOf(bcOccurrences.get(bc)/(double)sum));
				bw.newLine();
			}
		}

	}

}

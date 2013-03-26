package nextgen.core.programs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.lang.Math;

import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

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
 * and a molecular barcode added to the read name. Optional
 * second pass through data adjusts markings at higher depths
 * to account for expected duplicate assignments of barcodes
 * to different molecules.
 * @author mgarber
 *
 */
public class MarkDuplicatesFromNameUMI {
	private static final Log log = Log.getInstance(MarkDuplicatesFromNameUMI.class);
	static final String USAGE = "MarkDuplicatesFromNameUMI Finds likely PCR duplicates by looking at BOTH location and a molecular barcode placed within the read name. This program assumes the barcode is the last N characters: " + 
			"\n\t-in <Path to the sorted BAM or SAM alignment file (NO standard input is supported must input a file path)>" +
			"\n\t-out <Ouput BAM/SAM file or standard out if none is provided>" +
			"\n\t-out2 <Adjusted Ouput BAM/SAM file, no adjustment if none is provided>" +
			"\n\t-bcSize <Number of characters from the end of the read name that constitute the molecular barcode >" +
			"\n\t-bcHistogram <File name to write the histogram of barcode ussage or stdout if none is provided>" +
			"\n\t-bcCountsFile <Optional file name to record counts: number unique, number barcodes duplicated and total number of duplicates for each position>" +
			"\n\t-bcStatsFile<Optional file name to write to write stats on total dups by number unique bcs>" +
			"\n\t-calculateThreshold1<Optional flag, TRUE indicates calculate threshold1, default TRUE>" +
			"\n\t-initialThreshold1<Optional initial upper cut-off for number of barcodes duplicated, if not calculating threshold1, final value as well, default 30>" +
			"\n\t-calculateThreshold2<Optional flag, TRUE indicates calculate threshold2, default TRUE>" +
			"\n\t-initialThreshold2<Optional initial upper cut-off for number of barcodes duplicated, if not calculating threshold2, final value as well, default 30>" +
			"\n\t-allowedFracOfMean<Optional double between 0 and 1, model parameter,default 0.02>" +
			"\n\t-numToIncrement<Optional postive integer, model parameter, default 5>" +
			"\n\t-revisedBCCountsFile <Optional file name to record revised counts by position for positions with number of barcodes duplicated above minNumBCsToApplyModel>" +
			"\n\t-revisedBCStatsFile<Optional file name to write revised stats on total dups by adjusted number unique bcs when  number of barcodes duplicated above minNumBCsToApplyModel>" +
			"\n";
	
	

	public static void main (String [] args) throws IOException {
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "default");
		String alnFile = argMap.getInput();
		int bcSize = argMap.getInteger("bcSize");
		int numBCs = (int) Math.pow(4, bcSize);
		
		String out = argMap.getOutput();
		

		SAMFileReader reader = new SAMFileReader(new File(alnFile));
		BAMFileWriter writer = new BAMFileWriter(new File(out));
		SAMFileHeader header = reader.getFileHeader();
		writer.setSortOrder(SortOrder.coordinate, false);
		header.addProgramRecord(new SAMProgramRecord("nextgen.core.programs.MarkDuplicatesFromNameUMI --in tophat_to_genome/50ng_R2.bam -out markdup.two_pass/50ng_R2.mark.dups.pass1.bam -out2 markdup.two_pass/50ng_R2.mark.dups.pass2.bam-bcSize 4 -bcCountsFile -bcHistogram markdup.two_pass/50ng_R2.mark.dups.mark.dups.bam.bchist.txt markdup.two_pass/50ng_R2.counts.txt -bcStatsFile markdup.two_pass/50ng_R2.stats.txt -revisedBCCountsFile markdup.two_pass/50ng_R2.counts.revised.txt -revisedBCStatsFile markdup.two_pass/50ng_R2.stats.revised.txt"));
		writer.setHeader(header);
		
		MarkDuplicatesResults result;					
		String out2;
		SAMFileReader reader2;
		BAMFileWriter writer2;
		boolean calcThreshold1;
		if (argMap.containsKey("calculateThreshold1")){
			calcThreshold1= argMap.isFlagTrue("calculateThreshold1");
		} else {
			calcThreshold1 = true;
		}
		int	threshold1 = argMap.getInteger("initialThreshold1",30);
		
		boolean calcThreshold2;
		if (argMap.containsKey("calculateThreshold1")){
			calcThreshold2= argMap.isFlagTrue("calculateThreshold2");
		} else {
			calcThreshold2 = true;
		}
		int	threshold2 = argMap.getInteger("initialThreshold2",30);
		
		
		double allowedFracOfMean = argMap.getDouble("allowedFracOfMean", 0.02);
		int numToIncrement = argMap.getInteger("numToIncrement",5);
		if (argMap.containsKey("out2")){ 
			if (threshold1 < 2){
				threshold1 = 30;
				System.out.println("Value of initialThreshold1 must be at least 2, using default value 30.");
			}
			
			if (threshold2 < threshold1){
				threshold2 = threshold1;
				System.out.println("Value of initialThreshold2 must be at least that of initialThreshold1.");				
				System.out.println("Increasing initialThreshold2 to " + threshold2 + ".");		
			}
			
					
			if (allowedFracOfMean <= 0){
				System.out.println("Value of allowedFracOfMean must be postiive.  Using default value of 0.02.");
				allowedFracOfMean = 0.02;
			}
			
			if (numToIncrement < 1){
				System.out.println("Value of numToIncrement must be > 0, using default value of 5.");
				numToIncrement = 5;			
			}
		}
		

		boolean exportCounts = false;
		if(argMap.containsKey("bcCountsFile")) {
			exportCounts = true;
		}

		
		try {
			result = markDuplicates(bcSize, reader, writer, exportCounts, threshold1, threshold2, allowedFracOfMean, numToIncrement);
			
			if (calcThreshold1){
				result.setThreshold1AndMeanNumDup(numBCs);
			} else {
				result.calcMeanNumDup();
				System.out.println("threshold1: " + threshold1);
			}
			threshold1 = result.threshold1;
			if (calcThreshold2){
				result.setThreshold2(numBCs,numToIncrement);  // put if statement to allow for passing in still
			} else {
				result.threshold2 = Math.max(threshold1,threshold2);
				System.out.println("threshold2: " + threshold2);
				
			}
			
			//result.computeMeanNumDup();  //deprecated using regression instead
			
			threshold2 = result.threshold2;
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
			
			if(argMap.containsKey("bcCountsFile")) {
				BufferedWriter bw = new BufferedWriter(new FileWriter(argMap.getMandatory("bcCountsFile")));
				result.writeCounts(bw);
				bw.close();
			}
			
			if(argMap.containsKey("bcStatsFile")) {
				BufferedWriter bw = new BufferedWriter(new FileWriter(argMap.getMandatory("bcStatsFile")));
				result.writeBCStats(bw);
				bw.close();
			}
	
			
		}finally {
			reader.close();
			writer.close();
		}
		
		if(!argMap.containsKey("out2")){
		
			System.out.println("No out2 file.  Adjustments to initial markings are not performed.");
		
		} else if(threshold2 >= Math.pow(4, bcSize)) {
		
			System.out.println("Adjustments to initial markings are not performed as threshold2 equals or exceeds the number of possible barcodes.");
		
		} else {
			out2 = argMap.getMandatory("out2");
			reader2 = new SAMFileReader(new File(out));
			writer2 = new BAMFileWriter(new File(out2));
			writer2.setSortOrder(SortOrder.coordinate, false);
			SAMFileHeader header2 = reader2.getFileHeader();
			writer2.setHeader(header2);
			try {
				
				adjustNumFlagged(bcSize, reader2, writer2, result);
				System.out.println("REPORT for 2nd Pass: ");
				System.out.println("Number aligned marked as dups " + result.numAlignedMarked);
				System.out.println("Number unmarked alignments " + result.numUnmarked);
				System.out.println("Number of aligned first marked, then changed back to unmarked " + result.numRemovedMark);
				
				if(argMap.containsKey("revisedBCCountsFile")) {
					BufferedWriter bw = new BufferedWriter(new FileWriter(argMap.getMandatory("revisedBCCountsFile")));
					result.writeRevisedCounts(bw, threshold2);
					bw.close();
				}
					
				if(argMap.containsKey("revisedBCStatsFile")) {
					BufferedWriter bw = new BufferedWriter(new FileWriter(argMap.getMandatory("revisedBCStatsFile")));
					result.writeRevisedBCStats(bw, threshold2);
					bw.close();
				}
				
			} finally{
				reader2.close();
				writer2.close();
			}
		}
		
		/*
		try {
			BAMIndexer.createAndWriteIndex(new File(out), new File(out+".bai"), false);
		} finally {
			log.debug("Created index");
		}*/


	}

	public static MarkDuplicatesResults markDuplicates(int bcSize, SAMFileReader reader, BAMFileWriter writer, 
			boolean exportCounts, int threshold1, int threshold2, double allowedFracOfMean, int numToIncrement) {

		MarkDuplicatesResults result = new MarkDuplicatesResults(threshold1, threshold2, allowedFracOfMean);

		int currentStart = 0;
		IntervalTree<Set<String>> currentBarcodeTree = new IntervalTree<Set<String>>();
		HashMap<String, Integer>  currentBarcodeDuplicates = new HashMap<String, Integer>(); 
		int numBCs = (int) Math.pow(4, bcSize);
		int curAlnNum = 0;
		int lastStart = 0;
		
		
		result.initHistNumUniqueBC(numBCs);
		result.initHistNumDupBC(numBCs);
		result.initStatsByNumDistinctBCSeen(numBCs);
		result.initStatsByNumBCDup(numBCs);
		int start = 0;
		int end = 0;
		

		SAMRecordIterator sri = reader.iterator();
		while(sri.hasNext()) {
			SAMRecord samR = sri.next();
			String readName = samR.getReadName();
			String bc = readName.substring(readName.length() - bcSize);
			
			start = samR.getAlignmentStart();
			end   = samR.getAlignmentEnd();
			result.alnNum++;
			curAlnNum++;
			if(currentStart != start || !sri.hasNext()) { //Its OK if we do not mark the last alignment as duplicated ... we'll leave
				lastStart = currentStart;
				currentStart = start;
				currentBarcodeTree.clear();
				int currentNumDup = 0;
				int currentNumUnique = 0;

				for (String seenBC : currentBarcodeDuplicates.keySet()) {
					int numOfSeenBC = currentBarcodeDuplicates.get(seenBC);
					
			
					if( numOfSeenBC == 1) {
						result.numUnique++;
						currentNumUnique++;
					} else {
						result.duplicationStats.addValue(numOfSeenBC);
						result.numDup += numOfSeenBC;
						currentNumDup += numOfSeenBC;
					}
				}
				int currentNumBCDup = currentBarcodeDuplicates.size() - currentNumUnique;
				
				if (currentBarcodeDuplicates.size() != 0){
					if (exportCounts){
						result.addCountsByPosition(currentNumUnique, currentNumBCDup,currentNumDup);
					}
					result.histNumUniqueBC[currentBarcodeDuplicates.size()-1] += 1;	
					result.updateStatsByNumDistinctBCSeen(currentBarcodeDuplicates.size(), currentNumBCDup);
					
					if (currentBarcodeDuplicates.size() != currentNumUnique){
			
						result.updateStatsByNumBCDup(currentNumBCDup, currentNumDup);
						
					}
				}
				
				
				if (currentNumBCDup != 0){
					result.histNumDupBC[currentNumBCDup - 1] += 1;
					if (currentNumBCDup >= threshold2){
						result.addStartToEvalList(result.alnNum - curAlnNum, lastStart);
						result.addRevisedCountsByPosition(currentNumUnique, currentNumBCDup,currentNumDup);  //these initial counts  used in adjusting marking and corrected on 2nd pass
					}
				}

				curAlnNum = 0;
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

	
	
	private static void adjustNumFlagged(int bcSize, SAMFileReader reader2, BAMFileWriter writer2, MarkDuplicatesResults result){	

		int currentStart = 0;
		double dupScore = 0;
		boolean evalFlag = false;
		int startsToEvalIndex = 0;
		int firstAlnToEval = result.startsToEval.get(0)[0];
		int firstStartToEval = result.startsToEval.get(0)[1];
		int currentNumUnique = 0;
		int currentAdjustedNumBCDup = 0; 
		int currentNumDup = 0;
		result.numUnmarked = 0;
		result.numAlignedMarked = 0;
		result.numRemovedMark = 0;
		int alnNum = 0;
		// using depth of 10 as estimate, lower should be better estimate as all dups expected to be pcr dups, but lowest depth numbers are a little funky.
		double estStDevDupPerBC = result.statsByNumBCDup[9].getStandardDeviation()/Math.sqrt(10);  
		SAMRecordIterator sri = reader2.iterator();
		while(sri.hasNext()) {
			SAMRecord samR = sri.next();
			int start = samR.getAlignmentStart();
			boolean readFlagged = samR.getDuplicateReadFlag();
			alnNum++;
			if(currentStart != start || !sri.hasNext()) { //Its OK if we do not mark the last alignment as duplicated ... we'll leave
				currentStart = start;
				if (evalFlag){
					
					int e[] = {currentNumUnique, currentAdjustedNumBCDup,currentNumDup};
					result.revisedPositionCounts.set(startsToEvalIndex - 1, e);
					if (currentAdjustedNumBCDup > result.maxAdjustedNumBCDup) {
						result.maxAdjustedNumBCDup = currentAdjustedNumBCDup;
					}
				}
				
				if (alnNum == firstAlnToEval && start == firstStartToEval){

				 	currentAdjustedNumBCDup = result.revisedPositionCounts.get(startsToEvalIndex)[1];
				 	currentNumDup = result.revisedPositionCounts.get(startsToEvalIndex)[2];
				 	if (currentAdjustedNumBCDup > result.threshold2 && 
				 			currentNumDup > (result.meanNumDup*currentAdjustedNumBCDup + 3*Math.sqrt(currentAdjustedNumBCDup)*estStDevDupPerBC)) {
						evalFlag = true;
						currentNumUnique = result.revisedPositionCounts.get(startsToEvalIndex)[0];				
						dupScore = result.meanNumDup/2 - currentAdjustedNumBCDup*result.meanNumDup - currentNumUnique; 				
					} else {
						evalFlag = false;
					}
				 	startsToEvalIndex ++;						
					if (startsToEvalIndex < result.startsToEval.size()){
						firstAlnToEval = result.startsToEval.get(startsToEvalIndex)[0];
						firstStartToEval = result.startsToEval.get(startsToEvalIndex)[1];
					} else {
						System.out.println("Completed adjustments.");
					}
							
				} else {
					evalFlag = false;
				}
				
								
			}
			if (evalFlag){
				dupScore++;
				if (!readFlagged){
					
					result.numUnmarked++;

				
				} else if (readFlagged){
					if (dupScore >= result.meanNumDup){
						samR.setDuplicateReadFlag(false);  
						dupScore = dupScore - result.meanNumDup;
						result.numRemovedMark++;
						result.numUnmarked++;
						currentAdjustedNumBCDup++;
						
					
					} else {
						result.numAlignedMarked++;
						
					}
				}
			} else {
				if (readFlagged){
					result.numAlignedMarked++;			
				} else {
					result.numUnmarked++;
				}
			}
			
			writer2.addAlignment(samR);
					
			
		}
		result.tallyRevisedStats();

	}

	public static class MarkDuplicatesResults {

		public int numUnmarked = 0;
		public int numAlignedMarked = 0;
		public SummaryStatistics duplicationStats = new SummaryStatistics();
		public int numUnique = 0;
		public int numDup = 0;
		public int alnNum = 0;
		HashMap<String, Integer> bcOccurrences = new HashMap<String, Integer>();
		List<int[]> positionCounts = new ArrayList<int[]>(); // currentNumUnique, currentNumBCdup, currentNumDup
		List<int[]> revisedPositionCounts = new ArrayList<int[]>(); // currentNumUnique, currentNumBCdup, currentNumDup
		public int[] histNumUniqueBC;  // eliminate? tracking but no longer using...add option to export?
		public int[] histNumDupBC;
		public SummaryStatistics[] statsByNumBCDup;
		public SummaryStatistics[] statsByNumDistinctBCSeen;
		public int maxAdjustedNumBCDup = 0;
		public SummaryStatistics[] revisedStatsByAdjustedNumBCDup;
		public double meanNumDup;
		public int numRemovedMark = 0;
		public List<int[]> startsToEval = new ArrayList<int[]>(); //1st alignment #, start 
		
		SimpleRegression lfit = new SimpleRegression(false);  // better name
		
		public int threshold1;
		public int threshold2;
		public double allowedFracOfMean;
		
		public MarkDuplicatesResults(int threshold1, int threshold2, double allowedFracOfMean){
			
			
			this.threshold1 = threshold1;
			this.threshold2 = threshold2;
			this.allowedFracOfMean = allowedFracOfMean;
			
		}
		

		public void addBarcode(String bc) {
			if(!bcOccurrences.containsKey(bc)) {
				this.bcOccurrences.put(bc, 1);
			} else {
				this.bcOccurrences.put(bc, bcOccurrences.get(bc) + 1);
			}

		}
		
		public void addCountsByPosition(int currentNumUnique, int currentNumBCDup, int currentNumDup){
			
			int[] e = {currentNumUnique,currentNumBCDup,currentNumDup};
			positionCounts.add(e);				
		
		}
		
		public void addRevisedCountsByPosition(int currentNumUnique, int currentNumBCDup, int currentNumDup){
			
			int[] e = {currentNumUnique,currentNumBCDup,currentNumDup};
			revisedPositionCounts.add(e);				
		
		}
		
		public void addStartToEvalList(int alnNumber, int start){
			int[] e = {alnNumber, start};
			startsToEval.add(e);
		}
		
		public void initHistNumUniqueBC(int numBCs){
			histNumUniqueBC = new int[numBCs];
		}
		
		public void initHistNumDupBC(int numBCs){
			
			histNumDupBC = new int[numBCs];
		}
		
		
		public void initStatsByNumBCDup(int numBCs){
			
			statsByNumBCDup = new SummaryStatistics[numBCs];
			
			for (int i = 0; i < numBCs; i++){
				
				statsByNumBCDup[i] = new SummaryStatistics();
			
			}	
		}
		
		public void updateStatsByNumBCDup(int currentNumBCDup, int currentNumDup){
			
			statsByNumBCDup[currentNumBCDup - 1].addValue(currentNumDup);
			
		}
		
		public void setThreshold1AndMeanNumDup(int numBCs){  
					
					double cumsum = 0;
					double numBCDupRep = 0;
					double RMSE = 0;
					double initRMSEBound = 0;
					double rmseBound = 0;
					
					
					for  (int i = 0; i < threshold1; i++){
						if (statsByNumBCDup[i].getN() > 0){
							lfit.addData( i + 1,statsByNumBCDup[i].getMean());
							cumsum = cumsum  + statsByNumBCDup[i].getMean();
							numBCDupRep++;
						}
						
						
	
					}
					// minimum 2 values to fit line:
					while (numBCDupRep < 2){
						threshold1++;
						if (statsByNumBCDup[threshold1].getN() > 0){
							lfit.addData( threshold1,statsByNumBCDup[threshold1 - 1].getMean());
							cumsum = cumsum  + statsByNumBCDup[threshold1 - 1].getMean();
							numBCDupRep++;
							
						}
					}

					
					double initRMSE = Math.sqrt(lfit.getMeanSquareError());
					allowedFracOfMean = Math.max(allowedFracOfMean,Math.ceil(100*initRMSE*numBCDupRep/cumsum)/100);
					initRMSEBound = Math.max(initRMSE,1);  
					initRMSEBound = Math.ceil(10*initRMSEBound)/10;
					
					
					for (int i = threshold1; i < numBCs; i++){
						if (statsByNumBCDup[i].getN() > 0){
							lfit.addData( i + 1,statsByNumBCDup[i].getMean());
							RMSE = Math.sqrt(lfit.getMeanSquareError());
							cumsum = cumsum  + statsByNumBCDup[i].getMean();
							numBCDupRep++;
							rmseBound = Math.max(initRMSEBound, allowedFracOfMean*cumsum/numBCDupRep);  // confirm cumsum not 0
							rmseBound = Math.ceil(10*rmseBound)/10;
							
						
						
							if (RMSE < rmseBound){
								threshold1++;
							} else {
								lfit.removeData( i + 1,statsByNumBCDup[i].getMean());
								
								break;
							}
						}
					}
					System.out.println("New threshold1: " + threshold1);
					meanNumDup = lfit.getSlope();
					System.out.println("Slope: " + lfit.getSlope());
					return;
				}
		
		public void calcMeanNumDup(){  // needed when not calculating threshold1
			
			double numBCDupRep = 0;
			
			for  (int i = 0; i < threshold1; i++){
				if (statsByNumBCDup[i].getN() > 0){
					lfit.addData( i + 1,statsByNumBCDup[i].getMean());
					numBCDupRep++;
	
				}
			}
			
			// minimum 2 values to fit line:
			while (numBCDupRep < 2){
				threshold1++;
				if (statsByNumBCDup[threshold1].getN() > 0){
					lfit.addData( threshold1,statsByNumBCDup[threshold1 - 1].getMean());
					numBCDupRep++;
					
				}
			}

			meanNumDup = lfit.getSlope();
			System.out.println("Slope: " + lfit.getSlope());

			
		}
		public void setThreshold2(int numBCs, int numToIncrement){ 
			
			threshold2 = Math.max(threshold1, threshold2);
			if (threshold2 == numBCs){
				return;
			} else {
				
				int lastThreshold2 = threshold2;
				threshold2 = Math.min(threshold2 + numToIncrement, numBCs);
				
			
				while (threshold2 < numBCs){
					double sse = 0;
					double sumActual = 0;
					double sumPredict = 0;
					int numBCRep = 0;
					
					for (int i = lastThreshold2; i < threshold2; i++){
						 if (statsByNumBCDup[i].getN() > 0){
						 
							 sse = sse + Math.pow(lfit.predict(i + 1) - statsByNumBCDup[i].getMean(), 2);
							 sumActual = sumActual + statsByNumBCDup[i].getMean(); 
							 sumPredict = sumPredict + lfit.predict(i + 1);
							 numBCRep ++;
						 }
					}
					
					
					
					if (numBCRep > 0){	

						
						sse = Math.sqrt(sse/numBCRep);  //? correct to use n-1 in denom?  going with no, no bias, no constants selected						
						// Require both sufficiently large rmse and mean increase >= 1 to reject model as representative of new data. 
						if ((sse*(threshold2 - lastThreshold2)/sumActual < allowedFracOfMean) || (sumActual - sumPredict < allowedFracOfMean*sumActual)){
							lastThreshold2 = threshold2;
							threshold2 = Math.min(threshold2 + numToIncrement, numBCs);
							if (threshold2 == lastThreshold2){
								break;
							}
						} else{
							threshold2 = lastThreshold2;
							break;

						}
					
					} else {
						
						threshold2 = Math.min(threshold2 + numToIncrement, numBCs);
						
					}
				
				}
			}
				
			System.out.println("New threshold2: " + threshold2);
			return;
			
		}
					
		
		public void initStatsByNumDistinctBCSeen(int numBCs){
			
			statsByNumDistinctBCSeen = new SummaryStatistics[numBCs];
			
			for (int i = 0; i < numBCs; i++){
				
				statsByNumDistinctBCSeen [i] = new SummaryStatistics();
				
			}	
			
		}
		
		public void updateStatsByNumDistinctBCSeen(int currentNumDistinctBC, int currentNumDup){
			
			statsByNumDistinctBCSeen[currentNumDistinctBC - 1].addValue((double)currentNumDup);
				
		}

		
		public void computeMeanNumDup(){		// deprecated...using regression
			

			
			int numBCsDupBelowThreshold = 0;  // counts barcode-position pairs that have been duplicated at positions where the count 
												  // remains below (and including) threshold
			double numDupsBelowThreshold = 0;	  // total count of corresponding duplicates. 	

			for (int i = 0; i < threshold1; i++){
				numBCsDupBelowThreshold += (i+1)*histNumDupBC[i];
				if (histNumDupBC[i] > 0 ){
					numDupsBelowThreshold += histNumDupBC[i]*statsByNumBCDup[i].getMean();
				}				
			}
			
			meanNumDup = (double) numDupsBelowThreshold/numBCsDupBelowThreshold;
			System.out.println("Mean number of duplicates: " + meanNumDup);
			
			
		}
		
		
		
		
		public void tallyRevisedStats(){
			
			int statsLength = maxAdjustedNumBCDup - threshold2 + 1;
			
			if (statsLength <= 0){
				
				return;
			}
			
			
			revisedStatsByAdjustedNumBCDup = new SummaryStatistics[statsLength];
			

			
			for (int i = 0; i < statsLength; i++){
				
				revisedStatsByAdjustedNumBCDup[i] = new SummaryStatistics();
				
			}
			
			for (int i = 0; i < revisedPositionCounts.size(); i++){
				
				if (revisedPositionCounts.get(i)[1] > threshold2){
					revisedStatsByAdjustedNumBCDup[revisedPositionCounts.get(i)[1] - threshold2].addValue((double) revisedPositionCounts.get(i)[2]);
				}

				
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
		
		public void writeCounts(BufferedWriter bw) throws IOException{
			
			bw.write("NumberUniqueBarcodes\tNumberDuplicatedBarCodes\tNumberDuplicateReads");
			
			bw.newLine();
			for (int i = 0; i < positionCounts.size(); i++){				
				bw.write(Integer.toString(positionCounts.get(i)[0]));
				bw.write("\t");
				bw.write(Integer.toString(positionCounts.get(i)[1]));
				bw.write("\t");
				bw.write(Integer.toString(positionCounts.get(i)[2]));
				bw.newLine();								
			}
			System.out.println("Saved counts.");
		
		}
		
		public void writeRevisedCounts(BufferedWriter bw, int threshold2) throws IOException{
			
			
			bw.write("NumberUniqueBarcodes\tAdjustedNumberDuplicatedBarCodes\tNumberDuplicateReads");
			bw.newLine();
			for (int i = 0; i < revisedPositionCounts.size(); i++){	
					if (revisedPositionCounts.get(i)[1] > threshold2){
					bw.write(Integer.toString(revisedPositionCounts.get(i)[0]));
					bw.write("\t");
					bw.write(Integer.toString(revisedPositionCounts.get(i)[1]));
					bw.write("\t");
					bw.write(Integer.toString(revisedPositionCounts.get(i)[2]));
					bw.newLine();
				}
			}
			System.out.println("Saved revised counts.");
		
		}
		
		public void writeBCStats(BufferedWriter bw) throws IOException{
			
			bw.write("NumberDuplicatedBarcodes\tFrequencyNumDuplicated\tMeanNumDup\tStandardDeviationNumDup");
			bw.newLine();
			for (int i = 0; i < histNumDupBC.length; i++){				
				bw.write(Integer.toString(i+1));
				bw.write("\t");
				bw.write(Integer.toString(histNumDupBC[i]));
				bw.write("\t");
				bw.write(Double.toString(statsByNumBCDup[i].getMean()));
				bw.write("\t");
				bw.write(Double.toString(statsByNumBCDup[i].getStandardDeviation()));
				bw.newLine();								
			}
			System.out.println("Saved stats.");
		
		}
		
		public void writeRevisedBCStats(BufferedWriter bw, int threshold2) throws IOException{
			
			bw.write("AdjustedNumberDuplicatedBarcodes\tFrequencyAdjustedNumDuplicated\tMeanNumDup\tStandardDeviationNumDup");
			bw.newLine();
			for (int i = 0; i < revisedStatsByAdjustedNumBCDup.length; i++){				
				bw.write(Integer.toString(i+threshold2));
				bw.write("\t");
				bw.write(Integer.toString((int) revisedStatsByAdjustedNumBCDup[i].getN()));
				bw.write("\t");
				bw.write(Double.toString(revisedStatsByAdjustedNumBCDup[i].getMean()));
				bw.write("\t");
				bw.write(Double.toString(revisedStatsByAdjustedNumBCDup[i].getStandardDeviation()));
				bw.newLine();								
			}
			System.out.println("Saved revised stats.");
		}
	}
}

	




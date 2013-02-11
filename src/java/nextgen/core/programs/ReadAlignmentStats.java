package nextgen.core.programs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.logging.Logger;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.fastq.FastqSequence;

import net.sf.picard.util.Log;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.programs.MarkDuplicatesFromNameUMI.MarkDuplicatesResults;

/**
 * Utility to compute general alignment stats. 
 * Given two sample descriptors, one pointing to fastq, the other to sample bam files
 * This program reports basic data about Reads en each library, reads aligned, multimapper reads and PCR duplicates (as marked in the alignment file).
 * @author mgarber
 *
 */
public class ReadAlignmentStats {
	
	LinkedHashMap<String, String> sampleFastqs;
	LinkedHashMap<String, String> sampleBams;
	private int maxThreads;
	private HashMap<String, Integer> fastqCounts = new HashMap<String, Integer>();
	private HashMap<String, SAMStats> samStats = new HashMap<String, ReadAlignmentStats.SAMStats>();
	
	private static final Log log = Log.getInstance(ReadAlignmentStats.class);
	static final String USAGE = "Utility to compute general alignment stats. " +
				"\n\t-sampleFastq <A two  column tab-separated file, first column the sample name, second column the path to the fastq file if data is paired just specify the first pair fastq file>" +
				"\n\t-sampleAln <A two column tab-separated file, first column the sample name, second column the path to the alignment file>" +
				"\n\t-maxThreads <Maximum number of threads >" +
				"\n\t-out <Output stats file or standard output if none are provided>" +
				"\n";

	
	public ReadAlignmentStats(String fastqTableFile, String bamTableFile) throws IOException {
		sampleFastqs = readTable(fastqTableFile);
		sampleBams   = readTable(bamTableFile);
		
		checkCompatibility();
	}

	

	private LinkedHashMap<String, String> readTable(String tableFile) throws IOException {
		LinkedHashMap<String, String> table = new LinkedHashMap<String, String>();
		
		BufferedReader br = new BufferedReader(new FileReader(tableFile));
		String line  = null;
		try {
			int lineNum = 0;
			while( (line = br.readLine()) != null) {
				lineNum++;
				line = line.trim();
				if(line.startsWith("#") || line.length() == 0) {
					continue;
				}
				String [] data = line.split("\t");
				if(data.length < 2) {
					throw new IllegalArgumentException("Table file must have at least two columns, but line " + lineNum + ":  " + line + " has less than 2");
				}
				
				table.put(data[0], data[1]);
				
			}
			
		} finally {
			br.close();
		}
		return table;
	}


	private void checkCompatibility() throws IllegalArgumentException{
		if(sampleFastqs.size() != sampleBams.size()) {
			throw new IllegalArgumentException("Size of the sampleFastq and sampleAln are different, both files must have the same number of samples");
		}
		
		List<String> sortedFastqSamples = new ArrayList<String>(sampleFastqs.keySet()); 
		List<String> sortedBamSamples   = new ArrayList<String>(sampleBams.keySet());
		
		Collections.sort(sortedFastqSamples);
		Collections.sort(sortedBamSamples);
		
		for (int i = 0; i < sortedFastqSamples.size(); i++) {
			if(!sortedFastqSamples.get(i).equals(sortedBamSamples.get(i))) {
				throw new IllegalArgumentException("Samples must be the same for both sampleFastq and sampleAln, however the keys at position " + i + " (after sorting keys) were " + sortedFastqSamples.get(i) + " and " + sortedBamSamples.get(i) + ", respectively");
			}
		}
		
	}
	
	protected int getMaxThreads() {
		return maxThreads;
	}


	protected void setMaxThreads(int maxThreads) {
		this.maxThreads = maxThreads;
	}




	public static void main (String [] args) throws IOException {
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "default");
		String fastqTableFile = argMap.getMandatory("sampleFastq");
		String bamTableFile = argMap.getMandatory("sampleAln");
		
		ReadAlignmentStats stats = new ReadAlignmentStats(fastqTableFile, bamTableFile);
		stats.setMaxThreads(argMap.containsKey("maxThreads") ? argMap.getInteger("maxThreads") : 1);
		
		stats.computeReadNumber();
		
		stats.computeAlignmentStats();
		BufferedWriter bw = argMap.getOutputWriter();
		bw.write("sample\tread_num\ttotalAlignments\ttotalAlignedReads\ttotalDuplicates\ttotalWeightedAlignments");
		bw.newLine();

		for (String sample : stats.sampleBams.keySet()) {
			int reads = stats.fastqCounts.get(sample);
			SAMStats bamStats = stats.samStats.get(sample);
			log.debug(sample + " -- " + reads + " " + bamStats);
			bw.write("sample");
			bw.write("\t");
			bw.write(reads+"\t");
			bw.write(bamStats.toString());
			bw.newLine();
		}
		
		bw.flush();
		bw.close();
	}
	
	private void computeAlignmentStats() {
		List<Thread> threads= new ArrayList<Thread>();
		for (String sample : sampleBams.keySet()) {
			SAMStatsRunner samStats = new SAMStatsRunner(sample, sampleBams.get(sample));
			
			samStats.addListner(this);
			Thread fastqReadCounterThread = new Thread(samStats);
			threads.add(fastqReadCounterThread);
			fastqReadCounterThread.start();	
			log.debug("Started SAM stats thread for sample " + sample);
		}
		for(Thread t : threads) {
			try {
				t.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}



	public synchronized void doneCounting(String sample, String fastqFile, int num) {
		fastqCounts.put(sample, num);
	}
	
	private void computeReadNumber() {
		List<Thread> threads= new ArrayList<Thread>();
		for (String sample : sampleFastqs.keySet()) {
			FastqReadCounter fastqReadCounter = new FastqReadCounter(sample, sampleFastqs.get(sample));
			fastqReadCounter.addListner(this);
			Thread fastqReadCounterThread = new Thread(fastqReadCounter);
			fastqReadCounterThread.start();	
			log.debug("Started fastq counter thread for sample " + sample);
			threads.add(fastqReadCounterThread);
		}	
		
		for(Thread t : threads) {
			try {
				t.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public void doneAlnAnalysis(SAMStats stats) {
		log.debug("got stats for sample "+ stats.sample + " - " + stats );
		samStats.put(stats.sample, stats);
		
	}
	
	public static class FastqReadCounter implements Runnable {

		private String sample;
		private String fastqFile;
		private ReadAlignmentStats listener;

		public FastqReadCounter(String sample, String fastqFile) {
			this.sample = sample;
			this.fastqFile = fastqFile;
		}

		public void addListner(ReadAlignmentStats listner) {
			this.listener = listner;
			
		}
		
		
		public void run() {
			int num = 0;
			FastqParser fqp = new FastqParser();
			try {
				fqp.start(new File(fastqFile));
				while(fqp.hasNext()) {
					FastqSequence seq = fqp.next();
					num++;
				}
				fqp.close();
			}  catch (IOException e) {
				// TODO Auto-generated catch block
				throw new RuntimeException("Error in fastq counting thread parsing file " + fastqFile);
			}
			
			listener.doneCounting(sample, fastqFile, num);
		}
		
	}

	public static class SAMStatsRunner implements Runnable {

		private String sample;
		private String aln;
		private ReadAlignmentStats listener;

		public SAMStatsRunner(String sample, String samAln) {
			this.sample = sample;
			this.aln = samAln;
		}

		public void addListner(ReadAlignmentStats listner) {
			this.listener = listner;
			
		}

		public void run() {
			int num = 0;
			double weightedNum = 0;
			int numDuplicated = 0;
			HashSet<String> readNames = new HashSet<String>();
			log.debug("opening " + aln);
			SAMFileReader reader = new SAMFileReader(new File(aln));
			SAMRecordIterator sri = reader.iterator();
			while(sri.hasNext()) {
				SAMRecord samR = sri.next();
				num++;
				if (samR.getDuplicateReadFlag()) {
					numDuplicated++;
				}
				int multiplicity = samR.getIntegerAttribute("NH");
				weightedNum += weightedNum + 1/(double) multiplicity;
				readNames.add(samR.getReadName());
			}  
			reader.close();
			SAMStats stats = new SAMStats(sample, aln , num, weightedNum, numDuplicated, readNames.size());
			listener.doneAlnAnalysis(stats );
		}
		
	}
	
	public static  class SAMStats {
		public SAMStats(String sample, String aln, int num, double weightedNum,
				int numDuplicated, int totalReads) {
			this.totalAlignments = num;
			this.totalAlignedReads = totalReads;
			this.totalDuplicates = numDuplicated;
			this.totalWeigthedAlignments = weightedNum;
			this.sample = sample;
			this.alignmentFile = aln;
		}

		public String toString() {
			StringBuilder sb = new StringBuilder ();
			sb.append(totalAlignments).append("\t")
				.append(totalAlignedReads).append("\t")
				.append(totalDuplicates).append("\t")
				.append(totalWeigthedAlignments);
			return sb.toString();
		}
		int totalAlignments;
		double totalWeigthedAlignments;
		int totalAlignedReads;
		int totalDuplicates;
		String sample;
		String alignmentFile;
	}

}

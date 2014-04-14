package nextgen.sequentialbarcode;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.OGSJob;
import nextgen.core.pipeline.util.FastqUtils;
import nextgen.core.pipeline.util.OGSUtils;
import nextgen.sequentialbarcode.readlayout.AnySequence;
import nextgen.sequentialbarcode.readlayout.Barcode;
import nextgen.sequentialbarcode.readlayout.BarcodeSet;
import nextgen.sequentialbarcode.readlayout.FixedSequence;
import nextgen.sequentialbarcode.readlayout.ReadLayout;
import nextgen.sequentialbarcode.readlayout.ReadLayoutFactory;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.fastq.FastqSequence;

/**
 * Barcode identification and other analyses
 * @author prussell
 *
 */
public class BarcodeAnalysis {

	public static Logger logger = Logger.getLogger(BarcodeAnalysis.class.getName());
	private static Session drmaaSession;
	
	/**
	 * For RNA-DNA 3D barcoding method
	 * Identify barcodes in reads and write to a table
	 * @param fastq Fastq file
	 * @param readLength Read length
	 * @param oddBarcodeTableFile Table of odd barcodes (line format: barcode_ID   barcode_sequence)
	 * @param evenBarcodeTableFile Table of even barcodes (line format: barcode_ID   barcode_sequence)
	 * @param rpm RPM sequence
	 * @param maxMismatchBarcode Max number of mismatches when matching barcodes to reads
	 * @param maxMismatchRpm Max number of mismatches when finding RPM in reads
	 * @param enforceOddEven Require odd and even barcodes to alternate
	 * @param outFile Output table
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private static void findRead2BarcodesRnaDna3D(String fastq, int readLength, String oddBarcodeTableFile, String evenBarcodeTableFile, String rpm, int maxMismatchBarcode, int maxMismatchRpm, boolean enforceOddEven, String outFile) throws IOException {
		findRead2BarcodesRnaDna3D(fastq, readLength, oddBarcodeTableFile, evenBarcodeTableFile, rpm, maxMismatchBarcode, maxMismatchRpm, enforceOddEven, outFile, false);
	}
	
	/**
	 * For RNA-DNA 3D barcoding method
	 * Identify barcodes in reads and write to a table
	 * Split the fastq file into several smaller files and batch out the writing of the table
	 * @param barcodeAnalysisJar Jar file for this main method
	 * @param numSplitFastqs Number of smaller fastqs to divide into
	 * @param fastq Fastq file
	 * @param readLength Read length
	 * @param oddBarcodeTableFile Table of odd barcodes (line format: barcode_ID   barcode_sequence)
	 * @param evenBarcodeTableFile Table of even barcodes (line format: barcode_ID   barcode_sequence)
	 * @param rpm RPM sequence
	 * @param maxMismatchBarcode Max number of mismatches when matching barcodes to reads
	 * @param maxMismatchRpm Max number of mismatches when finding RPM in reads
	 * @param enforceOddEven Require odd and even barcodes to alternate
	 * @param outFile Output table
	 * @param verbose Verbose table output
	 * @throws IOException
	 * @throws DrmaaException
	 * @throws InterruptedException
	 */
	private static void divideFastqAndFindRead2BarcodesRnaDna3D(String barcodeAnalysisJar, int numSplitFastqs, String fastq, int readLength, String oddBarcodeTableFile, String evenBarcodeTableFile, String rpm, int maxMismatchBarcode, int maxMismatchRpm, boolean enforceOddEven, String outFile, boolean verbose, String email) throws IOException, DrmaaException, InterruptedException {
		Collection<String> splitFastqs = FastqUtils.divideFastqFile(fastq, numSplitFastqs);
		Map<String, String> splitTables = new TreeMap<String, String>();
		int i = 0;
		for(String fq : splitFastqs) {
			splitTables.put(fq, outFile + "." + i);
			i++;
		}
		Collection<Job> jobs = new ArrayList<Job>();
		for(String fq : splitTables.keySet()) {
			String cmmd = "java -jar -Xmx25g -Xms15g -Xmn10g " + barcodeAnalysisJar;
			cmmd += " -ob3d " + splitTables.get(fq);
			cmmd += " -f2 " + fq;
			cmmd += " -rl " + readLength;
			cmmd += " -ob " + oddBarcodeTableFile;
			cmmd += " -eb " + evenBarcodeTableFile;
			cmmd += " -rpm " + rpm;
			cmmd += " -mmb " + maxMismatchBarcode;
			cmmd += " -mmr " + maxMismatchRpm;
			cmmd += " -oe " + enforceOddEven;
			cmmd += " -v " + verbose;
			String jobName = "OGS_job_" + fq;
			OGSJob job = new OGSJob(drmaaSession, cmmd, true, jobName, email);
			job.submit();
			jobs.add(job);
		}
		JobUtils.waitForAll(jobs);
		FileWriter w = new FileWriter(outFile);
		for(String fq : splitTables.keySet()) {
			String table = splitTables.get(fq);
			FileReader r = new FileReader(table);
			BufferedReader b = new BufferedReader(r);
			while(b.ready()) {
				w.write(b.readLine() + "\n");
			}
			r.close();
			b.close();
			File f = new File(table);
			f.delete();
		}
		w.close();
		for(String fq : splitFastqs) {
			File f = new File(fq);
			f.delete();
		}
	}
	
	
	/**
	 * For RNA-DNA 3D barcoding method
	 * Identify barcodes in reads and write to a table
	 * @param fastq Fastq file
	 * @param readLength Read length
	 * @param oddBarcodeTableFile Table of odd barcodes (line format: barcode_ID   barcode_sequence)
	 * @param evenBarcodeTableFile Table of even barcodes (line format: barcode_ID   barcode_sequence)
	 * @param rpm RPM sequence
	 * @param maxMismatchBarcode Max number of mismatches when matching barcodes to reads
	 * @param maxMismatchRpm Max number of mismatches when finding RPM in reads
	 * @param enforceOddEven Require odd and even barcodes to alternate
	 * @param outFile Output table
	 * @param verbose Verbose table output
	 * @throws IOException
	 */
	private static void findRead2BarcodesRnaDna3D(String fastq, int readLength, String oddBarcodeTableFile, String evenBarcodeTableFile, String rpm, int maxMismatchBarcode, int maxMismatchRpm, boolean enforceOddEven, String outFile, boolean verbose) throws IOException {
		logger.info("");
		logger.info("Identifying read2 barcodes for RNA-DNA-3D and writing to table " + outFile + "...");
		Collection<Barcode> oddBarcodes = Barcode.createBarcodesFromTable(oddBarcodeTableFile, maxMismatchBarcode);
		Collection<Barcode> evenBarcodes = Barcode.createBarcodesFromTable(evenBarcodeTableFile, maxMismatchBarcode);
		FileWriter w = new FileWriter(outFile);
		ReadLayout layout = ReadLayoutFactory.getRead2LayoutRnaDna3D(evenBarcodes, oddBarcodes, 6, rpm, readLength, maxMismatchBarcode, maxMismatchRpm, enforceOddEven);
		FastqParser iter = new FastqParser();
		iter.start(new File(fastq));
		int numDone = 0;
		while(iter.hasNext()) {
			numDone++;
			if(numDone % 10000 == 0) {
				logger.info("Finished " + numDone + " reads.");
			}
			FastqSequence record = iter.next();
			if(record == null) {
				continue;
			}
			record.removeAtSymbolFromName();
			String seq = record.getSequence();
			String name = record.getName();
			String line = StringParser.firstField(name) + "\t";
			if(layout.getMatchedElements(seq) != null) {
				BarcodedFragment f = new BarcodedFragmentImpl(name, null, seq, null, layout, maxMismatchBarcode);
				BarcodeSequence barcodes = f.getBarcodes();
				if(verbose) line += barcodes.getNumBarcodes() + "\t";
				line += barcodes.toString() + "\t";
				if(verbose) line += seq + "\t";
				w.write(line + "\n");
				f = null;
				seq = null;
				name = null;
				line = null;
				record = null;
				continue;
			}
			seq = null;
			name = null;
			line = null;
			record = null;
		}
		w.close();
	}
	
	
	/**
	 * For RNA-DNA 3D barcoding method
	 * Count instances of each barcode and print totals
	 * @param fastq Fastq file
	 * @param readLength Read length
	 * @param oddBarcodeTableFile Table of odd barcodes (line format: barcode_ID   barcode_sequence)
	 * @param evenBarcodeTableFile Table of even barcodes (line format: barcode_ID   barcode_sequence)
	 * @param rpm RPM sequence
	 * @param maxMismatchBarcode Max number of mismatches when matching barcodes to reads
	 * @param maxMismatchRpm Max number of mismatches when finding RPM in reads
	 * @throws IOException
	 */
	private static void countRead2BarcodesRnaDna3D(String fastq, int readLength, String oddBarcodeTableFile, String evenBarcodeTableFile, String rpm, int maxMismatchBarcode, int maxMismatchRpm) throws IOException {
		logger.info("");
		logger.info("Counting read2 barcodes for RNA-DNA-3D...");
		Collection<Barcode> oddBarcodes = Barcode.createBarcodesFromTable(oddBarcodeTableFile, maxMismatchBarcode);
		Collection<Barcode> evenBarcodes = Barcode.createBarcodesFromTable(evenBarcodeTableFile, maxMismatchBarcode);
		Map<Barcode, Integer> barcodeCounts = new TreeMap<Barcode, Integer>();
		for(Barcode b : oddBarcodes) {
			barcodeCounts.put(b, Integer.valueOf(0));
		}
		for(Barcode b : evenBarcodes) {
			barcodeCounts.put(b, Integer.valueOf(0));
		}
		ReadLayout layout = ReadLayoutFactory.getRead2LayoutRnaDna3D(evenBarcodes, oddBarcodes, 6, rpm, readLength, maxMismatchBarcode, maxMismatchRpm, false);
		FastqParser iter = new FastqParser();
		iter.start(new File(fastq));
		int numDone = 0;
		while(iter.hasNext()) {
			FastqSequence record = iter.next();
			String seq = record.getSequence();
			String name = record.getName();
			if(layout.getMatchedElements(seq) != null) {
				BarcodedFragment f = new BarcodedFragmentImpl(name, null, seq, null, layout, maxMismatchBarcode);
				BarcodeSequence barcodes = f.getBarcodes();
				for(Barcode b : barcodes.getBarcodes()) {
					barcodeCounts.put(b, Integer.valueOf(barcodeCounts.get(b).intValue() + 1));
				}
			}
			if(numDone > 0 && numDone % 10000 == 0) {
				logger.info("");
				logger.info("Finished " + numDone + " reads.");
				for(Barcode b : barcodeCounts.keySet()) {
					System.out.println(b.getId() + "\t" + b.getSequence() + "\t" + barcodeCounts.get(b).intValue());
				}
			}
			numDone++;
		}
	}
	

	
	
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static void main(String[] args) throws IOException, DrmaaException, InterruptedException {
		
		drmaaSession = OGSUtils.getDrmaaSession();
		
		CommandLineParser p = new CommandLineParser();
		p.addBooleanArg("-d", "Debug logging on", false, false);
		p.addStringArg("-ob3d", "Output file for barcode identification for RNA-DNA-3D (requires -f2, -rl, -ob, -eb, -rpm, -mmb, -mmr)", false, null);
		p.addStringArg("-f2", "Read 2 fastq file", true);
		p.addIntArg("-rl", "Read length", true);
		p.addStringArg("-ob", "Odd barcode table file (format: barcode_id	barcode_seq)", false, null);
		p.addStringArg("-eb", "Even barcode table file (format: barcode_id	barcode_seq)", false, null);
		p.addStringArg("-rpm", "RPM sequence", false, null);
		p.addIntArg("-mmb", "Max mismatches in barcode", true);
		p.addIntArg("-mmr", "Max mismatches in RPM", true);
		p.addBooleanArg("-oe", "Enforce odd/even alternation for barcodes", false, false);
		p.addBooleanArg("-cb", "Count barcodes", false, false);
		p.addBooleanArg("-v", "Verbose output table for barcode identification", false, false);
		p.addBooleanArg("-bt", "Batch out writing of barcode table", false, false);
		p.addStringArg("-j", "This jar file for batched jobs", false, null);
		p.addIntArg("-nf", "Number of fastq files to divide into if batching", false, 20);
		p.addStringArg("-e", "Email address for OGS jobs", false, null);
		p.parse(args);
		if(p.getBooleanArg("-d")) {
			ReadLayout.logger.setLevel(Level.DEBUG);
			BarcodeSet.logger.setLevel(Level.DEBUG);
			ReadLayoutFactory.logger.setLevel(Level.DEBUG);
			BarcodeSequence.logger.setLevel(Level.DEBUG);
			Barcode.logger.setLevel(Level.DEBUG);
			BarcodedFragmentImpl.logger.setLevel(Level.DEBUG);
			FixedSequence.logger.setLevel(Level.DEBUG);
			BarcodedDNAFragment.logger.setLevel(Level.DEBUG);
			BarcodedRNAFragment.logger.setLevel(Level.DEBUG);
			AnySequence.logger.setLevel(Level.DEBUG);
			BarcodeAnalysis.logger.setLevel(Level.DEBUG);
		}
		
		String outRD3 = p.getStringArg("-ob3d");
		String email = p.getStringArg("-e");
		String read2fastq = p.getStringArg("-f2");
		int readLength = p.getIntArg("-rl");
		String oddBarcodeList = p.getStringArg("-ob");
		String evenBarcodeList = p.getStringArg("-eb");
		String rpm = p.getStringArg("-rpm");
		int maxMismatchBarcode = p.getIntArg("-mmb");
		int maxMismatchRpm = p.getIntArg("-mmr");
		boolean enforceOddEven = p.getBooleanArg("-oe");
		boolean countBarcodes = p.getBooleanArg("-cb");
		boolean verbose = p.getBooleanArg("-v");
		boolean batch = p.getBooleanArg("-bt");
		String jar = p.getStringArg("-j");
		int numFastq = p.getIntArg("-nf");
		
		if(numFastq < 1) {
			throw new IllegalArgumentException("Number of fastq files must be at least 1.");
		}
		
		if(countBarcodes) {
			countRead2BarcodesRnaDna3D(read2fastq, readLength, oddBarcodeList, evenBarcodeList, rpm, maxMismatchBarcode, maxMismatchRpm);
		}
		
		if(outRD3 != null) {
			if(batch) {
				if(jar == null) {
					throw new IllegalArgumentException("Must provide jar file");
				}
				divideFastqAndFindRead2BarcodesRnaDna3D(jar, numFastq, read2fastq, readLength, oddBarcodeList, evenBarcodeList, rpm, maxMismatchBarcode, maxMismatchRpm, enforceOddEven, outRD3, verbose, email);
			} else {
				findRead2BarcodesRnaDna3D(read2fastq, readLength, oddBarcodeList, evenBarcodeList, rpm, maxMismatchBarcode, maxMismatchRpm, enforceOddEven, outRD3, verbose);
			}
		}
		
		logger.info("");
		logger.info("All done.");
				
	}

}

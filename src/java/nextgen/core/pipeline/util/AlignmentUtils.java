package nextgen.core.pipeline.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;


import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.LSFJob;
import nextgen.core.job.OGSJob;
import nextgen.core.pipeline.Scheduler;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import broad.core.parser.StringParser;

/**
 * @author prussell
 *
 */
public class AlignmentUtils {
		
	private static Logger logger = Logger.getLogger(AlignmentUtils.class.getName());
	
	/**
	 * Check whether the string represents a valid sam file line
	 * @param line The string
	 * @return True if the line is a sam header line or a valid alignment line, false otherwise
	 */
	public static boolean isSamLine(String line) {
		StringParser p = new StringParser();
		return isSamLine(p, line);
	}

	/**
	 * Check whether the string represents a valid sam file line using an exisiting StringParser object
	 * @param p The StringParser object
	 * @param line The string
	 * @return True if the line is a sam header line or a valid alignment line, false otherwise
	 */		
	@SuppressWarnings("unused")
	public static boolean isSamLine(StringParser p, String line) {
		p.parse(line);
		if(p.getFieldCount() == 0) return false;
		if(p.getFieldCount() >= 11) {
			// Line might be an alignment line
			// Check if correct fields are ints
			try {
				int i1 = p.asInt(1);
				int i3 = p.asInt(3);
				int i4 = p.asInt(4);
				int i7 = p.asInt(7);
				int i8 = p.asInt(8);
				return true;
			} catch(NumberFormatException e) {
				// The value in a sam int position is not an int
				return false;
			}
		}
		// Check if this is a header line
		String first = p.asString(0);
		if(first.equals("@HD") || first.equals("@SQ") || first.equals("@RG") || first.equals("@PG") || first.equals("@CO")) return true;
		// Line is neither a header line nor an alignment line
		return false;
	}
	

	
	/**
	 * Count total number of mapped reads in a sam or bam file
	 * @param samtoolsExecutable Samtools executable file
	 * @param samFile The sam file
	 * @param logDir Directory to write log information to
	 * @param samFormat True if sam format, false if bam format
	 * @return The alignment count obtained by running the command "samtools view -c -F 4"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static int countAlignments(String samtoolsExecutable, String samFile, String logDir, boolean samFormat) throws IOException, InterruptedException {
		return countAlignments(samtoolsExecutable, samFile, logDir, false, samFormat);
	}
	
	/**
	 * Count total number of mapped or unmapped reads in a sam file
	 * @param samtoolsExecutable Samtools executable file
	 * @param samFile The sam file
	 * @param logDir Directory to write log information to
	 * @param countUnmapped True if counting unmapped reads, false if counting mapped reads
	 * @param samFormat True if sam format, false if bam format
	 * @return The alignment count obtained by running the command "samtools view -c -F 4"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static int countAlignments(String samtoolsExecutable, String samFile, String logDir, boolean countUnmapped, boolean samFormat) throws IOException, InterruptedException {

		File dir = new File(logDir);
		@SuppressWarnings("unused")
		boolean madeDir = dir.mkdir();
		
		// Use samtools to count alignments
		String cmmd = samtoolsExecutable + " view -c ";
		if(samFormat) cmmd += " -S ";
		if(!countUnmapped) cmmd += " -F 4 ";
		else cmmd += " -f 4 ";
		cmmd += samFile;
		
		// Capture output of samtools process and read into int
		Process p = Runtime.getRuntime().exec(cmmd);
		InputStream is = p.getInputStream();
		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);
		String line = br.readLine();
		br.close();
		int count = 0;
		try {
			count = Integer.parseInt(line);
		} catch(NumberFormatException e) {
			logger.error("Command did not return valid integer result:");
			logger.error(cmmd);
			logger.error("Result:");
			logger.error(line);
			throw e;
		}
		is.close();
		isr.close();
		
		String type = countUnmapped ? "unmapped" : "mapped";
		logger.info("There are " + count + " " + type + " reads in sam file " + samFile + ".");
		
		return count;
		
	}

	
	
	/**
	 * Make bowtie2 index of fasta file
	 * @param fastaFile The fasta file
	 * @param outBtIndexBase Output index file without .bt2 extension
	 * @param bowtie2BuildExecutable The bowtie2-build executable
	 * @param bsubOutputDir Output directory for bsub file
	 * @param scheduler Scheduler
	 * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method. 
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void makeBowtie2Index(String fastaFile, String outBtIndexBase, String bowtie2BuildExecutable, String bsubOutputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		logger.info("");
		logger.info("Writing bowtie2 index for file " + fastaFile + " to files " + outBtIndexBase);
		File bsubDir = new File(bsubOutputDir);
		@SuppressWarnings("unused")
		boolean madeDir = bsubDir.mkdir();
		
		// Check if index already exists
		String f1 = outBtIndexBase + ".1.bt2";
		String f2 = outBtIndexBase + ".2.bt2";
		String f3 = outBtIndexBase + ".3.bt2";
		String f4 = outBtIndexBase + ".4.bt2";
		String r1 = outBtIndexBase + ".rev.1.bt2";
		String r2 = outBtIndexBase + ".rev.2.bt2";
		File f1file = new File(f1);
		File f2file = new File(f2);
		File f3file = new File(f3);
		File f4file = new File(f4);
		File r1file = new File(r1);
		File r2file = new File(r2);
		if(f1file.exists() && f2file.exists() && f3file.exists() && f4file.exists() && r1file.exists() && r2file.exists()) {
			logger.warn("Bowtie2 index files already exist. Not writing new files.");
			return;
		}
		
		String cmmd = bowtie2BuildExecutable + " " + fastaFile + " " + outBtIndexBase;
		logger.info("Submitting command " + cmmd);
		String jobID = Long.valueOf(System.currentTimeMillis()).toString();
		String output = bsubOutputDir + "/make_bowtie_index_" + jobID + ".bsub";
		switch(scheduler) {
			case LSF:
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, output, "week", 4);
				lsfJob.submit();
				lsfJob.waitFor();
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "bowtie2_index");
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                ogsJob.waitFor();
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		logger.info("Done creating bowtie2 index for file " + fastaFile);
	}
	
	/**
	 * Run bowtie2 with unpaired reads and default max insert size
	 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
	 * @param options 
	 * @param reads Fastq file containing unpaired reads
	 * @param outSamFile Output sam file
	 * @param outUnalignedFastq Output fastq file for unaligned reads or null
	 * @param bowtie2Executable Bowtie2 executable
	 * @param bsubOutDir Output directory for bsub file
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @return The job object
	 * @throws DrmaaException 
	 */
	public static Job runBowtie2(String bowtie2IndexBase, Map<String, String> options, String reads, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return runBowtie2(bowtie2IndexBase, options, reads, null, outSamFile, outUnalignedFastq, bowtie2Executable, bsubOutDir, false, scheduler, drmaaSession);
	}

	
	/**
	 * Run bowtie2 with paired reads and default max insert size
	 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
	 * @param options Bowtie2 option flags and values
	 * @param read1Fastq Fastq file containing read1 if paired or single end reads if unpaired
	 * @param read2Fastq Fastq file containing read2 if paired or null if unpaired
	 * @param outSamFile Output sam file
	 * @param outUnalignedFastq Output fastq file for unaligned reads or null
	 * @param bowtie2Executable Bowtie2 executable
	 * @param bsubOutDir Output directory for bsub file
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @return The job object
	 * @throws DrmaaException 
	 */
	public static Job runBowtie2(String bowtie2IndexBase, Map<String, String> options, String read1Fastq, String read2Fastq, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return runBowtie2(bowtie2IndexBase, options, read1Fastq, read2Fastq, outSamFile, outUnalignedFastq, bowtie2Executable, bsubOutDir, true, scheduler, drmaaSession);
	}
	
	/**
	 * Run bowtie2 with paired reads and specified max insert size
	 * @param bowtie2IndexBase Index filename prefix (minus trailing .X.bt2)
	 * @param options Bowtie2 option flags and values
	 * @param read1Fastq Fastq file containing read1 if paired or single end reads if unpaired
	 * @param read2Fastq Fastq file containing read2 if paired or null if unpaired
	 * @param outSamFile Output sam file
	 * @param outUnalignedFastq Output fastq file for unaligned reads or null
	 * @param bowtie2Executable Bowtie2 executable
	 * @param bsubOutDir Output directory for bsub files
	 * @param readsPaired Whether the reads are paired
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @return The bsub job ID
	 * @throws DrmaaException 
	 */
	public static Job runBowtie2(String bowtie2IndexBase, Map<String, String> options, String read1Fastq, String read2Fastq, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir, boolean readsPaired, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {

		File bsubOutDirectory = new File(bsubOutDir);
		@SuppressWarnings("unused")
		boolean madeDir = bsubOutDirectory.mkdir();
		
		String executable = bowtie2Executable + " ";
		String index = "-x " + bowtie2IndexBase + " ";
		
		String reads;
		if(read2Fastq != null) reads = "-1 " + read1Fastq + " -2 " + read2Fastq + " ";
		else reads = "-U " + read1Fastq + " ";

		String optionString = "";
		for(String flag : options.keySet()) {
			optionString += flag + " " + options.get(flag) + " ";
		}
		
		if(outUnalignedFastq != null) {
			if(!readsPaired) optionString += "--un " + outUnalignedFastq + " ";
			else optionString += "--un-conc " + outUnalignedFastq + " ";
		}
		
		String sam = "-S " + outSamFile;
		
		String cmmd = executable + optionString + index + reads + sam;
		
		logger.info("");
		logger.info("Running bowtie2 command:");
		logger.info(cmmd);
		
		switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				String output = bsubOutDir + "/run_bowtie_" + jobID + ".bsub";
				logger.info("Writing bsub output to file " + output);
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, output, "week", 4);
				lsfJob.submit();
				logger.info("Job ID is " + jobID);
				return lsfJob;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "bowtie2");
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                return ogsJob;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
		}
		
	}
	
	/**
	 * Run tophat for single end reads against genome and count alignments
	 * @param tophatExecutable Tophat executable
	 * @param samtools Samtools executable
	 * @param leftFastqs Map of sample name to read1
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDirsPerSample Directory to write tophat output to, by sample name
	 * @param tophatOutputBamPerSample Bam file to write by sample name
	 * @param tophatBamFinalPath Final bam file for tophat output, by sample name
	 * @param genomeIndex Genome fasta index
	 * @return Map of sample name to bam file generated from this TopHat run
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static Map<String, String> runTophat(String tophatExecutable, String samtools, Map<String, String> leftFastqs, Map<String, String> tophatOptions, Map<String, String> tophatDirsPerSample, Map<String, String> tophatOutputBamPerSample, Map<String, String> tophatBamFinalPath, String genomeIndex) throws IOException, InterruptedException {
		return runTophat(tophatExecutable, samtools, leftFastqs, tophatOptions, tophatDirsPerSample, tophatOutputBamPerSample, tophatBamFinalPath, genomeIndex);
	}
	
	/**
	 * Run tophat against genome and count alignments
	 * @param tophatExecutable Tophat executable
	 * @param samtools Samtools executable
	 * @param leftFastqs Map of sample name to read1
	 * @param rightFastqs Map of sample name to read2 for samples with paired reads, or null if no samples are paired
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDirsPerSample Directory to write tophat output to, by sample name
	 * @param tophatOutputBamPerSample Bam file to write by sample name
	 * @param tophatBamFinalPath Final bam file for tophat output, by sample name
	 * @param genomeIndex Genome fasta index
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Map of sample name to bam file generated from this TopHat run
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static Map<String, String> runTophat(String tophatExecutable, String samtools, Map<String, String> leftFastqs, Map<String, String> rightFastqs, Map<String, String> tophatOptions, Map<String, String> tophatDirsPerSample, Map<String, String> tophatOutputBamPerSample, Map<String, String> tophatBamFinalPath, String genomeIndex, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return runTophat(tophatExecutable, samtools, leftFastqs, rightFastqs, tophatOptions, tophatDirsPerSample, tophatOutputBamPerSample, tophatBamFinalPath, genomeIndex, ".", scheduler, drmaaSession);
	}
	
	/**
	 * Run tophat against genome and count alignments
	 * @param tophatExecutable Tophat executable
	 * @param samtools Samtools executable
	 * @param sampleName Sample name
	 * @param leftFastq Left fastq
	 * @param rightFastq Right fastq
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDir Tophat directory
	 * @param tophatBamOutput Bam file to write
	 * @param finalTophatBam Final bam file; skip if already exists
	 * @param genomeIndex Genome fasta index
	 * @param outputDir Output directory
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Bam file generated from this TopHat run
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static String runTophat(String tophatExecutable, String samtools, String sampleName, String leftFastq, String rightFastq, Map<String, String> tophatOptions, String tophatDir, String tophatBamOutput, String finalTophatBam, String genomeIndex, String outputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return runTophat(tophatExecutable, samtools, sampleName, leftFastq, rightFastq, tophatOptions, tophatDir, tophatBamOutput, finalTophatBam, genomeIndex, "week", outputDir, scheduler, drmaaSession);
	}
	
	/**
	 * Run tophat against genome and count alignments
	 * @param tophatExecutable Tophat executable
	 * @param samtools Samtools executable
	 * @param leftFastqs Map of sample name to read1
	 * @param rightFastqs Map of sample name to read2 for samples with paired reads, or null if no samples are paired
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDirsPerSample Directory to write tophat output to, by sample name
	 * @param tophatOutputBamPerSample Bam file to write by sample name
	 * @param tophatBamFinalPath Final bam file for tophat output, by sample name
	 * @param genomeIndex Genome fasta index
	 * @param outputDir Output directory
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Map of sample name to bam file generated from this TopHat run
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static Map<String, String> runTophat(String tophatExecutable, String samtools, Map<String, String> leftFastqs, Map<String, String> rightFastqs, Map<String, String> tophatOptions, Map<String, String> tophatDirsPerSample, Map<String, String> tophatOutputBamPerSample, Map<String, String> tophatBamFinalPath, String genomeIndex, String outputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return runTophat(tophatExecutable, samtools, leftFastqs, rightFastqs, tophatOptions, tophatDirsPerSample, tophatOutputBamPerSample, tophatBamFinalPath, genomeIndex, "week", outputDir, scheduler, drmaaSession);
	}
	
	/**
	 * Run tophat against genome and count alignments
	 * @param tophatExecutable Tophat executable
	 * @param samtools Samtools executable
	 * @param leftFastqs Map of sample name to read1
	 * @param rightFastqs Map of sample name to read2 for samples with paired reads, or null if no samples are paired
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDirsPerSample Directory to write tophat output to, by sample name
	 * @param tophatOutputBamPerSample Bam file to write by sample name
	 * @param tophatBamFinalPath Final bam file for tophat output, by sample name
	 * @param genomeIndex Genome fasta index
	 * @param queueName LSF queue name
	 * @param outputDir Output directory
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Map of sample name to bam file generated from this TopHat run
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static Map<String, String> runTophat(String tophatExecutable, String samtools, Map<String, String> leftFastqs, Map<String, String> rightFastqs, Map<String, String> tophatOptions, Map<String, String> tophatDirsPerSample, Map<String, String> tophatOutputBamPerSample, Map<String, String> tophatBamFinalPath, String genomeIndex, String queueName, String outputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return runTophat(tophatExecutable, samtools, leftFastqs.keySet(), leftFastqs, rightFastqs, tophatOptions, tophatDirsPerSample, tophatOutputBamPerSample, tophatBamFinalPath, genomeIndex, queueName, outputDir, scheduler, drmaaSession);
	}
	
	/**
	 * Run tophat against genome and count alignments
	 * @param tophatExecutable Tophat executable
	 * @param samtools Samtools executable
	 * @param sampleName Sample name
	 * @param leftFastq Left fastq
	 * @param rightFastq Right fastq
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDir Tophat directory
	 * @param tophatOutputBam Bam file to write
	 * @param finalTophatBam Final bam file; skip if already exists
	 * @param genomeIndex Genome fasta index
	 * @param queueName LSF queue name
	 * @param outputDir Output directory
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Bam file generated from this TopHat run
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static String runTophat(String tophatExecutable, String samtools, String sampleName, String leftFastq, String rightFastq, Map<String, String> tophatOptions, String tophatDir, String tophatOutputBam, String finalTophatBam, String genomeIndex, String queueName, String outputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Collection<String> sampleNames = new ArrayList<String>();
		sampleNames.add(sampleName);
		Map<String, String> leftFastqs = new TreeMap<String, String>();
		Map<String, String> rightFastqs = new TreeMap<String, String>();
		Map<String, String> tophatDirsPerSample = new TreeMap<String, String>();
		Map<String, String> tophatBamFinalPath = new TreeMap<String, String>();
		Map<String, String> tophatOutputBamPerSample = new TreeMap<String, String>();
		tophatOutputBamPerSample.put(sampleName, tophatOutputBam);
		leftFastqs.put(sampleName, leftFastq);
		rightFastqs.put(sampleName, rightFastq);
		tophatDirsPerSample.put(sampleName, tophatDir);
		tophatBamFinalPath.put(sampleName, finalTophatBam);
		return runTophat(tophatExecutable, samtools, sampleNames, leftFastqs, rightFastqs, tophatOptions, tophatDirsPerSample, tophatOutputBamPerSample, tophatBamFinalPath, genomeIndex, queueName, outputDir, scheduler, drmaaSession).values().iterator().next();
	}
	
	/**
	 * Run tophat against genome and count alignments
	 * @param tophatExecutable Tophat executable
	 * @param samtools Samtools executable
	 * @param sampleNames Collection of sample names
	 * @param leftFastqs Map of sample name to read1
	 * @param rightFastqs Map of sample name to read2 for samples with paired reads, or null if no samples are paired
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDirsPerSample Directory to write tophat output to, by sample name
	 * @param tophatOutputBamPerSample Bam file to write by sample name
	 * @param tophatBamFinalPath Final bam file for tophat output by sample name. Tophat is not run if this file already exists.
	 * @param genomeIndex Genome fasta index
	 * @param queueName LSF queue name
	 * @param outputDir Output directory
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Map of sample name to bam file generated from this TopHat run
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static Map<String, String> runTophat(String tophatExecutable, String samtools, Collection<String> sampleNames, Map<String, String> leftFastqs, Map<String, String> rightFastqs, Map<String, String> tophatOptions, Map<String, String> tophatDirsPerSample, Map<String, String> tophatOutputBamPerSample, Map<String, String> tophatBamFinalPath, String genomeIndex, String queueName, String outputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		
		Map<String, String> rtrn = new TreeMap<String, String>();
		
		// Identify tophat version
		boolean version2 = tophatExecutable.substring(tophatExecutable.length()-1).equals("2");
		
		ArrayList<Job> tophatJobs = new ArrayList<Job>();
		// Run tophat
		for(String sample : sampleNames) {
			File outdir = new File(tophatDirsPerSample.get(sample));
			outdir.mkdir();
			File tmpBamFile = new File(tophatOutputBamPerSample.get(sample));
			// Check if tophat has already been run
			if(tmpBamFile.exists()) {
				logger.warn("Alignment file " +tophatOutputBamPerSample.get(sample) + " already exists. Not rerunning alignment.");
				continue;
			}
			File finalBamFile = new File(tophatBamFinalPath.get(sample));
			// Check if tophat has already been run
			if(finalBamFile.exists()) {
				logger.warn("Alignment file " +tophatBamFinalPath.get(sample) + " already exists. Not rerunning alignment.");
				continue;
			}
			// Make tophat command
			String optionsString = "";
			for(String flag : tophatOptions.keySet()) {
				optionsString += flag + " " + tophatOptions.get(flag) + " ";
			}
			String tophatCmmd = tophatExecutable + " ";
			tophatCmmd += optionsString + " ";
			tophatCmmd += "--output-dir " + outdir + " "; 
			tophatCmmd += genomeIndex + " ";
			tophatCmmd += leftFastqs.get(sample) + " ";
			if(rightFastqs != null) {
				if(rightFastqs.containsKey(sample)) tophatCmmd += rightFastqs.get(sample);
			}
			logger.info("Writing tophat output for sample " + sample + " to directory " + outdir + ".");
			logger.info("Running tophat command: " + tophatCmmd);
			
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, tophatCmmd, outdir + "/tophat_" + jobID + ".bsub", queueName, 16);
				tophatJobs.add(lsfJob);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				lsfJob.submit();
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, tophatCmmd, "tophat");
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                tophatJobs.add(ogsJob);
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		// Wait for tophat jobs to finish
		logger.info("Waiting for tophat jobs to finish...");
		JobUtils.waitForAll(tophatJobs);
		logger.info("All samples done aligning to genome.");
		
		// Move all tophat bam files to one directory
		logger.info("");
		logger.info("Moving all tophat bam files to directory " + outputDir + "...");
		for(String sample : sampleNames) {
			String cmmd = "mv " + tophatDirsPerSample.get(sample) + "/accepted_hits.bam " + tophatOutputBamPerSample.get(sample);
			Process p = Runtime.getRuntime().exec(cmmd, null);
			p.waitFor();
			// Update current bam files
			rtrn.put(sample, tophatOutputBamPerSample.get(sample));
		}
		
		// Count aligned and unaligned reads
		// If tophat version 1, unzip unmapped fastq files and rename
		logger.info("");
		logger.info("Counting mapped and unmapped reads...");
		String countFile = outputDir + "/mapped_unmapped_count.out";
		String pctFile = outputDir + "/mapped_unmapped_percentage.out";
		
		// Check if files already exist
		File cFile = new File(countFile);
		File pFile = new File(pctFile);
		if(cFile.exists() && pFile.exists()) {
			logger.warn("Files " + countFile +" and " + pctFile +" already exist. Not recomputing mapping counts.");
		} else {
		
			FileWriter w = new FileWriter(countFile);
			FileWriter wp = new FileWriter(pctFile);
			String header = "Sample\tMapped\tUnmapped\n";
			w.write(header);
			wp.write(header);
			for(String sample : sampleNames) {
				int mapped = countAlignments(samtools, tophatOutputBamPerSample.get(sample), tophatDirsPerSample.get(sample), false);
				int unmapped = 0;
				if(version2) unmapped = countAlignments(samtools, tophatDirsPerSample.get(sample) + "/unmapped.bam", tophatDirsPerSample.get(sample), true, false);
				else {
					// Tophat version 1 writes separate unmapped fastq files
					// Unzip files
					String cmmd1 = "gunzip " + tophatDirsPerSample.get(sample) + "/unmapped_left.fq.z";
					Process p1 = Runtime.getRuntime().exec(cmmd1);
					p1.waitFor();
					String cmmd2 = "mv " + tophatDirsPerSample.get(sample) + "/unmapped_left.fq " + tophatDirsPerSample.get(sample) + "/unmapped_1.fq";
					Process p2 = Runtime.getRuntime().exec(cmmd2);
					p2.waitFor();
	
					String cmmd3 = "gunzip " + tophatDirsPerSample.get(sample) + "/unmapped_right.fq.z";
					Process p3 = Runtime.getRuntime().exec(cmmd3);
					p3.waitFor();
					String cmmd4 = "mv " + tophatDirsPerSample.get(sample) + "/unmapped_right.fq " + tophatDirsPerSample.get(sample) + "/unmapped_2.fq";
					Process p4 = Runtime.getRuntime().exec(cmmd4);
					p4.waitFor();
	
					// Count unmapped reads
					int totalLines = 0;
					FileReader r = new FileReader(tophatDirsPerSample.get(sample) + "/unmapped_1.fq");
					BufferedReader b = new BufferedReader(r);
					while(b.ready()) {
						@SuppressWarnings("unused")
						String line = b.readLine();
						totalLines++;
					}
					b.close();
					if(rightFastqs != null) {
						if(rightFastqs.containsKey(sample)) {
							FileReader r2 = new FileReader(tophatDirsPerSample.get(sample) + "/unmapped_2.fq");
							BufferedReader b2 = new BufferedReader(r2);
							while(b2.ready()) {
								@SuppressWarnings("unused")
								String line = b2.readLine();
								totalLines++;
							}
							b2.close();
						}
					}
					unmapped = totalLines / 4;
				
				}
				int total = mapped + unmapped;
				w.write(sample + "\t" + mapped + "\t" + unmapped + "\n");
				wp.write(sample + "\t" + (double)mapped/(double)total + "\t" + (double)unmapped/(double)total + "\n");
			}
			w.close();
			wp.close();
		}
		
		logger.info("Wrote table of counts to file " + countFile);
		logger.info("Wrote table of percentages to file " + pctFile);
		
		return rtrn;
		
	}

	
}
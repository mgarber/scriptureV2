package broad.core.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.readFilters.FirstOfPairFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.readFilters.SecondOfPairFilter;
import nextgen.core.writers.PairedEndWriter;
import nextgen.core.writers.WigWriter;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.math.EmpiricalDistribution;
import broad.core.parser.CommandLineParser;
import broad.core.parser.ConfigFileParser;
import broad.core.parser.StringParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.countreads.FastqLibraryStats;
import broad.pda.countreads.LibraryCompositionByRnaClass;

/**
 * This class will use an input list of files and a config file to run an automated pipeline to align and process sequencing data
 * @author skadri 
 * @author prussell
 *
 */
public class PipelineAutomator {

	static Logger logger = Logger.getLogger(PipelineAutomator.class.getName());

	/*
	 * For a file with paired data in the fq list file, there will be 4 columns
	 * Name	left.fq	right.fq	condition
	 * Name is unique
	 * Condition represents a string describing the sample/expt. Replicates have the same string in the condition
	 */
	private static int PAIRED_COLUMNS = 4;
	private static int UNPAIRED_COLUMNS = PAIRED_COLUMNS - 1;
	
	ConfigFileParser configP;
	
	/**
	 * The sample names
	 */
	TreeSet<String> sampleNames;
	
	/*
	 * Mapping of name to FastQ file paths
	 */
	Map<String,String> leftFqs;
	Map<String,String> rightFqs;
	/*
	 * Mapping of name to whether there is paired data
	 */
	Map<String,Boolean> pairedData;
	
	/*
	 * Mappings of names to CURRENT fqs being used. 
	 * This is because if the user filters out rRNA, the name of the fqs to align changes.
	 */
	Map<String,String> currentLeftFqs;
	Map<String,String> currentRightFqs;

	/**
	 * The paths of the current bam files being used
	 * This is because we might run extra alignments and merge bam files
	 */
	Map<String,String> currentBamFiles;
		
	/**
	 * The directory containing the current bam files
	 */
	String currentBamDir;
	
	Map<String,String> nameToCondition;

	String queueName;
	
	/**
	 * The delimiter that separates the pair number from the main part of the read ID in fastq file
	 */
	private String fastqReadIdPairNumberDelimiter;
	
	/*
	 * LIST OF COMMANDS (in order for ease of coding)
	 */
	static String SPLIT_TRIM_BARCODES = "SPLIT_TRIM_BARCODES";
	static String CLIP_ADAPTERS = "TRIM_ADAPTERS";
	static String LIBRARY_STATS = "LIBRARY_STATS";
	static String RNA_CLASSES = "RNA_CLASSES";
	static String FILTER_RRNA = "FILTER_RRNA";
	static String ALIGN = "ALIGN";
	static String ALIGN_TO_TRANSCRIPTS = "ALIGN_TO_TRANSCRIPTS";
	static String RUN_DGE = "RUN_DGE";
	
	/**
	 * Output directories
	 */
	static String TOPHAT_DIRECTORY= "tophat_to_genome";
	static String NOVOALIGN_DIRECTORY = "novoalign_unmapped_to_genome";
	static String MERGED_TOPHAT_NOVOALIGN_DIRECTORY = "merged_alignments_tophat_novoalign";
	static String LIBRARY_STATS_DIRECTORY = "library_stats";
	static String FILTER_RRNA_DIRECTORY = "filter_rRNA";
	static String ALIGN_TO_TRANSCRIPTS_DIRECTORY = "bowtie_to_transcripts";


	
	/**
	 * @param inputListFile The input list of fastq files
	 * @param configFileName The config file
	 * @throws IOException
	 * @throws ParseException
	 * @throws MathException
	 * @throws InterruptedException
	 */
	public PipelineAutomator(String inputListFile,String configFileName) throws IOException, ParseException, InterruptedException {
		
		Globals.setHeadless(true);
		logger.info("Using Version R4.4");
		logger.debug("DEBUG ON");

		configP = new ConfigFileParser(configFileName);
		fastqReadIdPairNumberDelimiter = configP.basicOptions.getFastqReadIdPairDelimiter();
		
		//If this flag is false, then the expected input file formats are different
		boolean runBasic = false;
		/*
		 * Format of the inputListFile can depend on the first task in question
		 */
		// IF THE CONFIG FILE CONTAINS A BASIC COMMAND
		if(containsBasicCommand()){
			
			runBasic = true;
			
			//Initialize the FQ file lists
			sampleNames = new TreeSet<String>();
			leftFqs = new TreeMap<String,String>();
			rightFqs = new TreeMap<String,String>();
			currentLeftFqs = new TreeMap<String,String>();
			currentRightFqs = new TreeMap<String,String>();
			pairedData = new TreeMap<String,Boolean>();
			currentBamFiles = new TreeMap<String,String>();
			nameToCondition = new TreeMap<String,String>();
			queueName = configP.basicOptions.getQueueName();
			
			//Read the Fq list
			readFqList(inputListFile);

			// Split and trim barcodes
			if(configP.hasCommandFor(SPLIT_TRIM_BARCODES)){
				splitTrimBarcodes();
			}
			
			// Clip sequencing adapters
			if(configP.hasCommandFor(CLIP_ADAPTERS)) {
				String fastxDir = configP.basicOptions.getFastxDirectory();
				String adapter1 = configP.basicOptions.getSequencingAdapterRead1();
				String adapter2 = null;
				if(configP.basicOptions.hasCommandFor(ConfigFileParser.OPTION_SEQUENCING_ADAPTER_READ_2)) adapter2 = configP.basicOptions.getSequencingAdapterRead2();
				clipAdapters(fastxDir, adapter1, adapter2);
			}
			
			// Count reads
			// Quantify duplicates
			// Estimate library size
			if(configP.hasCommandFor(LIBRARY_STATS)){
				calculateLibraryStats();
			}
			
			// Characterize library composition by RNA class
			if(configP.hasCommandFor(RNA_CLASSES)){
				quantifyRNAClasses();
			}
			
			// Update current fastqs with reads not matching rRNA
			if(configP.hasCommandFor(FILTER_RRNA)){
				filterrRNA();
			}
			
			// Align to transcript sequences
			// Calculate median fragment size for each sequence
			if(configP.hasCommandFor(ALIGN_TO_TRANSCRIPTS)) {
				alignToTranscripts();
			}
			
			// Align to genome
			// Generate bam and tdf files
			// Generate fragment size distributions
			if(configP.hasCommandFor(ALIGN)){
				alignToGenome();
			}
			
		}
		
		if(configP.hasCommandFor(RUN_DGE)){
			
			String DGEInputFile= null;
			if(configP.hasCommandFor(ALIGN)){
				//Input file format
				//Name	Bam_file	condition
				//TODO: Not needed right now but provision for future flexibility
				DGEInputFile = inputListFile;
				//Make the DGE command using all the options
				
			}
			else{
				/*
				 * Create my own DGE input file using the output of above
				 */
				DGEInputFile = createDGEInputFile(inputListFile);
			}
			
			runDGE(runBasic,DGEInputFile);
			
		}
	}
	
	/**
	 * Returns true if the config file contains one or more basic commands
	 * @return
	 */
	private boolean containsBasicCommand(){
		return (configP.hasCommandFor(SPLIT_TRIM_BARCODES) || configP.hasCommandFor(LIBRARY_STATS) 
				|| configP.hasCommandFor(RNA_CLASSES) || configP.hasCommandFor(FILTER_RRNA) || configP.hasCommandFor(ALIGN) || configP.hasCommandFor(ALIGN_TO_TRANSCRIPTS));
	}
	
	/**
	 * 
	 * @param fileName
	 * @throws IOException 
	 */
	private void readFqList(String fileName) throws IOException{
		
		/*	Input file is a list
		*	Name		Left.fq		Right.fq	condition
		*	If unpaired, only Name	Left.fq	condition
		*	The NAME must be unique
		*/
		//TODO: If the name is not unique, assign a unique name
		logger.info("");
		logger.info("Reading list of fastq files from " + fileName + "...");
		
		FileReader r = new FileReader(fileName);
		BufferedReader b = new BufferedReader(r);
		if(!b.ready()) {
			throw new IllegalArgumentException("File " + fileName + " is empty.");
		}
		StringParser p = new StringParser();
		
		boolean first = true;
		boolean unpaired = true;
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			int cols = p.getFieldCount();
			if(cols == 0) continue;
			
			if(first) {
				if(cols < UNPAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}
				if(cols == PAIRED_COLUMNS) unpaired = false;
				first = false;
			}
			
			String sampleName = p.asString(0);
			sampleNames.add(sampleName);
			String leftReads = p.asString(1);
			leftFqs.put(sampleName, leftReads);
			
			if(unpaired) {
				if(cols != UNPAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}
				pairedData.put(sampleName, Boolean.valueOf(false));
				nameToCondition.put(sampleName, p.asString(2));
			} else {
				if(cols != PAIRED_COLUMNS) {
					throw new IllegalArgumentException("Illegal number of columns in " + fileName);
				}				
				rightFqs.put(sampleName, p.asString(2));
				pairedData.put(sampleName, Boolean.valueOf(true));
				nameToCondition.put(sampleName, p.asString(3));
			}
			
		}
		
		b.close();
		r.close();
		
		if(first) {
			throw new IllegalArgumentException("File " + fileName + " is invalid.");
		}
		
		currentLeftFqs = leftFqs;
		currentRightFqs = rightFqs;
		
		int pa = currentRightFqs.isEmpty() ? 0 : currentRightFqs.keySet().size();
		int u = currentLeftFqs.keySet().size() - pa;
		logger.info("There are " + pa + " sets of paired end reads and " + u + " sets of single end reads.");
		
	}
	
	/**
	 * TASK 1: SPLIT_TRIM_BARCODES
	 */
	private void splitTrimBarcodes(){
		// TODO: finish
	}
	
	/**
	 * Clip sequencing adapters
	 * @param fastxDir Directory containing fastx binaries
	 * @param adapter1 Sequencing adapter for read 1
	 * @param adapter2 Sequencing adapter for read 2
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private void clipAdapters(String fastxDir, String adapter1, String adapter2) throws IOException, InterruptedException {
		
		logger.info("");
		logger.info("Trimming sequencing adapters...");
		Map<String, String> outTmpFilesLeft = new TreeMap<String, String>();
		Map<String, String> outClippedFilesLeft = new TreeMap<String, String>();
		Map<String, String> outTmpFilesRight = new TreeMap<String, String>();
		Map<String, String> outClippedFilesRight = new TreeMap<String, String>();
		
		ArrayList<String> jobIDs = new ArrayList<String>();
		
		// Clip left fastq files
		for(String sampleName : sampleNames) {
			String inFile = currentLeftFqs.get(sampleName);
			String outTmpFile = inFile + ".clipped_tmp";
			String outClippedFile = inFile + ".clipped.fq";
			outTmpFilesLeft.put(sampleName, outTmpFile);
			outClippedFilesLeft.put(sampleName, outClippedFile);
			File finalClipped = new File(outClippedFile);
			if(finalClipped.exists()) {
				logger.warn("Clipped file " + outClippedFile + " already exists. Not rerunning fastx_clipper.");
				continue;
			}
			File tmpFile = new File(outTmpFile);
			if(!tmpFile.exists()) {
				// Use fastx program fastx_clipper
				String cmmd = fastxDir + "/fastx_clipper -a " + adapter1 + " -Q 33 -n -i " + inFile + " -o " ;
				if(pairedData.get(sampleName).booleanValue()) {
					cmmd += outTmpFile;
				} else {
					cmmd += outClippedFile;
				}
				logger.info("Running fastx command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				jobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, "fastx_clipper_" + jobID + ".bsub", "week", 4);
			} else {
				logger.warn("Temp clipped file " + outTmpFile + " already exists. Not rerunning fastx_clipper. Starting from temp file.");
			}
		}
		
		// Clip right fastq files
		for(String sampleName : sampleNames) {
			if(pairedData.get(sampleName).booleanValue()) {
				String inFile = currentRightFqs.get(sampleName);
				String outTmpFile = inFile + ".clipped_tmp";
				String outClippedFile = inFile + ".clipped.fq";
				outTmpFilesRight.put(sampleName, outTmpFile);
				outClippedFilesRight.put(sampleName, outClippedFile);
				File finalClipped = new File(outClippedFile);
				if(finalClipped.exists()) {
					logger.warn("Clipped file " + outClippedFile + " already exists. Not rerunning fastx_clipper.");
					continue;
				}
				File tmpFile = new File(outTmpFile);
				if(!tmpFile.exists()) {
					// Use fastx program fastx_clipper
					String cmmd = fastxDir + "/fastx_clipper -a " + adapter2 + " -Q 33 -n -i " + inFile + " -o " + outTmpFile;
					logger.info("Running fastx command: " + cmmd);
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					jobIDs.add(jobID);
					logger.info("LSF job ID is " + jobID + ".");
					PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, "fastx_clipper_" + jobID + ".bsub", "week", 4);
				} else {
					logger.warn("Temp clipped file " + outTmpFile + " already exists. Not rerunning fastx_clipper. Starting from temp file.");
				}
			}
		}
		
		logger.info("Waiting for fastx_clipper jobs to finish...");
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		
		logger.info("Done running fastx_clipper.");
		
		for(String sampleName : sampleNames) {
			if(!pairedData.get(sampleName).booleanValue()) {
				currentLeftFqs.put(sampleName, outClippedFilesLeft.get(sampleName));
				logger.info("Current left fq file for sample " + sampleName + " is " + currentLeftFqs.get(sampleName) + ".");
				continue;
			}
			String inFastq1 = outTmpFilesLeft.get(sampleName);
			String inFastq2 = outTmpFilesRight.get(sampleName);
			String outFastq1 = outClippedFilesLeft.get(sampleName);
			String outFastq2 = outClippedFilesRight.get(sampleName);
			File file1 = new File(outFastq1);
			File file2 = new File(outFastq1);
			if(file1.exists() && file2.exists()) {
				logger.warn("Using existing files " + outFastq1 + " and " + outFastq2 + ".");
				currentLeftFqs.put(sampleName, outFastq1);
				currentRightFqs.put(sampleName, outFastq2);
				logger.info("Current fastq files for sample " + sampleName + " are " + currentLeftFqs.get(sampleName) + " and " + currentRightFqs.get(sampleName) + ".");
				continue;
			}
			logger.info("Removing reads from files " + inFastq1 + " and " + inFastq2 + " that are not in both files and writing new files to " + outFastq1 + " and " + outFastq2 + "...");
			filterPairedFastqFilesMissingReads(inFastq1, inFastq2, outFastq1, outFastq2);
			logger.info("Done writing cleaned fastq files. To save storage, delete temporary files " + inFastq1 + " and " + inFastq2 + ".");
			currentLeftFqs.put(sampleName, outFastq1);
			currentRightFqs.put(sampleName, outFastq2);
			logger.info("Current fastq files for sample " + sampleName + " are " + currentLeftFqs.get(sampleName) + " and " + currentRightFqs.get(sampleName) + ".");
			
		}
		
		
		
		
	}
	
	/**
	 * Remove reads that are missing from file1 or file2 and write new fastq files
	 * @param inFastq1 Input fastq file read 1
	 * @param inFastq2 Input fastq file read 2
	 * @param outFastq1 Output fastq file read 1
	 * @param outFastq2 Output fastq file read 2
	 * @throws IOException 
	 */
	private void filterPairedFastqFilesMissingReads(String inFastq1, String inFastq2, String outFastq1, String outFastq2) throws IOException {
		
		StringParser sp = new StringParser();
		
		// Save all read names from infile 1 to a treeset
		TreeSet<String> inFastq1Ids = new TreeSet<String>();
		FileReader r1 = new FileReader(inFastq1);
		BufferedReader b1 = new BufferedReader(r1);
		int linesRead1 = 0;
		while(b1.ready()) {
			String line = b1.readLine();
			linesRead1++;
			if(linesRead1 % 4 == 1) {
				sp.parse(line, fastqReadIdPairNumberDelimiter);
				String id = sp.asString(0);
				// Avoid rare duplicated read IDs
				if(inFastq1Ids.contains(id)) { 
					inFastq1Ids.remove(id);
					logger.warn("Skipping read " + id + " because appears in file " + inFastq1 + " twice.");
					continue;
				}
				inFastq1Ids.add(id);
				id = null;
			}
			line = null;
		}
		r1.close();
		b1.close();
		
		// Save read names that appear in both files to a treeset
		// Write those reads from input file 2 to new file
		FileWriter w2 = new FileWriter(outFastq2);
		TreeSet<String> bothFastqIds = new TreeSet<String>();
		TreeSet<String> inFastq2Ids = new TreeSet<String>();
		FileReader r2 = new FileReader(inFastq2);
		BufferedReader b2 = new BufferedReader(r2);
		int linesRead2 = 0;
		while(b2.ready()) {
			String line = b2.readLine();
			linesRead2++;
			if(linesRead2 % 4 == 1) {
				sp.parse(line, fastqReadIdPairNumberDelimiter);
				String id = sp.asString(0);
				// Avoid rare duplicated read IDs
				if(inFastq2Ids.contains(id)) {
					logger.warn("Skipping read " + id + " because appears in file " + inFastq2 + " twice.");
					inFastq2Ids.remove(id);
					continue;
				}
				inFastq2Ids.add(id);
				if(inFastq1Ids.contains(id)) {
					bothFastqIds.add(id);
					String line2 = b2.readLine();
					String line3 = b2.readLine();
					String line4 = b2.readLine();
					linesRead2 += 3;
					w2.write(line + "\n" + line2 + "\n" + line3 + "\n" + line4 + "\n");			
					line2 = null;
					line3 = null;
					line4 = null;
				}
				id = null;
			}
			line = null;
		}
		r2.close();
		b2.close();
		w2.close();
		
		// Write reads from infile 1 that appear in both files to new file
		FileWriter w1 = new FileWriter(outFastq1);
		FileReader r3 = new FileReader(inFastq1);
		BufferedReader b3 = new BufferedReader(r3);
		int linesRead3 = 0;
		while(b3.ready()) {
			String line = b3.readLine();
			linesRead3++;
			if(linesRead3 % 4 == 1) {
				sp.parse(line, fastqReadIdPairNumberDelimiter);
				if(bothFastqIds.contains(sp.asString(0))) {
					String line2 = b3.readLine();
					String line3 = b3.readLine();
					String line4 = b3.readLine();
					linesRead3 += 3;
					w1.write(line + "\n" + line2 + "\n" + line3 + "\n" + line4 + "\n");			
					line2 = null;
					line3 = null;
					line4 = null;
				}
			}
			line = null;
		}
		r3.close();
		b3.close();
		w1.close();
		
		
		
		
	}
	
	/**
	 * TASK 2: CALCULATE LIBRARY STATS
	 * For each library, count total reads, unique reads, percent duplicates, and estimated library size
	 * @author prussell
	 * @throws IOException 
	 */
	private void calculateLibraryStats() throws IOException {
		
		logger.info("");
		logger.info("Calculating library stats...");
		
		// Make directory for library stats
		File dir = new File(LIBRARY_STATS_DIRECTORY);
		boolean madeDir = dir.mkdir();
		if(!dir.exists()) {
			throw new IOException("Could not create directory " + LIBRARY_STATS_DIRECTORY);
		}
		String output = LIBRARY_STATS_DIRECTORY + "/" + "library_stats.out";
		File outputFile = new File(output);
		if(outputFile.exists()) {
			logger.warn("Library stats file " + output + " exists. Not calculating new library stats.");
			return;
		}
		
		logger.info("Writing library stats to file " + output + ".");
		FileWriter w = new FileWriter(output);
		
		// Get stats for each sample
		boolean first = true;
		for(String sample : leftFqs.keySet()) {
			logger.info("Calculating library stats for sample " + sample + "...");
			if(first && pairedData.get(sample).booleanValue()) {
				w.write("Sample\tTotal_read_pairs\tUnique_read_pairs\tPercent_duplicated\tEst_library_size\n");
			}
			if(first && !pairedData.get(sample).booleanValue()) {
				w.write("Sample\tTotal_reads\tUnique_reads\tPercent_duplicated\tEst_library_size\n");
			}
			first = false;
			FastqLibraryStats d = new FastqLibraryStats(leftFqs.get(sample), pairedData.get(sample).booleanValue() ? rightFqs.get(sample) : null);
			w.write(sample + "\t" + d.getTotalReads() + "\t" + d.getNumUniqueReads() + "\t" + d.getPercentDuplicated() + "\t" + d.getEstimatedLibrarySize() + "\n");
		}
		w.close();
		logger.info("");
		logger.info("Done calculating library stats.");
	}
	
	/**
	 * TASK 3: CALCULATE RNA CLASSES
	 * Quantify percentage of reads originating from each RNA class
	 * @author prussell
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	private void quantifyRNAClasses() throws IOException, InterruptedException{
		
		logger.info("");
		logger.info("Quantifying RNA classes...");
		
		// Make directory for RNA classes
		String rnaClassDir = "rna_class_counts";
		File dir = new File(rnaClassDir);
		boolean madeDir = dir.mkdir();
		if(!dir.exists()) {
			throw new IOException("Could not create directory " + rnaClassDir);
		}
		
		// Get options from config file
		Map<String, String> classFiles = configP.basicOptions.getRNAClassFileNames();
		String genomeBowtieIndex = configP.basicOptions.getGenomeBowtieIndex();
		String bowtie2Executable = configP.basicOptions.getBowtie2ExecutablePath();
		String bowtie2BuildExecutable = configP.basicOptions.getBowtie2BuildExecutablePath();
		String samtoolsExecutable = configP.basicOptions.getSamtoolsPath();
		
		// Files to write tables to
		String countFileName = rnaClassDir + "/counts_by_class.out";
		String pctFileName = rnaClassDir + "/percentages_by_class.out";
		File countFile = new File(countFileName);
		File pctFile = new File(pctFileName);
		if(countFile.exists() && pctFile.exists()) {
			logger.warn("RNA class files " + countFileName + " and " + pctFileName + " already exist. Not rerunning RNA class counts.");
			logger.info("");
			logger.info("Done quantifying RNA classes.");
			return;
		}
		
		// Align and count reads for each library and RNA class
		LibraryCompositionByRnaClass lcrc = new LibraryCompositionByRnaClass(genomeBowtieIndex, classFiles, leftFqs, rightFqs, logger);
		Map<String, String> bowtie2options = configP.basicOptions.getBowtie2Options();
		Map<String, Integer> totalReadCounts = lcrc.getTotalReadCounts();
		Map<String, Map<String, Integer>> classCounts = lcrc.alignAndGetCounts(samtoolsExecutable, bowtie2Executable, bowtie2options, bowtie2BuildExecutable, rnaClassDir);
		
		logger.info("Writing table of counts to file " + countFileName);
		logger.info("Writing table of percentages to file " + pctFileName);
		FileWriter countWriter = new FileWriter(countFileName);
		FileWriter pctWriter = new FileWriter(pctFileName);
		
		// Get set of class names
		Iterator<String> iter = classCounts.keySet().iterator();
		Set<String> classNames = classCounts.get(iter.next()).keySet();
		String header = "Sample\t";
		for(String className : classNames) {
			header += className;
			header += "\t";
		}
		header += "\n";
		countWriter.write(header);
		pctWriter.write(header);
		
		// Get counts and make percentages
		for(String sample : classCounts.keySet()) {
			String countLine = sample + "\t";
			String pctLine = sample + "\t";
			for(String className : classNames) {
				countLine += classCounts.get(sample).get(className).toString();
				String pct = Double.valueOf((double)classCounts.get(sample).get(className).intValue() / (double)totalReadCounts.get(sample).intValue()).toString();
				pctLine += pct;
				countLine += "\t";
				pctLine += "\t";
			}
			countLine += "\n";
			pctLine += "\n";
			countWriter.write(countLine);
			pctWriter.write(pctLine);
		}

		countWriter.close();
		pctWriter.close();
		
		logger.info("");
		logger.info("Done quantifying RNA classes.");
		logger.info("Delete sam and fastq files to save disk space. If rerunning pipeline remove " + RNA_CLASSES + " option from config file.");

	}
	
	/**
	 * TASK 4: FILTER OUT RRNA
	 * Align reads to rRNA, retain reads that do not align, and replace current fastq files with files of non aligning reads
	 * @author prussell
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private void filterrRNA() throws IOException, InterruptedException{
		
		// Establish paths and software locations
		File outDirFile = new File(FILTER_RRNA_DIRECTORY);
		boolean madeDir = outDirFile.mkdir();
		String outIndex = FILTER_RRNA_DIRECTORY + "/rRNA";
		String rRnaFasta = configP.basicOptions.getRrnaSeqsForDepletion();
		String bowtieBuild = configP.basicOptions.getBowtie2BuildExecutablePath();
		String bowtie = configP.basicOptions.getBowtie2ExecutablePath();
		
		logger.info("");
		logger.info("Filtering ribosomal RNA by removing reads that map to sequences in " + rRnaFasta);
		
		// Make bowtie2 index for ribosomal RNA
		AlignmentUtils.makeBowtie2Index(rRnaFasta, outIndex, bowtieBuild, FILTER_RRNA_DIRECTORY);
		
		// Establish output file names
		Map<String, String> outRibosomalMap = new TreeMap<String, String>();
		Map<String, String> outFilteredUnpairedMap = new TreeMap<String, String>();
		Map<String, String> outFilteredPairedArgMap = new TreeMap<String, String>();
		Map<String, String> outFilteredPaired1Map = new TreeMap<String, String>();
		Map<String, String> outFilteredPaired2Map = new TreeMap<String, String>();
		ArrayList<String> jobIDs = new ArrayList<String>();
		
		// Align each sample to rRNA
		for(String sample : sampleNames) {

			String outRibosomal = FILTER_RRNA_DIRECTORY + "/" + sample + "_ribosomal_mappings.sam";
			outRibosomalMap.put(sample, outRibosomal);
			String outFilteredUnpaired = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA.fq";
			outFilteredUnpairedMap.put(sample, outFilteredUnpaired);
			String outFilteredPairedArg = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA_%.fq";
			outFilteredPairedArgMap.put(sample, outFilteredPairedArg);
			String outFilteredPaired1 = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA_1.fq";
			outFilteredPaired1Map.put(sample,outFilteredPaired1);
			String outFilteredPaired2 = FILTER_RRNA_DIRECTORY + "/" + sample + "_filtered_rRNA_2.fq";
			outFilteredPaired2Map.put(sample,outFilteredPaired2);
			
			boolean paired = pairedData.get(sample).booleanValue();
			
			// check if filtered fastq files exist
			if(!paired) {
				File filteredFile = new File(outFilteredUnpaired);
				if(filteredFile.exists()) {
					logger.warn("Filtered file for sample " + sample + " already exists: " + outFilteredUnpaired + ". Not creating new file.");
					continue;
				}
			} else {
				logger.info("Filtering sample " + sample + "...");
				File filteredFile1 = new File(outFilteredPaired1);
				File filteredFile2 = new File(outFilteredPaired2);
				if(filteredFile1.exists() && filteredFile2.exists()) {
					logger.warn("Filtered files for sample " + sample + " already exist: " + outFilteredPaired1 + " and " + outFilteredPaired2 + ". Not creating new files.");
					continue;
				}
				
			}
			
			// Align to rRNA and keep unmapped reads
			Map<String, String> bowtie2options = configP.basicOptions.getBowtie2Options();
			String jobID = AlignmentUtils.runBowtie2(outIndex, bowtie2options, currentLeftFqs.get(sample), paired ? currentRightFqs.get(sample) : null, outRibosomal, paired ? outFilteredPairedArg : outFilteredUnpaired, bowtie, FILTER_RRNA_DIRECTORY, paired);
			jobIDs.add(jobID);
			
		}
		
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		
		logger.info("");
		logger.info("Done aligning to rRNAs. Updating current fastq files to filtered files.");
		
		// Update current fastq files to the ribosome filtered files
		for(String sample : leftFqs.keySet()) {
			if(pairedData.get(sample).booleanValue()) {
				currentLeftFqs.put(sample, outFilteredPaired1Map.get(sample));
				currentRightFqs.put(sample, outFilteredPaired2Map.get(sample));
				logger.info("Current fastq files for sample " + sample + " are " + currentLeftFqs.get(sample) + " and " + currentRightFqs.get(sample) + ".");
			} else {
				currentLeftFqs.put(sample, outFilteredUnpairedMap.get(sample));
				logger.info("Current fastq file for sample " + sample + " is " + currentLeftFqs.get(sample) + ".");
			}
		}
		
		logger.info("");
		logger.info("Done filtering rRNA and updating current fastq files.");
		logger.info("Delete sam files to save disk space. Keep fastq files for future pipeline runs.");
	}
	
	private void alignToTranscripts() throws IOException, InterruptedException{
		
		// Establish paths and software locations
		File outDirFile = new File(ALIGN_TO_TRANSCRIPTS_DIRECTORY);
		boolean madeDir = outDirFile.mkdir();
		String outIndex = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/transcripts";
		String fasta = configP.basicOptions.getTranscriptSeqsForAlignment();
		String bowtieBuild = configP.basicOptions.getBowtie2BuildExecutablePath();
		String bowtie = configP.basicOptions.getBowtie2ExecutablePath();
		Map<String, String> bowtie2options = configP.basicOptions.getBowtie2Options();
		String samtools = configP.basicOptions.getSamtoolsPath();
		String picardJarDir = configP.basicOptions.getPicardDirectory();
		
		logger.info("");
		logger.info("Aligning to transcript sequences in " + fasta);
		
		// Get transcript names and sizes
		Map<String, Integer> sequenceSizes = new TreeMap<String, Integer>();
		FastaSequenceIO fsio = new FastaSequenceIO(fasta);
		List<Sequence> seqs = fsio.loadAll();
		for(Sequence seq : seqs) {
			String name = seq.getId();
			int len = seq.getLength();
			sequenceSizes.put(name, Integer.valueOf(len));
			logger.info("Got sequence " + name + "\tlength=" + len);
		}
		
		// Make bowtie2 index for transcripts
		AlignmentUtils.makeBowtie2Index(fasta, outIndex, bowtieBuild, ALIGN_TO_TRANSCRIPTS_DIRECTORY);
		
		// Establish output file names
		ArrayList<String> jobIDs = new ArrayList<String>();
		Map<String, String> samOutput = new TreeMap<String, String>();
		Map<String, String> unsortedBamOutput = new TreeMap<String, String>();
		Map<String, String> sortedBamOutput = new TreeMap<String, String>();
		Map<String, String> peBamOutput = new TreeMap<String, String>();
		
		// Align each sample
		for(String sample : sampleNames) {

			String sam = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/" + sample + "_transcript_mappings.sam";
			String unsortedBam = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/" + sample + "_transcript_mappings_unsorted.bam";
			String sortedBam = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/" + sample + "_transcript_mappings.bam";
			String peBam = sortedBam + PairedEndWriter.PAIRED_END_EXTENSION;
			samOutput.put(sample, sam);
			unsortedBamOutput.put(sample, unsortedBam);
			sortedBamOutput.put(sample, sortedBam);
			peBamOutput.put(sample, peBam);
			
			boolean paired = pairedData.get(sample).booleanValue();
			
			// Check if bam files exist
			File bamFile = new File(sortedBamOutput.get(sample));
			File pebamFile = new File(peBamOutput.get(sample));
			if(bamFile.exists()) {
				if(pairedData.get(sample).booleanValue() && pebamFile.exists()) {
					logger.warn("Bam file and paired end bam file for sample " + sample + " already exist. Not rerunning alignment.");
					continue;
				} 
				if(!pairedData.get(sample).booleanValue()) {
					logger.warn("Bam file for sample " + sample + " already exists. Not rerunning alignment.");
					continue;	
				}
			}
					
			// Delete old paired end bam file if one exists
			if(pebamFile.exists()) {
				boolean deleted = pebamFile.delete();
				if(!deleted) {
					logger.warn("Could not delete existing paired end bam file " + peBamOutput.get(sample) + ".");
				}
			}
			
			// Align to transcripts
			String jobID = AlignmentUtils.runBowtie2(outIndex, bowtie2options, currentLeftFqs.get(sample), paired ? currentRightFqs.get(sample) : null, sam, null, bowtie, ALIGN_TO_TRANSCRIPTS_DIRECTORY, paired);
			jobIDs.add(jobID);
			
		}
		
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		
		logger.info("");
		logger.info("Done aligning to transcripts. Converting sam to bam files.");
		Map<String, String> bsubDir = new TreeMap<String, String>();
		for(String sample : sampleNames) bsubDir.put(sample, ALIGN_TO_TRANSCRIPTS_DIRECTORY);
		samToBam(samOutput, unsortedBamOutput, sortedBamOutput, bsubDir, samtools);
		
		logger.info("");
		logger.info("Done converting sam to bam files. Delete sam files to save storage.");
		logger.info("Sorting bam files.");
		sortBamFiles(unsortedBamOutput, sortedBamOutput, bsubDir, sortedBamOutput, picardJarDir);
		// Delete unsorted bam files
		for(String sample : sampleNames) {
			File unsorted = new File(unsortedBamOutput.get(sample));
			boolean deleted = unsorted.delete();
			if(!deleted) logger.warn("Could not delete unsorted bam file " + unsortedBamOutput.get(sample) + ". Delete manually.");
		}
		
		logger.info("");
		logger.info("Done sorting bam files. Indexing sorted bam files.");
		indexBamFiles(sortedBamOutput, samtools);
				
		logger.info("");
		logger.info("Getting median fragment sizes per transcript.");
		String outputMedians = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/fragment_size_median";
		File outputMedianFile = new File(outputMedians);
		if(outputMedianFile.exists()) {
			logger.warn("File " + outputMedians + " already exists. Not recalculating median fragment sizes.");
		} else {
			// Calculate median of each sequence for each sample
			Map< String, Map<String, Double> > mediansBySample = new TreeMap<String, Map< String, Double>>();
			GenomicSpace gs = new GenomicSpace(sequenceSizes);
			for(String sample : sampleNames) {
				Map<String, Double> medianBySequence = new TreeMap<String, Double>();
				ScanStatisticDataAlignmentModel data = new ScanStatisticDataAlignmentModel(sortedBamOutput.get(sample), gs);
				data.addFilter(new ProperPairFilter());
				data.addFilter(new GenomicSpanFilter(configP.fragmentSizeOptions.getMaxGenomicSpan()));
				for(String seqName : sequenceSizes.keySet()) {
					Annotation seq = gs.getReferenceAnnotation(seqName);
					EmpiricalDistribution seqDistrib = data.getReadSizeDistribution(seq, gs, configP.fragmentSizeOptions.getMaxFragmentSize(), configP.fragmentSizeOptions.getNumBins());
					if(seqDistrib.getAllDataValues().isEmpty()) continue;
					double median = seqDistrib.getMedianOfAllDataValues();
					medianBySequence.put(seqName, Double.valueOf(median));
				}
				mediansBySample.put(sample, medianBySequence);
			}
			// Write file
			FileWriter w = new FileWriter(outputMedians);
			String header = "Sample\t";
			for(String seqName : sequenceSizes.keySet()) {
				header += seqName + "\t";
			}
			w.write(header + "\n");
			for(String sample : sampleNames) {
				String line = sample + "\t";
				for(String seqName : sequenceSizes.keySet()) {
					line += mediansBySample.get(sample).get(seqName) + "\t";
				}
				w.write(line + "\n");
			}
			w.close();
		}
		
		// Make paired end bam files
		logger.info("");
		logger.info("Making paired end bam files.");
		for(String sample : sampleNames) {
			if(!pairedData.get(sample).booleanValue()) {
				continue;
			}
			File pebamFile = new File(peBamOutput.get(sample));
			if(pebamFile.exists()) {
				logger.warn("Paired end bam file for sample " + sample + " already exists. Not regenerating file.");
				continue;
			}
			ScanStatisticDataAlignmentModel dam = new ScanStatisticDataAlignmentModel(sortedBamOutput.get(sample), new GenomicSpace(sequenceSizes));
		}
		
		// Make tdf of paired end bam files
		logger.info("");
		logger.info("Making tdf coverage files of paired end bam files.");
		logger.info("Indexing fasta file.");
		String indexedFasta = fasta + ".fai";
		File indexedFastaFile = new File(indexedFasta);
		if(indexedFastaFile.exists()) {
			logger.warn("Fasta index " + indexedFasta + " already exists. Not remaking fasta index.");
		} else {
			indexFastaFile(fasta, samtools, ALIGN_TO_TRANSCRIPTS_DIRECTORY);
			logger.info("Done indexing fasta file.");
		}
		logger.info("Making tdf files.");
		makeTdfs(peBamOutput, ALIGN_TO_TRANSCRIPTS_DIRECTORY, fasta, configP.basicOptions.getIgvtoolsPath());
		
		// Make wig and bigwig files of fragment ends
		logger.info("");
		logger.info("Making wig and bigwig files of fragment end points.");
		writeWigFragmentEnds(sortedBamOutput, ALIGN_TO_TRANSCRIPTS_DIRECTORY, fasta, null);
		logger.info("");
		logger.info("Done aligning to transcripts.");
		

	}

	
	/**
	 * TASK 5: ALIGN TO GENOME
	 * 1. Align to genome with tophat
	 * 2. If specified, use Novoalign to align unmapped reads and merge Tophat+Novoalign bam files
	 * 3. Report alignment statistics
	 * 4. Index bam files
	 * 5. Make tdf and cda files
	 * 6. Run Picard metrics
	 * @author prussell
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private void alignToGenome() throws IOException, InterruptedException{
		
		logger.info("");
		logger.info("Aligning to genome...");

		logger.info("");
		logger.info("Aligning current fastq files to genome using Tophat...");

		// Establish index files and executables
		String genomeIndex = configP.basicOptions.getGenomeBowtieIndex();
		String samtools = configP.basicOptions.getSamtoolsPath();
		String tophat = configP.basicOptions.getTophatPath();
		String picardJarDir = configP.basicOptions.getPicardDirectory();
		File tophatDirFile = new File(TOPHAT_DIRECTORY);
		tophatDirFile.mkdir();
		
		// Get tophat options
		Map<String, String> tophatOptions = configP.basicOptions.getTophatOptions();
		if(tophatOptions.containsKey("-o") || tophatOptions.containsKey("--output-dir")) {
			logger.warn("Overriding tophat output directory provided in config file. Creating directories for each sample.");
			tophatOptions.remove("-o");
			tophatOptions.remove("--output-dir");
		}
		
		// Establish output directories
		Map<String, String> tophatDirsPerSample = new TreeMap<String, String>();
		Map<String, String> tophatSubdirsPerSample = new TreeMap<String, String>();
		Map<String, String> tophatBamFinalPath = new TreeMap<String, String>();  // where the bam file will be moved to
		for(String sample : sampleNames) {
			tophatDirsPerSample.put(sample, TOPHAT_DIRECTORY + "/" + TOPHAT_DIRECTORY + "_" + sample);
			tophatSubdirsPerSample.put(sample, TOPHAT_DIRECTORY + "_" + sample);
			tophatBamFinalPath.put(sample, TOPHAT_DIRECTORY + "/" + sample + ".bam");
		}
		
		// Run tophat
		runTophat(tophat, samtools, tophatOptions, tophatDirsPerSample, tophatBamFinalPath, genomeIndex);
		logger.info("Done running tophat.");
		
		// *** Novoalign steps ***
		// Only run if novoalign path was provided in config file
		if(configP.basicOptions.hasNovoalignPath()) {
		
			logger.info("");
			logger.info("Entering steps to align unmapped reads with Novoalign...");
			
			// Make sure genome novoindex was provided in config file
			if(!configP.basicOptions.getAllOptionMaps().containsKey(ConfigFileParser.OPTION_GENOME_NOVOINDEX)) {
				throw new IllegalArgumentException("Novoalign index for genome is required. Specify in config file with option genome_novoindex.");
			}
			
			String novoalign = configP.basicOptions.getNovoalignPath();
			String novoIndex = configP.basicOptions.getGenomeNovoindex();
			
			// Convert tophat unmapped.bam files to fastq files
			logger.info("");
			logger.info("Getting unaligned reads in fastq format...");
			boolean version2 = tophat.substring(tophat.length()-1).equals("2");
			Map<String, String[]> unmappedFastq = unmappedToFastq(tophatDirsPerSample, picardJarDir, version2);
			logger.info("Got unaligned reads in fastq format.");
		
			// Run novoalign with unmapped reads
			logger.info("");
			logger.info("Aligning unmapped reads to genome with Novoalign...");
			
			// Establish file locations
			File novoDirFile = new File(NOVOALIGN_DIRECTORY);
			boolean madeNovDir = novoDirFile.mkdir();
			Map<String, String> novoDirsPerSample = new TreeMap<String, String>();
			Map<String, String> novoSubdirsPerSample = new TreeMap<String, String>();
			Map<String, String> novoSamOutput = new TreeMap<String, String>();
			Map<String, String> novoSamOutputNoHeader = new TreeMap<String, String>();
			Map<String, String> novoBamOutput = new TreeMap<String, String>();
			Map<String, String> novoSamOutputReheadered = new TreeMap<String, String>();
			Map<String, String> novoSortedBam = new TreeMap<String, String>();
			Map<String, String> novoBamFinalPath = new TreeMap<String, String>(); // the final location for bam file
			
			// Make directories for each sample
			for(String sample : sampleNames) {
				String dir = NOVOALIGN_DIRECTORY + "/" + NOVOALIGN_DIRECTORY + "_" + sample;
				File dirFile = new File(dir);
				boolean madeSampleDir = dirFile.mkdir();
				novoDirsPerSample.put(sample,dir);
				novoSubdirsPerSample.put(sample, NOVOALIGN_DIRECTORY + "_" + sample);
				novoSamOutput.put(sample, dir + "/" + novoalign + "_" + sample + ".sam");
				novoSamOutputReheadered.put(sample, dir + "/" + novoalign + "_" + sample + "_reheadered.sam");
				novoBamOutput.put(sample, dir + "/" + novoalign + "_" + sample + ".bam");
				novoSortedBam.put(sample, dir + "/" + novoalign + "_" + sample + ".sorted.bam");
				novoBamFinalPath.put(sample, NOVOALIGN_DIRECTORY + "/" + sample + ".bam");
				novoSamOutputNoHeader.put(sample, dir + "/" + novoalign + "_" + sample + "_noheader.sam");
			}
			
			// Run Novoalign
			runNovoalignOnUnmappedReads(tophatOptions, novoIndex, novoalign, novoDirsPerSample, novoSamOutput, novoBamFinalPath, unmappedFastq, version2);
			logger.info("Done running novoalign.");
			
			// Reheader novoalign sam files to match tophat sam header
			logger.info("");
			logger.info("Replacing headers in novoalign sam files with header from tophat alignments...");
			replaceNovoalignSamHeaders(tophatBamFinalPath, novoSamOutput, novoSamOutputNoHeader, novoSamOutputReheadered, novoBamFinalPath, samtools);
			logger.info("Done replacing sam headers.");
			
			// Convert novoalign sam files to bam format
			logger.info("");
			logger.info("Converting novoalign sam files to bam format...");
			samToBam(novoSamOutput, novoBamOutput, novoBamFinalPath, novoDirsPerSample, samtools);
			logger.info("All samples done converting to bam format.");
			
			// Sort the bam files
			logger.info("");
			logger.info("Sorting novoalign bam files...");
			sortBamFiles(novoBamOutput, novoSortedBam, novoDirsPerSample, novoBamFinalPath, picardJarDir);
			logger.info("All bam files sorted.");
			
			// Move all novoalign bam files to one directory
			logger.info("");
			logger.info("Moving all sorted novoalign bam files to directory " + NOVOALIGN_DIRECTORY + "...");
			for(String sample : sampleNames) {
				File finalBam = new File(novoBamFinalPath.get(sample));
				if(finalBam.exists()) {
					logger.warn("Alignment file " + finalBam + " already exists. Not replacing file.");
					continue;					
				}
				String cmmd = "mv " + novoSortedBam.get(sample) + " " + novoBamFinalPath.get(sample);
				Process p = Runtime.getRuntime().exec(cmmd, null);
				p.waitFor();
			}
	
			// Merge tophat and novoalign bam files
			logger.info("");
			logger.info("Merging tophat and novoalign bam files in directory " + MERGED_TOPHAT_NOVOALIGN_DIRECTORY	+ "...");
			mergeTophatNovoalign(tophatBamFinalPath, novoBamFinalPath, picardJarDir);
			logger.info("All bam files merged.");
			
			// Update current bam files to merged bams
			logger.info("");
			logger.info("Updating current bam files...");
			for(String sample : sampleNames) {
				currentBamFiles.put(sample, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam");
				logger.info("Current bam file for sample " + sample + " is " + currentBamFiles.get(sample));
			}
			
			// Update current bam directory
			currentBamDir = MERGED_TOPHAT_NOVOALIGN_DIRECTORY;
	
			// Count mapped and unmapped reads
			logger.info("");
			logger.info("Counting mapped and unmapped reads...");
			countMappingsMergedAlignments(samtools);
			logger.info("Done counting mappings.");
	
	
		} // *** Done with novoalign steps ***
		
		// Merge bam files
		logger.info("");
		logger.info("Merging bam files...");
		try {
			mergeBamFiles(picardJarDir);
		} catch (IllegalArgumentException e) {
			// Make sure current bam files are ordered the same way as reference genome
			logger.info("");
			logger.info("Reordering bam files in directory " + currentBamDir + "...");
			reorderCurrentBams(picardJarDir);
			logger.info("Done reordering bam files.");					
		}
		mergeBamFiles(picardJarDir);
		logger.info("Done merging bam files.");
		
		// Index current bam files
		logger.info("");
		logger.info("Indexing bam files...");
		indexCurrentBams(samtools);
		logger.info("All bam files indexed.");

		// Collect Picard metrics
		logger.info("");
		logger.info("Collecting Picard metrics for bam files in directory " + currentBamDir + "...");
		//collectPicardMetricsCurrentBams(picardJarDir);
		logger.info("All Picard metrics done.");
		
		// Make tdf files from current bam files
		logger.info("");
		logger.info("Making tdf files for bam files...");
		makeTdfs(currentBamFiles, currentBamDir, configP.basicOptions.getGenomeAssembly(), configP.basicOptions.getIgvtoolsPath());
		logger.info("All tdf files created.");
		
		// Make fragment size distributions
		if(configP.fragmentSizeOptions.hasOptionMakeFragmentSizeDistribution()) {
			logger.info("");
			logger.info("Making fragment size distributions for bam files...");
			makeFragmentSizeDistributionCurrentBams();
			logger.info("All fragment size distributions created.");
		}
		
		// Make wig and bigwig files of fragment ends
		if(configP.basicOptions.getWigToBigWigExecutable() != null && configP.basicOptions.getBedFileForFragmentEndWig() != null) {
			logger.info("");
			logger.info("Making wig and bigwig files of fragment end points.");
			writeWigFragmentEnds(currentBamFiles, currentBamDir, configP.basicOptions.getGenomeFasta(), configP.basicOptions.getBedFileForFragmentEndWig());
			logger.info("");
			logger.info("Done writing wig files.\n");
		}

		// Compute global transcriptome space stats
		if(configP.basicOptions.getAlignmentGlobalStatsJar() != null) {
			if(configP.basicOptions.getBedFileForTranscriptomeSpaceStats() != null || configP.basicOptions.getChrSizeFileForGenomicSpaceStats() != null) {
				writeAlignmentGlobalStats(configP.basicOptions.getAlignmentGlobalStatsJar(), configP.basicOptions.getBedFileForTranscriptomeSpaceStats(), configP.basicOptions.getChrSizeFileForGenomicSpaceStats());
			}
		}
	}
	

	/**
	 * Run tophat against genome
	 * @param tophat Tophat executable
	 * @param samtools Samtools executable
	 * @param tophatOptions Mapping of tophat flags to values
	 * @param tophatDirsPerSample Directory to write tophat output to, by sample name
	 * @param tophatBamFinalPath Location for final tophat output, by sample name
	 * @param genomeIndex Genome fasta index
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void runTophat(String tophat, String samtools, Map<String, String> tophatOptions, Map<String, String> tophatDirsPerSample, Map<String, String> tophatBamFinalPath, String genomeIndex) throws IOException, InterruptedException {
		
		// Identify tophat version
		boolean version2 = tophat.substring(tophat.length()-1).equals("2");
		
		ArrayList<String> tophatJobIDs = new ArrayList<String>();
		// Run tophat
		for(String sample : sampleNames) {
			File outdir = new File(tophatDirsPerSample.get(sample));
			outdir.mkdir();
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
			String tophatCmmd = tophat + " ";
			tophatCmmd += optionsString + " ";
			tophatCmmd += "--output-dir " + outdir + " "; 
			tophatCmmd += genomeIndex + " ";
			tophatCmmd += currentLeftFqs.get(sample) + " ";
			if(pairedData.get(sample).booleanValue()) tophatCmmd += currentRightFqs.get(sample);
			logger.info("Writing tophat output for sample " + sample + " to directory " + outdir + ".");
			logger.info("Running tophat command: " + tophatCmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			tophatJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, tophatCmmd, outdir + "/tophat_" + jobID + ".bsub", queueName, 16);
		}
		// Wait for tophat jobs to finish
		logger.info("Waiting for tophat jobs to finish...");
		PipelineUtils.waitForAllJobs(tophatJobIDs, Runtime.getRuntime());
		logger.info("All samples done aligning to genome.");
		
		// Move all tophat bam files to one directory
		logger.info("");
		logger.info("Moving all tophat bam files to directory " + TOPHAT_DIRECTORY + "...");
		for(String sample : sampleNames) {
			String cmmd = "mv " + tophatDirsPerSample.get(sample) + "/accepted_hits.bam " + tophatBamFinalPath.get(sample);
			Process p = Runtime.getRuntime().exec(cmmd, null);
			p.waitFor();
			// Update current bam files
			currentBamFiles.put(sample, tophatBamFinalPath.get(sample));
		}
		// Update current bam directory
		currentBamDir = TOPHAT_DIRECTORY;
		
		// Count aligned and unaligned reads
		// If tophat version 1, unzip unmapped fastq files and rename
		logger.info("");
		logger.info("Counting mapped and unmapped reads...");
		String countFile = TOPHAT_DIRECTORY + "/mapped_unmapped_count.out";
		String pctFile = TOPHAT_DIRECTORY + "/mapped_unmapped_percentage.out";
		
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
				int mapped = AlignmentUtils.countAlignments(samtools, TOPHAT_DIRECTORY + "/" + sample + ".bam", tophatDirsPerSample.get(sample), false);
				int unmapped = 0;
				if(version2) AlignmentUtils.countAlignments(samtools, tophatDirsPerSample.get(sample) + "/unmapped.bam", tophatDirsPerSample.get(sample), true, false);
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
						String line = b.readLine();
						totalLines++;
					}
					if(pairedData.get(sample)) {
						FileReader r2 = new FileReader(tophatDirsPerSample.get(sample) + "/unmapped_2.fq");
						BufferedReader b2 = new BufferedReader(r);
						while(b2.ready()) {
							String line = b2.readLine();
							totalLines++;
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

	}
	
	/**
	 * Convert unmapped reads to fastq format if necessary and get fastq file names
	 * @param tophatDirsPerSample Directories containing tophat output, by sample name
	 * @param picardJarDir Directory containing Picard jar files
	 * @param version2 Whether tophat2 was used
	 * @return Map associating sample name with unmapped fastq files for read1 and read2. Read2 file is null if tophat2 was used.
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private Map<String, String[]> unmappedToFastq(Map<String, String> tophatDirsPerSample, String picardJarDir, boolean version2) throws IOException, InterruptedException {
		Map<String, String[]> rtrn = new TreeMap<String, String[]>();
		if(version2) {
			// Convert the unmapped.bam files to fastq format and get file names
			Map<String, String> fastq1 = unmappedFastqTophat2(tophatDirsPerSample, picardJarDir);
			for(String sample : fastq1.keySet()) {
				String[] files = new String[2];
				files[0] = fastq1.get(sample);
				files[1] = null;
				rtrn.put(sample, files);
			}
			return rtrn;
		}
		for(String sample : sampleNames) {
			// Tophat version 1 writes unmapped files in fastq format
			// Just get the names
			String[] files = new String[2];
			files[0] = tophatDirsPerSample.get(sample) + "/unmapped_1.fq";
			files[1] = tophatDirsPerSample.get(sample) + "/unmapped_2.fq";
			rtrn.put(sample, files);
		}
		return rtrn;
	}
	
	
	/**
	 * Convert tophat unmapped.bam files to fastq format for further alignment
	 * @param tophatDirsPerSample Directories containing tophat2 output, by sample name
	 * @param picardJarDir Directory containing Picard jar files
	 * @throws IOException
	 * @throws InterruptedException
	 */	
	private Map<String, String> unmappedFastqTophat2(Map<String, String> tophatDirsPerSample, String picardJarDir) throws IOException, InterruptedException {
		ArrayList<String> convertJobIDs = new ArrayList<String>();
		// Store names of fastq files
		Map<String, String> unmappedFastq1 = new TreeMap<String,String>();
		for(String sample : sampleNames) {
			String dir = tophatDirsPerSample.get(sample);
			// Use Picard program SamToFastq
			String cmmd = "java -jar " + picardJarDir + "/SamToFastq.jar INPUT=" + dir + "/unmapped.bam VALIDATION_STRINGENCY=SILENT ";
			/* Tophat seems to discard pair information when writing unmapped reads to file unmapped.bam
			 * So proceed with unmapped reads as if they were single end
			 */
			String fastq1 = dir + "/unmapped.fq";
			unmappedFastq1.put(sample, fastq1);
			File fastq1file = new File(unmappedFastq1.get(sample));
			// Check if fastq file already exists
			if(fastq1file.exists()) {
				logger.warn("Fastq file " + unmappedFastq1.get(sample) + " already exists. Not rerunning format conversion.");
				continue;
			}
			// Complete command
			cmmd += " FASTQ=" + fastq1;
			logger.info("Running Picard command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			convertJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, dir + "/sam_to_fastq_" + jobID + ".bsub", "hour", 16);
		}
		logger.info("Waiting for SamToFastq jobs to finish...");
		PipelineUtils.waitForAllJobs(convertJobIDs, Runtime.getRuntime());
		return unmappedFastq1;

	}
	
	/**
	 * Use Novoalign to map reads that did not align with Tophat
	 * @param tophatOptions Map of tophat flag to option
	 * @param novoindex Executable to create novoindex
	 * @param novoalign Novoalign executable
	 * @param novoDirsPerSample Directory to write novoalign output to by sample name
	 * @param novoSamOutput Sam file to write novoalign output to by sample name
	 * @param novoBamFinalPath Final bam file for novoalign output, by sample name
	 * @param unmappedFastq Unmapped reads by sample name
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void runNovoalignOnUnmappedReads(Map<String, String> tophatOptions, String novoindex, String novoalign, Map<String, String> novoDirsPerSample, Map<String, String> novoSamOutput, Map<String, String> novoBamFinalPath, Map<String, String[]> unmappedFastq, boolean tophat2) throws IOException, InterruptedException {
		// Get novoalign options
		Map<String, String> novoalignOptions = configP.basicOptions.getNovoalignOptions();
		// Override certain flags if provided in config file
		if(novoalignOptions.containsKey("-F")) {
			logger.warn("Overriding novoalign option -F provided in config file. Using STDFQ.");
			tophatOptions.remove("-F");
		}
		if(novoalignOptions.containsKey("-f")) {
			logger.warn("Overriding novoalign option -f provided in config file. Using unmapped reads from Tophat.");
			tophatOptions.remove("-f");
		}
		if(novoalignOptions.containsKey("-o")) {
			logger.warn("Overriding novoalign option -o provided in config file. Using SAM.");
			tophatOptions.remove("-o");
		}
		if(novoalignOptions.containsKey("-d")) {
			logger.warn("Overriding novoalign option -d provided in config file. Using " + novoindex);
			tophatOptions.remove("-d");
		}
		// Make string of Novoalign options
		String novoOptionsString = " -d " + novoindex + " -F STDFQ -o SAM ";
		for(String flag : novoalignOptions.keySet()) {
			novoOptionsString += flag + " " + novoalignOptions.get(flag) + " ";
		}
		
		// Run novoalign
		String cmmdBase = novoalign + novoOptionsString;
		ArrayList<String> novoJobIDs = new ArrayList<String>();
		Map<String, String> novoBsubFiles = new TreeMap<String, String>();
		for(String sample : sampleNames) {
			File outdir = new File(novoDirsPerSample.get(sample));
			outdir.mkdir();
			File sam = new File(novoSamOutput.get(sample));
			File finalBam = new File(novoBamFinalPath.get(sample));
			// Check if Novoalign has already been run
			if(sam.exists()) {
				logger.warn("Alignment file " + sam + " already exists. Not rerunning alignment.");
				continue;
			}
			if(finalBam.exists()) {
				logger.warn("Alignment file " + finalBam + " already exists. Not rerunning alignment.");
				continue;					
			}
			String cmmd = cmmdBase + " -f " + unmappedFastq.get(sample)[0];
			if(!tophat2) cmmd += " " + unmappedFastq.get(sample)[1];
			logger.info("Writing novoalign output for sample " + sample + " to directory " + outdir + ".");
			logger.info("Running novoalign command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			novoJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			String bsubFile = outdir + "/novoalign_" + jobID + ".bsub";
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bsubFile, "week", 8);
			novoBsubFiles.put(sample,bsubFile);
		}
		// Wait for novoalign jobs to finish
		logger.info("Waiting for novoalign jobs to finish...");
		PipelineUtils.waitForAllJobs(novoJobIDs, Runtime.getRuntime());
		logger.info("All samples done aligning unmapped reads with novoalign.");
		// Parse bsub output to sam files
		logger.info("");
		logger.info("Parsing LSF output to sam files...");
		StringParser stringparse = new StringParser();
		for(String sample : sampleNames) {
			File sam = new File(novoSamOutput.get(sample));
			// Check if Novoalign has already been run
			if(sam.exists()) {
				logger.warn("Alignment file " + sam + " already exists. Not checking for or parsing bsub file.");
				continue;
			}
			File finalBam = new File(novoBamFinalPath.get(sample));
			if(finalBam.exists()) {
				logger.warn("Alignment file " + finalBam + " already exists. Not checking for or parsing bsub file.");
				continue;					
			}
			String bsub = novoBsubFiles.get(sample);
			String bsub_only = bsub + ".bsub_output_only";
			logger.info("Parsing sam lines from bsub file " + bsub + " to sam file " + sam + "...");
			logger.info("Saving bsub output (excluding sam lines) to file " + bsub_only + "...");
			FileReader fr = new FileReader(bsub);
			BufferedReader br = new BufferedReader(fr);
			FileWriter fw = new FileWriter(sam);
			FileWriter fwb = new FileWriter(bsub_only);
			while(br.ready()) {
				String line = br.readLine();
				if(AlignmentUtils.isSamLine(stringparse, line)) {
					fw.write(line + "\n");
				} else {
					fwb.write(line + "\n");
				}
			}
			fw.close();
			fwb.close();
			fr.close();
			logger.info("Deleting bsub file " + bsub);
			File bsubFile = new File(bsub);
			bsubFile.delete();
		}
		logger.info("Done creating sam files.");

	}
	
	/**
	 * Replace headers in novoalign sam files with header from tophat files
	 * @param tophatBamFinalPath Bam files produced by tophat, to get header
	 * @param novoSamOutput Sam files containing novoalign alignments
	 * @param novoSamOutputNoHeader Files to write novoalign alignments with no header
	 * @param novoSamOutputReheadered Files to write reheadered novoalign sam alignments
	 * @param novoBamFinalPath Final novoalign bam files; skip if already exists
	 * @param samtools Samtools executable
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void replaceNovoalignSamHeaders(Map<String, String> tophatBamFinalPath, Map<String, String> novoSamOutput, Map<String, String> novoSamOutputNoHeader, Map<String, String> novoSamOutputReheadered, Map<String, String> novoBamFinalPath, String samtools) throws IOException, InterruptedException {
		// First get tophat header and write to novoalign directory
		Iterator<String> tophatBamIter = tophatBamFinalPath.keySet().iterator();
		String firstTophatBam = tophatBamFinalPath.get(tophatBamIter.next());
		String getHeaderCmmd = samtools + " view -H -o ";
		String tmpHeader = NOVOALIGN_DIRECTORY + "/tophat_sam_header_sorted.sam";
		String tophatHeader = NOVOALIGN_DIRECTORY + "/tophat_sam_header.sam";
		getHeaderCmmd += tmpHeader + " ";
		getHeaderCmmd += firstTophatBam;
		String getHeaderJobID = Long.valueOf(System.currentTimeMillis()).toString();
		logger.info("");
		logger.info("Getting sam header from tophat alignments to reheader novoalign sam files");
		File tophatHeaderFile = new File(tophatHeader);
		if(tophatHeaderFile.exists()) {
			logger.warn("Header file " + tophatHeader + " already exists. Not remaking header.");
		} else {
			logger.info("Running command: " + getHeaderCmmd);
			logger.info("LSF job ID is " + getHeaderJobID + ".");
			PipelineUtils.bsubProcess(Runtime.getRuntime(), getHeaderJobID, getHeaderCmmd, NOVOALIGN_DIRECTORY + "/get_sam_header_" + getHeaderJobID + ".bsub", "hour", 1);
			logger.info("Waiting for samtools view to finish...");
			PipelineUtils.waitForJobs(getHeaderJobID, Runtime.getRuntime());
			// Change sort order to unsorted
			FileReader r = new FileReader(tmpHeader);
			BufferedReader b = new BufferedReader(r);
			FileWriter wr = new FileWriter(tophatHeader);
			while(b.ready()) {
				String line = b.readLine();
				wr.write(line.replaceAll("coordinate", "unsorted") + "\n");
			}
			r.close();
			b.close();
			wr.close();
			File f = new File(tmpHeader);
			f.delete();
		}
		logger.info("Done getting header.");
		
		// Now reheader each novoalign sam

		for(String sample : sampleNames) {
			File finalBam = new File(novoBamFinalPath.get(sample));
			if(finalBam.exists()) {
				logger.warn("Alignment file " + finalBam + " already exists. Not looking for sam file or replacing header.");
				continue;					
			}
			// Cat the new header and the sam file
			String novo = novoSamOutput.get(sample);
			String reheadered = novoSamOutputReheadered.get(sample);
			String noheader = novoSamOutputNoHeader.get(sample);
			// Get the novo sam file without header
			String cmmd = samtools + " view -S -o " + noheader + " " + novo;
			Process p = Runtime.getRuntime().exec(cmmd);
			p.waitFor();
			// Cat the files
			FileReader r1 = new FileReader(tophatHeader);
			BufferedReader b1 = new BufferedReader(r1);
			FileReader r2 = new FileReader(noheader);
			BufferedReader b2 = new BufferedReader(r2);
			FileWriter w1 = new FileWriter(reheadered);
			while(b1.ready()) w1.write(b1.readLine() + "\n");
			while(b2.ready()) w1.write(b2.readLine() + "\n");
			r1.close();
			r2.close();
			b1.close();
			b2.close();
			w1.close();
			
			// Replace old sam file with reheadered file
			String mvcmmd = "mv " + reheadered + " " + novo;
			Process mvp = Runtime.getRuntime().exec(mvcmmd);
			mvp.waitFor();
			// Remove the no header file
			File noheaderfile = new File(noheader);
			noheaderfile.delete();
		}
	}
	
	/**
	 * Convert sam files to bam files
	 * @param samFiles The sam files by sample name
	 * @param bamFiles The bam files to write, by sample name
	 * @param finalBamFiles Final bam files; skip if they already exist
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param samtools Samtools executables
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void samToBam(Map<String, String> samFiles, Map<String, String> bamFiles, Map<String, String> finalBamFiles, Map<String, String> bsubOutputDirs, String samtools) throws IOException, InterruptedException {
		ArrayList<String> cbJobIDs = new ArrayList<String>();
		for(String sample : sampleNames) {
			File bam = new File(bamFiles.get(sample));
			// Check if bam files already exist
			if(bam.exists()) {
				logger.warn("Bam file " + bam + " already exists. Not redoing format conversion from sam.");
				continue;
			}
			File finalBam = new File(finalBamFiles.get(sample));
			if(finalBam.exists()) {
				logger.warn("Alignment file " + finalBam + " already exists. Not looking for or converting sam file.");
				continue;					
			}
			// Use samtools view
			String cmmd = samtools + " view -Sb -o " + bamFiles.get(sample) + " " + samFiles.get(sample);
			logger.info("Running Samtools command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			cbJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bsubOutputDirs.get(sample) + "/sam_to_bam_" + jobID + ".bsub", "hour", 1);
		}
		// Wait for jobs to finish
		logger.info("Waiting for samtools view jobs to finish...");
		PipelineUtils.waitForAllJobs(cbJobIDs, Runtime.getRuntime());
	}
	
	/**
	 * Sort all bam files and write sorted files to specified locations
	 * @param unsortedBams The unsorted bam files to sort, by sample name
	 * @param sortedBams The sorted bam files to write, by sample name
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param finalBams Final bam files; skip if they already exist
	 * @param picardJarDir Directory containing Picard jar files
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void sortBamFiles(Map<String, String> unsortedBams, Map<String, String> sortedBams, Map<String, String> bsubOutputDirs, Map<String, String> finalBams, String picardJarDir) throws IOException, InterruptedException {
		ArrayList<String> sbJobIDs = new ArrayList<String>();
		for(String sample : sampleNames) {
			File finalBam = new File(finalBams.get(sample));
			if(finalBam.exists()) {
				logger.warn("Alignment file " + finalBam + " already exists. Not sorting.");
				continue;					
			}
			// Use Picard program SortSam
			String cmmd = "java -jar " + picardJarDir + "/SortSam.jar INPUT=" + unsortedBams.get(sample) + " OUTPUT=" + sortedBams.get(sample) + " SORT_ORDER=coordinate";
			logger.info("Running Picard command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			sbJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bsubOutputDirs.get(sample) + "/sort_bam_" + jobID + ".bsub", "hour", 4);
		}
		// Wait for jobs to finish
		logger.info("Waiting for SortSam jobs to finish...");
		PipelineUtils.waitForAllJobs(sbJobIDs, Runtime.getRuntime());
	}
	
	/**
	 * Index a fasta file using samtools faidx
	 * @param fastaFileName The fasta file to index
	 * @param samtoolsExecutable Path to samtools executable
	 * @param bsubOutDir Directory to write bsub output to
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void indexFastaFile(String fastaFileName, String samtoolsExecutable, String bsubOutDir) throws IOException, InterruptedException {
		String cmmd = samtoolsExecutable + " faidx " + fastaFileName;
		logger.info("Running samtools command: " + cmmd);
		String jobID = Long.valueOf(System.currentTimeMillis()).toString();
		logger.info("LSF job ID is " + jobID + ".");
		PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/index_fasta_" + jobID + ".bsub", "hour", 4);
		logger.info("Waiting for samtools faidx job to finish...");
		PipelineUtils.waitForJobs(jobID, Runtime.getRuntime());
	}
	
	/**
	 * Merge tophat and novoalign bam files
	 * @param tophatBamFinalPath Mapping of sample name to final tophat bam file
	 * @param novoBamFinalPath Mapping of sample name to final novoalign bam file
	 * @param picardJarDir Directory containing Picard jar files
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void mergeTophatNovoalign(Map<String, String> tophatBamFinalPath, Map<String, String> novoBamFinalPath, String picardJarDir) throws IOException, InterruptedException {
		File mergedDir = new File(MERGED_TOPHAT_NOVOALIGN_DIRECTORY);
		boolean madeMergedDir = mergedDir.mkdir();
		ArrayList<String> mergeJobIDs = new ArrayList<String>();
		for(String sample : sampleNames) {
			String tophatBam = tophatBamFinalPath.get(sample);
			String novoBam = novoBamFinalPath.get(sample); 
			String mergedBam = MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam";
			File mergedFile = new File(mergedBam);
			// Check if merged files already exist
			if(mergedFile.exists()) {
				logger.warn("Merged bam file " + mergedBam + " already exists. Not re-merging tophat and novoalign files.");
				continue;					
			}
			// Use Picard program MergeSamFiles
			String cmmd = "java -jar " + picardJarDir + "/MergeSamFiles.jar INPUT=" + tophatBam + " INPUT=" + novoBam + " OUTPUT=" + mergedBam;
			logger.info("Running Picard command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			mergeJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/merge_bams_" + jobID + ".bsub", "hour", 1);
		}
		// Wait for jobs to finish
		logger.info("Waiting for MergeSamFiles jobs to finish...");
		PipelineUtils.waitForAllJobs(mergeJobIDs, Runtime.getRuntime());
	}
	
	/**
	 * Count merged tophat and novoalign alignments
	 * @param samtools Samtools executable
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void countMappingsMergedAlignments(String samtools) throws IOException, InterruptedException {
		String mergedCountFile = MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/mapped_unmapped_count.out";
		String mergedPctFile = MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/mapped_unmapped_percentage.out";
		FileWriter mw = new FileWriter(mergedCountFile);
		FileWriter mwp = new FileWriter(mergedPctFile);
		String mergedHeader = "Sample\tMapped\tUnmapped\n";
		mw.write(mergedHeader);
		mwp.write(mergedHeader);
		for(String sample : sampleNames) {
			int mapped = AlignmentUtils.countAlignments(samtools, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam", MERGED_TOPHAT_NOVOALIGN_DIRECTORY, false, false);
			int unmapped = AlignmentUtils.countAlignments(samtools, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/" + sample + ".bam", MERGED_TOPHAT_NOVOALIGN_DIRECTORY, true, false);
			int total = mapped + unmapped;
			mw.write(sample + "\t" + mapped + "\t" + unmapped + "\n");
			mwp.write(sample + "\t" + (double)mapped/(double)total + "\t" + (double)unmapped/(double)total + "\n");
		}
		mw.close();
		mwp.close();
		logger.info("Wrote table of counts to file " + mergedCountFile);
		logger.info("Wrote table of percentages to file " + mergedPctFile);
	}
	
	/**
	 * Index current bam files and write bai files to current bam directory
	 * @param samtools Samtools executable
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void indexCurrentBams(String samtools) throws IOException, InterruptedException {
		indexBamFiles(currentBamFiles, samtools);
	}
	
	/**
	 * Write fragment end points to wig and bigwig files
	 * @param bamFiles Bam files by sample name
	 * @param bamDir Directory containing bam files
	 * @param refFasta Fasta file of sequences these bam files were aligned against
	 * @param geneBedFile Bed file of genes to count reads in or null if using genomic space
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void writeWigFragmentEnds(Map<String, String> bamFiles, String bamDir, String refFasta, String geneBedFile) throws IOException, InterruptedException {
		
		Map<String, String> read1wig = new TreeMap<String, String>();
		Map<String, String> read2wig = new TreeMap<String, String>();
		Map<String, String> read1bigwig = new TreeMap<String, String>();
		Map<String, String> read2bigwig = new TreeMap<String, String>();
		String wigToBigWig = configP.basicOptions.getWigToBigWigExecutable();
		ArrayList<String> bigwigJobIDs = new ArrayList<String>();
		
		// Chromosome size file to pass to UCSC program wigToBigWig
		String chrSizesForWigToBigWig = refFasta + ".sizes";
		
		// Chromosome size file to pass to wig writer if using genomic space
		// Null if using transcriptome space
		String chrSizesForWigWriter = null;
		if(geneBedFile == null) {
			chrSizesForWigWriter = chrSizesForWigToBigWig;
		}
		File chrSizeFile = new File(chrSizesForWigToBigWig);
		
		// Write chromosome size file
		if(!chrSizeFile.exists()) {
			FileWriter w = new FileWriter(chrSizesForWigToBigWig);
			logger.info("Writing chromosome sizes to file " + chrSizesForWigToBigWig);
			FastaSequenceIO fsio = new FastaSequenceIO(refFasta);
			Collection<Sequence> seqs = fsio.loadAll();
			for(Sequence seq : seqs) {
				w.write(seq.getId() + "\t" + seq.getLength() + "\n");
			}
			w.close();
		}
		
		
		for(String sampleName : sampleNames) {
			
			String bamFile = bamFiles.get(sampleName);
			read1wig.put(sampleName, bamFile + ".read1.wig");
			read2wig.put(sampleName, bamFile + ".read2.wig");
			read1bigwig.put(sampleName, bamFile + ".read1.bw");
			read2bigwig.put(sampleName, bamFile + ".read2.bw");
			
			// Write wig file for read1
			String wig1 = read1wig.get(sampleName);
			File read1wigFile = new File(wig1);
			if(read1wigFile.exists()) {
				logger.warn("Wig file " + wig1 + " already exists. Not remaking file.");
			} else {
				logger.info("Writing fragment ends of read 1 from bam file " + bamFile + " to wig file " + wig1 + ".");
				WigWriter read1ww = new WigWriter(bamFile, geneBedFile, chrSizesForWigWriter, true, false);
				read1ww.addReadFilter(new FirstOfPairFilter());
				read1ww.addReadFilter(new ProperPairFilter());
				read1ww.writeWig(wig1);
				logger.info("Done writing file " + wig1 + ".");
			}
			// Write bigwig file for read1
			String bigwig1 = read1bigwig.get(sampleName);
			File read1bigwigFile = new File(bigwig1);
			if(read1bigwigFile.exists()) {
				logger.warn("Bigwig file " + bigwig1 + " already exists. Not remaking file.");
			} else {
				String cmmd = wigToBigWig + " " + wig1 + " " + chrSizesForWigToBigWig + " " + bigwig1;
				logger.info("");
				logger.info("Making bigwig file for wig file " + wig1 + ".");
				logger.info("Running UCSC command " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				bigwigJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_" + jobID + ".bsub", "hour", 4);
			}
			
			if(pairedData.get(sampleName).booleanValue()) {
				// Write wig file for read2
				String wig2 = read2wig.get(sampleName);
				File read2wigFile = new File(wig2);
				if(read2wigFile.exists()) {
					logger.warn("Wig file " + wig2 + " already exists. Not remaking file.");
				} else {
					logger.info("Writing fragment ends of read 2 from bam file " + bamFile + " to wig file " + wig2 + ".");
					WigWriter read2ww = new WigWriter(bamFile, geneBedFile, chrSizesForWigWriter, true, false);
					read2ww.addReadFilter(new SecondOfPairFilter());
					read2ww.addReadFilter(new ProperPairFilter());
					read2ww.writeWig(wig2);
					logger.info("Done writing file " + wig2 + ".");
				}
				// Write bigwig file for read2
				String bigwig2 = read2bigwig.get(sampleName);
				File read2bigwigFile = new File(bigwig2);
				if(read2bigwigFile.exists()) {
					logger.warn("Bigwig file " + bigwig2 + " already exists. Not remaking file.");
				} else {
					String cmmd = wigToBigWig + " " + wig2 + " " + chrSizesForWigToBigWig + " " + bigwig2;
					logger.info("");
					logger.info("Making bigwig file for wig file " + wig2 + ".");
					logger.info("Running UCSC command " + cmmd);
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					bigwigJobIDs.add(jobID);
					logger.info("LSF job ID is " + jobID + ".");
					PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_" + jobID + ".bsub", "hour", 4);
				}
			}
		
			logger.info("Waiting for wigToBigWig jobs to finish...");
			PipelineUtils.waitForAllJobs(bigwigJobIDs, Runtime.getRuntime());
			
		}
		
	}
	
	/**
	 * Index a set of bam files
	 * @param bamFiles Map of sample name to bam file to index
	 * @param outDir Output directory for all indexed bam files
	 * @param samtools Samtools executable
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void indexBamFiles(Map<String, String> bamFiles, String samtools) throws IOException, InterruptedException {
		ArrayList<String> indexJobIDs = new ArrayList<String>();
		for(String sample : sampleNames) {
			String bam = bamFiles.get(sample);
			String index = bam + ".bai";
			File indexfile = new File(index);
			// Check if bai files exist
			if(indexfile.exists()) {
				logger.warn("Index " + index + " already exists. Not re-indexing bam file.");
				continue;
			}
			String cmmd = samtools + " index " + bam + " " + index;
			logger.info("Running samtools command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			indexJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, indexfile.getParent() + "/index_bam_" + jobID + ".bsub", "hour", 1);
		}
		logger.info("Waiting for samtools jobs to finish...");
		PipelineUtils.waitForAllJobs(indexJobIDs, Runtime.getRuntime());
	}
	
	/**
	 * Reorder current bam files to match reference genome
	 * @param picardDir Directory containing Picard jar files
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void reorderCurrentBams(String picardDir) throws IOException, InterruptedException {
		ArrayList<String> reorderJobIDs = new ArrayList<String>();
		Map<String, String> reordered = new TreeMap<String, String>();
		for(String sample : sampleNames) {
			String bam = currentBamFiles.get(sample);
			reordered.put(sample, bam + ".reordered");
			String cmmd = "java -jar " + picardDir + "/ReorderSam.jar I=" + bam + " O=" + reordered.get(sample) + " R=" + configP.basicOptions.getGenomeFasta();
			logger.info("Running picard command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			reorderJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, currentBamDir + "/reorder_bam_" + jobID + ".bsub", "week", 16);
		}
		logger.info("Waiting for picard jobs to finish...");
		PipelineUtils.waitForAllJobs(reorderJobIDs, Runtime.getRuntime());
		// Replace bam files with reordered files
		for(String sample : sampleNames) {
			String mvcmmd = "mv " + reordered.get(sample) + " " + currentBamFiles.get(sample);
			Process p = Runtime.getRuntime().exec(mvcmmd);
			p.waitFor();
		}
	}
	
	/**
	 * Write global stats for current bam files
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bedFile Bed file for transcriptome space stats. To skip, pass null.
	 * @param chrSizeFile Chromosome size file for genomic space stats. To skip, pass null.
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void writeAlignmentGlobalStats(String alignmentGlobalStatsJar, String bedFile, String chrSizeFile) throws IOException, InterruptedException {
		logger.info("Writing global stats for alignments...");
		ArrayList<String> jobIDs = new ArrayList<String>();
		if(bedFile != null) {
			Collection<String> tJobIDs = writeTranscriptomeSpaceStats(currentBamFiles, bedFile, alignmentGlobalStatsJar);
			jobIDs.addAll(tJobIDs);
		}
		if(chrSizeFile != null) {
			Collection<String> gJobIDs = writeGenomicSpaceStats(currentBamFiles, chrSizeFile, alignmentGlobalStatsJar);
			jobIDs.addAll(gJobIDs);
		}
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		logger.info("Done writing all global alignment stats.");
	}
	
	/**
	 * Write global transcriptome space stats for bam files
	 * @param bamFiles Bam files by sample name
	 * @param bedFile Bed annotation for transcriptome space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	private Collection<String> writeTranscriptomeSpaceStats(Map<String, String> bamFiles, String bedFile, String alignmentGlobalStatsJar) throws IOException, InterruptedException {
		logger.info("Calculating transcriptome space stats from annotation in " + bedFile);
		ArrayList<String> jobIDs = new ArrayList<String>();
		for(String sampleName : bamFiles.keySet()) {
			String bamFile = bamFiles.get(sampleName);
			logger.info("Calculating transcriptome space stats for sample " + sampleName + "...");
			String cmmd = "java -jar -Xmx30g -Xms20g -Xmn10g " + alignmentGlobalStatsJar + " -b " + bamFile + " -t " + bedFile;
			logger.info("Running command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			jobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, currentBamDir + "/compute_transcriptome_space_stats_" + jobID + ".bsub", "week", 32);		
		}
		return jobIDs;
	}

	/**
	 * Write global genomic space stats for bam files
	 * @param bamFiles Bam files by sample name
	 * @param chrSizeFile Chromosome size file for genomic space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	private Collection<String> writeGenomicSpaceStats(Map<String, String> bamFiles, String chrSizeFile, String alignmentGlobalStatsJar) throws IOException, InterruptedException {
		logger.info("Calculating genomic space stats from chromosome sizes in " + chrSizeFile);
		ArrayList<String> jobIDs = new ArrayList<String>();
		for(String sampleName : bamFiles.keySet()) {
			String bamFile = bamFiles.get(sampleName);
			logger.info("Calculating genomic space stats for sample " + sampleName + "...");
			String cmmd = "java -jar -Xmx30g -Xms20g -Xmn10g " + alignmentGlobalStatsJar + " -b " + bamFile + " -g " + chrSizeFile;
			logger.info("Running command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			jobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, currentBamDir + "/compute_genomic_space_stats_" + jobID + ".bsub", "week", 32);		
		}
		return jobIDs;
	}


	/**
	 * Merge specified samples into new bam files
	 * Add merged samples to sample name list and current bam files
	 * @param picardJarDir Directory containing Picard executables
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void mergeBamFiles(String picardJarDir) throws IOException, InterruptedException {
		Map<String, Collection<String>> setsToMerge = configP.basicOptions.getSamplesToMerge();
		if(setsToMerge.isEmpty()) {
			logger.info("No sets to merge.");
			return;
		}
		
		// Make sure all samples exist
		for(String mergedName : setsToMerge.keySet()) {
			for(String sampleName : setsToMerge.get(mergedName)) {
				if(!sampleNames.contains(sampleName)) {
					throw new IllegalArgumentException("Can't merge sample: " + sampleName + ". Sample does not exist.");
				}
			}
		}
		
		Map<String, String> mergedFiles = new TreeMap<String, String>();
		Map<String, Boolean> paired = new TreeMap<String, Boolean>();
		for(String mergedName : setsToMerge.keySet()) {
			String mergedBam = currentBamDir + "/" + mergedName + ".bam";
			mergedFiles.put(mergedName, mergedBam);
		}
		ArrayList<String> jobIDs = new ArrayList<String>();
		for(String mergedName : setsToMerge.keySet()) {
			
			// Check if all files to merge are same sequencing format
			Collection<String> samplesToMerge = setsToMerge.get(mergedName);
			boolean isPaired = pairedData.get(samplesToMerge.iterator().next()).booleanValue();
			for(String sample : samplesToMerge) {
				if(pairedData.get(sample).booleanValue() != isPaired) {
					throw new IllegalArgumentException("All samples to merge must be same format (paired or unpaired)");
				}
			}
			paired.put(mergedName, Boolean.valueOf(isPaired));
			
			File file = new File(mergedFiles.get(mergedName));
			if(file.exists()) {
				logger.info("Merged bam file " + file + " already exists. Not rerunning bam file merge.");
				continue;
			}
			logger.info("Creating merged bam file " + file + ".");
			String inputs = "";
			for(String sample : samplesToMerge) {
				inputs += "INPUT=" + currentBamFiles.get(sample) + " ";
			}
			String output = "OUTPUT=" + mergedFiles.get(mergedName);
			String cmmd = "java -jar " + picardJarDir + "/MergeSamFiles.jar " + inputs + " " + output + " ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true";
			logger.info("Running picard command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			jobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, currentBamDir + "/merge_bam_files_" + jobID + ".bsub", "week", 8);			
		}
		logger.info("Waiting for picard jobs to finish...");
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		
		// Update current bam files and sample names
		currentBamFiles.putAll(mergedFiles);
		sampleNames.addAll(mergedFiles.keySet());
		pairedData.putAll(paired);
	}
	
	/**
	 * Make tdf files for bam files
	 * @param bamFilesBySampleName Bam file name by sample name
	 * @param bamDirectory Directory containing bam files
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void makeTdfs(Map<String, String> bamFilesBySampleName, String bamDirectory, String assemblyFasta, String igvtoolsExecutable) throws IOException, InterruptedException {
		ArrayList<String> tdfJobIDs = new ArrayList<String>();
		Map<String, String> outTdf = new TreeMap<String, String>();
		for(String sample : sampleNames) {
			String bam = bamFilesBySampleName.get(sample);
			String tdf = bam + ".tdf";
			outTdf.put(sample, tdf);
			File tdffile = new File(tdf);
			// Check if tdf files exist
			if(tdffile.exists()) {
				logger.warn("Tdf file " + tdf + " already exists. Not remaking tdf file.");
				continue;
			}
			// Use igvtools count
			String cmmd = igvtoolsExecutable + " count -w 3 " + bam + " " + tdf + " " + assemblyFasta;
			logger.info("Running igvtools command: " + cmmd);
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			tdfJobIDs.add(jobID);
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bamDirectory + "/make_tdf_" + jobID + ".bsub", "hour", 1);
		}
		if(tdfJobIDs.isEmpty()) return;
		logger.info("Waiting for igvtools jobs to finish...");
		logger.info("Note: igvtools count always exits with code -1 even though it worked, so disregard failure notifications from LSF.");
		// Igvtools count always ends by crashing even though it worked, so catch the exception and check if files were really created
		try {
			PipelineUtils.waitForAllJobs(tdfJobIDs, Runtime.getRuntime());
		} catch(IllegalArgumentException e) {
			boolean ok = true;
			String errMsg = "";
			for(String sample : sampleNames) {
				String tdf = outTdf.get(sample);
				File tdffile = new File(tdf);
				if(!tdffile.exists()) {
					ok = false;
					errMsg += "Tdf file " + tdf + " does not exist.";
					break;
				}
			}
			if(!ok) {
				throw new IllegalArgumentException(errMsg);
			}
		}

	}
	
	/**
	 * Write fragment size distribution for all paired bam files
	 * @throws IOException
	 */
	private void makeFragmentSizeDistributionCurrentBams() throws IOException {
		
		String distFileName = currentBamDir + "/fragment_size_histogram";
		File distFile = new File(distFileName);
		String medianFileName = currentBamDir + "/fragment_size_median";
		File medianFile = new File(medianFileName);
		String indGeneFileName = currentBamDir + "/fragment_size_median_individual_genes";
		File indGeneFile = new File(indGeneFileName);
		if(distFile.exists() && medianFile.exists() && indGeneFile.exists()) {
			logger.warn("Fragment size distribution file " + distFileName + " and median files " + medianFileName + " and " + indGeneFileName + " already exist. Not remaking files.");
			return;
		}
		
		String annotation = configP.fragmentSizeOptions.getBedFile();
		TranscriptomeSpace coord = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(annotation)));
		int maxFragmentSize = configP.fragmentSizeOptions.getMaxFragmentSize();
		int maxGenomicSpan = configP.fragmentSizeOptions.getMaxGenomicSpan();
		int numBins = configP.fragmentSizeOptions.getNumBins();
		boolean properPairsOnly = configP.fragmentSizeOptions.getProperPairsOnly();
		Map<String, EmpiricalDistribution> distributions = new TreeMap<String, EmpiricalDistribution>();
		Map<String, ScanStatisticDataAlignmentModel> data = new TreeMap<String, ScanStatisticDataAlignmentModel>();
		
		for(String sample : sampleNames) {
			if(!pairedData.get(sample).booleanValue()) {
				logger.info("Not making fragment size histogram for sample " + sample + " because data is not paired.");
				continue;
			}
			logger.info("Making fragment size distribution for sample " + sample + ".");
			logger.info("Max fragment size: " + maxFragmentSize);
			logger.info("Max genomic span: " + maxGenomicSpan);
			logger.info("Number of bins: " + numBins);
			String bam = currentBamFiles.get(sample);
			ScanStatisticDataAlignmentModel d = new ScanStatisticDataAlignmentModel(bam,coord);
			d.addFilter(new GenomicSpanFilter(maxGenomicSpan));
			if(properPairsOnly) d.addFilter(new ProperPairFilter());
			data.put(sample, d);
			distributions.put(sample, d.getReadSizeDistribution(coord, maxFragmentSize, numBins));
		}
		if(distributions.isEmpty()) return;

		logger.info("Writing fragment size distributions to file " + distFileName);
		FileWriter w = new FileWriter(distFileName);
		
		String header = "bin\t";
		for(String sample : distributions.keySet()) {
			header += sample + "\t";
		}
		header += "\n";
		w.write(header);
		
		for(int i=0; i<numBins; i++) {
			String firstSample = distributions.keySet().iterator().next();
			String line = (int)distributions.get(firstSample).getBinStart(i) + "-" + (int)distributions.get(firstSample).getBinEnd(i) + "\t";
			for(String sample : distributions.keySet()) {
				line += distributions.get(sample).getHistogram(i) + "\t";
			}
			line += "\n";
			w.write(line);
		}
		
		w.close();
		
		logger.info("Writing median fragment sizes to file " + medianFileName);
		FileWriter w2 = new FileWriter(medianFileName);
		for(String sample : distributions.keySet()) {
			w2.write(sample + "\t" + distributions.get(sample).getMedianOfAllDataValues() + "\n");
		}
		w2.close();
		
		logger.info("Writing individual gene medians to file " + indGeneFileName);
		FileWriter w3 = new FileWriter(indGeneFileName);
		Collection<String> indGeneNames = configP.fragmentSizeOptions.getIndividualGeneNames();
		if(!indGeneNames.isEmpty()) {
			Map<String, Gene> genesByName = BEDFileParser.loadDataByName(new File(annotation));
			String header2 = "Gene\t";
			for(String sampleName : sampleNames) {
				header2 += sampleName + "\t";
			}
			w3.write(header2 + "\n");
			for(String geneName : indGeneNames) {
				logger.info(geneName);
				String line = geneName + "\t";
				for(String sampleName : sampleNames) {
					String value;
					try {
						value = Double.valueOf(data.get(sampleName).getReadSizeDistribution(genesByName.get(geneName), coord, maxFragmentSize, numBins).getMedianOfAllDataValues()).toString();
						line += value + "\t";
					} catch (IllegalStateException e) {
						line += "NA\t";
					}
				}
				w3.write(line + "\n");
			}
		}
		w3.close();
		
	}
	
	
	/**
	 * Collect several sets of Picard metrics
	 * @param picardExecutableDir Directory containing Picard jar files
	 * @param picardMetricsDir Directory to write metrics to
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void collectPicardMetricsCurrentBams(String picardExecutableDir) throws IOException, InterruptedException {
		String picardMetricsDir = currentBamDir + "/picard_metrics";
		logger.info("Writing all Picard metrics in directory " + picardMetricsDir + "....");
		File picardMetricsDirFile = new File(picardMetricsDir);
		boolean madePicardMetricsDir = picardMetricsDirFile.mkdir();
		ArrayList<String> pmJobIDs = new ArrayList<String>();
		// Inputs to picard metrics
		String refFlat = configP.basicOptions.getPicardRefFlat();
		String ribIntervals = configP.basicOptions.getPicardRibosomalIntervals();
		String genomeFasta = configP.basicOptions.getGenomeFasta();
		String strandSpecificity = configP.basicOptions.getPicardStrandSpecificity();
		
		for(String sample : sampleNames) {
			
			// Establish output files
			String asMetrics = picardMetricsDir + "/" + sample + "_alignmentSummaryMetrics.out";
			String isMetrics = picardMetricsDir + "/" + sample + "_insertSizeMetrics.out";
			String isHistogram = picardMetricsDir + "/" + sample + "_insertSizeMetrics.histogram";
			String rsMetrics = picardMetricsDir + "/" + sample + "_rnaSeqMetrics.out";
			String rsChart = picardMetricsDir + "/" + sample + "_rnaSeqMetrics_NormalizedPositionCoverage.pdf";
			File asFile = new File(asMetrics);
			File isFile = new File(isMetrics);
			File ishFile = new File(isHistogram);
			File rsFile = new File(rsMetrics);
			File rscFile = new File(rsChart);
			
			// Run alignment summary metrics
			if(asFile.exists()) {
				logger.warn("Alignment summary metrics file " + asMetrics + " already exists. Not rerunning alignment summary metrics.");
			} else {
				// Use Picard program CollectAlignmentSummaryMetrics
				String cmmd = "java -jar " + picardExecutableDir + "/CollectAlignmentSummaryMetrics.jar";
				cmmd += " INPUT=" + currentBamFiles.get(sample);
				cmmd += " OUTPUT=" + asMetrics;
				cmmd += " REFERENCE_SEQUENCE=" + genomeFasta;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				pmJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_alignment_summary_metrics_" + jobID + ".bsub", "hour", 4);
			}
			
			// Run insert size metrics
			if(isFile.exists() && ishFile.exists()) {
				logger.warn("Insert size metrics files " + isMetrics + " and " + isHistogram + " already exist. Not rerunning insert size metrics.");
			} else {
				// Use Picard program CollectInsertSizeMetrics
				String cmmd = "java -jar " + picardExecutableDir + "/CollectInsertSizeMetrics.jar";
				cmmd += " INPUT=" + currentBamFiles.get(sample);
				cmmd += " OUTPUT=" + isMetrics;
				cmmd += " REFERENCE_SEQUENCE=" + genomeFasta;
				cmmd += " HISTOGRAM_FILE=" + isHistogram;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				pmJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_insert_size_metrics_" + jobID + ".bsub", "hour", 4);
			}
			
			// Run RNA-seq metrics
			if(rsFile.exists() && rscFile.exists()) {
				logger.warn("RNA-seq metrics files " + rsMetrics + " and " + rsChart + " already exist. Not rerunning RNA-seq metrics.");
			} else {
				// Use Picard program CollectRnaSeqMetrics
				String cmmd = "java -jar ";
				// Use user-provided CollectRnaSeqMetrics executable if provided
				if(configP.basicOptions.getPicardCollectRnaSeqMetrics() != null) cmmd += configP.basicOptions.getPicardCollectRnaSeqMetrics();
				else cmmd += picardExecutableDir + "/CollectRnaSeqMetrics.jar";
				cmmd += " INPUT=" + currentBamFiles.get(sample);
				cmmd += " OUTPUT=" + rsMetrics;
				cmmd += " REFERENCE_SEQUENCE=" + genomeFasta;
				cmmd += " CHART_OUTPUT=" + rsChart;
				cmmd += " REF_FLAT=" + refFlat;
				if(ribIntervals != null) cmmd += " RIBOSOMAL_INTERVALS=" + ribIntervals;
				cmmd += " STRAND_SPECIFICITY=" + strandSpecificity;
				logger.info("Running Picard command: " + cmmd);
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				pmJobIDs.add(jobID);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_rnaseq_metrics_" + jobID + ".bsub", "hour", 4);
			}
		}
		
		// Wait for jobs to finish
		logger.info("Waiting for Picard metrics jobs to finish...");
		try {
			PipelineUtils.waitForAllJobs(pmJobIDs, Runtime.getRuntime());
		} catch(IllegalArgumentException e) {
			logger.info("");
			logger.warn("Caught exception; not all Picard metrics jobs completed successfully.");
			logger.info("");
		}
	}
	
	
	private String createDGEInputFile(String inputListFile) throws IOException{
		
		String newInputFile = inputListFile+".dgeInputFile";
		FileWriter w = new FileWriter(newInputFile);
		for(String name: currentBamFiles.keySet()){
			w.write(name+"\t"+currentBamFiles.get(name)+"\t"+nameToCondition.get(name)+"\n");
		}
		w.close();
		
		return newInputFile;
	}

	/**
	 * TASK 6: RUN DGE
	 * @author skadri
	 * @throws MathException 
	 * @throws ParseException 
	 * @throws IOException 
	 */
	private void runDGE(boolean runBasic,String DGEInputFile) throws IOException, ParseException{
		throw new UnsupportedOperationException("This needs to be fixed!");
		
		//ALIGNMENTS
		/*if("multiple".equalsIgnoreCase(configP.dgeOptions.getrunType())){
			
			 //Construct the DGE argument list
			
			ArrayList<String> arguments =  new ArrayList<String>();
			//3'DGE or 5'DGE
			//single/multiple
			//THis determines task
			arguments.add("-task");
			if("3DGE".equalsIgnoreCase(configP.dgeOptions.getdgeType())){
				arguments.add("score3PMultiple");
			}
			else if("5DGE".equalsIgnoreCase(configP.dgeOptions.getdgeType())){
				arguments.add("score5PMultiple");
			}
			else
				throw new IllegalArgumentException("Illegal DGE type (3DGE/5DGE): "+configP.dgeOptions.getdgeType());
			
			//Changing the output to a directory
			configP.dgeOptions.setOutputInDirectory();
			
			arguments.add("-alignments");
			arguments.add(DGEInputFile);
			for(String flag:configP.getRunSpecificDGEOptions()){
				arguments.add("-"+flag);
				arguments.add(configP.dgeOptions.getOption(flag));
			}
			logger.info("Parameters to DGE:" + arguments);
			DGE dge = new DGE(l2a(arguments));
		}
		else if("single".equalsIgnoreCase(configP.dgeOptions.getrunType())){
			DGEInput inp = readDGEInputFile(DGEInputFile);
			
			//DGE arguments
			String optionsString = "java -jar "+configP.dgeOptions.getdgePath()+" -task ";
			if("3DGE".equalsIgnoreCase(configP.dgeOptions.getdgeType())){
				optionsString += "score3P ";
			}
			else if("5DGE".equalsIgnoreCase(configP.dgeOptions.getdgeType())){
				optionsString += "score5P ";
			}
			else{
				throw new IllegalArgumentException("Illegal DGE type (3DGE/5DGE): "+configP.dgeOptions.getdgeType());
			}
			for(String flag:configP.getRunSpecificDGEOptions()){
				if(flag.equalsIgnoreCase("out"))
					;//dont add
				else{
					optionsString += "-"+flag+" ";
					optionsString += configP.dgeOptions.getOption(flag)+" ";
				}
			}
			
			//OUTPUT FILE NAMES ARE DIFF
			for(String sample:inp.getSamples()){
				//MULTIPLE BSUB ROUTINES FOR EACH DGE
				String thisCommand = optionsString;
				//ALIGNMENT
				thisCommand +="-alignment ";
				thisCommand += inp.getBamFileFor(sample)+" ";
				//OTHER FLAGS
				
				//Changing the output to a directory
				String outputName = inp.getBamFileFor(sample)+".dge/"+inp.getBamFileFor(sample)+".dge.exp";
				thisCommand += "-out ";
				thisCommand += outputName;
				
				//Submit the job
				//PipelineUtils.
			}
		}
		else{
			throw new IllegalArgumentException("Illegal Run type (single/multiple): "+configP.dgeOptions.getrunType());
		}*/
	}

	/**
	 * Reads the "Name		Bam_file	condition" list file
	 * Helper function to runDGE()
	 * @param inputListFile
	 * @throws FileNotFoundException 
	 */
	private static DGEInput readDGEInputFile(String inputListFile) throws FileNotFoundException {
		
		DGEInput dge = new DGEInput();
		Scanner reader = new Scanner(new File(inputListFile));
		while (reader.hasNextLine()) {
			String[] str = StringParser.getTokens(reader.nextLine());
			dge.addInput(str);
		}
		return dge;
	}


	/**
	 * Helper function: Converts list to array and returns the array
	 * @param list
	 * @return
	 */
	private static String[] l2a(List<String> list){
		String[] rtrn=new String[list.size()];
		int i=0;
		for(String val: list){rtrn[i++]=val;}
		return rtrn;
	}

	public static void main (String [] args) throws IOException, ParseException, InterruptedException{
		
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("\n*** Configurable pipeline for RNA-seq read processing and analysis ***");
		p.addStringArg("-r", "File containing list of fastq files. \n\t\tLine format: \n\t\t<sample_name> <left.fq> <right.fq> <condition> \n\t\tOR \n\t\t<sample_name> <unpaired.fq> <condition>", true);
		p.addStringArg("-c", "Config file", true);
		p.parse(args);
		String fastqList = p.getStringArg("-r");
		String configFile = p.getStringArg("-c");
		
		PipelineAutomator PA = new PipelineAutomator(fastqList,configFile);
		logger.info("Pipeline all done.");
	}

	public static class DGEInput {
		Map<String,String> nameToBam;
		Map<String,String> nameToCondition;
		
		public DGEInput(){
			nameToBam = new HashMap<String,String>();
			nameToCondition = new HashMap<String,String>();
		}
		
		/**
		 * 
		 * @param str
		 */
		public void addInput(String[] str){
			//Mandatory 3 columns are required.
			if(str.length<3){
				throw new IllegalArgumentException("Illegal number of columns in input DGE");
			}
			nameToBam.put(str[0], str[1]);
			nameToCondition.put(str[0], str[2]);
		}
		
		public Set<String> getSamples(){
			return nameToBam.keySet();
		}
		
		public String getBamFileFor(String name){
			return nameToBam.get(name);
		}
	}
	
	/**
	 * Bowtie command lines needed for this pipeline
	 * @author prussell
	 *
	 */
	public static class AlignmentUtils {
		
		static int MAX_INSERT_SIZE = 1000;
		
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
			boolean madeDir = dir.mkdir();
			
			// Use samtools to count alignments
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			String output = logDir + "/count_alignments_" + jobID + ".bsub";
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
			br.close();
			
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
		 * @throws IOException
		 * @throws InterruptedException
		 */
		public static void makeBowtie2Index(String fastaFile, String outBtIndexBase, String bowtie2BuildExecutable, String bsubOutputDir) throws IOException, InterruptedException {
			logger.info("");
			logger.info("Writing bowtie2 index for file " + fastaFile + " to files " + outBtIndexBase);
			File bsubDir = new File(bsubOutputDir);
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
			int prc = PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, output, "week", 4);
			PipelineUtils.waitForJobs(jobID, Runtime.getRuntime());
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
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie2(String bowtie2IndexBase, Map<String, String> options, String reads, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir) throws IOException, InterruptedException {
			return runBowtie2(bowtie2IndexBase, options, reads, null, outSamFile, outUnalignedFastq, bowtie2Executable, bsubOutDir, false);
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
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie2(String bowtie2IndexBase, Map<String, String> options, String read1Fastq, String read2Fastq, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir) throws IOException, InterruptedException {
			return runBowtie2(bowtie2IndexBase, options, read1Fastq, read2Fastq, outSamFile, outUnalignedFastq, bowtie2Executable, bsubOutDir, true);
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
		 * @throws IOException
		 * @throws InterruptedException
		 * @return The bsub job ID
		 */
		public static String runBowtie2(String bowtie2IndexBase, Map<String, String> options, String read1Fastq, String read2Fastq, String outSamFile, String outUnalignedFastq, String bowtie2Executable, String bsubOutDir, boolean readsPaired) throws IOException, InterruptedException {

			File bsubOutDirectory = new File(bsubOutDir);
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
			
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			String output = bsubOutDir + "/run_bowtie_" + jobID + ".bsub";
			logger.info("Writing bsub output to file " + output);
			int prc = PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, output, "week", 4);
			logger.info("Job ID is " + jobID);
			return jobID;
			
		}

		
	}
	
}


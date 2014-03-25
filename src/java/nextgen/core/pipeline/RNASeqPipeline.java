package nextgen.core.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
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
import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.LSFJob;
import nextgen.core.job.OGSJob;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.pipeline.util.AlignmentUtils;
import nextgen.core.pipeline.util.BamUtils;
import nextgen.core.pipeline.util.FastaUtils;
import nextgen.core.pipeline.util.FastqUtils;
import nextgen.core.pipeline.util.WigUtils;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.writers.PairedEndWriter;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import broad.core.error.ParseException;
import broad.core.math.EmpiricalDistribution;
import broad.core.parser.CommandLineParser;
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
public class RNASeqPipeline {

	private static Logger logger = Logger.getLogger(RNASeqPipeline.class.getName());
	private Session drmaaSession;


	/*
	 * For a file with paired data in the fq list file, there will be 4 columns
	 * Name	left.fq	right.fq	condition
	 * Name is unique
	 * Condition represents a string describing the sample/expt. Replicates have the same string in the condition
	 */
	private static int PAIRED_COLUMNS = 4;
	private static int UNPAIRED_COLUMNS = PAIRED_COLUMNS - 1;
	
	private ConfigFile configFile;
	
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
	 * Scheduler e.g. LSF or OGS
	 */
	Scheduler scheduler;
	
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
	 * @param drmaasession DRMAA session. There should only be one session at a time. Instantiate in the main method.
	 * @throws IOException
	 * @throws ParseException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public RNASeqPipeline(String inputListFile,String configFileName, Session drmaasession) throws IOException, ParseException, InterruptedException, DrmaaException {
		
		Globals.setHeadless(true);
		logger.info("Using Version R4.4");
		logger.debug("DEBUG ON");

		drmaaSession = drmaasession;
		configFile = getConfigFile(configFileName);
		fastqReadIdPairNumberDelimiter = configFile.getSingleValueField(sectionBasicOptions, optionFastqReadNumberDelimiter);
		scheduler = Scheduler.fromString(configFile.getSingleValueField(sectionScheduler, optionScheduler));
		
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
			queueName = configFile.getSingleValueField(sectionBasicOptions, optionQueueName);
			
			//Read the Fq list
			readFqList(inputListFile);

			// Split and trim barcodes
			if(configFile.hasOption(sectionCommands, optionSplitTrimBarcodes)){
				splitTrimBarcodes();
			}
			
			// Clip sequencing adapters
			if(configFile.hasOption(sectionCommands, optionTrimAdapters)) {
				if(!configFile.hasOption(sectionBasicOptions, optionFastxDirectory)) {
					throw new IllegalArgumentException("In order to clip adapters, must provide config file option " + optionFastxDirectory.getName());
				}
				if(!configFile.hasOption(sectionBasicOptions, optionRead1Adapter)) {
					throw new IllegalArgumentException("In order to clip adapters, must provide config file option " + optionRead1Adapter.getName());
				}
				if(!configFile.hasOption(sectionBasicOptions, optionFastqUtilsJar)) {
					throw new IllegalArgumentException("In order to clip adapters, must provide config file option " + optionFastqUtilsJar.getName());
				}
				String fastxDir = configFile.getSingleValueField(sectionBasicOptions, optionFastxDirectory);
				String adapter1 = configFile.getSingleValueField(sectionBasicOptions, optionRead1Adapter);
				String adapter2 = null;
				if(configFile.hasOption(sectionBasicOptions, optionRead2Adapter)) adapter2 = configFile.getSingleValueField(sectionBasicOptions, optionRead2Adapter);
				Map<String, ArrayList<String>> clippedFqs = FastqUtils.clipAdapters(fastxDir, currentLeftFqs, currentRightFqs, adapter1, adapter2, fastqReadIdPairNumberDelimiter, scheduler, drmaaSession, configFile.getSingleValueField(sectionBasicOptions, optionFastqUtilsJar));
				// Update current fastq files
				for(String sample : clippedFqs.keySet()) {
					currentLeftFqs.put(sample, clippedFqs.get(sample).get(0));
					if(clippedFqs.get(sample).size() > 1) {
						currentRightFqs.put(sample, clippedFqs.get(sample).get(1));
					}
				}
			}
			
			// Count reads
			// Quantify duplicates
			// Estimate library size
			if(configFile.hasOption(sectionCommands, optionComputeLibraryStats)){
				calculateLibraryStats();
			}
			
			// Characterize library composition by RNA class
			if(configFile.hasOption(sectionCommands, optionCountRnaClasses)){
				quantifyRNAClasses();
			}
			
			// Update current fastqs with reads not matching rRNA
			if(configFile.hasOption(sectionCommands, optionFilterRrna)){
				filterrRNA();
			}
			
			// Align to transcript sequences
			// Calculate median fragment size for each sequence
			if(configFile.hasOption(sectionCommands, optionAlignToTranscripts)) {
				alignToTranscripts();
			}
			
			// Align to genome
			// Generate bam and tdf files
			// Generate fragment size distributions
			if(configFile.hasOption(sectionCommands, optionAlign)){
				alignToGenome();
			}
			
		}
		
		if(configFile.hasOption(sectionCommands, optionRunDge)){
			
			String DGEInputFile= null;
			if(configFile.hasOption(sectionCommands, optionAlign)){
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
		return (configFile.hasOption(sectionCommands, optionSplitTrimBarcodes) || configFile.hasOption(sectionCommands, optionComputeLibraryStats) 
				|| configFile.hasOption(sectionCommands, optionCountRnaClasses) || configFile.hasOption(sectionCommands, optionFilterRrna) || configFile.hasOption(sectionCommands, optionAlign) || configFile.hasOption(sectionCommands, optionAlignToTranscripts));
	}
	
	/**
	 * 
	 * @param fileName
	 * @throws IOException 
	 */
	@SuppressWarnings("resource")
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
		@SuppressWarnings("unused")
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
	 * @throws DrmaaException 
	 */
	private void quantifyRNAClasses() throws IOException, InterruptedException, DrmaaException{
		
		logger.info("");
		logger.info("Quantifying RNA classes...");
		
		// Make directory for RNA classes
		String rnaClassDir = "rna_class_counts";
		File dir = new File(rnaClassDir);
		@SuppressWarnings("unused")
		boolean madeDir = dir.mkdir();
		if(!dir.exists()) {
			throw new IOException("Could not create directory " + rnaClassDir);
		}
		
		// Get options from config file
		if(!configFile.hasOption(sectionBasicOptions, optionRnaClassFastaFile)) {
			throw new IllegalArgumentException("In order to quantify RNA classes, config file must specify option " + optionRnaClassFastaFile.getName());
		}
		Map<String, String> classFiles = new TreeMap<String, String>();
		for(ConfigFileOptionValue value : configFile.getOptionValues(sectionBasicOptions, optionRnaClassFastaFile)) {
			classFiles.put(value.asString(1), value.asString(2));
		}
		
		if(!configFile.hasOption(sectionBasicOptions, optionGenomeBowtieIndex)) {
			throw new IllegalArgumentException("In order to quantify RNA classes, config file must specify option " + optionGenomeBowtieIndex.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionBowtie2Executable)) {
			throw new IllegalArgumentException("In order to quantify RNA classes, config file must specify option " + optionBowtie2Executable.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionBowtie2BuildExecutable)) {
			throw new IllegalArgumentException("In order to quantify RNA classes, config file must specify option " + optionBowtie2BuildExecutable.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionSamtoolsPath)) {
			throw new IllegalArgumentException("In order to quantify RNA classes, config file must specify option " + optionSamtoolsPath.getName());
		}
		
		
		String genomeBowtieIndex = configFile.getSingleValueField(sectionBasicOptions, optionGenomeBowtieIndex);
		String bowtie2Executable = configFile.getSingleValueField(sectionBasicOptions, optionBowtie2Executable);
		String bowtie2BuildExecutable = configFile.getSingleValueField(sectionBasicOptions, optionBowtie2BuildExecutable);
		String samtoolsExecutable = configFile.getSingleValueField(sectionBasicOptions, optionSamtoolsPath);
		
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
		LibraryCompositionByRnaClass lcrc = new LibraryCompositionByRnaClass(genomeBowtieIndex, classFiles, currentLeftFqs, currentRightFqs, logger);
		Map<String, String> bowtie2options = new TreeMap<String, String>();
		for(ConfigFileOptionValue value : configFile.getOptionValues(sectionBasicOptions, optionBowtie2Option)) {
			bowtie2options.put(value.asString(1), value.getLastFields(2));
		}
		Map<String, Integer> totalReadCounts = lcrc.getTotalReadCounts();
		Map<String, Map<String, Integer>> classCounts = lcrc.alignAndGetCounts(samtoolsExecutable, bowtie2Executable, bowtie2options, bowtie2BuildExecutable, rnaClassDir, scheduler, drmaaSession);
		
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
		logger.info("Delete sam and fastq files to save disk space. If rerunning pipeline remove " + optionCountRnaClasses.getName() + " option from config file.");

	}
	
	/**
	 * TASK 4: FILTER OUT RRNA
	 * Align reads to rRNA, retain reads that do not align, and replace current fastq files with files of non aligning reads
	 * @author prussell
	 * @throws InterruptedException 
	 * @throws IOException 
	 * @throws DrmaaException 
	 */
	private void filterrRNA() throws IOException, InterruptedException, DrmaaException{
		
		// Establish paths and software locations
		File outDirFile = new File(FILTER_RRNA_DIRECTORY);
		@SuppressWarnings("unused")
		boolean madeDir = outDirFile.mkdir();
		String outIndex = FILTER_RRNA_DIRECTORY + "/rRNA";
		
		if(!configFile.hasOption(sectionBasicOptions, optionRrnaSequences)) {
			throw new IllegalArgumentException("In order to filter ribosomal RNA, must provide config file option " + optionRrnaSequences.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionBowtie2BuildExecutable)) {
			throw new IllegalArgumentException("In order to filter ribosomal RNA, must provide config file option " + optionBowtie2BuildExecutable.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionBowtie2Executable)) {
			throw new IllegalArgumentException("In order to filter ribosomal RNA, must provide config file option " + optionBowtie2Executable.getName());
		}

		
		String rRnaFasta = configFile.getSingleValueField(sectionBasicOptions, optionRrnaSequences);
		String bowtieBuild = configFile.getSingleValueField(sectionBasicOptions, optionBowtie2BuildExecutable);
		String bowtie = configFile.getSingleValueField(sectionBasicOptions, optionBowtie2Executable);
		
		logger.info("");
		logger.info("Filtering ribosomal RNA by removing reads that map to sequences in " + rRnaFasta);
		
		// Make bowtie2 index for ribosomal RNA
		AlignmentUtils.makeBowtie2Index(rRnaFasta, outIndex, bowtieBuild, FILTER_RRNA_DIRECTORY, scheduler, drmaaSession);
		
		// Establish output file names
		Map<String, String> outRibosomalMap = new TreeMap<String, String>();
		Map<String, String> outFilteredUnpairedMap = new TreeMap<String, String>();
		Map<String, String> outFilteredPairedArgMap = new TreeMap<String, String>();
		Map<String, String> outFilteredPaired1Map = new TreeMap<String, String>();
		Map<String, String> outFilteredPaired2Map = new TreeMap<String, String>();
		Collection<Job> jobs = new ArrayList<Job>();
		
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
			Map<String, String> bowtie2options = new TreeMap<String, String>();
			for(ConfigFileOptionValue value : configFile.getOptionValues(sectionBasicOptions, optionBowtie2Option)) {
				bowtie2options.put(value.asString(1), value.getLastFields(2));
			}
			Job job = AlignmentUtils.runBowtie2(outIndex, bowtie2options, currentLeftFqs.get(sample), paired ? currentRightFqs.get(sample) : null, outRibosomal, paired ? outFilteredPairedArg : outFilteredUnpaired, bowtie, FILTER_RRNA_DIRECTORY, paired, scheduler, drmaaSession);
			jobs.add(job);
			
		}
		
		JobUtils.waitForAll(jobs);
		
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
	
	private void alignToTranscripts() throws IOException, InterruptedException, DrmaaException{
		
		
		// Establish paths and software locations
		File outDirFile = new File(ALIGN_TO_TRANSCRIPTS_DIRECTORY);
		@SuppressWarnings("unused")
		boolean madeDir = outDirFile.mkdir();
		String outIndex = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/transcripts";
		
		if(!configFile.hasOption(sectionBasicOptions, optionTranscriptFasta)) {
			throw new IllegalArgumentException("In order to align to transcripts, must provide config file option " + optionTranscriptFasta.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionBowtie2BuildExecutable)) {
			throw new IllegalArgumentException("In order to align to transcripts, must provide config file option " + optionBowtie2BuildExecutable.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionBowtie2Executable)) {
			throw new IllegalArgumentException("In order to align to transcripts, must provide config file option " + optionBowtie2Executable.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionSamtoolsPath)) {
			throw new IllegalArgumentException("In order to align to transcripts, must provide config file option " + optionSamtoolsPath.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionPicardDir)) {
			throw new IllegalArgumentException("In order to align to transcripts, must provide config file option " + optionPicardDir.getName());
		}

		
		String fasta = configFile.getSingleValueField(sectionBasicOptions, optionTranscriptFasta);
		String bowtieBuild = configFile.getSingleValueField(sectionBasicOptions, optionBowtie2BuildExecutable);
		String bowtie = configFile.getSingleValueField(sectionBasicOptions, optionBowtie2Executable);
		Map<String, String> bowtie2options = new TreeMap<String, String>();
		for(ConfigFileOptionValue value : configFile.getOptionValues(sectionBasicOptions, optionBowtie2Option)) {
			bowtie2options.put(value.asString(1), value.getLastFields(2));
		}
		String samtools = configFile.getSingleValueField(sectionBasicOptions, optionSamtoolsPath);
		String picardJarDir = configFile.getSingleValueField(sectionBasicOptions, optionPicardDir);
		
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
		AlignmentUtils.makeBowtie2Index(fasta, outIndex, bowtieBuild, ALIGN_TO_TRANSCRIPTS_DIRECTORY, scheduler, drmaaSession);
		
		// Establish output file names
		ArrayList<Job> jobs = new ArrayList<Job>();
		Map<String, String> samOutput = new TreeMap<String, String>();
		Map<String, String> unsortedBamOutput = new TreeMap<String, String>();
		Map<String, String> sortedBamOutput = new TreeMap<String, String>();
		Map<String, String> peBamOutput = new TreeMap<String, String>();
		for(String sample : sampleNames) {
			String sam = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/" + sample + "_transcript_mappings.sam";
			String unsortedBam = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/" + sample + "_transcript_mappings_unsorted.bam";
			String sortedBam = ALIGN_TO_TRANSCRIPTS_DIRECTORY + "/" + sample + "_transcript_mappings.bam";
			String peBam = sortedBam + PairedEndWriter.PAIRED_END_EXTENSION;
			samOutput.put(sample, sam);
			unsortedBamOutput.put(sample, unsortedBam);
			sortedBamOutput.put(sample, sortedBam);
			peBamOutput.put(sample, peBam);
		}
		
		// Align each sample
		for(String sample : sampleNames) {
			boolean paired = pairedData.get(sample).booleanValue();
			String sam = samOutput.get(sample);
			
			// Check if bam files exist
			File bamFile = new File(sortedBamOutput.get(sample));
			File pebamFile = new File(peBamOutput.get(sample));
			if(bamFile.exists()) {
				logger.warn("Bam file for sample " + sample + " already exists. Not rerunning alignment.");
				continue;	
			}
					
			// Delete old paired end bam file if one exists
			if(pebamFile.exists()) {
				boolean deleted = pebamFile.delete();
				if(!deleted) {
					logger.warn("Could not delete existing paired end bam file " + peBamOutput.get(sample) + ".");
				}
			}
			
			// Align to transcripts
			Job job = AlignmentUtils.runBowtie2(outIndex, bowtie2options, currentLeftFqs.get(sample), paired ? currentRightFqs.get(sample) : null, sam, null, bowtie, ALIGN_TO_TRANSCRIPTS_DIRECTORY, paired, scheduler, drmaaSession);
			jobs.add(job);
			
		}
		
		JobUtils.waitForAll(jobs);
		
		logger.info("");
		logger.info("Done aligning to transcripts. Converting sam to bam files.");
		Map<String, String> bsubDir = new TreeMap<String, String>();
		for(String sample : sampleNames) bsubDir.put(sample, ALIGN_TO_TRANSCRIPTS_DIRECTORY);
		BamUtils.samToBam(samOutput, unsortedBamOutput, sortedBamOutput, bsubDir, samtools, scheduler, drmaaSession);
		
		logger.info("");
		logger.info("Done converting sam to bam files. Delete sam files to save storage.");
		logger.info("Sorting bam files.");
		BamUtils.sortBamFiles(unsortedBamOutput, sortedBamOutput, bsubDir, sortedBamOutput, picardJarDir, scheduler, drmaaSession);
		// Delete unsorted bam files
		for(String sample : sampleNames) {
			File unsorted = new File(unsortedBamOutput.get(sample));
			boolean deleted = unsorted.delete();
			if(!deleted) logger.warn("Could not delete unsorted bam file " + unsortedBamOutput.get(sample) + ". Delete manually.");
		}
		
		logger.info("");
		logger.info("Done sorting bam files. Indexing sorted bam files.");
		BamUtils.indexBamFiles(sortedBamOutput, samtools, scheduler, drmaaSession);
				
		if(configFile.hasSection(sectionFragmentSizeDistribution)) {
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
					data.addFilter(new GenomicSpanFilter(configFile.getSingleValue(sectionFragmentSizeDistribution, optionFragmentSizeDistMaxSize).asInt(1)));
					for(String seqName : sequenceSizes.keySet()) {
						Annotation seq = gs.getReferenceAnnotation(seqName);
						try {
							double median = data.getMedianReadSize(seq, gs, configFile.getSingleValue(sectionFragmentSizeDistribution, optionFragmentSizeDistMaxSize).asInt(1), configFile.getSingleValue(sectionFragmentSizeDistribution, optionFragmentSizeDistNumBins).asInt(1));
							medianBySequence.put(seqName, Double.valueOf(median));
						} catch (IllegalArgumentException e) {
							continue;
						}
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
		}
		
		// Make paired end bam files
		/*logger.info("");
		logger.info("Making paired end bam files.");
		Collection<String> bamFilesToTranslate = new TreeSet<String>();
		for(String sample : sampleNames) {
			if(!pairedData.get(sample).booleanValue()) {
				continue;
			}
			File pebamFile = new File(peBamOutput.get(sample));
			if(pebamFile.exists()) {
				logger.warn("Paired end bam file for sample " + sample + " already exists. Not regenerating file.");
				continue;
			}
			bamFilesToTranslate.add(sortedBamOutput.get(sample));
		}
		BamUtils.createPairedEndBamFiles(bamFilesToTranslate, configFile.getSingleValueField(sectionBasicOptions, optionTranscriptionRead), ".", configFile.getSingleValueField(sectionBasicOptions, optionPairedEndWriterJar), scheduler, drmaaSession);
		*/
		
		// Make tdf of paired end bam files
		logger.info("");
		logger.info("Making tdf coverage files of paired end bam files.");
		logger.info("Indexing fasta file.");
		String indexedFasta = fasta + ".fai";
		File indexedFastaFile = new File(indexedFasta);
		if(indexedFastaFile.exists()) {
			logger.warn("Fasta index " + indexedFasta + " already exists. Not remaking fasta index.");
		} else {
			FastaUtils.indexFastaFile(fasta, samtools, ALIGN_TO_TRANSCRIPTS_DIRECTORY, scheduler, drmaaSession);
			logger.info("Done indexing fasta file.");
		}
		logger.info("Making tdf files.");
		
		if(!configFile.hasOption(sectionBasicOptions, optionIgvToolsExecutable)) {
			throw new IllegalArgumentException("In order to make tdf, must provide config file option " + optionIgvToolsExecutable.getName());
		}
		
		BamUtils.makeTdfs(peBamOutput, ALIGN_TO_TRANSCRIPTS_DIRECTORY, fasta, configFile.getSingleValueField(sectionBasicOptions, optionIgvToolsExecutable), scheduler, drmaaSession);
		
		// Make wig and bigwig files of fragment ends
		logger.info("");
		logger.info("Making wig and bigwig files of fragment end points.");
		writeWigFragmentEndsAndMidpoints(sortedBamOutput, ALIGN_TO_TRANSCRIPTS_DIRECTORY, fasta, null, configFile.getSingleValueField(sectionBasicOptions, optionWigWriterJar));
		WigUtils.writeWigPositionCount(sortedBamOutput, ALIGN_TO_TRANSCRIPTS_DIRECTORY, null, fasta, configFile.getSingleValueField(sectionBasicOptions, optionWigToBigWigExecutable), configFile.getSingleValueField(sectionBasicOptions, optionWigWriterJar), scheduler, drmaaSession);
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
	 * @throws DrmaaException 
	 */
	private void alignToGenome() throws IOException, InterruptedException, DrmaaException{
		
		logger.info("");
		logger.info("Aligning to genome...");

		logger.info("");
		logger.info("Aligning current fastq files to genome using Tophat...");

		if(!configFile.hasOption(sectionBasicOptions, optionGenomeBowtieIndex)) {
			throw new IllegalArgumentException("In order to align to genome, config file must specify option " + optionGenomeBowtieIndex.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionSamtoolsPath)) {
			throw new IllegalArgumentException("In order to align to genome, config file must specify option " + optionSamtoolsPath.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionTophatExecutable)) {
			throw new IllegalArgumentException("In order to align to genome, config file must specify option " + optionTophatExecutable.getName());
		}
		if(!configFile.hasOption(sectionBasicOptions, optionPicardDir)) {
			throw new IllegalArgumentException("In order to align to genome, config file must specify option " + optionPicardDir.getName());
		}

		
		// Establish index files and executables
		String genomeIndex = configFile.getSingleValueField(sectionBasicOptions, optionGenomeBowtieIndex);
		String samtools = configFile.getSingleValueField(sectionBasicOptions, optionSamtoolsPath);
		String tophat = configFile.getSingleValueField(sectionBasicOptions, optionTophatExecutable);
		String picardJarDir = configFile.getSingleValueField(sectionBasicOptions, optionPicardDir);
		File tophatDirFile = new File(TOPHAT_DIRECTORY);
		tophatDirFile.mkdir();
		
		// Get tophat options
		Map<String, String> tophatOptions = new TreeMap<String, String>();
		for(ConfigFileOptionValue value : configFile.getOptionValues(sectionBasicOptions, optionTophatOption)) {
			tophatOptions.put(value.asString(1), value.getLastFields(2));
		}
		if(tophatOptions.containsKey("-o") || tophatOptions.containsKey("--output-dir")) {
			logger.warn("Overriding tophat output directory provided in config file. Creating directories for each sample.");
			tophatOptions.remove("-o");
			tophatOptions.remove("--output-dir");
		}
		
		// Establish output directories
		Map<String, String> tophatDirsPerSample = new TreeMap<String, String>();
		Map<String, String> tophatSubdirsPerSample = new TreeMap<String, String>();
		Map<String, String> tophatBamUnsorted = new TreeMap<String, String>();
		Map<String, String> tophatBamSorted = new TreeMap<String, String>();  // where the bam file will be moved to
		for(String sample : sampleNames) {
			tophatDirsPerSample.put(sample, TOPHAT_DIRECTORY + "/" + TOPHAT_DIRECTORY + "_" + sample);
			tophatSubdirsPerSample.put(sample, TOPHAT_DIRECTORY + "_" + sample);
			tophatBamSorted.put(sample, TOPHAT_DIRECTORY + "/" + sample + ".bam");
			tophatBamUnsorted.put(sample, TOPHAT_DIRECTORY + "/" + sample + ".unsorted.bam");
		}
		
		// Run tophat
		Map<String, String> unsortedRun = AlignmentUtils.runTophat(tophat, samtools, sampleNames, leftFqs, rightFqs, tophatOptions, tophatDirsPerSample, tophatBamUnsorted, tophatBamSorted, genomeIndex, queueName, TOPHAT_DIRECTORY, scheduler, drmaaSession);
		// Sort bam files
		BamUtils.sortBamFiles(unsortedRun, tophatBamSorted, tophatDirsPerSample, tophatBamSorted, picardJarDir, scheduler, drmaaSession);
		// Delete unsorted bam files
		for(String sample : sampleNames) {
			File unsorted = new File(unsortedRun.get(sample));
			boolean deleted = unsorted.delete();
			if(!deleted) logger.warn("Could not delete unsorted bam file " + tophatBamUnsorted.get(sample) + ". Delete manually.");
		}
		currentBamFiles.putAll(tophatBamSorted);
		// Update current bam directory
		currentBamDir = TOPHAT_DIRECTORY;
		logger.info("Done running tophat.");
		
		// *** Novoalign steps ***
		// Only run if novoalign path was provided in config file
		if(configFile.hasOption(sectionBasicOptions, optionNovoalignExecutable)) {
		
			logger.info("");
			logger.info("Entering steps to align unmapped reads with Novoalign...");
			
			// Make sure genome novoindex was provided in config file
			if(!configFile.hasOption(sectionBasicOptions, optionGenomeNovoindex)) {
				throw new IllegalArgumentException("Novoalign index for genome is required. Specify in config file with option genome_novoindex.");
			}
			
			String novoalign = configFile.getSingleValueField(sectionBasicOptions, optionNovoalignExecutable);
			String novoIndex = configFile.getSingleValueField(sectionBasicOptions, optionGenomeNovoindex);
			
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
			@SuppressWarnings("unused")
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
				@SuppressWarnings("unused")
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
			replaceNovoalignSamHeaders(tophatBamSorted, novoSamOutput, novoSamOutputNoHeader, novoSamOutputReheadered, novoBamFinalPath, samtools);
			logger.info("Done replacing sam headers.");
			
			// Convert novoalign sam files to bam format
			logger.info("");
			logger.info("Converting novoalign sam files to bam format...");
			BamUtils.samToBam(novoSamOutput, novoBamOutput, novoBamFinalPath, novoDirsPerSample, samtools, scheduler, drmaaSession);
			logger.info("All samples done converting to bam format.");
			
			// Sort the bam files
			logger.info("");
			logger.info("Sorting novoalign bam files...");
			BamUtils.sortBamFiles(novoBamOutput, novoSortedBam, novoDirsPerSample, novoBamFinalPath, picardJarDir, scheduler, drmaaSession);
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
			mergeTophatNovoalign(tophatBamSorted, novoBamFinalPath, picardJarDir);
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
		if(configFile.hasOption(sectionBasicOptions, optionIgvToolsExecutable)) {
			logger.info("");
			logger.info("Making tdf files for bam files...");
			BamUtils.makeTdfs(currentBamFiles, currentBamDir, configFile.getSingleValueField(sectionBasicOptions, optionGenomeAssemblyName), configFile.getSingleValueField(sectionBasicOptions, optionIgvToolsExecutable), scheduler, drmaaSession);
			logger.info("All tdf files created.");
		}
		
		// Make fragment size distributions
		if(configFile.hasSection(sectionFragmentSizeDistribution)) {
			logger.info("");
			logger.info("Making fragment size distributions for bam files...");
			makeFragmentSizeDistributionCurrentBams();
			logger.info("All fragment size distributions created.");
		}
		
		// Compute global transcriptome space stats
		if(configFile.hasOption(sectionBasicOptions, optionAlignGlobalStatsJar)) {
			if(configFile.hasOption(sectionBasicOptions, optionTranscriptomeSpaceStatsBedFile) || configFile.hasOption(sectionBasicOptions, optionGenomicSpaceStatsSizeFile)) {
				writeAlignmentGlobalStats(configFile.getSingleValueField(sectionBasicOptions, optionAlignGlobalStatsJar), configFile.getSingleValueField(sectionBasicOptions, optionTranscriptomeSpaceStatsBedFile), configFile.getSingleValueField(sectionBasicOptions, optionGenomicSpaceStatsSizeFile));
			}
		}

		// Make wig and bigwig files of fragment ends and midpoints
		if(configFile.hasOption(sectionBasicOptions, optionWigToBigWigExecutable) && configFile.hasOption(sectionBasicOptions, optionBedFileForWig)) {
			logger.info("");
			logger.info("Making wig and bigwig files of fragment end points.");
			writeWigFragmentEndsAndMidpoints(currentBamFiles, currentBamDir, configFile.getSingleValueField(sectionBasicOptions, optionGenomeFasta), configFile.getSingleValueField(sectionBasicOptions, optionBedFileForWig), configFile.getSingleValueField(sectionBasicOptions, optionWigWriterJar));
			logger.info("");
			logger.info("Done writing wig files.\n");
		}
		
		// Make wig files of position fragment counts and normalized position counts
		if(configFile.hasOption(sectionBasicOptions, optionWigWriterJar) && configFile.hasOption(sectionBasicOptions, optionBedFileForWig)) {
			logger.info("");
			logger.info("Making wig files of position fragment counts and position counts normalized to transcript average coverage.");
			writeWigPositionCount(currentBamFiles, currentBamDir, configFile.getSingleValueField(sectionBasicOptions, optionBedFileForWig), configFile.getSingleValueField(sectionBasicOptions, optionGenomeFasta));
			logger.info("");
			logger.info("Done writing wig files.\n");
		}
		
		
		
	}
	

	/**
	 * Convert unmapped reads to fastq format if necessary and get fastq file names
	 * @param tophatDirsPerSample Directories containing tophat output, by sample name
	 * @param picardJarDir Directory containing Picard jar files
	 * @param version2 Whether tophat2 was used
	 * @return Map associating sample name with unmapped fastq files for read1 and read2. Read2 file is null if tophat2 was used.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private Map<String, String[]> unmappedToFastq(Map<String, String> tophatDirsPerSample, String picardJarDir, boolean version2) throws IOException, InterruptedException, DrmaaException {
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
	 * @throws DrmaaException 
	 */	
	private Map<String, String> unmappedFastqTophat2(Map<String, String> tophatDirsPerSample, String picardJarDir) throws IOException, InterruptedException, DrmaaException {
		ArrayList<Job> convertJobs = new ArrayList<Job>();
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
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, dir + "/sam_to_fastq_" + jobID + ".bsub", "hour", 16);
				job.submit();
				convertJobs.add(job);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "unmapped_fastq", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                convertJobs.add(ogsJob);
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		logger.info("Waiting for SamToFastq jobs to finish...");
		JobUtils.waitForAll(convertJobs);
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
	 * @throws DrmaaException 
	 */
	private void runNovoalignOnUnmappedReads(Map<String, String> tophatOptions, String novoindex, String novoalign, Map<String, String> novoDirsPerSample, Map<String, String> novoSamOutput, Map<String, String> novoBamFinalPath, Map<String, String[]> unmappedFastq, boolean tophat2) throws IOException, InterruptedException, DrmaaException {
		// Get novoalign options
		Map<String, String> novoalignOptions = new TreeMap<String, String>();
		for(ConfigFileOptionValue value : configFile.getOptionValues(sectionBasicOptions, optionNovoalignOption)) {
			novoalignOptions.put(value.asString(1), value.getLastFields(2));
		}
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
		ArrayList<Job> novoJobs = new ArrayList<Job>();
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
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				String bsubFile = outdir + "/novoalign_" + jobID + ".bsub";
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubFile, "week", 8);
				job.submit();
				novoJobs.add(job);
				novoBsubFiles.put(sample,bsubFile);
				break;
            case OGS:
            	throw new IllegalStateException("Not properly implemented for OGS. Change the code that expects bsub output files from LSF.");
                /*if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd);
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                ogsJob.submit();
                novoJobs.add(ogsJob);
                break;*/
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		// Wait for novoalign jobs to finish
		logger.info("Waiting for novoalign jobs to finish...");
		JobUtils.waitForAll(novoJobs);
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
	 * @throws DrmaaException 
	 */
	private void replaceNovoalignSamHeaders(Map<String, String> tophatBamFinalPath, Map<String, String> novoSamOutput, Map<String, String> novoSamOutputNoHeader, Map<String, String> novoSamOutputReheadered, Map<String, String> novoBamFinalPath, String samtools) throws IOException, InterruptedException, DrmaaException {
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
			switch(scheduler) {
			case LSF:
				logger.info("LSF job ID is " + getHeaderJobID + ".");
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), getHeaderJobID, getHeaderCmmd, NOVOALIGN_DIRECTORY + "/get_sam_header_" + getHeaderJobID + ".bsub", "hour", 1);
				lsfJob.submit();
				logger.info("Waiting for samtools view to finish...");
				lsfJob.waitFor();
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, getHeaderCmmd, "replace_sam_header", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
				logger.info("Waiting for samtools view to finish...");
				ogsJob.waitFor();
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
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
	 * Merge tophat and novoalign bam files
	 * @param tophatBamFinalPath Mapping of sample name to final tophat bam file
	 * @param novoBamFinalPath Mapping of sample name to final novoalign bam file
	 * @param picardJarDir Directory containing Picard jar files
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private void mergeTophatNovoalign(Map<String, String> tophatBamFinalPath, Map<String, String> novoBamFinalPath, String picardJarDir) throws IOException, InterruptedException, DrmaaException {
		File mergedDir = new File(MERGED_TOPHAT_NOVOALIGN_DIRECTORY);
		@SuppressWarnings("unused")
		boolean madeMergedDir = mergedDir.mkdir();
		ArrayList<Job> mergeJobs = new ArrayList<Job>();
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
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, MERGED_TOPHAT_NOVOALIGN_DIRECTORY + "/merge_bams_" + jobID + ".bsub", "hour", 1);
				job.submit();
				mergeJobs.add(job);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "merge_tophat_novoalign", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                mergeJobs.add(ogsJob);
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		// Wait for jobs to finish
		logger.info("Waiting for MergeSamFiles jobs to finish...");
		JobUtils.waitForAll(mergeJobs);
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
	 * @throws DrmaaException 
	 */
	private void indexCurrentBams(String samtools) throws IOException, InterruptedException, DrmaaException {
		BamUtils.indexBamFiles(currentBamFiles, samtools, scheduler, drmaaSession);
	}
	
	
	/**
	 * Write wig files of position count normalized by average coverage over gene
	 * @param bamFiles Bam files by sample name
	 * @param bamDir Bam directory
	 * @param geneBedFile Bed file of genes to use
	 * @param refFasta Reference fasta file
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private void writeWigPositionCount(Map<String, String> bamFiles, String bamDir, String geneBedFile, String refFasta) throws IOException, InterruptedException, DrmaaException {
		if(!configFile.hasOption(sectionBasicOptions, optionWigToBigWigExecutable)) {
			throw new IllegalArgumentException("In order to write wig file, must specify " + optionWigToBigWigExecutable.getName() + " in config file.");
		}
		if(!configFile.hasOption(sectionBasicOptions, optionWigToBigWigExecutable)) {
			throw new IllegalArgumentException("In order to write wig file, must specify " + optionWigWriterJar.getName() + " in config file.");
		}
		String wigToBigWig = configFile.getSingleValueField(sectionBasicOptions, optionWigToBigWigExecutable);
		String wigWriter = configFile.getSingleValueField(sectionBasicOptions, optionWigWriterJar);
		
		WigUtils.writeWigPositionCount(bamFiles, bamDir, geneBedFile, refFasta, wigToBigWig, wigWriter, scheduler, drmaaSession);
		
	}
	
	/**
	 * Write fragment end points and midpoints to wig and bigwig files
	 * @param bamFiles Bam files by sample name
	 * @param bamDir Directory containing bam files
	 * @param refFasta Fasta file of sequences these bam files were aligned against
	 * @param geneBedFile Bed file of genes to count reads in or null if using genomic space
	 * @param wigWriterJar WigWriter jar file
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private void writeWigFragmentEndsAndMidpoints(Map<String, String> bamFiles, String bamDir, String refFasta, String geneBedFile, String wigWriterJar) throws IOException, InterruptedException, DrmaaException {
		if(!configFile.hasOption(sectionBasicOptions, optionWigToBigWigExecutable)) {
			throw new IllegalArgumentException("In order to write wig file, must specify " + optionWigToBigWigExecutable.getName() + " in config file.");
		}
		String wigToBigWig = configFile.getSingleValueField(sectionBasicOptions, optionWigToBigWigExecutable);
		WigUtils.writeWigFragmentEndsAndMidpoints(bamFiles, pairedData, bamDir, refFasta, geneBedFile, wigWriterJar, wigToBigWig, scheduler, drmaaSession);
	}
	
	
	
	/**
	 * Reorder current bam files to match reference genome
	 * @param picardDir Directory containing Picard jar files
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private void reorderCurrentBams(String picardDir) throws IOException, InterruptedException, DrmaaException {
		ArrayList<Job> reorderJobs = new ArrayList<Job>();
		Map<String, String> reordered = new TreeMap<String, String>();
		for(String sample : sampleNames) {
			String bam = currentBamFiles.get(sample);
			reordered.put(sample, bam + ".reordered");
			String cmmd = "java -jar " + picardDir + "/ReorderSam.jar I=" + bam + " O=" + reordered.get(sample) + " R=" + configFile.getSingleValue(sectionBasicOptions, optionGenomeFasta);
			logger.info("Running picard command: " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, currentBamDir + "/reorder_bam_" + jobID + ".bsub", "week", 16);
				job.submit();
				reorderJobs.add(job);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "reorder_bam", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                reorderJobs.add(ogsJob);
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		logger.info("Waiting for picard jobs to finish...");
		JobUtils.waitForAll(reorderJobs);
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
	 * @throws DrmaaException 
	 */
	private void writeAlignmentGlobalStats(String alignmentGlobalStatsJar, String bedFile, String chrSizeFile) throws IOException, InterruptedException, DrmaaException {
		logger.info("Writing global stats for alignments...");
		ArrayList<Job> jobs = new ArrayList<Job>();
		if(bedFile != null) {
			Collection<Job> tJobs = BamUtils.writeTranscriptomeSpaceStats(currentBamFiles, bedFile, alignmentGlobalStatsJar, currentBamDir, scheduler, drmaaSession);
			jobs.addAll(tJobs);
		}
		if(chrSizeFile != null) {
			Collection<Job> gJobs = BamUtils.writeGenomicSpaceStats(currentBamFiles, chrSizeFile, alignmentGlobalStatsJar, currentBamDir, scheduler, drmaaSession);
			jobs.addAll(gJobs);
		}
		JobUtils.waitForAll(jobs);
		logger.info("Done writing all global alignment stats.");
	}
	



	/**
	 * Merge specified samples into new bam files
	 * Add merged samples to sample name list and current bam files
	 * @param picardJarDir Directory containing Picard executables
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private void mergeBamFiles(String picardJarDir) throws IOException, InterruptedException, DrmaaException {
		Map<String, Collection<String>> setsToMerge = new TreeMap<String, Collection<String>>();
		if(configFile.hasOption(sectionBasicOptions, optionSamplesToMerge)) {
			for(ConfigFileOptionValue value : configFile.getOptionValues(sectionBasicOptions, optionSamplesToMerge)) {
				String newName = value.asString(1);
				Collection<String> oldNames = new TreeSet<String>();
				for(int i = 2; i < value.getActualNumValues(); i++) {
					oldNames.add(value.asString(i));
				}
				setsToMerge.put(newName, oldNames);
			}
		}
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
		ArrayList<Job> jobs = new ArrayList<Job>();
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
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, currentBamDir + "/merge_bam_files_" + jobID + ".bsub", "week", 8);	
				job.submit();
				jobs.add(job);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "merge_bam_files", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                jobs.add(ogsJob);
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
			}
		}
		logger.info("Waiting for picard jobs to finish...");
		JobUtils.waitForAll(jobs);
		
		// Update current bam files and sample names
		currentBamFiles.putAll(mergedFiles);
		sampleNames.addAll(mergedFiles.keySet());
		pairedData.putAll(paired);
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
		
		String annotation = configFile.getSingleValueField(sectionFragmentSizeDistribution, optionFragmentSizeDistBedAnnotation);
		TranscriptomeSpace coord = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(annotation)));
		int maxFragmentSize = configFile.getSingleValue(sectionFragmentSizeDistribution, optionFragmentSizeDistMaxSize).asInt(1);
		int maxGenomicSpan = configFile.getSingleValue(sectionFragmentSizeDistribution, optionFragmentSizeDistMaxGenomicSpan).asInt(1);
		int numBins = configFile.getSingleValue(sectionFragmentSizeDistribution, optionFragmentSizeDistNumBins).asInt(1);
		boolean properPairsOnly = configFile.hasOption(sectionFragmentSizeDistribution, optionFragmentSizeDistProperPairsOnly);
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
		Collection<String> indGeneNames = ConfigFile.valuesAsStrings(configFile.getOptionValues(sectionFragmentSizeDistribution, optionFragmentSizeDistIndividualGene), 1);
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
	 * @throws DrmaaException 
	 */
	@SuppressWarnings("unused")
	private void collectPicardMetricsCurrentBams(String picardExecutableDir) throws IOException, InterruptedException, DrmaaException {
		String picardMetricsDir = currentBamDir + "/picard_metrics";
		logger.info("Writing all Picard metrics in directory " + picardMetricsDir + "....");
		File picardMetricsDirFile = new File(picardMetricsDir);
		boolean madePicardMetricsDir = picardMetricsDirFile.mkdir();
		ArrayList<Job> pmJobs = new ArrayList<Job>();
		// Inputs to picard metrics
		
		String refFlat = configFile.getSingleValueField(sectionBasicOptions, optionPicardRefFlat);
		String ribIntervals = configFile.getSingleValueField(sectionBasicOptions, optionPicardRibosomalIntervals);
		String genomeFasta = configFile.getSingleValueField(sectionBasicOptions, optionGenomeFasta);
		String strandSpecificity = configFile.getSingleValueField(sectionBasicOptions, optionPicardStrandSpecificity);
		
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
				switch(scheduler) {
				case LSF:
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					logger.info("LSF job ID is " + jobID + ".");
					// Submit job
					LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_alignment_summary_metrics_" + jobID + ".bsub", "hour", 4);
					job.submit();
					pmJobs.add(job);
					break;
                case OGS:
                    if(drmaaSession == null) {
                            throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                    }
                    OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "picard_metrics", null);
                    ogsJob.submit();
                    logger.info("OGS job ID is " + ogsJob.getID() + ".");
                    pmJobs.add(ogsJob);
                    break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
				}
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
				switch(scheduler) {
				case LSF:
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					logger.info("LSF job ID is " + jobID + ".");
					// Submit job
					LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_insert_size_metrics_" + jobID + ".bsub", "hour", 4);
					pmJobs.add(lsfJob);
					lsfJob.submit();
					break;
                case OGS:
                    if(drmaaSession == null) {
                            throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                    }
                    OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "picard_metrics", null);
                    ogsJob.submit();
                    logger.info("OGS job ID is " + ogsJob.getID() + ".");
                    pmJobs.add(ogsJob);
                    break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
				}
			}
			
			// Run RNA-seq metrics
			if(rsFile.exists() && rscFile.exists()) {
				logger.warn("RNA-seq metrics files " + rsMetrics + " and " + rsChart + " already exist. Not rerunning RNA-seq metrics.");
			} else {
				// Use Picard program CollectRnaSeqMetrics
				String cmmd = "java -jar ";
				// Use user-provided CollectRnaSeqMetrics executable if provided
				if(configFile.hasOption(sectionBasicOptions, optionPicardCollectRnaSeqMetricsJar)) cmmd += configFile.getSingleValue(sectionBasicOptions, optionPicardCollectRnaSeqMetricsJar);
				else cmmd += picardExecutableDir + "/CollectRnaSeqMetrics.jar";
				cmmd += " INPUT=" + currentBamFiles.get(sample);
				cmmd += " OUTPUT=" + rsMetrics;
				cmmd += " REFERENCE_SEQUENCE=" + genomeFasta;
				cmmd += " CHART_OUTPUT=" + rsChart;
				cmmd += " REF_FLAT=" + refFlat;
				if(ribIntervals != null) cmmd += " RIBOSOMAL_INTERVALS=" + ribIntervals;
				cmmd += " STRAND_SPECIFICITY=" + strandSpecificity;
				logger.info("Running Picard command: " + cmmd);
				switch(scheduler) {
				case LSF:
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					logger.info("LSF job ID is " + jobID + ".");
					// Submit job
					LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, picardMetricsDir + "/picard_rnaseq_metrics_" + jobID + ".bsub", "hour", 4);
					job.submit();
					pmJobs.add(job);
					break;
                case OGS:
                    if(drmaaSession == null) {
                            throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                    }
                    OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "picard_metrics", null);
                    ogsJob.submit();
                    logger.info("OGS job ID is " + ogsJob.getID() + ".");
                    pmJobs.add(ogsJob);
                    break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
				}
			}
		}
		
		// Wait for jobs to finish
		logger.info("Waiting for Picard metrics jobs to finish...");
		try {
			JobUtils.waitForAll(pmJobs);
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

	public static void main (String [] args) throws IOException, ParseException, InterruptedException, DrmaaException{
		
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("\n*** Configurable pipeline for RNA-seq read processing and analysis ***");
		p.addStringArg("-r", "File containing list of fastq files. \n\t\tLine format: \n\t\t<sample_name> <left.fq> <right.fq> <condition> \n\t\tOR \n\t\t<sample_name> <unpaired.fq> <condition>", true);
		p.addStringArg("-c", "Config file", true);
		p.parse(args);
		String fastqList = p.getStringArg("-r");
		String configFile = p.getStringArg("-c");
		
		ConfigFile c = getConfigFile(configFile);
		Scheduler s = Scheduler.fromString(c.getSingleValueField(sectionScheduler, optionScheduler));
		Session drmaaSession = s.equals(Scheduler.OGS) ? OGSUtils.getDrmaaSession() : null;

		
		RNASeqPipeline PA = new RNASeqPipeline(fastqList,configFile, drmaaSession);
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
	
	
	
	/*
	 * Config file options and sections 
	 */
	
	private static ConfigFileSection sectionScheduler = new ConfigFileSection("scheduler", true);
	private static ConfigFileSection sectionCommands = new ConfigFileSection("commands", true);
	private static ConfigFileSection sectionSpecies = new ConfigFileSection("species", false);
	private static ConfigFileSection sectionBasicOptions = new ConfigFileSection("basic_options", true);
	private static ConfigFileSection sectionFragmentSizeDistribution = new ConfigFileSection("fragment_size_distribution", false);
	
	private static ConfigFileOption optionScheduler = new ConfigFileOption("scheduler", 2, false, false, true);
	private static ConfigFileOption optionTranscriptionRead = new ConfigFileOption("transcription_read", 2, false, false, true);
	//private static ConfigFileOption optionPairedEndWriterJar = new ConfigFileOption("paired_end_writer_jar", 2, false, false, true);
	private static ConfigFileOption optionSplitTrimBarcodes = new ConfigFileOption("SPLIT_TRIM_BARCODES", 1, false, false, false);
	private static ConfigFileOption optionTrimAdapters = new ConfigFileOption("TRIM_ADAPTERS", 1, false, false, false);
	private static ConfigFileOption optionComputeLibraryStats = new ConfigFileOption("LIBRARY_STATS", 1, false, false, false);
	private static ConfigFileOption optionCountRnaClasses = new ConfigFileOption("RNA_CLASSES", 1, false, false, false);
	private static ConfigFileOption optionFilterRrna = new ConfigFileOption("FILTER_RRNA", 1, false, false, false);
	private static ConfigFileOption optionAlign = new ConfigFileOption("ALIGN", 1, false, false, false);
	private static ConfigFileOption optionAlignToTranscripts = new ConfigFileOption("ALIGN_TO_TRANSCRIPTS", 1, false, false, false);
	private static ConfigFileOption optionRunDge = new ConfigFileOption("RUN_DGE", 1, false, false, false);
	private static ConfigFileOption optionGenomeFasta = new ConfigFileOption("genome_fasta", 2, false, false, true);
	private static ConfigFileOption optionQueueName = new ConfigFileOption("queue_name", 2, false, false, false, "week");
	private static ConfigFileOption optionGenomeBowtieIndex = new ConfigFileOption("genome_bowtie", 2, false, false, true);
	private static ConfigFileOption optionGenomeAssemblyName = new ConfigFileOption("genome_assembly", 2, false, false, true);
	private static ConfigFileOption optionGtfAnnotation = new ConfigFileOption("annotation", 2, false, false, false);
	private static ConfigFileOption optionRnaClassFastaFile = new ConfigFileOption("rna_classes", 3, false, true, false);
	private static ConfigFileOption optionRrnaSequences = new ConfigFileOption("filter_rrna", 2, false, false, false);
	private static ConfigFileOption optionTranscriptFasta = new ConfigFileOption("align_to_transcripts", 2, false, false, false);
	private static ConfigFileOption optionSamplesToMerge = new ConfigFileOption("merge_samples", 10, true, true, false);
	private static ConfigFileOption optionTophatExecutable = new ConfigFileOption("tophat_path", 2, false, false, true);
	private static ConfigFileOption optionBedFileForWig = new ConfigFileOption("bed_file_for_wig", 2, false, false, false);
	private static ConfigFileOption optionWigToBigWigExecutable = new ConfigFileOption("wig_to_bigwig_path", 2, false, false, false);
	private static ConfigFileOption optionWigWriterJar = new ConfigFileOption("wig_writer_jar", 2, false, false, true);
	private static ConfigFileOption optionChrSizeFile = new ConfigFileOption("chr_size_file", 2, false, false, true);
	private static ConfigFileOption optionAlignGlobalStatsJar = new ConfigFileOption("alignment_global_stats_jar_file", 2, false, false, false);
	private static ConfigFileOption optionFastqUtilsJar = new ConfigFileOption("fastq_utils_jar_file", 2, false, false, false);
	private static ConfigFileOption optionTophatOption = new ConfigFileOption("tophat_options", 3, true, true, false);
	private static ConfigFileOption optionBowtie2Option = new ConfigFileOption("bowtie2_options", 3, true, true, false);
	private static ConfigFileOption optionNovoalignOption = new ConfigFileOption("novoalign_options", 3, true, true, false);
	private static ConfigFileOption optionFragmentSizeDistOption = new ConfigFileOption("fragment_size_dist_options", 3, true, true, false);
	private static ConfigFileOption optionBowtie2Executable = new ConfigFileOption("bowtie2_path", 2, false, false, true);
	private static ConfigFileOption optionBowtie2BuildExecutable = new ConfigFileOption("bowtie2_build_path", 2, false, false, true);
	private static ConfigFileOption optionSamtoolsPath = new ConfigFileOption("samtools_path", 2, false, false, true);
	private static ConfigFileOption optionIgvToolsExecutable = new ConfigFileOption("igvtools_path", 2, false, false, true);
	private static ConfigFileOption optionFastqReadNumberDelimiter = new ConfigFileOption("fastq_read_number_special_delimiter", 2, false, false, false, null);
	private static ConfigFileOption optionPicardDir = new ConfigFileOption("picard_directory", 2, false, false, true);
	private static ConfigFileOption optionFastxDirectory = new ConfigFileOption("fastx_directory", 2, false, false, false);
	private static ConfigFileOption optionRead1Adapter = new ConfigFileOption("sequencing_adapter_read1", 2, false, false, false);
	private static ConfigFileOption optionRead2Adapter = new ConfigFileOption("sequencing_adapter_read2", 2, false, false, false);
	private static ConfigFileOption optionNovoalignExecutable = new ConfigFileOption("novoalign_path", 2, false, false, false);
	private static ConfigFileOption optionGenomeNovoindex = new ConfigFileOption("genome_novoindex", 2, false, false, false);
	private static ConfigFileOption optionPicardRefFlat = new ConfigFileOption("picard_ref_flat", 2, false, false, false);
	private static ConfigFileOption optionPicardRibosomalIntervals = new ConfigFileOption("picard_ribosomal_intervals", 2, false, false, false);
	private static ConfigFileOption optionPicardStrandSpecificity = new ConfigFileOption("picard_strand_specificity", 2, false, false, false);
	private static ConfigFileOption optionPicardCollectRnaSeqMetricsJar = new ConfigFileOption("picard_collect_rnaseq_metrics", 2, false, false, false);
	private static ConfigFileOption optionTranscriptomeSpaceStatsBedFile = new ConfigFileOption("bed_transcriptome_space_stats", 2, false, false, false);
	private static ConfigFileOption optionGenomicSpaceStatsSizeFile = new ConfigFileOption("size_file_genomic_space_stats", 2, false, false, false);
	private static ConfigFileOption optionSpeciesName = new ConfigFileOption("species", 2, false, false, false, "mouse");
	
	private static ConfigFileOption optionFragmentSizeDistBedAnnotation = new ConfigFileOption("fragment_size_dist_bed_annotation", 2, false, false, true);
	private static ConfigFileOption optionFragmentSizeDistMaxSize = new ConfigFileOption("fragment_size_dist_max_fragment_length", 2, false, false, false, "1000");
	private static ConfigFileOption optionFragmentSizeDistMaxGenomicSpan = new ConfigFileOption("fragment_size_dist_max_genomic_span", 2, false, false, false, "300000");
	private static ConfigFileOption optionFragmentSizeDistNumBins = new ConfigFileOption("fragment_size_dist_num_bins", 2, false, false, false, "100");
	private static ConfigFileOption optionFragmentSizeDistProperPairsOnly = new ConfigFileOption("fragment_size_dist_proper_pairs_only", 1, false, false, false);
	private static ConfigFileOption optionFragmentSizeDistIndividualGene = new ConfigFileOption("fragment_size_dist_individual_gene", 2, false, true, false);

	
	private static ConfigFile getConfigFile(String fileName) throws IOException {
		
		sectionScheduler.addAllowableOption(optionScheduler);
		
		sectionCommands.addAllowableOption(optionAlign);
		sectionCommands.addAllowableOption(optionAlignToTranscripts);
		sectionCommands.addAllowableOption(optionComputeLibraryStats);
		sectionCommands.addAllowableOption(optionCountRnaClasses);
		sectionCommands.addAllowableOption(optionFilterRrna);
		sectionCommands.addAllowableOption(optionRunDge);
		sectionCommands.addAllowableOption(optionSplitTrimBarcodes);
		sectionCommands.addAllowableOption(optionTrimAdapters);
		
		sectionSpecies.addAllowableOption(optionSpeciesName);
		
		sectionBasicOptions.addAllowableOption(optionGenomeFasta);
		sectionBasicOptions.addAllowableOption(optionQueueName);
		sectionBasicOptions.addAllowableOption(optionGenomeBowtieIndex);
		sectionBasicOptions.addAllowableOption(optionGenomeAssemblyName);
		sectionBasicOptions.addAllowableOption(optionGtfAnnotation);
		sectionBasicOptions.addAllowableOption(optionRnaClassFastaFile);
		sectionBasicOptions.addAllowableOption(optionRrnaSequences);
		sectionBasicOptions.addAllowableOption(optionTranscriptFasta);
		sectionBasicOptions.addAllowableOption(optionSamplesToMerge);
		sectionBasicOptions.addAllowableOption(optionTophatExecutable);
		sectionBasicOptions.addAllowableOption(optionBedFileForWig);
		sectionBasicOptions.addAllowableOption(optionWigToBigWigExecutable);
		sectionBasicOptions.addAllowableOption(optionWigWriterJar);
		sectionBasicOptions.addAllowableOption(optionChrSizeFile);
		sectionBasicOptions.addAllowableOption(optionAlignGlobalStatsJar);
		sectionBasicOptions.addAllowableOption(optionFastqUtilsJar);
		sectionBasicOptions.addAllowableOption(optionTophatOption);
		sectionBasicOptions.addAllowableOption(optionBowtie2Option);
		sectionBasicOptions.addAllowableOption(optionNovoalignOption);
		sectionBasicOptions.addAllowableOption(optionFragmentSizeDistOption);
		sectionBasicOptions.addAllowableOption(optionBowtie2Executable);
		sectionBasicOptions.addAllowableOption(optionBowtie2BuildExecutable);
		sectionBasicOptions.addAllowableOption(optionSamtoolsPath);
		sectionBasicOptions.addAllowableOption(optionIgvToolsExecutable);
		sectionBasicOptions.addAllowableOption(optionFastqReadNumberDelimiter);
		sectionBasicOptions.addAllowableOption(optionPicardDir);
		sectionBasicOptions.addAllowableOption(optionFastxDirectory);
		sectionBasicOptions.addAllowableOption(optionRead1Adapter);
		sectionBasicOptions.addAllowableOption(optionRead2Adapter);
		sectionBasicOptions.addAllowableOption(optionNovoalignExecutable);
		sectionBasicOptions.addAllowableOption(optionGenomeNovoindex);
		sectionBasicOptions.addAllowableOption(optionPicardRefFlat);
		sectionBasicOptions.addAllowableOption(optionPicardRibosomalIntervals);
		sectionBasicOptions.addAllowableOption(optionPicardStrandSpecificity);
		sectionBasicOptions.addAllowableOption(optionPicardCollectRnaSeqMetricsJar);
		sectionBasicOptions.addAllowableOption(optionTranscriptomeSpaceStatsBedFile);
		sectionBasicOptions.addAllowableOption(optionGenomicSpaceStatsSizeFile);
		sectionBasicOptions.addAllowableOption(optionTranscriptionRead);
		//sectionBasicOptions.addAllowableOption(optionPairedEndWriterJar);
		
		sectionFragmentSizeDistribution.addAllowableOption(optionFragmentSizeDistBedAnnotation);
		sectionFragmentSizeDistribution.addAllowableOption(optionFragmentSizeDistMaxSize);
		sectionFragmentSizeDistribution.addAllowableOption(optionFragmentSizeDistMaxGenomicSpan);
		sectionFragmentSizeDistribution.addAllowableOption(optionFragmentSizeDistNumBins);
		sectionFragmentSizeDistribution.addAllowableOption(optionFragmentSizeDistProperPairsOnly);
		sectionFragmentSizeDistribution.addAllowableOption(optionFragmentSizeDistIndividualGene);
		
		Collection<ConfigFileSection> sections = new ArrayList<ConfigFileSection>();
		sections.add(sectionScheduler);
		sections.add(sectionCommands);
		sections.add(sectionSpecies);
		sections.add(sectionBasicOptions);
		sections.add(sectionFragmentSizeDistribution);
		
		ConfigFile cf = new ConfigFile(sections, fileName);
		return cf;

	}

	
}


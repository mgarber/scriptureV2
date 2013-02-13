package broad.core.parser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

/**
 * This class will parse the config file provided to it.
 * Helper class to broad.core.util.PipelineAutomator
 * @author skadri
 *
 */
public class ConfigFileParser {

	static Logger logger = Logger.getLogger(ConfigFileParser.class.getName());

	private List<String> commands;
	private String species;
	public static int START_OF_STRING = 1;
	private static String DEFAULT_SPECIES = "mouse";
	private static String DEFAULT_QUEUE_NAME = "week";
	public BasicOptions basicOptions  = new BasicOptions();
	public FragmentSizeDistributionOptions fragmentSizeOptions = new FragmentSizeDistributionOptions();
	public DGEOptions dgeOptions;
	
	// Basic options
	/**
	 * Config file option for path to genome fasta file
	 */
	public static String OPTION_GENOME_FASTA = "genome_fasta";
	public static String OPTION_QUEUE_NAME = "queue_name";
	
	/**
	 * Config file option for prefix of genome bowtie index files
	 */
	public static String OPTION_GENOME_BOWTIE = "genome_bowtie";
	
	/**
	 * Config file option for genome assembly name e.g. mm9
	 */
	public static String OPTION_GENOME_ASSEMBLY = "genome_assembly";
	
	/**
	 * Config file option for genome annotation in gtf format for use by tophat
	 */
	public static String OPTION_ANNOTATION = "annotation";
	
	/**
	 * Config file option to provide an RNA class fasta file e.g. rna_classes mRNA mRNA.fa
	 */
	public static String OPTION_RNA_CLASS = "rna_classes";
	
	/**
	 * Config file option to provide a fasta file of rRNA sequences for depletion
	 */
	public static String OPTION_FILTER_RRNA = "filter_rrna";
	
	/**
	 * Config file option to align to a fasta file of transcript sequences
	 */
	public static String OPTION_ALIGN_TO_TRANSCRIPTS = "align_to_transcripts";
	
	/**
	 * Merge a list of samples. Keep separate files too.
	 */
	public static String OPTION_MERGE_SAMPLES = "merge_samples";
	
	/**
	 * Config file option for tophat executable file
	 */
	public static String OPTION_TOPHAT_PATH = "tophat_path";
	
	/**
	 * Config file option to write fragment endpoints to wig and bigwig file
	 * Writes counts for positions overlapping genes in this bed file
	 */
	public static String OPTION_BED_FILE_FOR_FRAGMENT_END_WIG = "bed_file_fragment_ends";
	
	/**
	 * Path to UCSC WigToBigWig executable
	 */
	public static String OPTION_WIG_TO_BIGWIG_PATH = "wig_to_bigwig_path";
	
	/**
	 * Chromosome size file
	 */
	public static String OPTION_CHR_SIZE_FILE = "chr_size_file";
	
	/**
	 * Config file option to provide options to tophat e.g. tophat_options --max-multihits 1
	 */
	public static String OPTION_TOPHAT_OPTIONS = "tophat_options";
	
	/**
	 * Config file option to provide options to bowtie2
	 */
	public static String OPTION_BOWTIE2_OPTIONS = "bowtie2_options";

	
	/**
	 * Config file option to provide options to novoalign e.g. novoalign_options -R 5
	 */
	public static String OPTION_NOVOALIGN_OPTIONS = "novoalign_options";
	
	/**
	 * Config file option to provide parameters to fragment size distribution calculation
	 */
	public static String OPTION_FRAGMENT_SIZE_DIST_OPTIONS = "fragment_size_dist_options";
	
	/**
	 * Config file option for bowtie2 executable file
	 */
	public static String OPTION_BOWTIE2_PATH = "bowtie2_path";
	
	/**
	 * Config file option for bowtie2-build executable file
	 */
	public static String OPTION_BOWTIE2_BUILD_PATH = "bowtie2_build_path";
	
	/**
	 * Config file option for samtools executable
	 */
	public static String OPTION_SAMTOOLS_PATH = "samtools_path";
	
	/**
	 * Config file option for igvtools executable
	 */
	public static String OPTION_IGVTOOLS_PATH = "igvtools_path";
	
	/*
	 * The delimiter that separates the main part of read ID from the pair number in fastq files
	 */
	public static String OPTION_FASTQ_READ_PAIR_NUMBER_DELIMITER = "fastq_read_number_delimiter";
	
	/*
	 * Default for the delimiter that separates the main part of read ID from the pair number in fastq files
	 */
	public static String DEFAULT_FASTQ_READ_PAIR_NUMBER_DELIMITER = "\\s++";

	/**
	 * Config file option for directory containing picard executables
	 */
	public static String OPTION_PICARD_DIRECTORY = "picard_directory";
	
	/**
	 * Directory containing fastx binaries
	 */
	public static String OPTION_FASTX_DIRECTORY = "fastx_directory";
	
	/**
	 * Sequencing adapter sequence for read 1
	 */
	public static String OPTION_SEQUENCING_ADAPTER_READ_1 = "sequencing_adapter_read1";
	
	/**
	 * Sequencing adapter sequence for read 2
	 */
	public static String OPTION_SEQUENCING_ADAPTER_READ_2 = "sequencing_adapter_read2";
	
	/**
	 * Config file option for novoalign executable
	 */
	public static String OPTION_NOVOALIGN_PATH = "novoalign_path";
	
	/**
	 * Config file option for novoindex of genome
	 */
	public static String OPTION_GENOME_NOVOINDEX = "genome_novoindex";
	
	/**
	 * Config file option to provide ref_flat annotation to Picard
	 */
	public static String OPTION_PICARD_REF_FLAT = "picard_ref_flat";
	
	/**
	 * Config file option to provide ribosomal intervals to Picard
	 */
	public static String OPTION_PICARD_RIBOSOMAL_INTERVALS = "picard_ribosomal_intervals";
	
	/**
	 * Config file option to provide strand specificity argument to Picard
	 */
	public static String OPTION_PICARD_STRAND_SPECIFICITY = "picard_strand_specificity";
	
	/**
	 * Config file option to provide an alternative jar file for Picard program CollectRnaSeqMetrics
	 */
	public static String OPTION_PICARD_COLLECT_RNASEQ_METRICS = "picard_collect_rnaseq_metrics";

	
	/**
	 * Constructor will read the input config file and store the commands and parameters in the respective 
	 * @param configFileName
	 * @throws IOException 
	 */
	public ConfigFileParser(String configFileName) throws IOException {
				
		Globals.setHeadless(true);
		System.out.println("Using Version R4.4");
		logger.debug("DEBUG ON");
		
		Scanner reader = new Scanner(new File(configFileName));
		commands = new ArrayList<String>();
		String nextLine;
		while (reader.hasNextLine()) {
			
			nextLine = reader.nextLine();
			if(startOfASection(nextLine)){
				/*
				 * READ ALL COMMANDS
				 */
				if(isCommandSection(nextLine)){
					logger.info("Entering the commands section");
					
					//Read the next line
					nextLine = reader.nextLine();
					
					//While the next section doesnt start or end of file
					while(!startOfNextSection(nextLine,"commands") && (reader.hasNextLine())){						
						addCommand(nextLine);
						nextLine = reader.nextLine();
					}
				}
				/*
				 * READ SPECIES
				 */
				if(isSpeciesSection(nextLine)){
					logger.info("Entering the species section");
					//Read the next line
					nextLine = reader.nextLine();

					//While the next section doesnt start or end of file
					while(!startOfNextSection(nextLine,"species") && (reader.hasNextLine())){
						addSpecies(nextLine);
						nextLine = reader.nextLine();
					}
				}
				/*
				 * READ THE BASIC OPTIONS FOR THE SPECIES IN QUESTION
				 */
				if(isBasicOptionsForSpecies(nextLine)){
					logger.info("Entering the basic options section for species : "+species);

					StringParser parser = new StringParser();
					/*
					 * Add all the options for the species
					 */
					//While the next section doesnt start or end of file
					while(!startOfNextSection(nextLine,"basic") && (reader.hasNextLine())){
						//Read the next line
						nextLine = reader.nextLine();
						
						/*
						 * BASIC OPTIONS
						 */
						if(!isComment(nextLine)){
							parser.parse(nextLine);
							if(parser.getFieldCount()>1) {
								if(parser.getFieldCount()==2 && this.isTophatPath(nextLine)) {
									basicOptions.setTophatPath(parser.asString(1));
									continue;
								}
								if(parser.getFieldCount()==2 && this.isBowtie2Path(nextLine)) {
									basicOptions.bt2Options.setBowtie2Path(parser.asString(1));
									basicOptions.addOption(parser.asString(0), parser.asString(1));
									continue;
								}
								if(parser.getFieldCount()==2 && this.isNovoalignPath(nextLine)) {
									basicOptions.setNovoalignPath(parser.asString(1));
									continue;
								}
								if(parser.getFieldCount()==2 && this.isQueueName(nextLine)) {
									basicOptions.setQueueName(parser.asString(1));
									continue;
								}
								if(this.isTophatOption(nextLine)){
									String value = "";
									if(parser.getFieldCount() > 2) {
										value += parser.asString(2);
									}
									basicOptions.addTophatOption(parser.asString(1), value);
									continue;
								}
								if(this.isBowtie2Option(nextLine)){
									String value = "";
									if(parser.getFieldCount() > 2) {
										value += parser.asString(2);
									}
									basicOptions.addBowtie2Option(parser.asString(1), value);
									continue;
								}
								if(this.isNovoalignOption(nextLine)) {
									if(parser.getFieldCount() == 2) basicOptions.addNovoalignOption(parser.asString(1), "");
									if(parser.getFieldCount() == 3) basicOptions.addNovoalignOption(parser.asString(1), parser.asString(2));
									if(parser.getFieldCount() == 4) basicOptions.addNovoalignOption(parser.asString(1), parser.asString(2) + " " + parser.asString(3));
									continue;
								}
								if(parser.getFieldCount() >= 3 && isRNAClass(nextLine)){
									basicOptions.addRNAClass(parser.asString(1), parser.asString(2));
									continue;
								}
								if(this.isFragmentSizeDistOption(nextLine)) {
									if(parser.getFieldCount() == 3 && parser.asString(1).equals(FragmentSizeDistributionOptions.OPTION_MAX_GENOMIC_SPAN)) {fragmentSizeOptions.setMaxGenomicSpan(parser.asInt(2)); continue;}
									if(parser.getFieldCount() == 2 && parser.asString(1).equals(FragmentSizeDistributionOptions.OPTION_INCLUDE_IMPROPER_PAIRS)) {fragmentSizeOptions.setProperPairsOnly(false); continue;}
									if(parser.getFieldCount() == 3 && parser.asString(1).equals(FragmentSizeDistributionOptions.OPTION_BED_FILE)) {fragmentSizeOptions.setBedFile(parser.asString(2)); continue;}
									if(parser.getFieldCount() == 3 && parser.asString(1).equals(FragmentSizeDistributionOptions.OPTION_MAX_FRAGMENT_SIZE)) {fragmentSizeOptions.setMaxFragmentSize(parser.asInt(2)); continue;}
									if(parser.getFieldCount() == 3 && parser.asString(1).equals(FragmentSizeDistributionOptions.OPTION_NUM_BINS)) {fragmentSizeOptions.setNumBins(parser.asInt(2)); continue;}
									if(parser.getFieldCount() == 3 && parser.asString(1).equals(FragmentSizeDistributionOptions.OPTION_IND_GENE)) {fragmentSizeOptions.addIndividualGene(parser.asString(2)); continue;}
									throw new IllegalArgumentException("Not a valid option for fragment size distribution: " + parser.asString(1));
								}
								if(parser.getFieldCount()==2) {
									basicOptions.addOption(parser.asString(0), parser.asString(1));
									continue;
								}
							}
						}
					}
				}
				
				/*
				 * READ THE OPTIONS FOR DGE FOR THE SPECIES IN QUESTION
				 */
				if(isDGEOptionsForSpecies(nextLine)){
					logger.info("Entering the DGE options section for species : "+species);

					dgeOptions = new DGEOptions();
					StringParser parser = new StringParser();
					/*
					 * Add all the options for the species
					 */
					logger.info(nextLine);
					//While the next section doesnt start or end of file
					while(!startOfNextSection(nextLine,"dge") && (reader.hasNextLine()) ){
						//Read the next line
						nextLine = reader.nextLine();
						/*
						 * DGE OPTIONS
						 */
						if(!isComment(nextLine)){
							parser.parse(nextLine);
							if(parser.getFieldCount()>1){
								if(parser.getFieldCount()==2){
									//DGE type
									if(parser.asString(0).contains("run3DGE")){
										if(Boolean.parseBoolean(parser.asString(1)))
											dgeOptions.setdgeType("3DGE");
										else
											dgeOptions.setdgeType("5DGE");
									}
									//RUN TYPE
									else if(parser.asString(0).equalsIgnoreCase("runs")){
										dgeOptions.setRunType(parser.asString(1));
									}
									//DGE PATH
									else if(parser.asString(0).equalsIgnoreCase("dge_jar_path")){
										dgeOptions.setdgePath(parser.asString(1));
									}
									else{
										dgeOptions.addOption(parser.asString(0), parser.asString(1));
									}
								}
								else{
									throw new IllegalArgumentException("Incorrect number of columns for DGE options");
								}
							}
						}
					}
				}
			}
		}
		
	/*(	for(String c:dgeOptions.options.keySet()){
			logger.info(c+" : "+dgeOptions.options.get(c));
		}*/
		
		if(species.isEmpty()){
			species = DEFAULT_SPECIES;
		}
		logger.info("Config File Read");

/*		for(String c:basicOptions.getAllOptionMaps().keySet()){
			logger.info(basicOptions.getAllOptionMaps().get(c));
		}
		for(String c:this.commands){
			logger.info(c);
		}
*/
	}
	
	/**
	 * Returns true is the specified line is a comment
	 * @param nextLine
	 * @return
	 */
	private boolean isComment(String nextLine){
		return(nextLine.startsWith("#"));
	}
	/**
	 * Returns true if the input line is the start of a new section which is not the passed section
	 * @param nextLine
	 * @return
	 */
	private boolean startOfNextSection(String nextLine,String currentSection) {
		return (nextLine.trim().length() > 0 && nextLine.startsWith(":") && !(nextLine.contains(currentSection)&& nextLine.contains(species)));
	}
	
	/**
	 * Returns true if the input line is the start of a new section 
	 * @param nextLine
	 * @return
	 */
	private boolean startOfASection(String nextLine) {
		return (nextLine.trim().length() > 0 && nextLine.startsWith(":") );
	}

	/**
	 * Returns true if the input line indicates the start of the command section
	 * @param nextLine
	 * @return
	 */
	private boolean isCommandSection(String nextLine) {
		return nextLine.trim().length() > 0 && "commands".equalsIgnoreCase(nextLine.substring(START_OF_STRING));
	}

	/**
	 * Returns true if the input line indicates the start of the species section
	 * @param nextLine
	 * @return
	 */
	private boolean isSpeciesSection(String nextLine) {
		return nextLine.trim().length() > 0 && "species".equalsIgnoreCase(nextLine.substring(START_OF_STRING));
	}
	
	/**
	 * Returns true if the input line indicates the start of the options for species section
	 * @param nextLine
	 * @return
	 */
	private boolean isBasicOptionsForSpecies(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(species) && nextLine.contains("basic");
	}
	
	/**
	 * Returns true if the input line indicates the start of the DGE options for species section
	 * @param nextLine
	 * @return
	 */
	private boolean isDGEOptionsForSpecies(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(species) && nextLine.contains("dge");
	}

	/**
	 * Adds the command in the current line if it is not prefixed with a "#"
	 * @param nextLine
	 */
	public void addCommand(String nextLine) {
		if (nextLine.trim().length() > 0 && !isComment(nextLine)){
			this.commands.add(nextLine);
		}
	}
	
	/**
	 * Adds the species
	 * @param nextLine
	 */
	public void addSpecies(String nextLine) {
		if (nextLine.trim().length() > 0 ){
			this.species = (nextLine.trim());
		}
	}

	/**
	 * Returns true if the specified command was included as runnable in the file
	 * @param command
	 * @return
	 */
	public boolean hasCommandFor(String command){
		return (this.commands.contains(command));
	}
	
	/**
	 * Returns the command at the specified index
	 * @param index
	 * @return
	 */
	public String getCommand(int index){
		return (this.commands.get(index));
	}

	/**
	 * Returns all commands
	 * @return
	 */
	public List<String> getAllCommands(){
		return (this.commands);
	}
	
	/**
	 * Returns all run specific DGE options
	 * @return
	 */
	public Set<String> getRunSpecificDGEOptions(){
		return (this.dgeOptions.runSpecificKeySet());
	}
	
	/**
	 * Returns true if the input line indicates that it contains the option for a rna class file
	 * @param nextLine
	 * @return
	 */
	private boolean isRNAClass(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_RNA_CLASS);
	}

	/**
	 * Returns true if the input line indicates that it contains the option queue name
	 * @param nextLine
	 * @return
	 */
	private boolean isQueueName(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_QUEUE_NAME);
	}
	
	/**
	 * Returns true if the input line indicates that it contains a tophat option 
	 * @param nextLine
	 * @return
	 */
	private boolean isTophatOption(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_TOPHAT_OPTIONS);
	}
		
	/**
	 * Returns true if the input line indicates that it contains a bowtie2 option 
	 * @param nextLine
	 * @return
	 */
	private boolean isBowtie2Option(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_BOWTIE2_OPTIONS);
	}

	
	/**
	 * Returns true if the input line indicates that it contains a novoalign option 
	 * @param nextLine
	 * @return
	 */
	private boolean isNovoalignOption(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_NOVOALIGN_OPTIONS);
	}

	/**
	 * Returns true if the input line indicates that it contains an option for calculating the fragment size distribution
	 * @param nextLine
	 * @return
	 */
	private boolean isFragmentSizeDistOption(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_FRAGMENT_SIZE_DIST_OPTIONS);
	}
	
	/**
	 * Returns true if the input line indicates that it contains an option for merging samples
	 * @param nextLine
	 * @return
	 */
	private boolean isMergeSamplesOption(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_MERGE_SAMPLES);
	}
	
	
	
	
	/**
	 * Returns true if the input line indicates that it contains a tophat path
	 * @param nextLine
	 * @return
	 */
	private boolean isTophatPath(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_TOPHAT_PATH);
	}

	/**
	 * Returns true if the input line indicates that it contains a bowtie2 path
	 * @param nextLine
	 * @return
	 */
	private boolean isBowtie2Path(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_BOWTIE2_PATH);
	}

	
	/**
	 * Returns true if the input line indicates that it contains a novoalign path
	 * @param nextLine
	 * @return
	 */
	private boolean isNovoalignPath(String nextLine) {
		return nextLine.trim().length() > 0 && nextLine.contains(OPTION_NOVOALIGN_PATH);
	}

	
	/**
	 * This is a helper class is a data structure mapping an RNA class to its fasta file path
	 * @author skadri
	 *
	 */
	public static class RNAClass extends HashMap<String, String> {

		public RNAClass() {
			super();
		}
		
		public RNAClass(HashMap<String,String> classFile){
			super(classFile);
		}
		
		public RNAClass(String rnaClass,String fileName){
			super();
			super.put(rnaClass, fileName);
		}
		
		public String get(String key) {
			String result = super.get(key);
			return result == null ? "" : result;
		}
	}
	
	/**
	 * This is a helper class is a data structure that maps a tophat flag to its value
	 * @author skadri
	 *
	 */
	public static class TophatOptions {
		
		HashMap<String, String> options; 
		String tophatPath;
		
		public TophatOptions(){
			options = new HashMap<String, String>();
		}
		
		public TophatOptions(HashMap<String,String> option){
			options = option;
		}
		
		public TophatOptions(String flag,String value){
			options = new HashMap<String, String>();
			options.put(flag, value);
		}
		
		public void setTophatPath(String path){
			tophatPath = path;	
		}
		
		public String getTophatPath(){
			return tophatPath;
		}
		
		public String get(String key) {
			String result = options.get(key);
			return result == null ? "" : result;
		}
		
		public HashMap<String,String> getAllOptions(){
			return options;
		}
		
		public void addOption(String flag,String value){
			options.put(flag,value);
		}
	}
	
	/**
	 * Sets of samples to merge into one bam file (leaving separate bam files too)
	 * @author prussell
	 *
	 */
	public static class MergeSamples {
		
		Collection<Collection<String>> sampleNameSets;
		
		public MergeSamples() {
			sampleNameSets = new ArrayList<Collection<String>>();
		}
		
		public void addSet(Collection<String> sampleNames) {
			sampleNameSets.add(sampleNames);
		}
		
		private static String makeSampleName(Collection<String> sampleNames) {
			String rtrn = "Merged";
			for(String name : sampleNames) {
				rtrn += "_" + name;
			}
			return rtrn;
		}
		
		public Map<String, Collection<String>> getSetsToMerge() {
			Map<String, Collection<String>> rtrn = new TreeMap<String, Collection<String>>();
			for(Collection<String> set : sampleNameSets) {
				rtrn.put(makeSampleName(set), set);
			}
			return rtrn;
		}
		
	}
	
	
	/**
	 * This is a helper class is a data structure that maps a bowtie2 flag to its value
	 * @author skadri
	 *
	 */
	public static class Bowtie2Options {
		
		HashMap<String, String> options; 
		String bowtie2Path;
		
		public Bowtie2Options(){
			options = new HashMap<String, String>();
		}
		
		public Bowtie2Options(HashMap<String,String> option){
			options = option;
		}
		
		public Bowtie2Options(String flag,String value){
			options = new HashMap<String, String>();
			options.put(flag, value);
		}
		
		public void setBowtie2Path(String path){
			bowtie2Path = path;	
		}
		
		public String getBowtie2Path(){
			return bowtie2Path;
		}
		
		public String get(String key) {
			String result = options.get(key);
			return result == null ? "" : result;
		}
		
		public HashMap<String,String> getAllOptions(){
			return options;
		}
		
		public void addOption(String flag,String value){
			options.put(flag,value);
		}
	}

	
	
	public static class FragmentSizeDistributionOptions {
				
		protected static String OPTION_MAX_GENOMIC_SPAN = "max_genomic_span";
		protected static String OPTION_INCLUDE_IMPROPER_PAIRS = "include_improper_pairs";
		protected static String OPTION_BED_FILE = "bed_file";
		protected static String OPTION_MAX_FRAGMENT_SIZE = "max_fragment_size";
		protected static String OPTION_NUM_BINS = "num_bins";
		protected static String OPTION_IND_GENE = "ind_gene";
		
		private final int DEFAULT_MAX_GENOMIC_SPAN = 300000;
		private final int DEFAULT_MAX_FRAGMENT_SIZE = 300;
		private final boolean DEFAULT_PROPER_PAIRS_ONLY = true;
		private final int DEFAULT_NUM_BINS = 10;
		
		private int maxGenomicSpan;
		private int maxFragmentSize;
		private String bedFile;
		private boolean properPairsOnly;
		private int numBins;
		private Collection<String> indGeneNames;
		
		public FragmentSizeDistributionOptions() {
			this.maxGenomicSpan = this.DEFAULT_MAX_GENOMIC_SPAN;
			this.maxFragmentSize = this.DEFAULT_MAX_FRAGMENT_SIZE;
			this.properPairsOnly = this.DEFAULT_PROPER_PAIRS_ONLY;
			this.numBins = this.DEFAULT_NUM_BINS;
			this.bedFile = null;
			this.indGeneNames = new TreeSet<String>();
		}
		
		public void setMaxGenomicSpan(int span) {
			this.maxGenomicSpan = span;
		}
		public void setMaxFragmentSize(int size) {
			this.maxFragmentSize = size;
		}
		public void setBedFile(String bed) {
			this.bedFile = bed;
		}
		public void setProperPairsOnly(boolean p) {
			this.properPairsOnly = p;
		}
		public void setNumBins(int bins) {
			this.numBins = bins;
		}
		public void addIndividualGene(String geneName) {
			this.indGeneNames.add(geneName);
		}
		
		public int getMaxGenomicSpan() {
			return this.maxGenomicSpan;
		}
		public int getMaxFragmentSize() {
			return this.maxFragmentSize;
		}
		
		public boolean hasOptionMakeFragmentSizeDistribution() {
			return this.bedFile != null;
		}
		
		public String getBedFile() {
			if(this.bedFile == null) {
				throw new IllegalArgumentException("Must specify a bed file for fragment size distribution with option " + ConfigFileParser.OPTION_FRAGMENT_SIZE_DIST_OPTIONS + "\t" + OPTION_BED_FILE);
			}
			return this.bedFile;
		}
		public boolean getProperPairsOnly() {
			return this.properPairsOnly;
		}
		public int getNumBins() {
			return this.numBins;
		}
		public Collection<String> getIndividualGeneNames() {
			return this.indGeneNames;
		}
		
		
	}

	/**
	 * This is a helper class is a data structure that maps a novoalign flag to its value
	 * @author prussell
	 *
	 */
	public static class NovoalignOptions {
		
		HashMap<String, String> options; 
		String novoalignPath;
		
		/**
		 * Instantiate object
		 */
		public NovoalignOptions(){
			this.options = new HashMap<String, String>();
		}
		
		/**
		 * Instantiate object with existing set of options
		 * @param option Map of flags to values
		 */
		public NovoalignOptions(HashMap<String,String> option){
			this.options = option;
		}
		
		/**
		 * Instantiate object with one existing option
		 * @param flag The flag
		 * @param value The value
		 */
		public NovoalignOptions(String flag,String value){
			this.options = new HashMap<String, String>();
			this.options.put(flag, value);
		}
		
		/**
		 * Set path to novoalign executable
		 * @param path The path
		 */
		public void setNovoalignPath(String path){
			this.novoalignPath = path;	
		}
		
		/**
		 * Get the path to novoalign executable
		 * @return The path to novoalign executable
		 */
		public String getNovoalignPath(){
			return this.novoalignPath;
		}
		
		/**
		 * Get all flags and options
		 * @return Map of flags to values
		 */
		public HashMap<String,String> getAllOptions(){
			return this.options;
		}
		
		/**
		 * Add a flag and value
		 * @param flag The flag
		 * @param value The value
		 */
		public void addOption(String flag,String value){
			this.options.put(flag,value);
		}
	}

	
	/**
	 * This is a helper class that maps the values of species-specific basic options to their value
	 * @author skadri
	 *
	 */
	public static class BasicOptions {
		
		Map<String, String> optionMap;
		RNAClass classFileNames;
		TophatOptions thOptions;
		Bowtie2Options bt2Options;
		NovoalignOptions novoalignOptions;
		MergeSamples mergeSamples;
		String queueName;
		
		public Map<String, String> getRNAClassFileNames() {
			return this.classFileNames;
		}
		
		public BasicOptions(){
			classFileNames = new RNAClass();
			thOptions = new TophatOptions();
			bt2Options = new Bowtie2Options();
			novoalignOptions = new NovoalignOptions();
			optionMap = new HashMap<String,String>();
			mergeSamples = new MergeSamples();
			queueName = "";
		}

		public void addOption(String flag, String value){
			optionMap.put(flag, value);
		}
		
		public String getOptionValueFor(String type){
			String result = optionMap.get(type);
			if(result == null) {
				throw new IllegalArgumentException("Please supply value for flag: "+ type);
			}
			return result;
		}
		
		/**
		 * Whether the flag has been specified
		 * @param flag The flag
		 * @return Whether a value has been provided
		 */
		public boolean hasCommandFor(String flag) {
			return optionMap.containsKey(flag);
		}
		
		/**
		 * Get path to genome fasta file
		 * @return The genome fasta file path specified in config file
		 */
		public String getGenomeFasta() {
			try {
				return this.optionMap.get(OPTION_GENOME_FASTA);
			} catch(NullPointerException e) {
				throw new NullPointerException("Genome fasta file must be specified with option " + OPTION_GENOME_FASTA);
			}
		}
		
		/**
		 * Get name of genome assembly
		 * @return The genome assembly name specified in config file
		 */
		public String getGenomeAssembly() {
			try {
				return this.optionMap.get(OPTION_GENOME_ASSEMBLY);
			} catch(NullPointerException e) {
				throw new NullPointerException("Genome assembly name must be specified with option " + OPTION_GENOME_ASSEMBLY);
			}
		}

		/**
		 * Get the delimiter for read pair number in fastq file
		 */
		public String getFastqReadIdPairDelimiter() {
			if(this.optionMap.containsKey(OPTION_FASTQ_READ_PAIR_NUMBER_DELIMITER)) {
				return this.optionMap.get(OPTION_FASTQ_READ_PAIR_NUMBER_DELIMITER);
			}
			return DEFAULT_FASTQ_READ_PAIR_NUMBER_DELIMITER;
		}

		/**
		 * Get path to genome bowtie index
		 * @return The path to bowtie index specified in config file
		 */
		public String getGenomeBowtieIndex() {
			try {
				return this.optionMap.get(OPTION_GENOME_BOWTIE);
			} catch(NullPointerException e) {
				throw new NullPointerException("Genome bowtie index base must be specified with option " + OPTION_GENOME_BOWTIE);
			}
		}
		
		/**
		 * Get path to genome novoindex
		 * @return The path to novoindex specified in config file
		 */
		public String getGenomeNovoindex() {
			try {
				return this.optionMap.get(OPTION_GENOME_NOVOINDEX);
			} catch(NullPointerException e) {
				throw new NullPointerException("Genome novoindex base must be specified with option " + OPTION_GENOME_NOVOINDEX);
			}
		}

		/**
		 * Get path to ref_flat annotation
		 * @return The path to ref_flat annotation specified in config file
		 */
		public String getPicardRefFlat() {
			try {
				return this.optionMap.get(OPTION_PICARD_REF_FLAT);
			} catch(NullPointerException e) {
				throw new NullPointerException("Picard ref flat file must be specified with option " + OPTION_PICARD_REF_FLAT);
			}
		}
		
		/**
		 * Get strand specificity argument for picard
		 * @return The strand specificity argument specified in config file
		 */
		public String getPicardStrandSpecificity() {
			try {
				return this.optionMap.get(OPTION_PICARD_STRAND_SPECIFICITY);
			} catch(NullPointerException e) {
				throw new NullPointerException("Picard strand specificity argument must be specified with option " + OPTION_PICARD_STRAND_SPECIFICITY);
			}
		}

		/**
		 * Get path to ribosomal intervals for picard
		 * @return The path to ribosomal intervals file specified in config file, or null if not provided
		 */
		public String getPicardRibosomalIntervals() {
			try {
				return this.optionMap.get(OPTION_PICARD_RIBOSOMAL_INTERVALS);
			} catch(NullPointerException e) {
				return null;
			}
		}

		/**
		 * Get path to alternative jar file for CollectRnaSeqMetrics
		 * @return The path to alternative jar file, or null if not provided
		 */
		public String getPicardCollectRnaSeqMetrics() {
			try {
				return this.optionMap.get(OPTION_PICARD_COLLECT_RNASEQ_METRICS);
			} catch(NullPointerException e) {
				return null;
			}
		}

		/**
		 * Get path to igvtools executable file
		 * @return The path to igvtools executable specified in config file
		 */
		public String getIgvtoolsPath() {
			try {
				return this.optionMap.get(OPTION_IGVTOOLS_PATH);
			} catch(NullPointerException e) {
				throw new NullPointerException("Igvtools path must be specified with option " + OPTION_IGVTOOLS_PATH);
			}
		}
		
		/**
		 * Whether a novoalign executable was specified in config file
		 * @return True if novoalign executable was specified, false otherwise
		 */
		public boolean hasNovoalignPath() {
			try {
				if(this.novoalignOptions.novoalignPath.length() > 0) return true;
			} catch(Exception e) {}
			return false;
		}

		/**
		 * Get path to novoalign executable
		 * @return The path to novoalign executable file specified in config file
		 */
		public String getNovoalignPath() {
			try {
				return this.novoalignOptions.novoalignPath;
			} catch(NullPointerException e) {
				throw new NullPointerException("Novoalign path must be specified with option " + OPTION_NOVOALIGN_PATH);
			}
		}

		public String getQueueName() {
			if(this.queueName.isEmpty())
				return DEFAULT_QUEUE_NAME;
			else
				return this.queueName;
		}
		
		/**
		 * Get path to bowtie2 executable file
		 * @return The path to bowtie2 executable specified in config file
		 */
		public String getBowtie2ExecutablePath() {
			try {
				return this.optionMap.get(OPTION_BOWTIE2_PATH);
			} catch(NullPointerException e) {
				throw new NullPointerException("Bowtie2 executable file must be specified with option " + OPTION_BOWTIE2_PATH);
			}
		}
		
		/**
		 * Get path to samtools executable
		 * @return The path to samtools executable specified in config file
		 */
		public String getSamtoolsPath() {
			try {
				return this.optionMap.get(OPTION_SAMTOOLS_PATH);
			} catch(NullPointerException e) {
				throw new NullPointerException("Samtools executable file must be specified with option " + OPTION_SAMTOOLS_PATH);
			}
		}

		/**
		 * Get directory containing picard executables
		 * @return The path to picard directory specified in config file
		 */
		public String getPicardDirectory() {
			try {
				return this.optionMap.get(OPTION_PICARD_DIRECTORY);
			} catch(NullPointerException e) {
				throw new NullPointerException("Directory containing Picard executables must be specified with option " + OPTION_PICARD_DIRECTORY);
			}
		}

		/**
		 * Get directory containing fastx executables
		 * @return The path to fastx directory specified in config file
		 */
		public String getFastxDirectory() {
			try {
				return this.optionMap.get(OPTION_FASTX_DIRECTORY);
			} catch(NullPointerException e) {
				throw new NullPointerException("Directory containing fastx executables must be specified with option " + OPTION_FASTX_DIRECTORY);
			}
		}

		/**
		 * Get sequencing adapter sequence for read 1
		 * @return The adapter sequence
		 */
		public String getSequencingAdapterRead1() {
			try {
				return this.optionMap.get(OPTION_SEQUENCING_ADAPTER_READ_1);
			} catch(NullPointerException e) {
				throw new NullPointerException("Sequencing adapter for read 1 must be specified with option " + OPTION_SEQUENCING_ADAPTER_READ_1);
			}
		}

		/**
		 * Get sequencing adapter sequence for read 2
		 * @return The adapter sequence
		 */
		public String getSequencingAdapterRead2() {
			try {
				return this.optionMap.get(OPTION_SEQUENCING_ADAPTER_READ_2);
			} catch(NullPointerException e) {
				throw new NullPointerException("Sequencing adapter for read 2 must be specified with option " + OPTION_SEQUENCING_ADAPTER_READ_2);
			}
		}

		
		
		/**
		 * Get bed file that determines which genes to count when writing fragment end positions to wig and bigwig files
		 * @return The bed file path
		 */
		public String getBedFileForFragmentEndWig() {
			try {
				return this.optionMap.get(OPTION_BED_FILE_FOR_FRAGMENT_END_WIG);
			} catch(NullPointerException e) {
				throw new NullPointerException("Bed file to use when writing fragment endpoints to wig file must be specified with option " + OPTION_BED_FILE_FOR_FRAGMENT_END_WIG);
			}
		}
		
		/**
		 * Get UCSC WigToBigWig executable file
		 * @return The WigToBigWig file
		 */
		public String getWigToBigWigExecutable() {
			try {
				return this.optionMap.get(OPTION_WIG_TO_BIGWIG_PATH);
			} catch(NullPointerException e) {
				throw new NullPointerException("WigToBigWig executable file must be specified with option " + OPTION_WIG_TO_BIGWIG_PATH);
			}
		}
		
		/**
		 * Get chromosome size file
		 * @return The chromosome size file
		 */
		public String getChrSizeFile() {
			try {
				return this.optionMap.get(OPTION_CHR_SIZE_FILE);
			} catch(NullPointerException e) {
				throw new NullPointerException("Chromosome size file must be specified with option " + OPTION_CHR_SIZE_FILE);
			}
		}
		
		
		
		/**
		 * Get the path to bowtie2-build executable
		 * @return The path to bowtie2-build executable specified in config file
		 */
		public String getBowtie2BuildExecutablePath() {
			try {
				return optionMap.get(OPTION_BOWTIE2_BUILD_PATH);
			} catch(NullPointerException e) {
				throw new NullPointerException("Bowtie2-build path must be specified with option " + OPTION_BOWTIE2_BUILD_PATH);
			}
		}
		
		/**
		 * Add an RNA class fasta file
		 * @param rnaClass The class name
		 * @param fileName The fasta file
		 */
		public void addRNAClass(String rnaClass,String fileName){
			this.classFileNames.put(rnaClass,fileName);
		}
		
		/**
		 * Get fasta file for an RNA class
		 * @param rnaClass The class name
		 * @return The fasta file
		 */
		public String getRNAClassFile(String rnaClass){
			try {
				return this.classFileNames.get(rnaClass);
			} catch(NullPointerException e) {
				throw new NullPointerException("RNA class file " + rnaClass + " was not specified in config file.");
			}
		}
		
		/**
		 * Specify a tophat option
		 * @param flag The flag
		 * @param value The value
		 */
		public void addTophatOption(String flag,String value){
			this.thOptions.addOption(flag,value);
		}
		
		public void addMergeSamples(String entireLine) {
			StringParser p = new StringParser();
			p.parse(entireLine);
			if(!p.asString(0).equals(OPTION_MERGE_SAMPLES) || p.getFieldCount() < 3) {
				throw new IllegalArgumentException("Not a valid line containing samples to merge");
			}
			Collection<String> set = new TreeSet<String>();
			for(int i=1; i<p.getFieldCount(); i++) {
				set.add(p.asString(i));
			}
			mergeSamples.addSet(set);
		}
		
		/**
		 * Specify a bowtie2 option
		 * @param flag The flag
		 * @param value The value
		 */
		public void addBowtie2Option(String flag,String value){
			this.bt2Options.addOption(flag,value);
		}
		
		/**
		 * Specify a novoalign option
		 * @param flag The flag
		 * @param value The value including spaces if necessary
		 */
		public void addNovoalignOption(String flag,String value){
			this.novoalignOptions.addOption(flag,value);
		}

		/**
		 * Set path to tophat executable
		 * @param path The path
		 */
		public void setTophatPath(String path){
			this.thOptions.setTophatPath(path);
		}
		
		/**
		 * Set path to novoalign executable
		 * @param path The path
		 */
		public void setNovoalignPath(String path){
			this.novoalignOptions.setNovoalignPath(path);
		}

		public void setQueueName(String name){
			this.queueName = name;
		}
		
		/**
		 * Get path to tophat executable
		 * @return The path
		 */
		public String getTophatPath(){
			try {
				return this.thOptions.getTophatPath();
			} catch(NullPointerException e) {
				throw new NullPointerException("Tophat path must be specified with option " + OPTION_TOPHAT_PATH);
			}
		}
				
		/**
		 * Get all tophat flags and values
		 * @return Map of flags to values
		 */
		public Map<String, String> getTophatOptions() {
			return this.thOptions.getAllOptions();
		}
		
		/**
		 * Get all sets of samples to merge, by merged sample name
		 * Also will keep individual samples
		 * @return Sets of samples to merge by new name of merged sample
		 */
		public Map<String, Collection<String>> getSamplesToMerge() {
			return this.mergeSamples.getSetsToMerge();
		}
		
		/**
		 * Get all bowtie2 flags and values
		 * @return Map of flags to values
		 */
		public Map<String, String> getBowtie2Options() {
			return this.bt2Options.getAllOptions();
		}
		
		/**
		 * Get all novoalign flags and values
		 * @return Map of flags to values
		 */
		public Map<String, String> getNovoalignOptions() {
			return this.novoalignOptions.getAllOptions();
		}

		/**
		 * Get ribosomal RNA sequences to use for in silico rRNA depletion
		 * @return The fasta file of sequences
		 */
		public String getRrnaSeqsForDepletion() {
			try {
				return this.optionMap.get(OPTION_FILTER_RRNA);
			} catch(NullPointerException e) {
				throw new NullPointerException("Fasta file of rRNA sequences for depletion must be specified with option " + OPTION_FILTER_RRNA);
			}
		}
		
		/**
		 * Get fasta file of transcript sequences to align to
		 * @return The fasta file of sequences
		 */
		public String getTranscriptSeqsForAlignment() {
			try {
				return this.optionMap.get(OPTION_ALIGN_TO_TRANSCRIPTS);
			} catch(NullPointerException e) {
				throw new NullPointerException("Fasta file of transcript sequences to align to must be specified with option " + OPTION_ALIGN_TO_TRANSCRIPTS);
			}
		}

		
		/**
		 * Get all options in one map
		 * @return Map of all flags for all programs to their values
		 */
		public Map<String,String> getAllOptionMaps(){
			Map<String,String> m = new HashMap<String,String>();
			m.putAll(optionMap);
			m.putAll(classFileNames);
			m.putAll(novoalignOptions.getAllOptions());
			m.putAll(thOptions.getAllOptions());
			m.putAll(bt2Options.getAllOptions());
			return m;
		}
	}
	
	/**
	 * This is a helper class is a data structure that maps a DGE flag or option to its value or path
	 * @author skadri
	 *
	 */
	public static class DGEOptions {
		
		Map<String, String> options;
		//3' or 5'
		String dgeType;
		//single or multiple
		String runType;
		String dgeJarPath;
		Set<String> multipleSpecificOptions = new HashSet<String>(Arrays.asList
				(new String[] {"out","normalizedOutput","scoreFullGene","conditions"}));
		
		public DGEOptions(HashMap<String,String> option){
			options = new HashMap<String, String>(option);
		}
		
		public DGEOptions(){
			options = new HashMap<String, String>();
		}
		
		public DGEOptions(String flag,String value){
			options = new HashMap<String, String>();
			options.put(flag,value);
		}
		
		public void setRunType(String type){
			runType = type;
		}
		
		public String getdgeType(){
			return dgeType;
		}
		
		public String getrunType(){
			return runType;
		}
		
		public void setdgeType(String type){
			dgeType = type;
		}
		
		public void setdgePath(String path){
			dgeJarPath = path;
		}
		
		public String getdgePath(){
			return dgeJarPath;
		}

		public String getOption(String key) {
			String result = options.get(key);
			return result == null ? "" : result;
		}
		
		public Set<String> keySet() {
			return options.keySet();
		}
		
		public void setOutputInDirectory(){
			String outputName = options.get("out");
			File resultDirFile = new File(outputName+".dge.results/");
			resultDirFile.mkdir();
			String newOutputName = outputName+".dge.results/"+outputName;
			options.put("out", newOutputName);
		}
		
		public String getOutput(){
			return options.get("out");
		}
		/**
		 * If run type is multiple, returns all options
		 * If run type is single, returns all except those options in multipleSpecificOptions
		 * @return
		 */
		public Set<String> runSpecificKeySet(){
			Set<String> ops = options.keySet();
			
			if(runType.equalsIgnoreCase("single"))
				ops.removeAll(multipleSpecificOptions);
			
			return ops;
		}
		
		public Collection<String> values() {
			return options.values();
		}
		
		public void addOption(String flag,String value){
			options.put(flag,value);
		}
	}

}

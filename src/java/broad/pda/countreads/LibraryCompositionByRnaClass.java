/**
 * 
 */
package broad.pda.countreads;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.LSFJob;
import nextgen.core.pipeline.Scheduler;
import nextgen.core.pipeline.util.AlignmentUtils;

import org.apache.log4j.Logger;

import broad.core.parser.StringParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 *
 */
public class LibraryCompositionByRnaClass {

	/**
	 * The name used for the genome when counting read mappings
	 */
	public static String GENOME_CLASS_NAME = "Genome";
	
	/**
	 * The name used for unmapped reads when counting read mappings
	 */
	public static String MYSTERY_READS_CLASS_NAME = "Other";
	
	private String genomeBowtieIndex;
	private Map<String, String> classSequenceFiles;
	private Set<String> classNames;
	private Set<String> sampleNames;
	private Map<String, String> read1FastqFiles;
	private Map<String, String> read2FastqFiles;
	private Map<String, Integer> totalReadCounts;
	private boolean readsPaired;
	private Logger logger;
	
	/**
	 * Constructor establishes locations of files and counts total reads per sample
	 * @param genomeBowtieIndexBase Genome fasta file
	 * @param rnaClasses Map associating RNA class name with fasta file of sequences
	 * @param read1fastqs Fastq file containing read1 if paired or single end reads if unpaired
	 * @param read2fastqs Fastq file containing read2 if paired or null if unpaired
	 * @param log Logger object
	 * @throws IOException 
	 */
	public LibraryCompositionByRnaClass(String genomeBowtieIndexBase, Map<String,String> rnaClasses, Map<String,String> read1fastqs, Map<String,String> read2fastqs, Logger log) throws IOException {
		
		this.genomeBowtieIndex = genomeBowtieIndexBase;
		this.classSequenceFiles = rnaClasses;
		this.classNames = this.classSequenceFiles.keySet();
		this.read1FastqFiles = read1fastqs;
		this.sampleNames = this.read1FastqFiles.keySet();
		this.read2FastqFiles = read2fastqs;
		this.readsPaired = true;
		if(this.read2FastqFiles == null) this.readsPaired = false;
		this.logger = log;
		
		// Count total number of reads for each sample
		this.totalReadCounts = new TreeMap<String, Integer>();
		for(String sample : this.read1FastqFiles.keySet()) {
			logger.info("Counting total reads for sample " + sample + "...");
			FileReader r = new FileReader(this.read1FastqFiles.get(sample));
			BufferedReader b = new BufferedReader(r);
			int numLines = 0;
			while(b.ready()) {
				String line = b.readLine();
				numLines ++;
			}
			this.totalReadCounts.put(sample, this.readsPaired ? Integer.valueOf(numLines/2) : Integer.valueOf(numLines/4));
			if(!this.readsPaired) this.logger.info("Counted " + this.totalReadCounts.get(sample) + " total reads for sample " + sample + ".");
			else this.logger.info("Counted " + this.totalReadCounts.get(sample) + " total reads for sample " + sample + " (" + this.totalReadCounts.get(sample).intValue()/2 + " pairs).");
		}
		
	}
	
	/**
	 * Get the total number of reads per sample
	 * @return Map associating sample name with read count; regardless of whether reads are paired, the count is total single reads
	 */
	public Map<String, Integer> getTotalReadCounts() {
		return this.totalReadCounts;
	}
	
	/**
	 * Whether the reads are paired
	 * @return True if paired, false if unpaired
	 */
	public boolean readsPaired() {
		return this.readsPaired;
	}
	
	/**
	 * Combine all RNA classes into a single fasta file with the class name as the sequence name
	 * @param outFasta Output combined fasta file
	 * @throws IOException
	 */
	private void makeSingleFasta(String outFasta) throws IOException {
		
		List<Sequence> allSequences = new ArrayList<Sequence>();
		
		// Read sequences for each class
		for(String className : this.classSequenceFiles.keySet()) {
			FastaSequenceIO fsio = new FastaSequenceIO(new File(this.classSequenceFiles.get(className)));
			List<Sequence> thisClassSeqs = fsio.loadAll();
			// Rename sequences to the class name
			for(Sequence seq : thisClassSeqs) {
				seq.setId(className);
			}
			allSequences.addAll(thisClassSeqs);
		}
		
		// Write all sequences to fasta file
		FastaSequenceIO out = new FastaSequenceIO(new File(outFasta));
		out.write(allSequences);
	}
	
	

	/**
	 * Count alignments to RNA classes in a sam file
	 * @param alignmentFile The sam file
	 * @return Map associating each RNA class with the number of alignments in the sam file
	 * @throws IOException
	 */
	private Map<String, Integer> countClasses(String alignmentFile) throws IOException {
		
		// Count classes in a sam file
		// Ignore lines starting with @
		
		FileReader r = new FileReader(alignmentFile);
		BufferedReader b = new BufferedReader(r);
		StringParser p = new StringParser();
		
		Map<String, Integer> rtrn = new TreeMap<String, Integer>();
		for(String className : this.classNames) rtrn.put(className, Integer.valueOf(0));
		
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			if(p.getFieldCount() < 3) continue;
			String alignedClass = p.asString(2);
			for(String className : this.classNames) {
				if(alignedClass.equals(className)) {
					Integer newCount = Integer.valueOf(rtrn.get(alignedClass).intValue() + 1);
					rtrn.put(alignedClass, newCount);
				}
			}
		}
		
		b.close();
		
		return rtrn;
		
	}
	
	
	
	/**
	 * For each sample, get counts of reads mapping to each RNA class as well as genome and unmapped
	 * @param samtoolsExecutable Samtools executable file
	 * @param bowtie2Executable Bowtie2 executable file
	 * @param bowtie2options 
	 * @param bowtie2BuildExecutable Bowtie2-build executable file
	 * @param logDir Output directory for logs and alignments
	 * @param scheduler Scheduler
	 * @return Map from sample name to class name to count
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public Map<String, Map<String, Integer>> alignAndGetCounts(String samtoolsExecutable, String bowtie2Executable, Map<String, String> bowtie2options, String bowtie2BuildExecutable, String logDir, Scheduler scheduler) throws IOException, InterruptedException {
		
		File dir = new File(logDir);
		boolean madeDir = dir.mkdir();
		
		// Make single fasta file for transcriptome
		String singleFasta = logDir + "/combined_RNA_classes.fa";
		makeSingleFasta(singleFasta);		
		
		// Make bowtie2 index for transcriptome
		String btBase = logDir + "/combined_RNA_classes";
		String bt1 = btBase + ".1.bt2";
		String bt2 = btBase + ".2.bt2";
		String bt3 = btBase + ".3.bt2";
		String bt4 = btBase + ".4.bt2";
		String rev1 = btBase + ".rev.1.bt2";
		String rev2 = btBase + ".rev.2.bt2";
		boolean writeIndex = false;
		if(!new File(bt1).exists()) writeIndex = true;
		if(!new File(bt2).exists()) writeIndex = true;
		if(!new File(bt3).exists()) writeIndex = true;
		if(!new File(bt4).exists()) writeIndex = true;
		if(!new File(rev1).exists()) writeIndex = true;
		if(!new File(rev2).exists()) writeIndex = true;
		
		if(writeIndex) AlignmentUtils.makeBowtie2Index(singleFasta, btBase, bowtie2BuildExecutable, logDir, scheduler);
		else {
			logger.warn("Bowtie2 index files " + btBase + ".*.bt2 already exist. Not remaking index.");
		}
		
		// Establish file names
		Map<String, String> samToTranscriptome = new TreeMap<String, String>();
		Map<String, String> fastqNotTranscriptomeUnpaired = new TreeMap<String, String>();
		Map<String, String> fastqNotTranscriptomePairedBase = new TreeMap<String, String>();
		Map<String, String> fastqNotTranscriptomePaired1 = new TreeMap<String, String>();
		Map<String, String> fastqNotTranscriptomePaired2 = new TreeMap<String, String>();
		Map<String, String> samToGenome = new TreeMap<String, String>();
		Map<String, String> fastqNotGenomeUnpaired = new TreeMap<String, String>();
		Map<String, String> fastqNotGenomePairedBase = new TreeMap<String, String>();
		Map<String, String> fastqNotGenomePaired1 = new TreeMap<String, String>();
		Map<String, String> fastqNotGenomePaired2 = new TreeMap<String, String>();
		for(String sample : this.sampleNames) {
			samToTranscriptome.put(sample, logDir + "/" + sample + "_to_transcriptome.sam");
			fastqNotTranscriptomeUnpaired.put(sample, logDir + "/" + sample + "_not_aligned_to_transcriptome.fq");
			fastqNotTranscriptomePairedBase.put(sample, logDir + "/" + sample + "_not_aligned_to_transcriptome_%.fq");
			fastqNotTranscriptomePaired1.put(sample, logDir + "/" + sample + "_not_aligned_to_transcriptome_1.fq");
			fastqNotTranscriptomePaired2.put(sample, logDir + "/" + sample + "_not_aligned_to_transcriptome_2.fq");
			samToGenome.put(sample, logDir + "/" + sample + "_non_transcriptome_to_genome.sam");
			fastqNotGenomeUnpaired.put(sample, logDir + "/" + sample + "_mystery_reads.fq");
			fastqNotGenomePairedBase.put(sample, logDir + "/" + sample + "_mystery_reads_%.fq");
			fastqNotGenomePaired1.put(sample, logDir + "/" + sample + "_mystery_reads_1.fq");
			fastqNotGenomePaired2.put(sample, logDir + "/" + sample + "_mystery_reads_2.fq");
		}
		
		// Align each sample to transcriptome
		ArrayList<Job> transcriptomeJobs = new ArrayList<Job>();
		for(String sample : this.sampleNames) {
			File outSamFile = new File(samToTranscriptome.get(sample));
			if(!this.readsPaired) {
				File outUnalignedFile = new File(fastqNotTranscriptomeUnpaired.get(sample));
				// If files already exist do not rerun
				if(outSamFile.exists() && outUnalignedFile.exists()) {
					this.logger.info("WARNING: sam file and fastq file for sample " + sample + " already exist. Not rerunning alignment to transcriptome.");
					continue;
				}
				Job job = AlignmentUtils.runBowtie2(btBase, bowtie2options, this.read1FastqFiles.get(sample), samToTranscriptome.get(sample), fastqNotTranscriptomeUnpaired.get(sample), bowtie2Executable, logDir, scheduler);
				transcriptomeJobs.add(job);
			} else {
				File outUnalignedFile1 = new File(fastqNotTranscriptomePaired1.get(sample));
				File outUnalignedFile2 = new File(fastqNotTranscriptomePaired2.get(sample));
				// If files already exist do not rerun
				if(outSamFile.exists() && outUnalignedFile1.exists() && outUnalignedFile2.exists()) {
					this.logger.info("WARNING: sam file and fastq files for sample " + sample + " already exist. Not rerunning alignment to transcriptome.");
					continue;
				}
				Job job = AlignmentUtils.runBowtie2(btBase, bowtie2options, this.read1FastqFiles.get(sample), this.read2FastqFiles.get(sample), samToTranscriptome.get(sample), fastqNotTranscriptomePairedBase.get(sample), bowtie2Executable, logDir, scheduler);
				transcriptomeJobs.add(job);
			}		
		}
		// Wait for bowtie2 jobs to finish
		this.logger.info("Waiting for jobs to finish...");
		JobUtils.waitForAll(transcriptomeJobs);
		this.logger.info("Done aligning to transcriptome.");
				
		// Align unmapped reads for each sample to genome
		this.logger.info("");
		this.logger.info("Aligning unmapped reads to genome...");
		ArrayList<Job> genomeJobs = new ArrayList<Job>();
		for(String sample : this.sampleNames) {
			File outSamFile = new File(samToGenome.get(sample));
			if(!this.readsPaired) {
				File outUnalignedFile = new File(fastqNotGenomeUnpaired.get(sample));
				// If files already exist do not rerun
				if(outSamFile.exists() && outUnalignedFile.exists()) {
					this.logger.info("WARNING: sam file and fastq file for sample " + sample + " already exist. Not rerunning alignment to genome.");
					continue;
				}
				Job job = AlignmentUtils.runBowtie2(this.genomeBowtieIndex, bowtie2options, fastqNotTranscriptomeUnpaired.get(sample), samToGenome.get(sample), fastqNotGenomeUnpaired.get(sample), bowtie2Executable, logDir, scheduler);
				genomeJobs.add(job);
			} else {
				File outUnalignedFile1 = new File(fastqNotGenomePaired1.get(sample));
				File outUnalignedFile2 = new File(fastqNotGenomePaired2.get(sample));
				// If files already exist do not rerun
				if(outSamFile.exists() && outUnalignedFile1.exists() && outUnalignedFile2.exists()) {
					this.logger.info("WARNING: sam file and fastq files for sample " + sample + " already exist. Not rerunning alignment to genome.");
					continue;
				}
				Job job = AlignmentUtils.runBowtie2(this.genomeBowtieIndex, bowtie2options, fastqNotTranscriptomePaired1.get(sample), fastqNotTranscriptomePaired2.get(sample), samToGenome.get(sample), fastqNotGenomePairedBase.get(sample), bowtie2Executable, logDir, scheduler);
				genomeJobs.add(job);
			}		
		}
		
		// Wait for bowtie2 jobs to finish
		this.logger.info("Waiting for jobs to finish...");
		JobUtils.waitForAll(genomeJobs);
		this.logger.info("Done aligning unmapped reads to genome.");
		
		// Count classes for each sample
		// Count genome alignments for each sample
		// Calculate number of unaligned reads
		this.logger.info("");
		this.logger.info("Counting mappings to RNA classes and genome...");
		Map<String, Map<String, Integer>> rtrn = new TreeMap<String, Map<String, Integer>>();
		for(String sample : this.sampleNames) {
			Map<String, Integer> classCounts = countClasses(samToTranscriptome.get(sample));
			int transcriptomeCount = 0;
			for(String className : classCounts.keySet()) transcriptomeCount += classCounts.get(className).intValue();
			int genomeCount = AlignmentUtils.countAlignments(samtoolsExecutable, samToGenome.get(sample), logDir, true);
			int mysteryCount = Math.max(0, this.totalReadCounts.get(sample).intValue() - transcriptomeCount - genomeCount);
			classCounts.put(GENOME_CLASS_NAME, Integer.valueOf(genomeCount));
			classCounts.put(MYSTERY_READS_CLASS_NAME, Integer.valueOf(mysteryCount));
			rtrn.put(sample, classCounts);
		}
		this.logger.info("");
		this.logger.info("Done counting classes.");
		
		return rtrn;
	
	}
	

}

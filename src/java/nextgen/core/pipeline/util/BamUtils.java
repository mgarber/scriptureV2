package nextgen.core.pipeline.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import nextgen.core.annotation.Gene;
import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.LSFJob;
import nextgen.core.pipeline.Scheduler;
import nextgen.core.readFilters.FirstOfPairFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.readFilters.SecondOfPairFilter;
import nextgen.core.writers.WigWriter;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class BamUtils {
	
	private static Logger logger = Logger.getLogger(BamUtils.class.getName());
	
	/**
	 * Write global genomic space stats for bam file and write bsub output to working directory
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param chrSizeFile Chromosome size file for genomic space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static String writeGenomicSpaceStats(String bamFile, String sampleName, String chrSizeFile, String alignmentGlobalStatsJar) throws IOException, InterruptedException {
		return writeGenomicSpaceStats(bamFile, chrSizeFile, alignmentGlobalStatsJar, ".");
	}
	
	/**
	 * Write global genomic space stats for bam files and write bsub output to working directory
	 * @param bamFiles Bam files by sample name
	 * @param chrSizeFile Chromosome size file for genomic space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param scheduler Scheduler
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static Collection<Job> writeGenomicSpaceStats(Map<String, String> bamFiles, String chrSizeFile, String alignmentGlobalStatsJar, Scheduler scheduler) throws IOException, InterruptedException {
		return writeGenomicSpaceStats(bamFiles, chrSizeFile, alignmentGlobalStatsJar, ".", scheduler);
	}
	
	/**
	 * Write global genomic space stats for bam files
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param chrSizeFile Chromosome size file for genomic space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutputDir Directory to write LSF output files to
	 * @param scheduler Scheduler
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static Job writeGenomicSpaceStats(String bamFile, String sampleName, String chrSizeFile, String alignmentGlobalStatsJar, String bsubOutputDir, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		return writeGenomicSpaceStats(bamFiles, chrSizeFile, alignmentGlobalStatsJar, bsubOutputDir, scheduler).iterator().next();
	}

	
	/**
	 * Write global genomic space stats for bam files
	 * @param bamFiles Bam files by sample name
	 * @param chrSizeFile Chromosome size file for genomic space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutputDir Directory to write LSF output files to
	 * @param scheduler Scheduler
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static Collection<Job> writeGenomicSpaceStats(Map<String, String> bamFiles, String chrSizeFile, String alignmentGlobalStatsJar, String bsubOutputDir, Scheduler scheduler) throws IOException, InterruptedException {
		logger.info("Calculating genomic space stats from chromosome sizes in " + chrSizeFile);
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String sampleName : bamFiles.keySet()) {
			String bamFile = bamFiles.get(sampleName);
			logger.info("Calculating genomic space stats for sample " + sampleName + "...");
			String cmmd = "java -jar -Xmx30g -Xms20g -Xmn10g " + alignmentGlobalStatsJar + " -b " + bamFile + " -g " + chrSizeFile;
			logger.info("Running command: " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutputDir + "/compute_genomic_space_stats_" + jobID + ".bsub", "week", 32);		
				jobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		return jobs;
	}
	
	/**
	 * Write global transcriptome space stats for bam files and write bsub output in the working directory
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param bedFile Bed annotation for transcriptome space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static String writeTranscriptomeSpaceStats(String bamFile, String sampleName, String bedFile, String alignmentGlobalStatsJar) throws IOException, InterruptedException {
		return writeTranscriptomeSpaceStats(bamFile, bedFile, alignmentGlobalStatsJar, ".");
	}

	
	/**
	 * Write global transcriptome space stats for bam files and write bsub output in the working directory
	 * @param bamFiles Bam files by sample name
	 * @param bedFile Bed annotation for transcriptome space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param scheduler Scheduler
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static Collection<Job> writeTranscriptomeSpaceStats(Map<String, String> bamFiles, String bedFile, String alignmentGlobalStatsJar, Scheduler scheduler) throws IOException, InterruptedException {
		return writeTranscriptomeSpaceStats(bamFiles, bedFile, alignmentGlobalStatsJar, ".", scheduler);
	}
	
	/**
	 * Write global transcriptome space stats for bam files
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param bedFile Bed annotation for transcriptome space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutDir Directory to write bsub output to
	 * @param scheduler Scheduler
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static Job writeTranscriptomeSpaceStats(String bamFile, String sampleName, String bedFile, String alignmentGlobalStatsJar, String bsubOutDir, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		return writeTranscriptomeSpaceStats(bamFiles, bedFile, alignmentGlobalStatsJar, bsubOutDir, scheduler).iterator().next();
	}
	
	/**
	 * Write global transcriptome space stats for bam files
	 * @param bamFiles Bam files by sample name
	 * @param bedFile Bed annotation for transcriptome space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutDir Directory to write bsub output to
	 * @param scheduler Scheduler
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public static Collection<Job> writeTranscriptomeSpaceStats(Map<String, String> bamFiles, String bedFile, String alignmentGlobalStatsJar, String bsubOutDir, Scheduler scheduler) throws IOException, InterruptedException {
		logger.info("Calculating transcriptome space stats from annotation in " + bedFile);
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String sampleName : bamFiles.keySet()) {
			String bamFile = bamFiles.get(sampleName);
			logger.info("Calculating transcriptome space stats for sample " + sampleName + "...");
			String cmmd = "java -jar -Xmx30g -Xms20g -Xmn10g " + alignmentGlobalStatsJar + " -b " + bamFile + " -t " + bedFile;
			logger.info("Running command: " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/compute_transcriptome_space_stats_" + jobID + ".bsub", "week", 32);	
				jobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		return jobs;
	}

	/**
	 * Make tdf files for bam file and write bsub output to working directory
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void makeTdfs(String bamFile, String sampleName, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler) throws IOException, InterruptedException {
		makeTdfs(bamFile, sampleName, ".", assemblyFasta, igvtoolsExecutable, scheduler);
	}
	
	/**
	 * Make tdf files for bam files and write bsub output to working directory
	 * @param bamFilesBySampleName Bam file name by sample name
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void makeTdfs(Map<String, String> bamFilesBySampleName, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler) throws IOException, InterruptedException {
		makeTdfs(bamFilesBySampleName, ".", assemblyFasta, igvtoolsExecutable, scheduler);
	}
	
	/**
	 * Make tdf files for bam file
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param bsubOutDir Directory to write bsub output to
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void makeTdfs(String bamFile, String sampleName, String bsubOutDir, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		makeTdfs(bamFiles, bsubOutDir, assemblyFasta, igvtoolsExecutable, scheduler);
	}
	
	/**
	 * Make tdf files for bam files
	 * @param bamFilesBySampleName Bam file name by sample name
	 * @param bsubOutDir Directory to write bsub output to
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void makeTdfs(Map<String, String> bamFilesBySampleName, String bsubOutDir, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler) throws IOException, InterruptedException {
		ArrayList<Job> tdfJobs = new ArrayList<Job>();
		Map<String, String> outTdf = new TreeMap<String, String>();
		for(String sample : bamFilesBySampleName.keySet()) {
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
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/make_tdf_" + jobID + ".bsub", "hour", 1);
				job.submit();
				tdfJobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		if(tdfJobs.isEmpty()) return;
		logger.info("Waiting for igvtools jobs to finish...");
		logger.info("Note: igvtools count always exits with code -1 even though it worked, so disregard failure notifications from LSF.");
		// Igvtools count always ends by crashing even though it worked, so catch the exception and check if files were really created
		try {
			JobUtils.waitForAll(tdfJobs);
		} catch(IllegalArgumentException e) {
			boolean ok = true;
			String errMsg = "";
			for(String sample : bamFilesBySampleName.keySet()) {
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
	 * Convert sam file to bam file
	 * @param sampleName Sample name
	 * @param samFile Original sam file
	 * @param bamFile Bam file to write
	 * @param finalBamFile Final bam file; skip if already exists
	 * @param samtools Samtools executables
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void samToBam(String sampleName, String samFile, String bamFile, String finalBamFile, String samtools, Scheduler scheduler) throws IOException, InterruptedException {
		samToBam(sampleName, samFile, bamFile, finalBamFile, ".", samtools, scheduler);
	}
	
	/**
	 * Convert sam files to bam files
	 * @param samFiles The sam files by sample name
	 * @param bamFiles The bam files to write, by sample name
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param samtools Samtools executables
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void samToBam(Map<String, String> samFiles, Map<String, String> bamFiles, Map<String, String> bsubOutputDirs, String samtools, Scheduler scheduler) throws IOException, InterruptedException {
		samToBam(samFiles, bamFiles, null, bsubOutputDirs, samtools, scheduler);
	}
	
	/**
	 * Convert sam file to bam file
	 * @param sampleName Sample name
	 * @param samFile Original sam file
	 * @param bamFile Bam file to write
	 * @param finalBamFile Final bam file; skip if already exists
	 * @param bsubOutputDir Directory to write bsub output to
	 * @param samtools Samtools executables
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void samToBam(String sampleName, String samFile, String bamFile, String finalBamFile, String bsubOutputDir, String samtools, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		Map<String, String> samFiles = new TreeMap<String, String>();
		samFiles.put(sampleName, samFile);
		Map<String, String> finalBamFiles = new TreeMap<String, String>();
		finalBamFiles.put(sampleName, finalBamFile);
		Map<String, String> bsubOutputDirs = new TreeMap<String, String>();
		bsubOutputDirs.put(sampleName, bsubOutputDir);
		samToBam(samFiles, bamFiles, finalBamFiles, bsubOutputDirs, samtools, scheduler);
	}
	
	
	/**
	 * Convert sam files to bam files
	 * @param samFiles The sam files by sample name
	 * @param bamFiles The bam files to write, by sample name
	 * @param finalBamFiles Final bam files; skip if they already exist
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param samtools Samtools executables
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void samToBam(Map<String, String> samFiles, Map<String, String> bamFiles, Map<String, String> finalBamFiles, Map<String, String> bsubOutputDirs, String samtools, Scheduler scheduler) throws IOException, InterruptedException {
		ArrayList<Job> cbJobs = new ArrayList<Job>();
		for(String sample : samFiles.keySet()) {
			File bam = new File(bamFiles.get(sample));
			// Check if bam files already exist
			if(bam.exists()) {
				logger.warn("Bam file " + bam + " already exists. Not redoing format conversion from sam.");
				continue;
			}
			if(finalBamFiles != null) {
				if(finalBamFiles.containsKey(sample)) {
					File finalBam = new File(finalBamFiles.get(sample));
					if(finalBam.exists()) {
						logger.warn("Alignment file " + finalBam + " already exists. Not looking for or converting sam file.");
						continue;					
					}
				}
			}
			// Use samtools view
			String cmmd = samtools + " view -Sb -o " + bamFiles.get(sample) + " " + samFiles.get(sample);
			logger.info("Running Samtools command: " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutputDirs.get(sample) + "/sam_to_bam_" + jobID + ".bsub", "hour", 1);
				cbJobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		// Wait for jobs to finish
		logger.info("Waiting for samtools view jobs to finish...");
		JobUtils.waitForAll(cbJobs);
	}

	/**
	 * Sort bam file and write sorted file to specified location and write bsub output to working directory
	 * @param sampleName Sample name
	 * @param unsortedBam Unsorted bam file
	 * @param sortedBam Sorted bam file to write
	 * @param finalBam Final bam file; skip if already exists
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void sortBamFile(String sampleName, String unsortedBam, String sortedBam, String finalBam, String picardJarDir, Scheduler scheduler) throws IOException, InterruptedException {
		sortBamFile(sampleName, unsortedBam, sortedBam, ".", finalBam, picardJarDir, scheduler);
	}
	
	/**
	 * Sort all bam files and write sorted files to specified locations, and write bsub output to working directory
	 * @param unsortedBams The unsorted bam files to sort, by sample name
	 * @param sortedBams The sorted bam files to write, by sample name
	 * @param finalBams Final bam files; skip if they already exist
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void sortBamFiles(Map<String, String> unsortedBams, Map<String, String> sortedBams, Map<String, String> finalBams, String picardJarDir, Scheduler scheduler) throws IOException, InterruptedException {
		sortBamFiles(unsortedBams, sortedBams, ".", finalBams, picardJarDir, scheduler);
	}
	
	/**
	 * Sort all bam files and write sorted files to specified locations
	 * @param unsortedBams The unsorted bam files to sort, by sample name
	 * @param sortedBams The sorted bam files to write, by sample name
	 * @param bsubOutputDir Directory to write all bsub output files to
	 * @param finalBams Final bam files; skip if they already exist
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void sortBamFiles(Map<String, String> unsortedBams, Map<String, String> sortedBams, String bsubOutputDir, Map<String, String> finalBams, String picardJarDir, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> outDirs = new TreeMap<String, String>();
		for(String sample : unsortedBams.keySet()) {
			outDirs.put(sample, bsubOutputDir);
		}
		sortBamFiles(unsortedBams, sortedBams, outDirs, finalBams, picardJarDir, scheduler);
	}
	
	/**
	 * Sort bam file and write sorted file to specified location
	 * @param sampleName Sample name
	 * @param unsortedBam Unsorted bam file
	 * @param sortedBam Sorted bam file to write
	 * @param bsubOutputDir Directory to write bsub output to
	 * @param finalBam Final bam file; skip if already exists
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void sortBamFile(String sampleName, String unsortedBam, String sortedBam, String bsubOutputDir, String finalBam, String picardJarDir, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> unsortedBams = new TreeMap<String, String>();
		unsortedBams.put(sampleName, unsortedBam);
		Map<String, String> sortedBams = new TreeMap<String, String>();
		sortedBams.put(sampleName, sortedBam);
		Map<String, String> finalBamFiles = new TreeMap<String, String>();
		finalBamFiles.put(sampleName, finalBam);
		Map<String, String> bsubOutputDirs = new TreeMap<String, String>();
		bsubOutputDirs.put(sampleName, bsubOutputDir);
		sortBamFiles(unsortedBams, sortedBams, bsubOutputDirs, finalBamFiles, picardJarDir, scheduler);
	}
	
	
	/**
	 * Sort all bam files and write sorted files to specified locations
	 * @param unsortedBams The unsorted bam files to sort, by sample name
	 * @param sortedBams The sorted bam files to write, by sample name
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param finalBams Final bam files; skip if they already exist
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void sortBamFiles(Map<String, String> unsortedBams, Map<String, String> sortedBams, Map<String, String> bsubOutputDirs, Map<String, String> finalBams, String picardJarDir, Scheduler scheduler) throws IOException, InterruptedException {
		ArrayList<Job> sbJobs = new ArrayList<Job>();
		for(String sample : unsortedBams.keySet()) {
			File finalBam = new File(finalBams.get(sample));
			if(finalBam.exists()) {
				logger.warn("Alignment file " + finalBam + " already exists. Not sorting.");
				continue;					
			}
			// Use Picard program SortSam
			String cmmd = "java -jar " + picardJarDir + "/SortSam.jar INPUT=" + unsortedBams.get(sample) + " OUTPUT=" + sortedBams.get(sample) + " SORT_ORDER=coordinate";
			logger.info("Running Picard command: " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutputDirs.get(sample) + "/sort_bam_" + jobID + ".bsub", "hour", 4);
				job.submit();
				sbJobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		// Wait for jobs to finish
		logger.info("Waiting for SortSam jobs to finish...");
		JobUtils.waitForAll(sbJobs);
	}

	/**
	 * Write wig file of raw position count and position count normalized by average coverage over gene
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param bamDir Bam directory
	 * @param geneBedFile Bed file of genes to use
	 * @param refFasta Reference fasta file
	 * @param wigToBigWigExecutable WigToBigWig executable file
	 * @param wigWriterJar WigWriter jar file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void writeWigPositionCount(String bamFile, String sampleName, String bamDir, String geneBedFile, String refFasta, String wigToBigWigExecutable, String wigWriterJar, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		writeWigPositionCount(bamFiles, bamDir, geneBedFile, refFasta, wigToBigWigExecutable, wigWriterJar, scheduler);
	}
	
	/**
	 * Write wig files of raw position count and position count normalized by average coverage over gene
	 * @param bamFiles Bam files by sample name
	 * @param bamDir Bam directory
	 * @param geneBedFile Bed file of genes to use
	 * @param refFasta Reference fasta file
	 * @param wigToBigWigExecutable WigToBigWig executable file
	 * @param wigWriterJar WigWriter jar file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void writeWigPositionCount(Map<String, String> bamFiles, String bamDir, String geneBedFile, String refFasta, String wigToBigWigExecutable, String wigWriterJar, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneBedFile));
		Collection<String> chrNames = genes.keySet();
		Map<String, Map<String, String>> normalizedWigFiles = new TreeMap<String, Map<String, String>>();
		Map<String, String> fullNormalizedWigFiles = new TreeMap<String, String>();
		Map<String, String> fullNormalizedBigwigFiles = new TreeMap<String, String>();
		Map<String, Map<String, String>> unnormalizedWigFiles = new TreeMap<String, Map<String, String>>();
		Map<String, String> fullUnnormalizedWigFiles = new TreeMap<String, String>();
		Map<String, String> fullUnnormalizedBigwigFiles = new TreeMap<String, String>();
		ArrayList<Job> wigJobs = new ArrayList<Job>();
		for(String sample : bamFiles.keySet()) {
			String normalizedWigFile = bamDir + "/" + sample + ".normalized.wig";
			String normalizedBigwigFile = bamDir + "/" + sample + ".normalized.bw";
			fullNormalizedWigFiles.put(sample, normalizedWigFile);
			fullNormalizedBigwigFiles.put(sample, normalizedBigwigFile);
			File normalizedBw = new File(fullNormalizedBigwigFiles.get(sample));
			String unnormalizedWigFile = bamDir + "/" + sample + ".wig";
			String unnormalizedBigwigFile = bamDir + "/" + sample + ".bw";
			fullUnnormalizedWigFiles.put(sample, unnormalizedWigFile);
			fullUnnormalizedBigwigFiles.put(sample, unnormalizedBigwigFile);
			File unnormalizedBw = new File(fullUnnormalizedBigwigFiles.get(sample));
			if(normalizedBw.exists() && unnormalizedBw.exists()) {
				logger.warn("Bigwig files " + normalizedBw + " and " + unnormalizedBw + " already exist. Not remaking wigs or bigwigs.");
				continue;
			}
			File wfn = new File(fullNormalizedWigFiles.get(sample));
			File wfu = new File(fullUnnormalizedWigFiles.get(sample));
			if(wfn.exists() && wfu.exists()) {
				logger.warn("Wig files " + wfn + " and " + wfu + " already exist. Not remaking files.");
				continue;
			}
			String bamFile = bamFiles.get(sample);
			
			// Create normalized wig files
			Map<String, String> normalizedWigFilesByChr = new TreeMap<String, String>();
			for(String chr : chrNames) {
				String prefix = bamDir + "/" + sample + "_" + chr + ".normalized";
				String normalizedFile = prefix + ".wig";
				normalizedWigFilesByChr.put(chr, normalizedFile);
				File nf = new File(normalizedFile);
				if(nf.exists()) {
					logger.warn("Wig file " + normalizedFile + " already exists. Not remaking file.");
				} else {
					logger.info("Writing wig file " + normalizedFile + "...");
					String cmmd = "java -jar -Xmx15g -Xms10g -Xmn5g " + wigWriterJar + " -b " + bamFile + " -g " + geneBedFile + " -n true -o " + prefix + " -chr " + chr;
					logger.info("Running command: " + cmmd);
					switch(scheduler) {
					case LSF:
						String jobID = Long.valueOf(System.currentTimeMillis()).toString();
						logger.info("LSF job ID is " + jobID + ".");
						// Submit job
						LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/write_wig_normalized_" + sample + "_" + chr + "_" + jobID + ".bsub", "week", 16);
						job.submit();
						wigJobs.add(job);
						break;
					default:
						throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
					}
				}
			}
			normalizedWigFiles.put(sample, normalizedWigFilesByChr);
			
			// Create unnormalized wig files
			Map<String, String> unnormalizedWigFilesByChr = new TreeMap<String, String>();
			for(String chr : chrNames) {
				String prefix = bamDir + "/" + sample + "_" + chr;
				String unnormalizedFile = prefix + ".wig";
				unnormalizedWigFilesByChr.put(chr, unnormalizedFile);
				File uf = new File(unnormalizedFile);
				if(uf.exists()) {
					logger.warn("Wig file " + unnormalizedFile + " already exists. Not remaking file.");
				} else {
					logger.info("Writing wig file " + unnormalizedFile + "...");
					String cmmd = "java -jar -Xmx15g -Xms10g -Xmn5g " + wigWriterJar + " -b " + bamFile + " -g " + geneBedFile + " -n false -o " + prefix + " -chr " + chr;
					logger.info("Running command: " + cmmd);
					switch(scheduler) {
					case LSF:
						String jobID = Long.valueOf(System.currentTimeMillis()).toString();
						logger.info("LSF job ID is " + jobID + ".");
						// Submit job
						LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/write_wig_unnormalized_" + sample + "_" + chr + "_" + jobID + ".bsub", "week", 16);
						job.submit();
						wigJobs.add(job);
						break;
					default:
						throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
					}
				}
			}
			unnormalizedWigFiles.put(sample, unnormalizedWigFilesByChr);
		}
		logger.info("");
		logger.info("Waiting for wig writer jobs to finish...");
		JobUtils.waitForAll(wigJobs);
		logger.info("");
		logger.info("Combining normalized chromosome wig files...");
		for(String sample : normalizedWigFiles.keySet()) {
			logger.info(sample);
			String wigFile = fullNormalizedWigFiles.get(sample);
			FileWriter w = new FileWriter(wigFile);
			for(String chrFile : normalizedWigFiles.get(sample).values()) {
				FileReader r = new FileReader(chrFile);
				BufferedReader b = new BufferedReader(r);
				while(b.ready()) {
					w.write(b.readLine() + "\n");
				}
				r.close();
				b.close();
			}
			w.close();
		}
		logger.info("Done combining normalized chromosome wig files. Delete individual chromosome files to save storage.");
		logger.info("");
		logger.info("Combining unnormalized chromosome wig files...");
		for(String sample : unnormalizedWigFiles.keySet()) {
			logger.info(sample);
			String wigFile = fullUnnormalizedWigFiles.get(sample);
			FileWriter w = new FileWriter(wigFile);
			for(String chrFile : unnormalizedWigFiles.get(sample).values()) {
				FileReader r = new FileReader(chrFile);
				BufferedReader b = new BufferedReader(r);
				while(b.ready()) {
					w.write(b.readLine() + "\n");
				}
				r.close();
				b.close();
			}
			w.close();
		}
		logger.info("Done combining unnormalized chromosome wig files. Delete individual chromosome files to save storage.");
		
		logger.info("");
		logger.info("Making normalized and unnormalized bigwig files...");
		String chrSizeFile = FastaUtils.writeSizeFile(refFasta);
		ArrayList<Job> bigwigJobs = new ArrayList<Job>();
		// Submit jobs for normalized files
		for(String sample : fullNormalizedBigwigFiles.keySet()) {
			String wig = fullNormalizedWigFiles.get(sample);
			String bigwig = fullNormalizedBigwigFiles.get(sample);
			File b = new File(bigwig);
			if(b.exists()) {
				logger.warn("Bigwig file " +  bigwig + " already exists. Not remaking file.");
				continue;
			}
			String cmmd = wigToBigWigExecutable + " " + wig + " " + chrSizeFile + " " + bigwig;
			logger.info("");
			logger.info("Making bigwig file for wig file " + wig + ".");
			logger.info("Running UCSC command " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_normalized_" + jobID + ".bsub", "hour", 4);
				job.submit();
				bigwigJobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
			}
		}
		// Submit jobs for unnormalized files
		for(String sample : fullUnnormalizedBigwigFiles.keySet()) {
			String wig = fullUnnormalizedWigFiles.get(sample);
			String bigwig = fullUnnormalizedBigwigFiles.get(sample);
			File b = new File(bigwig);
			if(b.exists()) {
				logger.warn("Bigwig file " +  bigwig + " already exists. Not remaking file.");
				continue;
			}
			String cmmd = wigToBigWigExecutable + " " + wig + " " + chrSizeFile + " " + bigwig;
			logger.info("");
			logger.info("Making bigwig file for wig file " + wig + ".");
			logger.info("Running UCSC command " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_unnormalized_" + jobID + ".bsub", "hour", 4);
				job.submit();
				bigwigJobs.add(job);
				break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
			}
		}
		logger.info("Waiting for wigToBigWig jobs to finish...");
		JobUtils.waitForAll(bigwigJobs);
	}
	
	/**
	 * Write fragment end points and midpoints to wig and bigwig files
	 * @param sampleName Sample name
	 * @param bamFile Bam file
	 * @param pairedData Whether the sample has paired end reads
	 * @param bamDir Directory containing bam files
	 * @param refFasta Fasta file of sequences these bam files were aligned against
	 * @param geneBedFile Bed file of genes to count reads in or null if using genomic space
	 * @param wigToBigWigExecutable WigToBigWig executable file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void writeWigFragmentEndsAndMidpoints(String sampleName, String bamFile, boolean pairedData, String bamDir, String refFasta, String geneBedFile, String wigToBigWigExecutable, Scheduler scheduler) throws IOException, InterruptedException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		Map<String, Boolean> paired = new TreeMap<String, Boolean>();
		paired.put(sampleName, Boolean.valueOf(pairedData));
		writeWigFragmentEndsAndMidpoints(bamFiles, paired, bamDir, refFasta, geneBedFile, wigToBigWigExecutable, scheduler);
	}
	
	
	/**
	 * Write fragment end points and midpoints to wig and bigwig files
	 * @param bamFiles Bam files by sample name
	 * @param pairedData Whether each sample has paired end reads, by sample name
	 * @param bamDir Directory containing bam files
	 * @param refFasta Fasta file of sequences these bam files were aligned against
	 * @param geneBedFile Bed file of genes to count reads in or null if using genomic space
	 * @param wigToBigWigExecutable WigToBigWig executable file
	 * @param scheduler Scheduler
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void writeWigFragmentEndsAndMidpoints(Map<String, String> bamFiles, Map<String, Boolean> pairedData, String bamDir, String refFasta, String geneBedFile, String wigToBigWigExecutable, Scheduler scheduler) throws IOException, InterruptedException {
		
		Map<String, String> read1endWig = new TreeMap<String, String>();
		Map<String, String> read2endWig = new TreeMap<String, String>();
		Map<String, String> read1endBigwig = new TreeMap<String, String>();
		Map<String, String> read2endBigwig = new TreeMap<String, String>();
		Map<String, String> midpointWig = new TreeMap<String, String>();
		Map<String, String> midpointBigwig = new TreeMap<String, String>();
		ArrayList<Job> bigwigJobs = new ArrayList<Job>();
		
		// Chromosome size file to pass to wig writer if using genomic space
		// Null if using transcriptome space
		String chrSizesForWigToBigWig = FastaUtils.writeSizeFile(refFasta);
		String chrSizesForWigWriter = null;
		if(geneBedFile == null) {
			chrSizesForWigWriter = chrSizesForWigToBigWig;
		}

		for(String sampleName : bamFiles.keySet()) {
			
			String bamFile = bamFiles.get(sampleName);
			read1endWig.put(sampleName, bamFile + ".read1.wig");
			read2endWig.put(sampleName, bamFile + ".read2.wig");
			read1endBigwig.put(sampleName, bamFile + ".read1.bw");
			read2endBigwig.put(sampleName, bamFile + ".read2.bw");
			midpointWig.put(sampleName, bamFile + ".midpoint.wig");
			midpointBigwig.put(sampleName, bamFile + ".midpoint.bw");
			
			// Write wig file for read1
			String wig1 = read1endWig.get(sampleName);
			File read1wigFile = new File(wig1);
			String bigwig1 = read1endBigwig.get(sampleName);
			File read1bigwigFile = new File(bigwig1);
			if(read1wigFile.exists() || read1bigwigFile.exists()) {
				logger.warn("Read 1 wig file or bigwig file for sample " + sampleName + " already exists. Not remaking wig file.");
			} else {
				logger.info("Writing fragment ends of read 1 from bam file " + bamFile + " to wig file " + wig1 + ".");
				WigWriter read1ww = new WigWriter(bamFile, geneBedFile, chrSizesForWigWriter, WigWriter.BEGINNING_POSITION_DESCRIPTION, false, false);
				read1ww.addReadFilter(new FirstOfPairFilter());
				read1ww.addReadFilter(new ProperPairFilter());
				String prefix1 = wig1.replaceAll(".wig", "");
				read1ww.writeFullWig(prefix1);
				logger.info("Done writing file " + wig1 + ".");
			}
			// Write bigwig file for read1
			if(read1bigwigFile.exists()) {
				logger.warn("Bigwig file " + bigwig1 + " already exists. Not remaking file.");
			} else {
				String cmmd = wigToBigWigExecutable + " " + wig1 + " " + chrSizesForWigToBigWig + " " + bigwig1;
				logger.info("");
				logger.info("Making bigwig file for wig file " + wig1 + ".");
				logger.info("Running UCSC command " + cmmd);
				switch(scheduler) {
				case LSF:
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					logger.info("LSF job ID is " + jobID + ".");
					LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_" + jobID + ".bsub", "hour", 4);
					bigwigJobs.add(job);
					break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
				}
			}
			
			// Write midpoint wig file
			String midpointWigFileName = midpointWig.get(sampleName);
			File midpointWigFile = new File(midpointWigFileName);
			String midpointBigwigFileName = midpointBigwig.get(sampleName);
			File midpointBigwigFile = new File(midpointBigwigFileName);
			if(midpointWigFile.exists() || midpointBigwigFile.exists()) {
				logger.warn("Fragment midpoint wig file or bigwig file for sample " + sampleName + " already exists. Not remaking wig file.");
			} else {
				logger.info("Writing fragment midpoints from bam file " + bamFile + " to wig file " + midpointWigFile + ".");
				WigWriter ww = new WigWriter(bamFile, geneBedFile, chrSizesForWigWriter, WigWriter.MIDPOINT_POSITION_DESCRIPTION, true, false);
				ww.addReadFilter(new ProperPairFilter());
				String prefix = midpointWigFileName.replaceAll(".wig", "");
				ww.writeFullWig(prefix);
				logger.info("Done writing file " + midpointWigFileName + ".");
			}
			// Write midpoint bigwig file
			if(midpointBigwigFile.exists()) {
				logger.warn("Bigwig file " + midpointBigwigFileName + " already exists. Not remaking file.");
			} else {
				String cmmd = wigToBigWigExecutable + " " + midpointWigFileName + " " + chrSizesForWigToBigWig + " " + midpointBigwigFileName;
				logger.info("");
				logger.info("Making bigwig file for wig file " + midpointWigFileName + ".");
				logger.info("Running UCSC command " + cmmd);
				switch(scheduler) {
				case LSF:
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					logger.info("LSF job ID is " + jobID + ".");
					LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_" + jobID + ".bsub", "hour", 4);
					job.submit();
					bigwigJobs.add(job);
					break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
				}
			}

			
			if(pairedData.get(sampleName).booleanValue()) {
				// Write wig file for read2
				String wig2 = read2endWig.get(sampleName);
				File read2wigFile = new File(wig2);
				String bigwig2 = read2endBigwig.get(sampleName);
				File read2bigwigFile = new File(bigwig2);
				if(read2wigFile.exists() || read2bigwigFile.exists()) {
					logger.warn("Read 2 wig file or bigwig file for sample " + sampleName + " already exists. Not remaking wig file.");
				} else {
					logger.info("Writing fragment ends of read 2 from bam file " + bamFile + " to wig file " + wig2 + ".");
					WigWriter read2ww = new WigWriter(bamFile, geneBedFile, chrSizesForWigWriter, WigWriter.BEGINNING_POSITION_DESCRIPTION, false, false);
					read2ww.addReadFilter(new SecondOfPairFilter());
					read2ww.addReadFilter(new ProperPairFilter());
					String prefix2 = wig2.replaceAll(".wig", "");
					read2ww.writeFullWig(prefix2);
					logger.info("Done writing file " + wig2 + ".");
				}
				// Write bigwig file for read2
				if(read2bigwigFile.exists()) {
					logger.warn("Bigwig file " + bigwig2 + " already exists. Not remaking file.");
				} else {
					String cmmd = wigToBigWigExecutable + " " + wig2 + " " + chrSizesForWigToBigWig + " " + bigwig2;
					logger.info("");
					logger.info("Making bigwig file for wig file " + wig2 + ".");
					logger.info("Running UCSC command " + cmmd);
					switch(scheduler) {
					case LSF:
						String jobID = Long.valueOf(System.currentTimeMillis()).toString();
						logger.info("LSF job ID is " + jobID + ".");
						LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_" + jobID + ".bsub", "hour", 4);
						job.submit();
						bigwigJobs.add(job);
						break;
					default:
						throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
					}
				}
			}
		
			
		}

		logger.info("");
		logger.info("Waiting for wigToBigWig jobs to finish...");
		JobUtils.waitForAll(bigwigJobs);
		
	}

	
}

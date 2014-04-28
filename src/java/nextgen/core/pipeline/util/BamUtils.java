package nextgen.core.pipeline.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.picard.sam.SortSam;
import net.sf.samtools.SAMFileReader;
import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.LSFJob;
import nextgen.core.job.OGSJob;
import nextgen.core.pipeline.Scheduler;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;


/**
 * @author prussell
 *
 */
public class BamUtils {
	
	private static Logger logger = Logger.getLogger(BamUtils.class.getName());
	
	/**
	 * Merge bam files
	 * @param inputBams Names of bam files to merge
	 * @param outputMergedBam Output merged bam file to write
	 * @param scheduler Scheduler
	 * @param drmaaSession DRMAA session if using OGS or null otherwise
	 * @param picardJarDir Directory containing Picard jar files
	 * @param assumeSorted Assume input files are sorted
	 * @throws DrmaaException
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void mergeBams(Collection<String> inputBams, String outputMergedBam, Scheduler scheduler, Session drmaaSession, String picardJarDir, boolean assumeSorted) throws DrmaaException, IOException, InterruptedException {
		logger.info("Creating merged bam file " + outputMergedBam + ".");
		String inputs = "";
		for(String bam : inputBams) {
			inputs += "INPUT=" + bam + " ";
		}
		String output = "OUTPUT=" + outputMergedBam;
		String cmmd = "java -jar " + picardJarDir + "/MergeSamFiles.jar " + inputs + " " + output + " ASSUME_SORTED=" + assumeSorted + " MERGE_SEQUENCE_DICTIONARIES=true";
		logger.info("Running picard command: " + cmmd);
		switch(scheduler) {
		case LSF:
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, "merge_bam_files_" + jobID + ".bsub", "week", 8);	
			job.submit();
			job.waitFor();
			break;
        case OGS:
            if(drmaaSession == null) {
                    throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
            }
            OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "merge_bam_files", null);
            ogsJob.submit();
            logger.info("OGS job ID is " + ogsJob.getID() + ".");
            ogsJob.waitFor();
            break;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
		}

	}
	
	/**
	 * Sort a bam file in coordinate order and write to a new file
	 * @param input Input bam file
	 * @param output Output sorted bam file
	 */
	public static void sortBam(String input, String output) {
		sortBam(input, output, "coordinate");
	}	
	
	/**
	 * Sort a bam file and write to a new file
	 * @param input Input bam file
	 * @param output Output sorted bam file
	 * @param sortOrder Sort order ("coordinate" or "queryname")
	 */
	public static void sortBam(String input, String output, String sortOrder) {
		
		logger.info("Sorting file " + input + ". Writing sorted bam to file " + output + ".");
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		String[] a = new String[3];
		a[0] = "INPUT=" + input;
		a[1] = "OUTPUT=" + output;
		a[2] = "SORT_ORDER=" + sortOrder;
		SortSam.main(a);
		reader.close();

	}

	
	/**
	 * Index a bam file and write the index to a file
	 * @param input Input bam file
	 * @param output Output bam index file
	 */
	public static void indexBam(String input, String output) {
		
		logger.info("Indexing file " + input + ". Writing index to file " + output + ".");
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		String[] a = new String[2];
		a[0] = "INPUT=" + input;
		a[1] = "OUTPUT=" + output;
		BuildBamIndex.main(a);
		reader.close();

	}
	
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
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static Collection<Job> writeGenomicSpaceStats(Map<String, String> bamFiles, String chrSizeFile, String alignmentGlobalStatsJar, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return writeGenomicSpaceStats(bamFiles, chrSizeFile, alignmentGlobalStatsJar, ".", scheduler, drmaaSession);
	}
	
	/**
	 * Write global genomic space stats for bam files
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param chrSizeFile Chromosome size file for genomic space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutputDir Directory to write LSF output files to
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Set of jobs
	 * @throws IOException
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static Job writeGenomicSpaceStats(String bamFile, String sampleName, String chrSizeFile, String alignmentGlobalStatsJar, String bsubOutputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		return writeGenomicSpaceStats(bamFiles, chrSizeFile, alignmentGlobalStatsJar, bsubOutputDir, scheduler, drmaaSession).iterator().next();
	}

	
	/**
	 * Write global genomic space stats for bam files
	 * @param bamFiles Bam files by sample name
	 * @param chrSizeFile Chromosome size file for genomic space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutputDir Directory to write LSF output files to
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Set of jobs
	 * @throws IOException
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static Collection<Job> writeGenomicSpaceStats(Map<String, String> bamFiles, String chrSizeFile, String alignmentGlobalStatsJar, String bsubOutputDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
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
				job.submit();
				jobs.add(job);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "genomic_space_stats", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                jobs.add(ogsJob);
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
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static Collection<Job> writeTranscriptomeSpaceStats(Map<String, String> bamFiles, String bedFile, String alignmentGlobalStatsJar, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		return writeTranscriptomeSpaceStats(bamFiles, bedFile, alignmentGlobalStatsJar, ".", scheduler, drmaaSession);
	}
	
	/**
	 * Write global transcriptome space stats for bam files
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param bedFile Bed annotation for transcriptome space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutDir Directory to write bsub output to
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static Job writeTranscriptomeSpaceStats(String bamFile, String sampleName, String bedFile, String alignmentGlobalStatsJar, String bsubOutDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		return writeTranscriptomeSpaceStats(bamFiles, bedFile, alignmentGlobalStatsJar, bsubOutDir, scheduler, drmaaSession).iterator().next();
	}
	
	/**
	 * Write global transcriptome space stats for bam files
	 * @param bamFiles Bam files by sample name
	 * @param bedFile Bed annotation for transcriptome space
	 * @param alignmentGlobalStatsJar Jar file for alignment global stats
	 * @param bsubOutDir Directory to write bsub output to
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @return Set of LSF job IDs
	 * @throws IOException
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static Collection<Job> writeTranscriptomeSpaceStats(Map<String, String> bamFiles, String bedFile, String alignmentGlobalStatsJar, String bsubOutDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
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
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/compute_transcriptome_space_stats_" + jobID + ".bsub", "week", 32);	
				jobs.add(lsfJob);
				lsfJob.submit();
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "transcriptome_space_stats", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                jobs.add(ogsJob);
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
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void makeTdfs(String bamFile, String sampleName, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		makeTdfs(bamFile, sampleName, ".", assemblyFasta, igvtoolsExecutable, scheduler, drmaaSession);
	}
	
	/**
	 * Make tdf files for bam files and write bsub output to working directory
	 * @param bamFilesBySampleName Bam file name by sample name
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void makeTdfs(Map<String, String> bamFilesBySampleName, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		makeTdfs(bamFilesBySampleName, ".", assemblyFasta, igvtoolsExecutable, scheduler, drmaaSession);
	}
	
	/**
	 * Make tdf files for bam file
	 * @param bamFile Bam file
	 * @param sampleName Sample name
	 * @param bsubOutDir Directory to write bsub output to
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void makeTdfs(String bamFile, String sampleName, String bsubOutDir, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		makeTdfs(bamFiles, bsubOutDir, assemblyFasta, igvtoolsExecutable, scheduler, drmaaSession);
	}
	
	/**
	 * Create paired end bam files for a set of bam files
	 * @param bamFiles Bam files
	 * @param transcriptionRead Read in direction of transcription e.g. "first" or "second" or "unstranded"
	 * @param bsubOutDir Directory to write bsub output to
	 * @param pairedEndWriterJar PairedEndWriter jar file
	 * @param scheduler Scheduler
	 * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws DrmaaException
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void createPairedEndBamFiles(Collection<String> bamFiles, String transcriptionRead, String bsubOutDir, String pairedEndWriterJar, Scheduler scheduler, Session drmaaSession) throws DrmaaException, IOException, InterruptedException {
		logger.info("Creating paired end bam files for " + bamFiles.size() + " bam files.");
		ArrayList<Job> jobs = new ArrayList<Job>();
		for(String bamFile : bamFiles) {
			logger.info("Creating paired end bam file for " + bamFile + "...");
			String cmmd = "java -jar -Xmx30g -Xms20g -Xmn10g " + pairedEndWriterJar + " -b " + bamFile + " -t " + transcriptionRead;
			logger.info("Running command: " + cmmd);
			switch(scheduler) {
			case LSF:
				String jobID = Long.valueOf(System.currentTimeMillis()).toString();
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/create_paired_end_bam_" + jobID + ".bsub", "week", 32);		
				job.submit();
				jobs.add(job);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "paired_end_writer", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                jobs.add(ogsJob);
                break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
		}
		// Wait for jobs to finish
		logger.info("Waiting for PairedEndWriter jobs to finish...");
		JobUtils.waitForAll(jobs);
	}
	
	/**
	 * Make tdf files for bam files
	 * @param bamFilesBySampleName Bam file name by sample name
	 * @param bsubOutDir Directory to write bsub output to
	 * @param assemblyFasta Fasta file of assembly
	 * @param igvtoolsExecutable Igvtools executable file
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void makeTdfs(Map<String, String> bamFilesBySampleName, String bsubOutDir, String assemblyFasta, String igvtoolsExecutable, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
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
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/make_tdf_" + jobID + ".bsub", "hour", 1);
				lsfJob.submit();
				tdfJobs.add(lsfJob);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "make_tdf", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                tdfJobs.add(ogsJob);
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
		} catch(Exception e) {
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
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void samToBam(String sampleName, String samFile, String bamFile, String finalBamFile, String samtools, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		samToBam(sampleName, samFile, bamFile, finalBamFile, ".", samtools, scheduler, drmaaSession);
	}
	
	/**
	 * Convert sam files to bam files
	 * @param samFiles The sam files by sample name
	 * @param bamFiles The bam files to write, by sample name
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param samtools Samtools executables
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void samToBam(Map<String, String> samFiles, Map<String, String> bamFiles, Map<String, String> bsubOutputDirs, String samtools, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		samToBam(samFiles, bamFiles, null, bsubOutputDirs, samtools, scheduler, drmaaSession);
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
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void samToBam(String sampleName, String samFile, String bamFile, String finalBamFile, String bsubOutputDir, String samtools, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		Map<String, String> samFiles = new TreeMap<String, String>();
		samFiles.put(sampleName, samFile);
		Map<String, String> finalBamFiles = new TreeMap<String, String>();
		finalBamFiles.put(sampleName, finalBamFile);
		Map<String, String> bsubOutputDirs = new TreeMap<String, String>();
		bsubOutputDirs.put(sampleName, bsubOutputDir);
		samToBam(samFiles, bamFiles, finalBamFiles, bsubOutputDirs, samtools, scheduler, drmaaSession);
	}
	
	
	/**
	 * Convert sam files to bam files
	 * @param samFiles The sam files by sample name
	 * @param bamFiles The bam files to write, by sample name
	 * @param finalBamFiles Final bam files; skip if they already exist
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param samtools Samtools executables
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void samToBam(Map<String, String> samFiles, Map<String, String> bamFiles, Map<String, String> finalBamFiles, Map<String, String> bsubOutputDirs, String samtools, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
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
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutputDirs.get(sample) + "/sam_to_bam_" + jobID + ".bsub", "hour", 1);
				lsfJob.submit();
				cbJobs.add(lsfJob);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "sam_to_bam", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                cbJobs.add(ogsJob);
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
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void sortBamFile(String sampleName, String unsortedBam, String sortedBam, String finalBam, String picardJarDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		sortBamFile(sampleName, unsortedBam, sortedBam, ".", finalBam, picardJarDir, scheduler, drmaaSession);
	}
	
	/**
	 * Sort all bam files and write sorted files to specified locations, and write bsub output to working directory
	 * @param unsortedBams The unsorted bam files to sort, by sample name
	 * @param sortedBams The sorted bam files to write, by sample name
	 * @param finalBams Final bam files; skip if they already exist
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void sortBamFiles(Map<String, String> unsortedBams, Map<String, String> sortedBams, Map<String, String> finalBams, String picardJarDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		sortBamFiles(unsortedBams, sortedBams, ".", finalBams, picardJarDir, scheduler, drmaaSession);
	}
	
	/**
	 * Sort all bam files and write sorted files to specified locations
	 * @param unsortedBams The unsorted bam files to sort, by sample name
	 * @param sortedBams The sorted bam files to write, by sample name
	 * @param bsubOutputDir Directory to write all bsub output files to
	 * @param finalBams Final bam files; skip if they already exist
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void sortBamFiles(Map<String, String> unsortedBams, Map<String, String> sortedBams, String bsubOutputDir, Map<String, String> finalBams, String picardJarDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> outDirs = new TreeMap<String, String>();
		for(String sample : unsortedBams.keySet()) {
			outDirs.put(sample, bsubOutputDir);
		}
		sortBamFiles(unsortedBams, sortedBams, outDirs, finalBams, picardJarDir, scheduler, drmaaSession);
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
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void sortBamFile(String sampleName, String unsortedBam, String sortedBam, String bsubOutputDir, String finalBam, String picardJarDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> unsortedBams = new TreeMap<String, String>();
		unsortedBams.put(sampleName, unsortedBam);
		Map<String, String> sortedBams = new TreeMap<String, String>();
		sortedBams.put(sampleName, sortedBam);
		Map<String, String> finalBamFiles = new TreeMap<String, String>();
		finalBamFiles.put(sampleName, finalBam);
		Map<String, String> bsubOutputDirs = new TreeMap<String, String>();
		bsubOutputDirs.put(sampleName, bsubOutputDir);
		sortBamFiles(unsortedBams, sortedBams, bsubOutputDirs, finalBamFiles, picardJarDir, scheduler, drmaaSession);
	}
	
	
	/**
	 * Sort all bam files and write sorted files to specified locations
	 * @param unsortedBams The unsorted bam files to sort, by sample name
	 * @param sortedBams The sorted bam files to write, by sample name
	 * @param bsubOutputDirs The directories to write bsub output to, by sample name
	 * @param finalBams Final bam files; skip if they already exist
	 * @param picardJarDir Directory containing Picard jar files
	 * @param scheduler Scheduler
     * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void sortBamFiles(Map<String, String> unsortedBams, Map<String, String> sortedBams, Map<String, String> bsubOutputDirs, Map<String, String> finalBams, String picardJarDir, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
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
				LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutputDirs.get(sample) + "/sort_bam_" + jobID + ".bsub", "hour", 4);
				lsfJob.submit();
				sbJobs.add(lsfJob);
				break;
            case OGS:
                if(drmaaSession == null) {
                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
                }
                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "sort_bam", null);
                ogsJob.submit();
                logger.info("OGS job ID is " + ogsJob.getID() + ".");
                sbJobs.add(ogsJob);
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
	 * Index a set of bam files
	 * @param bamFiles Map of sample name to bam file to index
	 * @param samtools Samtools executable
	 * @param scheduler Scheduler
	 * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void indexBamFiles(Map<String, String> bamFiles, String samtools, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		ArrayList<Job> indexJobs = new ArrayList<Job>();
		for(String sample : bamFiles.keySet()) {
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
			switch(scheduler) {
				case LSF:
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					logger.info("LSF job ID is " + jobID + ".");
					// Submit job
					LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, indexfile.getParent() + "/index_bam_" + jobID + ".bsub", "hour", 1);
					indexJobs.add(job);
					job.submit();
					break;
	            case OGS:
	                if(drmaaSession == null) {
	                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	                }
	                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "index_bam", null);
	                ogsJob.submit();
	                logger.info("OGS job ID is " + ogsJob.getID() + ".");
	                indexJobs.add(ogsJob);
	                break;
	
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
			}
			
		}
		logger.info("Waiting for samtools jobs to finish...");
		JobUtils.waitForAll(indexJobs);
	}

	
}

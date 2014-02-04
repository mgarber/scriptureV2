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
import java.util.TreeSet;

import nextgen.core.annotation.Gene;
import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.LSFJob;
import nextgen.core.job.OGSJob;
import nextgen.core.pipeline.Scheduler;
import nextgen.core.readFilters.FirstOfPairFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.readFilters.SecondOfPairFilter;
import nextgen.core.writers.WigWriter;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import broad.core.sequence.FastaSequenceIO;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class WigUtils {

	private static Logger logger = Logger.getLogger(WigUtils.class.getName());
	
	/**
	 * Write fragment end points and midpoints to wig and bigwig files
	 * @param bamFiles Bam files by sample name
	 * @param pairedData Whether each sample has paired end reads, by sample name
	 * @param bamDir Directory containing bam files
	 * @param refFasta Fasta file of sequences these bam files were aligned against
	 * @param geneBedFile Bed file of genes to count reads in or null if using genomic space
	 * @param wigWriterJar WigWriter jar file
	 * @param wigToBigWigExecutable WigToBigWig executable file
	 * @param scheduler Scheduler
	 * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void writeWigFragmentEndsAndMidpoints(Map<String, String> bamFiles, Map<String, Boolean> pairedData, String bamDir, String refFasta, String geneBedFile, String wigWriterJar, String wigToBigWigExecutable, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		
		if(geneBedFile != null) {
			if(geneBedFile.equals("null")) {
				geneBedFile = null;
			}
		}
		
		Map<String, String> read1endWig = new TreeMap<String, String>();
		Map<String, String> read2endWig = new TreeMap<String, String>();
		Map<String, String> read1endBigwig = new TreeMap<String, String>();
		Map<String, String> read2endBigwig = new TreeMap<String, String>();
		Map<String, String> midpointWig = new TreeMap<String, String>();
		Map<String, String> midpointBigwig = new TreeMap<String, String>();
		ArrayList<Job> wigJobs = new ArrayList<Job>();
		ArrayList<Job> bigwigJobs = new ArrayList<Job>();
		
		// Chromosome size file to pass to wig writer if using genomic space
		// Null if using transcriptome space
		String chrSizesForWigToBigWig = FastaUtils.writeSizeFile(refFasta);
		String chrSizesForWigWriter = null;
		if(geneBedFile == null) {
			chrSizesForWigWriter = chrSizesForWigToBigWig;
		}
	
		// Make file names
		for(String sampleName : bamFiles.keySet()) {
			
			String bamFile = bamFiles.get(sampleName);
			read1endWig.put(sampleName, bamFile + ".read1.wig");
			read2endWig.put(sampleName, bamFile + ".read2.wig");
			read1endBigwig.put(sampleName, bamFile + ".read1.bw");
			read2endBigwig.put(sampleName, bamFile + ".read2.bw");
			midpointWig.put(sampleName, bamFile + ".midpoint.wig");
			midpointBigwig.put(sampleName, bamFile + ".midpoint.bw");
			
		}
		
		// Write wig files
		for(String sampleName : bamFiles.keySet()) {
			
			String bamFile = bamFiles.get(sampleName);
			
			// Write wig file for read1
			String wig1 = read1endWig.get(sampleName);
			File read1wigFile = new File(wig1);
			String bigwig1 = read1endBigwig.get(sampleName);
			File read1bigwigFile = new File(bigwig1);
			String prefix1 = wig1.replaceAll(".wig", "");
			if(read1wigFile.exists() || read1bigwigFile.exists()) {
				logger.warn("Read 1 wig file or bigwig file for sample " + sampleName + " already exists. Not remaking wig file.");
			} else {
				logger.info("Writing fragment ends of read 1 from bam file " + bamFile + " to wig file " + wig1 + ".");
				// Write wig file for read1
				if(read1wigFile.exists() || read1bigwigFile.exists()) {
					logger.warn("Read 1 wig file or bigwig file for sample " + sampleName + " already exists. Not remaking wig file.");
				} else {
					String cmmd = "java -Xmx30g -Xms29g -Xmn28g -jar " + wigWriterJar + " -b " + bamFile + " -g " + geneBedFile + " -o " + prefix1 + " -c " + chrSizesForWigWriter + " -sp beginning -pp true -r1 true -pe false"; 
					logger.info("");
					logger.info("Writing fragment ends of read 1 from bam file " + bamFile + " to wig file " + wig1 + ".");
					logger.info("Running WigWriter command " + cmmd);
					switch(scheduler) {
					case LSF:
						String jobID = Long.valueOf(System.currentTimeMillis()).toString();
						logger.info("LSF job ID is " + jobID + ".");
						LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_fragment_ends_" + jobID + ".bsub", "week", 32);
						lsfJob.submit();
						wigJobs.add(lsfJob);
						break;
		            case OGS:
		                if(drmaaSession == null) {
		                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
		                }
		                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "wig_fragment_ends");
		                ogsJob.submit();
		                logger.info("OGS job ID is " + ogsJob.getID() + ".");
		                wigJobs.add(ogsJob);
		                break;
					default:
						throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
					}
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
				String prefix = midpointWigFileName.replaceAll(".wig", "");
				String cmmd = "java -Xmx30g -Xms29g -Xmn28g -jar " + wigWriterJar + " -b " + bamFile + " -g " + geneBedFile + " -o " + prefix + " -c " + chrSizesForWigWriter + " -sp midpoint -pp true -pe true"; 
				logger.info("");
				logger.info("Writing fragment midpoints from bam file " + bamFile + " to wig file " + midpointWigFileName + ".");
				logger.info("Running WigWriter command " + cmmd);
				switch(scheduler) {
				case LSF:
					String jobID = Long.valueOf(System.currentTimeMillis()).toString();
					logger.info("LSF job ID is " + jobID + ".");
					LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_fragment_midpoints_" + jobID + ".bsub", "week", 32);
					lsfJob.submit();
					wigJobs.add(lsfJob);
					break;
	            case OGS:
	                if(drmaaSession == null) {
	                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	                }
	                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "wig_fragment_midpoints");
	                ogsJob.submit();
	                logger.info("OGS job ID is " + ogsJob.getID() + ".");
	                wigJobs.add(ogsJob);
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
					String prefix2 = wig2.replaceAll(".wig", "");
					String cmmd = "java -Xmx30g -Xms29g -Xmn28g -jar " + wigWriterJar + " -b " + bamFile + " -g " + geneBedFile + " -o " + prefix2 + " -c " + chrSizesForWigWriter + " -sp beginning -pp true -r2 true -pe false"; 
					logger.info("");
					logger.info("Writing fragment ends of read 2 from bam file " + bamFile + " to wig file " + wig2 + ".");
					logger.info("Running WigWriter command " + cmmd);
					switch(scheduler) {
					case LSF:
						String jobID = Long.valueOf(System.currentTimeMillis()).toString();
						logger.info("LSF job ID is " + jobID + ".");
						LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_fragment_ends_" + jobID + ".bsub", "week", 32);
						lsfJob.submit();
						wigJobs.add(lsfJob);
						break;
		            case OGS:
		                if(drmaaSession == null) {
		                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
		                }
		                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "wig_fragment_ends");
		                ogsJob.submit();
		                logger.info("OGS job ID is " + ogsJob.getID() + ".");
		                wigJobs.add(ogsJob);
		                break;
					default:
						throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
					}
				}
			}
		
			
		}
		logger.info("");
		logger.info("Waiting for WigWriter jobs to finish...");
		JobUtils.waitForAll(wigJobs);
	
		
		
		// Write bigwig files
		for(String sampleName : bamFiles.keySet()) {
			String wig1 = read1endWig.get(sampleName);
			String bigwig1 = read1endBigwig.get(sampleName);
			File read1bigwigFile = new File(bigwig1);
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
					LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bamDir + "/wig_to_bigwig_" + jobID + ".bsub", "hour", 4);
					lsfJob.submit();
					bigwigJobs.add(lsfJob);
					break;
	            case OGS:
	                if(drmaaSession == null) {
	                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	                }
	                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "bigwig_fragment_ends");
	                ogsJob.submit();
	                logger.info("OGS job ID is " + ogsJob.getID() + ".");
	                bigwigJobs.add(ogsJob);
	                break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
				}
			}
			
			// Write midpoint bigwig file
			String midpointWigFileName = midpointWig.get(sampleName);
			String midpointBigwigFileName = midpointBigwig.get(sampleName);
			File midpointBigwigFile = new File(midpointBigwigFileName);
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
	            case OGS:
	                if(drmaaSession == null) {
	                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	                }
	                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "bigwig_fragment_ends");
	                ogsJob.submit();
	                logger.info("OGS job ID is " + ogsJob.getID() + ".");
	                bigwigJobs.add(ogsJob);
	                break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
				}
			}
			
			// Write bigwig file for read2
			String wig2 = read2endWig.get(sampleName);
			String bigwig2 = read2endBigwig.get(sampleName);
			File read2bigwigFile = new File(bigwig2);
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
	            case OGS:
	                if(drmaaSession == null) {
	                        throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	                }
	                OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "bigwig_fragment_ends");
	                ogsJob.submit();
	                logger.info("OGS job ID is " + ogsJob.getID() + ".");
	                bigwigJobs.add(ogsJob);
	                break;
				default:
					throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " not supported.");
				}
			}
	
		}
		
		
	
		logger.info("");
		logger.info("Waiting for wigToBigWig jobs to finish...");
		JobUtils.waitForAll(bigwigJobs);
		
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
	 * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void writeWigPositionCount(String bamFile, String sampleName, String bamDir, String geneBedFile, String refFasta, String wigToBigWigExecutable, String wigWriterJar, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		WigUtils.writeWigPositionCount(bamFiles, bamDir, geneBedFile, refFasta, wigToBigWigExecutable, wigWriterJar, scheduler, drmaaSession);
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
	 * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void writeWigPositionCount(Map<String, String> bamFiles, String bamDir, String geneBedFile, String refFasta, String wigToBigWigExecutable, String wigWriterJar, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Collection<String> chrNames = new TreeSet<String>();
		String chrSizeFile = null;
		if(geneBedFile != null) {
			Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneBedFile));
			chrNames.addAll(genes.keySet());
		} else if(refFasta != null) {
			chrNames.addAll(FastaSequenceIO.getSequenceNames(refFasta));
			chrSizeFile = FastaSequenceIO.createSizeFile(refFasta);
		}
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
					String cmmd = "java -jar -Xmx30g -Xms29g -Xmn28g " + wigWriterJar + " -b " + bamFile + " -g " + geneBedFile + " -n true -o " + prefix + " -chr " + chr;
					if(chrSizeFile != null) {
						cmmd += " -c " + chrSizeFile;
					}
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
	                case OGS:
	                    if(drmaaSession == null) {
	                            throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	                    }
	                    OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "wig_position_count");
	                    ogsJob.submit();
	                    logger.info("OGS job ID is " + ogsJob.getID() + ".");
	                    wigJobs.add(ogsJob);
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
					String cmmd = "java -jar -Xmx30g -Xms29g -Xmn28g " + wigWriterJar + " -b " + bamFile + " -g " + geneBedFile + " -n false -o " + prefix + " -chr " + chr;
					if(chrSizeFile != null) {
						cmmd += " -c " + chrSizeFile;
					}
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
	                case OGS:
	                    if(drmaaSession == null) {
	                            throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	                    }
	                    OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "wig_position_count");
	                    ogsJob.submit();
	                    logger.info("OGS job ID is " + ogsJob.getID() + ".");
	                    wigJobs.add(ogsJob);
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
		String sizeFile = FastaUtils.writeSizeFile(refFasta);
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
			String cmmd = wigToBigWigExecutable + " " + wig + " " + sizeFile + " " + bigwig;
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
	        case OGS:
	            if(drmaaSession == null) {
	                    throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	            }
	            OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "wig_position_count");
	            ogsJob.submit();
	            logger.info("OGS job ID is " + ogsJob.getID() + ".");
	            bigwigJobs.add(ogsJob);
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
			String cmmd = wigToBigWigExecutable + " " + wig + " " + sizeFile + " " + bigwig;
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
	        case OGS:
	            if(drmaaSession == null) {
	                    throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
	            }
	            OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "wig_position_count");
	            ogsJob.submit();
	            logger.info("OGS job ID is " + ogsJob.getID() + ".");
	            bigwigJobs.add(ogsJob);
	            break;
			default:
				throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
			}
		}
		logger.info("Waiting for wigToBigWig jobs to finish...");
		JobUtils.waitForAll(bigwigJobs);
		logger.info("Done writing wig files.");
	}

	/**
	 * Write fragment end points and midpoints to wig and bigwig files
	 * @param sampleName Sample name
	 * @param bamFile Bam file
	 * @param pairedData Whether the sample has paired end reads
	 * @param bamDir Directory containing bam files
	 * @param refFasta Fasta file of sequences these bam files were aligned against
	 * @param geneBedFile Bed file of genes to count reads in or null if using genomic space
	 * @param wigWriterJar WigWriter jar file
	 * @param wigToBigWigExecutable WigToBigWig executable file
	 * @param scheduler Scheduler
	 * @param drmaaSession Active DRMAA session. Pass null if not using OGS. There should only be one active session at a time. Session should have been created in the main method of the class calling this method.
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public static void writeWigFragmentEndsAndMidpoints(String sampleName, String bamFile, boolean pairedData, String bamDir, String refFasta, String geneBedFile, String wigWriterJar, String wigToBigWigExecutable, Scheduler scheduler, Session drmaaSession) throws IOException, InterruptedException, DrmaaException {
		Map<String, String> bamFiles = new TreeMap<String, String>();
		bamFiles.put(sampleName, bamFile);
		Map<String, Boolean> paired = new TreeMap<String, Boolean>();
		paired.put(sampleName, Boolean.valueOf(pairedData));
		writeWigFragmentEndsAndMidpoints(bamFiles, paired, bamDir, refFasta, geneBedFile, wigWriterJar, wigToBigWigExecutable, scheduler, drmaaSession);
	}
	
	
	
}

/**
 * Call protein binding sites on RNA with Protect-seq data from multiple control and multiple signal samples
 */
package broad.pda.seq.protection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.math.MathUtil;
import broad.core.math.Statistics;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.util.PipelineUtils;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.analysis.PeakCaller;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.utils.AlignmentUtils;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 *
 */
public class MultiSampleScanPeakCaller implements PeakCaller {
	
	private TranscriptomeSpace coord;
	protected Map<String, Collection<Gene>> genes;
	private Map<Gene, Map<Annotation, Double>> tStatisticWindowScores;
	private Map<SampleData, Map<Gene, Map<Annotation, Double>>> singleSampleWindowEnrichmentOverGene;
	private Map<SampleData, Map<Gene, Collection<Annotation>>> singleSampleScanPeaks;
	protected GenomeSpaceSampleData expressionData;
	protected ArrayList<SampleData> controlSamples;
	protected ArrayList<SampleData> signalSamples;
	protected ArrayList<SampleData> allSamples;
	protected static Logger logger = Logger.getLogger(MultiSampleScanPeakCaller.class.getName());
	protected int windowSize;
	protected int stepSize;
	private static int DEFAULT_WINDOW_SIZE = 30;
	private static int DEFAULT_STEP_SIZE = 1;
	private static double DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF = 0.001;
	private static double DEFAULT_PEAK_WINDOW_COUNT_CUTOFF = 10;
	private static double DEFAULT_TRIM_PEAK_QUANTILE = 0.6;
	private static int DEFAULT_BATCH_MEM_REQUEST = 8;
	protected int numControls;
	protected int numSignals;
	protected int numSamples;
	protected Random random;
	private Collection<SamplePermutation> sampleIdentityPermutations;
	private double peakWindowScanPvalCutoff;
	private double peakWindowCountCutoff;
	private double trimQuantile;
	private boolean permutationScoring;
	private String sampleFile;
	private String bedAnnotationFile;
	private String sizeFile;
	private static int RGB_RED_WITH_GENE = 106;
	private static int RGB_GREEN_WITH_GENE = 7;
	private static int RGB_BLUE_WITH_GENE = 205;
	private static int RGB_RED_AGAINST_GENE = 218;
	private static int RGB_GREEN_AGAINST_GENE = 2;
	private static int RGB_BLUE_AGAINST_GENE = 38;
	private static int RGB_RED_UNKNOWN = 62;
	private static int RGB_GREEN_UNKNOWN = 63;
	private static int RGB_BLUE_UNKNOWN = 73;
	
	
	protected MultiSampleScanPeakCaller(MultiSampleScanPeakCaller other) throws IOException {
		this(other.sampleFile, other.bedAnnotationFile, other.sizeFile, other.windowSize, other.stepSize, other.peakWindowScanPvalCutoff, other.peakWindowCountCutoff, other.trimQuantile);
	}
	
	/**
	 * Instantiate with default parameters
	 * @param sampleListFile File containing sample list
	 * @param bedFile Bed gene annotation
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private MultiSampleScanPeakCaller(String sampleListFile, String bedFile, String chrSizeFile) throws IOException {
		this(sampleListFile, bedFile, chrSizeFile, DEFAULT_WINDOW_SIZE, DEFAULT_STEP_SIZE, DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF, DEFAULT_PEAK_WINDOW_COUNT_CUTOFF, DEFAULT_TRIM_PEAK_QUANTILE);
	}
	
	
	/**
	 * Specify parameters
	 * @param sampleListFile File containing sample list
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @param window Window size
	 * @param step Step size
	 * @param maxPermutations Max number of permutations
	 * @param peakScanPvalCutoff P value cutoff for scan statistic
	 * @param trimPeakQuantile Quantile for trim max contiguous algorithm
	 * @param doPermutationScoring Whether to do permutation test
	 * @throws IOException
	 */
	private MultiSampleScanPeakCaller(String sampleListFile, String bedFile, String chrSizeFile, int window, int step, double peakScanPvalCutoff, double peakCountCutoff, double trimPeakQuantile) throws IOException {
		
		// Set basic parameters
		sampleFile = sampleListFile;
		bedAnnotationFile = bedFile;
		sizeFile = chrSizeFile;
		windowSize = window;
		stepSize = step;
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		coord = new TranscriptomeSpace(genes);
		random = new Random();
		peakWindowScanPvalCutoff = peakScanPvalCutoff;
		peakWindowCountCutoff = peakCountCutoff;
		trimQuantile = trimPeakQuantile;
		
		// Read sample information
		SampleFileParser p = new SampleFileParser(sampleListFile, chrSizeFile);
		controlSamples = p.getControlDatasets();
		numControls = controlSamples.size();
		signalSamples = p.getSignalDatasets();
		expressionData = p.getExpressionData();
		numSignals = signalSamples.size();
		allSamples = new ArrayList<SampleData>();
		allSamples.addAll(controlSamples);
		allSamples.addAll(signalSamples);
		numSamples = allSamples.size();

		// Initialize score maps
		tStatisticWindowScores = new TreeMap<Gene, Map<Annotation, Double>>();
		singleSampleWindowEnrichmentOverGene = new HashMap<SampleData, Map<Gene, Map<Annotation, Double>>>();
		singleSampleScanPeaks = new HashMap<SampleData, Map<Gene, Collection<Annotation>>>();
		for(SampleData sample : allSamples) {
			if(!singleSampleWindowEnrichmentOverGene.containsKey(sample)) {
				Map<Gene, Map<Annotation, Double>> m = new TreeMap<Gene, Map<Annotation, Double>>();
				singleSampleWindowEnrichmentOverGene.put(sample, m);
			}
			if(!singleSampleScanPeaks.containsKey(sample)) {
				Map<Gene, Collection<Annotation>> m = new TreeMap<Gene, Collection<Annotation>>();
				singleSampleScanPeaks.put(sample, m);
			}
		}
		
	}
	
	/**
	 * Whether the gene is expressed in all control samples at the given significance level
	 * @param gene The gene
	 * @return Whether the gene is expressed by these criteria
	 */
	public boolean isExpressedInAllControlSamples(Gene gene) {
		for(SampleData control : controlSamples) {
			if(!control.isExpressed(gene)) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Whether the gene is significantly expressed in the special expression sample
	 * @param gene The gene
	 * @return True iff the gene is expressed in the expression sample
	 */
	public boolean isExpressed(Gene gene) {
		return expressionData.isExpressed(gene);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs) throws IOException, InterruptedException {
		batchWriteSingleSampleScanPeaksAllSamples(commandArgs, null, DEFAULT_BATCH_MEM_REQUEST);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs, String chrListFile) throws IOException, InterruptedException {
		batchWriteSingleSampleScanPeaksAllSamples(commandArgs, chrListFile, DEFAULT_BATCH_MEM_REQUEST);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs, int memRequestGb) throws IOException, InterruptedException {
		batchWriteSingleSampleScanPeaksAllSamples(commandArgs, null, memRequestGb);
	}
	
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs, String chrListFile, int memRequestGb) throws IOException, InterruptedException {
		
		logger.info("\nBatching out peak calling by sample and chromosome...\n");
		
		int xmx = (int)Math.floor(0.9 * memRequestGb);
		int xms = (int)Math.floor(0.7 * memRequestGb);
		int xmn = (int)Math.floor(0.5 * memRequestGb);
		
		TreeSet<String> chrs = new TreeSet<String>();
		if(chrListFile == null) {
			chrs.addAll(genes.keySet());
		} else {
			FileReader r = new FileReader(chrListFile);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			while(b.ready()) {
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() == 0) continue;
				if(!genes.keySet().contains(line)) {
					throw new IllegalArgumentException("Chromosome name " + line + " not recognized.");
				}
				chrs.add(line);
			}
			r.close();
			b.close();
		}
		
		String jar = commandLineBatchJar(commandArgs);
		ArrayList<String> jobIDs = new ArrayList<String>();
		String outDir = commandLineOutDir(commandArgs);
		File o = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = o.mkdir();
		if(!o.exists()) {
			throw new IOException("Could not create directory " + outDir);
		}
		
		Map<String, String> cmmds = new TreeMap<String, String>();
		
		for(SampleData sample : allSamples) {
			for(String chr : chrs) {
				String[] batchedCmmdArgs = BatchedMultiSampleScanPeakCaller.extendSuperArgsForSampleAndChr(commandArgs, sample.getSampleName(), chr);
				String args = "";
				for(int i=0; i < batchedCmmdArgs.length; i++) {
					args += batchedCmmdArgs[i] + " ";
				}
				String cmmd = "java -jar -Xmx" + xmx + "g -Xms" + xms + "g -Xmn" + xmn + "g " + jar + " " + args;
				logger.info("Running command: " + cmmd);
				String jobID = sample.getSampleName() + "_" + chr + "_" + Long.valueOf(System.currentTimeMillis()).toString();
				jobIDs.add(jobID);
				cmmds.put(jobID, cmmd);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				PipelineUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, outDir + "/" + jobID + ".bsub", "week", memRequestGb);

			}
		}

		boolean allJobsReservedHeapSpace = false;
		while(!allJobsReservedHeapSpace) {
			logger.info("Waiting for all jobs to start...");
			PipelineUtils.waitForAllJobsToStart(jobIDs, Runtime.getRuntime());
			logger.info("All jobs have started.");
			Thread.sleep(5000);
			ArrayList<String> jobsThatFailedHeapSpace = new ArrayList<String>();
			for(String jobID : jobIDs) {
				String out = outDir + "/" + jobID + ".bsub";
				File outFile = new File(out);
				if(outFile.exists()) {
					if(PipelineUtils.jobFailedCouldNotReserveHeapSpace(out)) {
						jobsThatFailedHeapSpace.add(jobID);
					}
				}
			}
			if(jobsThatFailedHeapSpace.isEmpty()) {
				allJobsReservedHeapSpace = true;
				continue;
			}
			for(String jobID : jobsThatFailedHeapSpace) {
				String cmmd = cmmds.get(jobID);
				logger.info("Resubmitting command because heap space reservation failed: " + cmmd);
				String newJobID = jobID + "_" + Long.valueOf(System.currentTimeMillis()).toString();
				jobIDs.add(newJobID);
				cmmds.put(newJobID, cmmd);
				logger.info("LSF job ID is " + newJobID + ".");
				PipelineUtils.bsubProcess(Runtime.getRuntime(), newJobID, cmmd, outDir + "/" + newJobID + ".bsub", "week", 32);
				jobIDs.remove(jobID);
				jobIDs.add(newJobID);
			}
		}
		
		logger.info("Waiting for jobs to finish...");
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		
		logger.info("\nAll jobs finished.\n");
		
	}
	
	/**
	 * Get name of bed file to write for peaks
	 * @param sample Sample
	 * @param outDir Output directory name or null if current directory
	 * @param chrName Chromosome name or null if all chromosomes
	 * @return File name
	 */
	protected String getPeakBedFileName(SampleData sample, String outDir, String chrName) {
		String rtrn = "";
		if(outDir != null) {
			rtrn += outDir + "/";
		}
		rtrn += sample.getSampleName() + "_scan_peaks_" + windowSize + "_" + stepSize + "_" + peakWindowScanPvalCutoff + "_" + trimQuantile;
		if(chrName != null) {
			rtrn += "_" + chrName;
		}
		rtrn += ".bed";
		return rtrn;
	}
	
	/**
	 * Write scan peaks for all samples to separate bed files
	 * @throws IOException
	 */
	private void writeSingleSampleScanPeaksAllSamples(String outDir) throws IOException {
		logger.info("Writing single sample scan peaks for each sample...");
		File o = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = o.mkdir();
		if(!o.exists()) {
			throw new IOException("Could not create directory " + outDir);
		}
		for(SampleData signal : signalSamples) {
			String outfile = getPeakBedFileName(signal, outDir, null);
			writeSingleSampleScanPeaks(signal, outfile);
		}		
		for(SampleData control : controlSamples) {
			String outfile = getPeakBedFileName(control, outDir, null);
			writeSingleSampleScanPeaks(control, outfile);
		}
		logger.info("Done writing single sample scan peaks for all samples.");
	}
	
	/**
	 * Write all single sample scan peaks for the sample to file
	 * @param sample The sample
	 * @param outFile Output bed file
	 * @param r Red value for bed file color
	 * @param g Green value for bed file color
	 * @param b Blue value for bed file color
	 * @throws IOException 
	 */
	private void writeSingleSampleScanPeaks(SampleData sample, String outFile) throws IOException {
		writeSingleSampleScanPeaks(sample, outFile, null);
	}
	
	/**
	 * Write all single sample scan peaks for the sample and chromosome to file
	 * @param sample The sample
	 * @param outFile Output bed file
	 * @param r Red value for bed file color
	 * @param g Green value for bed file color
	 * @param b Blue value for bed file color
	 * @param chrName Only write peaks for this chromosome
	 * @throws IOException 
	 */
	protected void writeSingleSampleScanPeaks(SampleData sample, String outFile, String chrName) throws IOException {
		logger.info("Writing single sample scan peaks for sample " + sample.getSampleName() + " to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		for(String chr : genes.keySet()) {
			if(chrName != null) {
				if(!chr.equals(chrName)) {
					continue;
				}
			}
			for(Gene gene : genes.get(chr)) {
				Collection<Annotation> peaks = getSingleSampleScanPeaks(sample, gene);
				for(Annotation window : peaks) {
					int r = RGB_RED_UNKNOWN;
					int g = RGB_GREEN_UNKNOWN;
					int b = RGB_BLUE_UNKNOWN;
					if(window.getOrientation().equals(Strand.UNKNOWN)) {
						String name = window.getName();
						name += "_STRAND_UNKNOWN";
						window.setName(name);						
					}
					if(window.getOrientation().equals(gene.getOrientation()) && !window.getOrientation().equals(Strand.UNKNOWN)) {
						r = RGB_RED_WITH_GENE;
						g = RGB_GREEN_WITH_GENE;
						b = RGB_BLUE_WITH_GENE;
					}
					if(!window.getOrientation().equals(gene.getOrientation()) && !window.getOrientation().equals(Strand.UNKNOWN)) {
						String name = window.getName();
						name += "_STRAND_AGAINST_GENE";
						window.setName(name);
						r = RGB_RED_AGAINST_GENE;
						g = RGB_GREEN_AGAINST_GENE;
						b = RGB_BLUE_AGAINST_GENE;
					}
					w.write(window.toBED(r, g, b) + "\n");
				}
			}
		}
		w.close();
		logger.info("Done writing scan peaks for sample " + sample.getSampleName() + ".");
	}
	
	/**
	 * Get scan peaks for the sample and the gene
	 * @param sample The sample
	 * @param gene The gene
	 * @return Significant scan peaks
	 * @throws IOException 
	 */
	public Collection<Annotation> getSingleSampleScanPeaks(SampleData sample, Gene gene) throws IOException {
		if(singleSampleScanPeaks.get(sample).containsKey(gene)) {
			return singleSampleScanPeaks.get(sample).get(gene);
		}
		identifySingleSampleScanPeaks(sample, gene);
		return singleSampleScanPeaks.get(sample).get(gene);
	}
	
	/**
	 * Identify significant peaks within the sample and the gene
	 * @param sample The sample
	 * @param gene The gene
	 * @throws IOException
	 */
	private void identifySingleSampleScanPeaks(SampleData sample, Gene gene) throws IOException {
		
		TreeSet<Annotation> finalPeaks = new TreeSet<Annotation>();
		
		TranscriptomeSpaceAlignmentModel data = sample.getData();
		
		// If gene is not expressed, skip
		if(!isExpressed(gene)) {
			logger.info("Gene " + gene.getName() + " (" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + ") not expressed in expression dataset.");
			singleSampleScanPeaks.get(sample).put(gene, finalPeaks);
			return;
		}
		
		double geneCount = data.getCount(gene);
		int geneSize = gene.getSize();
		double geneAvgCoverage = geneCount / geneSize;
		
		logger.info("Finding scan peaks for sample " + sample.getSampleName() + " and gene " + gene.getName() + " (" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + ")");
		
		// Get fixed size windows with sufficient count and significant scan statistic
		TreeSet<Annotation> scanSignificantWindows = new TreeSet<Annotation>();
		Map<Annotation, ScanStatisticScore> windowScores = sample.getWindowScores(gene);
		for(Annotation window : windowScores.keySet()) {
			ScanStatisticScore score = windowScores.get(window);
			double count = score.getCount();
			if(count < peakWindowCountCutoff) {
				continue;
			}
			double pval = score.getScanPvalue();
			if(pval < peakWindowScanPvalCutoff) {
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\t" + gene.getName());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\t" + window.toBED());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tglobal_length=" + score.getGlobalLength());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tglobal_count=" + score.getTotal());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tglobal_lambda=" + score.getGlobalLambda());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\twindow_size=" + score.getCoordinateSpace().getSize(window));
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\twindow_count=" + score.getCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tpval=" + score.getScanPvalue());
				ScanStatisticScore fragmentLengthFilterScore = sample.scoreWindowWithFragmentLengthFilter(gene, window);
				double pval2 = fragmentLengthFilterScore.getScanPvalue();
				if(pval2 < peakWindowScanPvalCutoff) {
					logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_length=" + fragmentLengthFilterScore.getGlobalLength());
					logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_count=" + fragmentLengthFilterScore.getTotal());
					logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_lambda=" + fragmentLengthFilterScore.getGlobalLambda());
					logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_size=" + fragmentLengthFilterScore.getCoordinateSpace().getSize(window));
					logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_count=" + fragmentLengthFilterScore.getCount());
					logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tpval=" + fragmentLengthFilterScore.getScanPvalue());				
					scanSignificantWindows.add(window);
				} else {
					logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_length=" + fragmentLengthFilterScore.getGlobalLength());
					logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_count=" + fragmentLengthFilterScore.getTotal());
					logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_lambda=" + fragmentLengthFilterScore.getGlobalLambda());
					logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_size=" + fragmentLengthFilterScore.getCoordinateSpace().getSize(window));
					logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_count=" + fragmentLengthFilterScore.getCount());
					logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tpval=" + fragmentLengthFilterScore.getScanPvalue());				
				}
			}
		}
		
		// Merge overlapping windows
		Collection<Annotation> mergedWindows = AnnotationUtils.mergeOverlappingBlocks(scanSignificantWindows);
		for(Annotation window : mergedWindows) {
			logger.debug("MERGED_WINDOW\t" + window.toBED());
		}
		
		
		// Trim each window
		TreeSet<Annotation> trimmedMergedWindows = new TreeSet<Annotation>();
		for(Annotation window : mergedWindows) {
			List<Double> coverageData = data.getPositionCountList(new Gene(window));
			Annotation trimmed = SampleData.trimMaxContiguous(window, coverageData, trimQuantile);
			trimmedMergedWindows.add(trimmed);
			logger.debug("MERGED_TRIMMED_WINDOW\t" + window.toBED());
		}
		
		// Filter by scan statistic again
		for(Annotation window : trimmedMergedWindows) {
			ScanStatisticScore score = sample.scoreWindow(gene, window);
			double p = score.getScanPvalue();
			if(p < peakWindowScanPvalCutoff) {
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\t" + gene.getName());
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\t" + window.toBED());
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tglobal_length=" + score.getGlobalLength());
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tglobal_count=" + score.getTotal());
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tglobal_lambda=" + score.getGlobalLambda());
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\twindow_size=" + score.getCoordinateSpace().getSize(window));
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\twindow_count=" + score.getCount());
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tpval=" + score.getScanPvalue());
				finalPeaks.add(window);
			}
		}
		
		// Add finishing touches to peaks
		for(Annotation window : finalPeaks) {
			
			// Name peaks
			window.setName(gene.getName() + ":" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd());
			
			// Set peak score to enrichment
			double windowCount = data.getCount(window);
			int windowSize1 = window.getSize();
			double windowAvgCoverage = windowCount / windowSize1;
			double enrichment = windowAvgCoverage / geneAvgCoverage;
			window.setScore(enrichment);
			
			// Assign orientation to peaks
			window.setOrientation(AlignmentUtils.assignOrientationToWindow(sample.getOriginalBamFile(), window, sample.firstReadTranscriptionStrand(), 0.9));
		
			logger.debug("FINAL_PEAK\t" + gene.getName());
			logger.debug("FINAL_PEAK\t" + window.toBED());
			logger.debug("FINAL_PEAK\tname=" + window.getName());
			logger.debug("FINAL_PEAK\twindow_count=" + windowCount);
			logger.debug("FINAL_PEAK\twindow_size=" + windowSize1);
			logger.debug("FINAL_PEAK\twindow_avg_coverage=" + windowAvgCoverage);
			logger.debug("FINAL_PEAK\tgene_count=" + geneCount);
			logger.debug("FINAL_PEAK\tgene_size=" + geneSize);
			logger.debug("FINAL_PEAK\tgene_avg_coverage=" + geneAvgCoverage);
			logger.debug("FINAL_PEAK\tenrichment_over_transcript=" + enrichment);
			logger.debug("FINAL_PEAK\tscore=" + window.getScore());
			logger.debug("FINAL_PEAK\torientation=" + window.getOrientation().toString());
			
		}
		
		singleSampleScanPeaks.get(sample).put(gene, finalPeaks);
	}
	
	private void computeSingleSampleWindowEnrichmentsOverGenes() {
		logger.info("Computing window enrichments for each sample...");
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				if(!isExpressed(gene)) {
					continue;
				}
				computeSingleSampleWindowEnrichmentOverGene(gene);
			}
		}		
		//writeSingleSampleWindowScoresToFileIfNeeded();
		logger.info("Done computing single sample window enrichments.");
	}
	
	/**
	 * Score all genes
	 * @throws IOException 
	 */
	@SuppressWarnings("unused")
	private void scoreGenesTStatisticScore() throws IOException {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		computeSingleSampleWindowEnrichmentsOverGenes();
		for(String chr : genes.keySet()) {
			logger.info("Scoring genes on chromosome " + chr);
			for(Gene gene : genes.get(chr)) {
				if(!isExpressed(gene)) {
					logger.info(gene.getName() + " is not expressed. Skipping.");
					continue;
				}
				logger.info("Scoring gene " + gene.getName());
				computeTStatisticWindowScores(gene);
			}
		}
	}
	
	
/*	*//**
	 * For samples that didn't read window scores from file, write to files
	 * @throws IOException
	 *//*
	private void writeSingleSampleWindowScoresToFileIfNeeded() throws IOException {
		for(SampleData sample : allSamples) {
			if(!sample.gotWindowScoresFromFile()) {
				sample.writeWindowScoresToFile(this);
			}
		}
		logger.info("Done writing window score files.");
	}
*/	
	/**
	 * For each sample compute the enrichment of each window over the gene average
	 * Cache the window enrichments
	 * @param gene The gene
	 */
	private void computeSingleSampleWindowEnrichmentOverGene(Gene gene) {
		logger.info("Getting enrichment for each window and each sample...");
		for(SampleData sample : allSamples) {
			Map<Annotation,Double> sampleWindowEnrichments = new TreeMap<Annotation,Double>();
			double geneAvgCoverage = sample.getGeneAverageCoverage(gene);
			logger.info("Sample " + sample.getSampleName() + " Gene average coverage = " + geneAvgCoverage);
			Map<Annotation, ScanStatisticScore> scores = sample.getWindowScores(gene);
			
			if(gene.getSize() < windowSize) {
				logger.info(gene.getName() + " is smaller than window size. Not computing single sample window enrichments.");
				singleSampleWindowEnrichmentOverGene.get(sample).put(gene, sampleWindowEnrichments);
				continue;
			}
			for(Annotation window : scores.keySet()) {
				double windowAvgCoverage = scores.get(window).getAverageCoverage();
				double enrichment = windowAvgCoverage / geneAvgCoverage;
				//logger.info(sample.getSampleName() + "\t" + gene.getName() + "\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + "\tavg_coverage=" + windowAvgCoverage + "\twindow_enrichment=" + enrichment);
				sampleWindowEnrichments.put(window, Double.valueOf(enrichment));
			}
			singleSampleWindowEnrichmentOverGene.get(sample).put(gene, sampleWindowEnrichments);
		}
	}
	
	/**
	 * Compute t statistic score for each window of gene and store scores
	 * @param gene The gene
	 */
	private void computeTStatisticWindowScores(Gene gene) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		Map<Annotation, Double> tStatisticScoresThisGene = new TreeMap<Annotation, Double>();
		if(gene.getSize() < windowSize) {
			logger.info(gene.getName() + " is smaller than window size. Not computing t statistic window scores.");
			tStatisticWindowScores.put(gene, tStatisticScoresThisGene);
			return;
		}		
		SampleData tmp = allSamples.iterator().next();
		Collection<Annotation> windows = singleSampleWindowEnrichmentOverGene.get(tmp).get(gene).keySet();
		for(Annotation window : windows) {
			double score = tStatisticWindowScore(gene, window, controlSamples, signalSamples);
			tStatisticScoresThisGene.put(window, Double.valueOf(score));
			String e = "";
			for(SampleData control : controlSamples) {
				e += control.getSampleName() + ":" + control.getWindowScores(gene).get(window).getCount() + ":" + singleSampleWindowEnrichmentOverGene.get(control).get(gene).get(window).toString() + "\t";
			}
			for(SampleData signal : controlSamples) {
				e += signal.getSampleName() + ":" + signal.getWindowScores(gene).get(window).getCount() + ":" + singleSampleWindowEnrichmentOverGene.get(signal).get(gene).get(window).toString() + "\t";
			}
			logger.info(gene.getName() + "\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + "\t" + score + "\t" + e);
		}
		tStatisticWindowScores.put(gene, tStatisticScoresThisGene);
	}
	
	/**
	 * The t statistic score for a single window
	 * T statistic between control and signal samples
	 * Statistic is positive iff mean of signal enrichments is greater than mean of control enrichments
	 * @param gene The parent gene
	 * @param window The window
	 * @param controls Control samples
	 * @param signals Signal samples
	 * @return The score for the window
	 */
	private double tStatisticWindowScore(Gene gene, Annotation window, Collection<SampleData> controls, Collection<SampleData> signals) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		List<Double> controlEnrichments = new ArrayList<Double>();
		List<Double> signalEnrichments = new ArrayList<Double>();
		for(SampleData control : controls) {
			controlEnrichments.add(singleSampleWindowEnrichmentOverGene.get(control).get(gene).get(window));
		}
		for(SampleData signal : signals) {
			signalEnrichments.add(singleSampleWindowEnrichmentOverGene.get(signal).get(gene).get(window));
		}
		return Statistics.tstat(signalEnrichments, controlEnrichments);
	}
	
	
	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#scoreWindows(java.util.Collection)
	 */
	@Override
	public void scoreWindows(Collection<Annotation> windows) {
		throw new UnsupportedOperationException("TODO");
	}
	
	/**
	 * Nominal P value calculated relative to null distribution of t statistic score
	 * Where null distribution is determined by permuting the contol/signal labels of samples
	 * @param gene The gene the window belongs to
	 * @param window The window
	 * @return The nominal P value for the score of the window
	 */
	@SuppressWarnings("unused")
	private double empiricalNominalPvalTStatisticScore(Gene gene, Annotation window) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		double windowScore = tStatisticWindowScores.get(gene).get(window).doubleValue();
		double numLess = 0;
		double numMore = 0;
		for(SamplePermutation perm : sampleIdentityPermutations) {
			double score = tStatisticWindowScore(gene, window, perm.getControls(), perm.getSignals());
			if(score > windowScore) numMore++;
			else numLess++;
		}
		return numMore / (numLess + numMore);
	}

	/**
	 * FDR for t statistic score of window
	 * @param gene The gene the window belongs to
	 * @param window The window
	 * @return The corrected P value
	 */
	@SuppressWarnings({ "static-method", "unused" })
	private double tStatisticWindowScoreFDR(Gene gene, Annotation window) {
		throw new UnsupportedOperationException("TODO");
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#filterWindows(java.util.Collection)
	 */
	@Override
	public Collection<Annotation> filterWindows(Collection<Annotation> windows) throws IOException {
		throw new UnsupportedOperationException("TODO");
 	}

	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#mergePeaks(java.util.Collection)
	 */
	@Override
	public Collection<Annotation> mergePeaks(Collection<Annotation> peaks) {
		throw new UnsupportedOperationException("TODO");
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#mergePeaks(java.util.Collection)
	 */	
	@Override
	public Annotation trimPeak(Annotation peak) {
		throw new UnsupportedOperationException("TODO");
	}


	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#writeResults(java.util.Collection, java.lang.String)
	 */
	@Override
	public void writeResults(Collection<Annotation> windows, String out) throws IOException {
		throw new UnsupportedOperationException("TODO");
	}

	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#writeResult(java.util.Collection, java.io.FileWriter)
	 */
	@Override
	public void writeResult(Collection<Annotation> windows, FileWriter writer) throws IOException {
		throw new UnsupportedOperationException("TODO");
	}

	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#setCoordinateSpace(nextgen.core.coordinatesystem.CoordinateSpace)
	 */
	@Override
	public void setCoordinateSpace(CoordinateSpace space) {
		coord = (TranscriptomeSpace) space;
	}

	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#getCoordinateSpace()
	 */
	@Override
	public CoordinateSpace getCoordinateSpace() {
		return coord;
	}
	
	
	/**
	 * Get a systematic iterator over all possible sample identity permutations
	 * @return An iterator over all permutations of sample identity (control or signal)
	 */
	@SuppressWarnings("unused")
	private Iterator<SamplePermutation> getIterAllSampleIdentityPermutations() {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		return new SystematicSamplePermutationIterator();
	}

	/**
	 * Get an iterator over a fixed number of random sample identity permutations
	 * @param numPermutations The number of permutations to get
	 * @return An iterator over sample identity permutations (control or signal)
	 */
	@SuppressWarnings("unused")
	private Iterator<SamplePermutation> getIterRandomSampleIdentityPermutations(int numPermutations) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		return new RandomSamplePermutationIterator(numPermutations);
	}
	
	/**
	 * Get an iterator over a set of sample identity permutations
	 * Either a fixed number of random permutations or all possible permutations, whichever is less
	 * @param maxNumPermutations The max number of permutations to get
	 * @return An iterator over at most the max number of sample identity permutations (control or signal)
	 */
	@SuppressWarnings("unused")
	private Iterator<SamplePermutation> getIterRandomOrAllSampleIdentityPermutations(int maxNumPermutations) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		SystematicSamplePermutationIterator systematicIter = new SystematicSamplePermutationIterator();
		if(systematicIter.getTotalNumPermutations() <= maxNumPermutations) {
			return systematicIter;
		}
		return new RandomSamplePermutationIterator(maxNumPermutations);
	}
	
	/**
	 * Get a random permutation of control and signal identities
	 * @return The random permutation
	 */
	protected SamplePermutation getOneRandomSamplePermutation() {
		ArrayList<Integer> controlPositions = new ArrayList<Integer>();
		int numControlsAssigned = 0;
		while(numControlsAssigned < numControls) {
			Integer randPos = Integer.valueOf(random.nextInt(numSamples));
			if(controlPositions.contains(randPos)) {
				continue;
			}
			controlPositions.add(randPos);
			numControlsAssigned++;
		}
		return new SamplePermutation(controlPositions);
	}

	/**
	 * Set logger levels for all samples
	 * @param level Level
	 */
	protected void setLoggerLevel(Level level) {
		logger.setLevel(level);
		for(SampleData sample : signalSamples) {
			sample.getLogger().setLevel(level);
		}
		for(SampleData sample : controlSamples) {
			sample.getLogger().setLevel(level);
		}
		expressionData.getLogger().setLevel(level);
	}
	
	private static CommandLineParser getCommandLineParser(String[] commandArgs) {
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-l", "Sample list file", true);
		p.addStringArg("-b", "Bed file of genes", true);
		p.addIntArg("-w", "Window size", false, DEFAULT_WINDOW_SIZE);
		p.addIntArg("-s", "Step size", false, DEFAULT_STEP_SIZE);
		p.addDoubleArg("-sp", "Scan P value cutoff for peak within gene", false, DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF);
		p.addDoubleArg("-cp", "Window count cutoff for peak", false, DEFAULT_PEAK_WINDOW_COUNT_CUTOFF);
		p.addBooleanArg("-batch", "Batch out peak writing by sample name and chromosome", false, false);
		p.addStringArg("-cl", "Chromosome list file for batched run", false, null);
		p.addIntArg("-m", "Memory request for batched processes", false, DEFAULT_BATCH_MEM_REQUEST);
		p.addStringArg("-bj", "Batched peak caller jar file", false, null);
		p.addStringArg("-o", "Output directory", false, null);
		p.addDoubleArg("-q", "Quantile for peak trimming by trim max contiguous algorithm", false, DEFAULT_TRIM_PEAK_QUANTILE);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.parse(commandArgs);
		return p;
	}
	
	protected static MultiSampleScanPeakCaller createFromCommandArgs(String[] commandArgs) throws IOException {
		CommandLineParser p = getCommandLineParser(commandArgs);
		String sampleListFile = p.getStringArg("-l");
		String bedFile = p.getStringArg("-b");
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-s");
		double scanPvalCutoff = p.getDoubleArg("-sp");
		double trimQuantile = p.getDoubleArg("-q");
		String chrSizeFile = p.getStringArg("-c");
		double windowCountCutoff = p.getDoubleArg("-cp");
		
		return new MultiSampleScanPeakCaller(sampleListFile, bedFile, chrSizeFile, windowSize, stepSize, scanPvalCutoff, windowCountCutoff, trimQuantile);
	}
	
	private static String commandLineBatchChrList(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getStringArg("-cl");
	}
	
	private static int commandLineBatchMemRequest(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getIntArg("-m");
	}
	
	protected static boolean commandLineHasDebugFlag(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getBooleanArg("-d");
	}
	
	private static boolean commandLineHasBatchFlag(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getBooleanArg("-batch");
	}
	
	protected static String commandLineOutDir(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getStringArg("-o");		
	}
	
	private static String commandLineBatchJar(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		String jar = p.getStringArg("-bj");
		if(jar == null) {
			throw new IllegalArgumentException("Must provide batch peak caller jar file with option -bj.");
		}
		return jar;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {

		MultiSampleScanPeakCaller m = createFromCommandArgs(args);
		
		if(commandLineHasDebugFlag(args)) {
			m.setLoggerLevel(Level.DEBUG);
		}
		
		if(commandLineHasBatchFlag(args)) {
			m.batchWriteSingleSampleScanPeaksAllSamples(args, commandLineBatchChrList(args), commandLineBatchMemRequest(args));
		} else {
			m.writeSingleSampleScanPeaksAllSamples(commandLineOutDir(args));
		}
		
		logger.info("");
		logger.info("All done.");
		
	}

	
	/**
	 * Helper class to parse file containing list of sample types and bam file names
	 * @author prussell
	 *
	 */
	private class SampleFileParser {
		
		private String EXPRESSION_LABEL = "Expression";
		private String EXPRESSION_PVAL_CUTOFF_LABEL = "Expression_pval_cutoff";
		private String EXPRESSION_AVG_COVERAGE_CUTOFF_LABEL = "Expression_avg_coverage_cutoff";
		private String CONTROL_LABEL = "Control";
		private String SIGNAL_LABEL = "Signal";
		private String FIRST_READ_TRANSCRIPTION_STRAND_LABEL = "first_read_transcription_strand";
		private String SECOND_READ_TRANSCRIPTION_STRAND_LABEL = "second_read_transcription_strand";
		private GenomeSpaceSampleData expressionSampleData;
		private ArrayList<SampleData> controlData;
		private ArrayList<SampleData> signalData;
		private String chrSizes;

		
		public SampleFileParser(String sampleFileName, String chrSizeFile) throws IOException {
			controlData = new ArrayList<SampleData>();
			signalData = new ArrayList<SampleData>();
			expressionSampleData = null;
			chrSizes = chrSizeFile;
			parseFile(sampleFileName);
		}
		
		/**
		 * Get the control datasets
		 * @return Set of control datasets
		 */
		public ArrayList<SampleData> getControlDatasets() {
			return controlData;
		}
		
		/**
		 * Get the signal datasets
		 * @return Set of signal datasets
		 */
		public ArrayList<SampleData> getSignalDatasets() {
			return signalData;
		}
		
		/**
		 * Get the expression dataset
		 * @return Expression dataset
		 */
		public GenomeSpaceSampleData getExpressionData() {
			return expressionSampleData;
		}
		
		/**
		 * Parse the sample file and populate data sets
		 * @throws IOException
		 */
		private void parseFile(String sampleFileName) throws IOException {
			boolean foundExpressionData = false;
			FileReader r = new FileReader(sampleFileName);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			
			if(!b.ready()) {
				crashWithHelpMessage();
			}
			
			String pvalLine = b.readLine();
			s.parse(pvalLine);
			boolean expByScanPval = false;
			if(s.getFieldCount() != 2) {
				crashWithHelpMessage();
			}
			if(s.asString(0).equals(EXPRESSION_AVG_COVERAGE_CUTOFF_LABEL)) {
				expByScanPval = false;
			} else if (s.asString(0).equals(EXPRESSION_PVAL_CUTOFF_LABEL)) {
				expByScanPval = true;
			} else {
				crashWithHelpMessage();
			}
			double cutoff = s.asDouble(1);
			
			if(!b.ready()) {
				crashWithHelpMessage();
			}
			
			while(b.ready()) {
				
				String line = b.readLine();
				s.parse(line);
				
				if(s.getFieldCount() == 0) continue;
				
				String label = s.asString(0);
				String bamFile = s.asString(1);

				if(label.equals(EXPRESSION_LABEL)) {
					if(s.getFieldCount() > 2) crashWithHelpMessage();
					if(foundExpressionData) {
						crashWithHelpMessage();
					}
					logger.info("Creating sample data object for gene expression from bam file " + bamFile);
					GenomeSpaceSampleData sample = new GenomeSpaceSampleData(bamFile, chrSizes, genes, windowSize, stepSize, cutoff);
					expressionSampleData = sample;
					foundExpressionData = true;
					continue;
				}
				
				if(label.equals(CONTROL_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					logger.info("Creating sample data object for bam file " + bamFile);
					boolean firstReadTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					SampleData sample = new SampleData(bamFile, firstReadTranscriptionStrand, genes, windowSize, stepSize, cutoff, expByScanPval);
					controlData.add(sample);
					continue;
				}
				
				if(label.equals(SIGNAL_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					logger.info("Creating sample data object for bam file " + bamFile);
					boolean firstReadTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					SampleData sample = new SampleData(bamFile, firstReadTranscriptionStrand, genes, windowSize, stepSize, cutoff, expByScanPval);
					signalData.add(sample);
					continue;
				}
				
				crashWithHelpMessage();
				
			}
			
			r.close();
			b.close();
			
			if(!foundExpressionData) {
				crashWithHelpMessage();
			}
			
		}
		
		/**
		 * Crash and print help message if sample file is invalid
		 */
		private void crashWithHelpMessage() {
			logger.error("");
			logger.error("**********");
			logger.error("");
			logger.error("Sample file not valid.");
			logger.error("");
			logger.error("First line must be:");
			logger.error(EXPRESSION_PVAL_CUTOFF_LABEL + "\t<pval_cutoff>");
			logger.error("-OR-");
			logger.error(EXPRESSION_AVG_COVERAGE_CUTOFF_LABEL + "\t<avg_depth_cutoff>");
			logger.error("");
			logger.error("Exactly one line must be of the form:");
			logger.error(EXPRESSION_LABEL + "\t<bam_file_name>");
			logger.error("");
			logger.error("Each additional line must be of the form:");
			logger.error(CONTROL_LABEL + "\t<bam_file_name>\t" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("- or -");
			logger.error(SIGNAL_LABEL + "\t<bam_file_name>" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("");
			logger.error("**********");
			logger.error("");
			throw new IllegalArgumentException("Sample file not valid.");
		}
		
		
	}

	
	/**
	 * Iterator over randomly generated sample identity permutations
	 * @author prussell
	 *
	 */
	private class RandomSamplePermutationIterator implements Iterator<SamplePermutation> {
		
		private int numPermutations;
		private int nextPosition;
		
		/**
		 * Construct with number of permutations to generate
		 * @param numToGenerate
		 */
		public RandomSamplePermutationIterator(int numToGenerate) {
			numPermutations = numToGenerate;
			nextPosition = 0;
		}
		
		/**
		 * Get the number of permutations to generate
		 * @return Number of permutations to generate
		 */
		@SuppressWarnings("unused")
		public int getNumToGenerate() {
			return numPermutations;
		}

		@Override
		public boolean hasNext() {
			return nextPosition < numPermutations;
		}

		@Override
		public SamplePermutation next() {
			nextPosition++;
			return getOneRandomSamplePermutation();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("TODO");
		}
		
	}
	
	/**
	 * Systematic iterator over all possible sample identity permutations
	 * @author prussell
	 *
	 */
	private class SystematicSamplePermutationIterator implements Iterator<SamplePermutation> {
		
		private int n;
		private int r;
		private int[] currentOneBasedPositions;
		private int[] nextOneBasedPositions;
		private int[] lastOneBasedPositions;
		
		/**
		 * Initialize first set of "control" positions to be 0,1,...,numControls
		 * Set last set to n-r,n-r+1,...,n-1
		 */
		public SystematicSamplePermutationIterator() {
			n = numSamples;
			r = numControls;
			nextOneBasedPositions = new int[r];
			for(int i = 0; i < r; i++) {
				nextOneBasedPositions[i] = algorithmI(i);
				lastOneBasedPositions[i] = algorithmI(n - r + i);
			}
		}
		
		/**
		 * Translate zero based position to one based position from combination generation algorithm
		 * @param currentPosition Zero based position
		 * @return Position in combination generation algorithm
		 */
		private int algorithmI(int currentPosition) {
			return currentPosition + 1;
		}
		
		/**
		 * Translate one based position from combination generation algorithm to zero based position
		 * @param algorithmI Position in combination generation algorithm
		 * @return Zero based position
		 */
		private int position(int algorithmI) {
			return algorithmI - 1;
		}
		
		/**
		 * Get the total number of all possible permutations
		 * @return The total number of permutations
		 */
		public long getTotalNumPermutations() {
			return MathUtil.binomialCoefficient(numSamples, numControls);
		}

		@Override
		public boolean hasNext() {
			return currentOneBasedPositions != lastOneBasedPositions;
		}

		@Override
		public SamplePermutation next() {
			
			currentOneBasedPositions = nextOneBasedPositions;
			// Increment next positions according to lexicographic order
			int k = 0;
			for(int i = 0; i < r; i++) {
				if(nextOneBasedPositions[i] < n - r + algorithmI(i)) {
					k = i;
				}
			}
			nextOneBasedPositions[k]++;
			for(int i = k + 1; i < r; i++) {
				nextOneBasedPositions[i] = nextOneBasedPositions[i-1]+1;
			}
			
			ArrayList<Integer> rtrnPositions = new ArrayList<Integer>();
			for(int i = 0; i < r; i++) {
				rtrnPositions.add(new Integer(position(currentOneBasedPositions[i])));
			}
			
			return new SamplePermutation(rtrnPositions);
			
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("TODO");
		}
		
	}
	
	/**
	 * A permutation of sample identities (control or signal)
	 * @author prussell
	 *
	 */
	private class SamplePermutation {
		
		private Collection<SampleData> controls;
		private Collection<SampleData> signals;
		private ArrayList<Integer> controlPositionsInRealSampleList;
		private ArrayList<Integer> signalPositionsInRealSampleList;
		
		/**
		 * Construct with positions of "control" and "signal" samples from the real sample list
		 * @param labeledControlsPositionsInRealSampleList Positions in real sample list of samples to label as controls
		 */
		public SamplePermutation(ArrayList<Integer> labeledControlsPositionsInRealSampleList) {
			controlPositionsInRealSampleList = labeledControlsPositionsInRealSampleList;
			signalPositionsInRealSampleList = getSignalPositionsInRealSampleList();
			controls = new ArrayList<SampleData>();
			for(int i=0; i < controlPositionsInRealSampleList.size(); i++) {
				controls.add(allSamples.get(controlPositionsInRealSampleList.get(i).intValue()));
			}
			signals = new ArrayList<SampleData>();
			for(int i=0; i < signalPositionsInRealSampleList.size(); i++) {
				signals.add(allSamples.get(signalPositionsInRealSampleList.get(i).intValue()));
			}
		}
		
		/**
		 * Get positions of "signal" labeled samples in real sample list
		 * Based on the samples labeled "control"
		 * @return Positions of "signal" labeled samples in real sample list
		 */
		private ArrayList<Integer> getSignalPositionsInRealSampleList() {
			ArrayList<Integer> rtrn = new ArrayList<Integer>();
			for(int i=0; i<numSamples; i++) {
				if(!controlPositionsInRealSampleList.contains(Integer.valueOf(i))) {
					rtrn.add(Integer.valueOf(i));
				}
			}
			return rtrn;
		}
		
		/**
		 * Get the samples labeled as controls in this permutation
		 * @return The controls in this permutation
		 */
		public Collection<SampleData> getControls() {
			return controls;
		}
		
		/**
		 * Get the samples labeled as signals in this permutation
		 * @return The signals in this permutation
		 */
		public Collection<SampleData> getSignals() {
			return signals;
		}
		
		
	}


}

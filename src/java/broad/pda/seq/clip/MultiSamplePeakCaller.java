/**
 * Call protein binding sites on RNA with Protect-seq data from multiple control and multiple signal samples
 */
package broad.pda.seq.clip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import broad.core.math.MathUtil;
import broad.core.math.Statistics;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.alignment.Alignment;
import nextgen.core.analysis.PeakCaller;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.job.Job;
import nextgen.core.job.JobUtils;
import nextgen.core.job.LSFJob;
import nextgen.core.job.OGSJob;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.MultiScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.pipeline.OGSUtils;
import nextgen.core.pipeline.Scheduler;
import nextgen.core.utils.AlignmentUtils;
import nextgen.core.utils.AnnotationUtils;
import nextgen.core.feature.GeneWindow;

/**
 * This class finds significant peaks using any/all of the scores included in MultiScore class.
 * Note: t-statistic analysis not currently performed, but methods to implement are included.
 * @author shari
 */
public class MultiSamplePeakCaller implements PeakCaller {
	
	private TranscriptomeSpace coord;
	protected Map<String, Collection<Gene>> genes;
	private Map<Gene, Map<Annotation, Double>> tStatisticWindowScores;
	private Map<SampleData, Map<Gene, Map<Annotation, Double>>> singleSampleWindowEnrichmentOverGene;
	private Map<SampleData,Map<Gene, Map<String, Collection<Gene>>>> singleSampleScanPeaks;
	protected Map<SampleData, Map<Gene, Map<Annotation, MultiScore>>> windowScores;
	protected Map<SampleData, Map<Gene, Map<Annotation, MultiScore>>> finalWindowScores;
	private Map<SampleData,WindowProcessor<MultiScore>> processors;
	protected GenomeSpaceSampleData expressionData;
	protected ArrayList<SampleData> controlSamples;
	protected ArrayList<SampleData> signalSamples;
	protected ArrayList<SampleData> allSamples;
	private SampleData ctrl;
	protected static Logger logger = Logger.getLogger(MultiSamplePeakCaller.class.getName());
	protected int windowSize;
	protected int stepSize;
	protected ArrayList<String> scoresToUse;
	private static String COMMENT_CHAR = "#";
	private static int DEFAULT_WINDOW_SIZE = 20;
	private static int DEFAULT_STEP_SIZE = 1;
	private static double DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF = 0.001;
	private static double DEFAULT_PEAK_WINDOW_COUNT_CUTOFF = 10;
	private static double DEFAULT_TRIM_PEAK_QUANTILE = 0.6;
	private static int DEFAULT_BATCH_MEM_REQUEST = 8;
	private static double DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF = 0.01;
	private static boolean DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND = false;
	private static double DEFAULT_PEAK_MAX_PCT_DUPLICATES = 0.5;
	private static boolean DEFAULT_FILTER_BY_STRAND = true;
	private static boolean DEFAULT_EXTRA_FIELDS = false;
	private static String BINOMIAL_KEY = "binomial";
	private static String SCAN_STAT_KEY = "scan_statistic";
	private static String[] DEFAULT_SCORES = { BINOMIAL_KEY, SCAN_STAT_KEY };
	private boolean filterByStrand;
	private boolean extraFields;
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
	protected String sizeFile;
	private boolean firstReadTranscriptionStrand;
	private double peakCutoffMaxReplicatePct;
	private Scheduler scheduler;
	private static Session drmaaSession;
	private static int RGB_RED_WITH_GENE = 106;
	private static int RGB_GREEN_WITH_GENE = 7;
	private static int RGB_BLUE_WITH_GENE = 205;
	private static int RGB_RED_AGAINST_GENE = 218;
	private static int RGB_GREEN_AGAINST_GENE = 2;
	private static int RGB_BLUE_AGAINST_GENE = 38;
	private static int RGB_RED_UNKNOWN = 62;
	private static int RGB_GREEN_UNKNOWN = 63;
	private static int RGB_BLUE_UNKNOWN = 73;
	private Map<SampleData, FileWriter> expressionFilterRejectWriters;
	private static String REJECT_FILE_PREFIX_EXPRESSION_FILTER = "expression_filter";
	private Map<SampleData, FileWriter> windowCountRejectFileWriters;
	private static String REJECT_FILE_PREFIX_WINDOW_COUNT = "window_count_filter";
	private Map<SampleData,Map<String, FileWriter>> windowScanPvalAllFragmentsRejectFileWriters;
	private static String REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_ALL_FRAGMENTS = "window_scan_pval_all_fragments_filter";
	private Map<SampleData,Map<String, FileWriter>> windowScanPvalWithFragmentLengthFilterRejectFileWriters;
	private static String REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_WITH_FRAGMENT_LENGTH = "fragment_length_filter";
	private Map<SampleData,Map<String, FileWriter>> peakScanPvalRejectWriters;
	private static String REJECT_FILE_PREFIX_PEAK_SCAN_PVAL = "peak_scan_pval_filter";
	private Map<SampleData,Map<String, FileWriter>> duplicateRejectWriters;
	private static String REJECT_FILE_PREFIX_DUPLICATE = "duplication_filter";
	private Map<SampleData,Map<String, FileWriter>> strandRejectWriters;
	private static String REJECT_FILE_PREFIX_STRAND = "strand_filter";
	private Map<String,ArrayList<FileWriter>> allRejectFileWriters;
	protected static String FILTER_REJECT_DIR = "filter_rejects";

	// Constructors
	/**
	 * Instantiate from another MultiSamplePeakCaller
	 * @param other 
	 * @throws IOException
	 * @throws DrmaaException 
	 */
	protected MultiSamplePeakCaller(MultiSamplePeakCaller other) throws IOException, DrmaaException {
		this(other.sampleFile, other.bedAnnotationFile, other.sizeFile, other.windowSize, other.stepSize);
		copyParameters(other);
	}
	
	/**
	 * Instantiate with default parameters
	 * @param sampleListFile File containing sample list
	 * @param bedFile Bed gene annotation
	 * @throws IOException
	 * @throws DrmaaException 
	 */
	@SuppressWarnings("unused")
	private MultiSamplePeakCaller(String sampleListFile, String bedFile, String chrSizeFile) throws IOException, DrmaaException {
		this(sampleListFile, bedFile, chrSizeFile, DEFAULT_WINDOW_SIZE, DEFAULT_STEP_SIZE);
	}
	
	/**
	 * Instantiate using file of sample information
	 * @param sampleListFile File containing sample list
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @param window Window size
	 * @param step Step size
	 * @throws IOException
	 * @throws DrmaaException 
	 */
	public MultiSamplePeakCaller(String sampleListFile, String bedFile, String chrSizeFile, int window, int step) throws IOException, DrmaaException {
		
		sampleFile = sampleListFile;
		bedAnnotationFile = bedFile;
		sizeFile = chrSizeFile;
		windowSize = window;
		stepSize = step;
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		coord = new TranscriptomeSpace(genes);
		random = new Random();
		scoresToUse = new ArrayList<String>(Arrays.asList(DEFAULT_SCORES));
		
		if(sampleFile != null) {
			initializeSamplesFromSampleListFile();
			initializeScoreMaps();
			initializeProcessors();
		}
		
		trimQuantile = DEFAULT_TRIM_PEAK_QUANTILE;
		peakWindowCountCutoff = DEFAULT_PEAK_WINDOW_COUNT_CUTOFF;
		peakWindowScanPvalCutoff = DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF;
		firstReadTranscriptionStrand = DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND;
		peakCutoffMaxReplicatePct = DEFAULT_PEAK_MAX_PCT_DUPLICATES;
		filterByStrand = DEFAULT_FILTER_BY_STRAND;
		extraFields = DEFAULT_EXTRA_FIELDS;
		
	}
	
	/**
	 * Instantiate with a single sample and an expression sample and default paramters
	 * @param expressionBamFile Bam file for expression sample
	 * @param signalBamFile Bam file for signal sample
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @throws IOException
	 * @throws DrmaaException 
	 */
	public MultiSamplePeakCaller(String expressionBamFile, String signalBamFile, String bedFile, String chrSizeFile) throws IOException, DrmaaException {
		this(expressionBamFile, signalBamFile, bedFile, chrSizeFile, DEFAULT_WINDOW_SIZE, DEFAULT_STEP_SIZE);
	}
	
	/**
	 * Instantiate with a single signal sample and an expression sample
	 * @param expressionBamFile Bam file for expression sample
	 * @param signalBamFile Bam file for signal sample
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @param window Window size
	 * @param step Step size
	 * @throws IOException
	 * @throws DrmaaException 
	 */
	public MultiSamplePeakCaller(String expressionBamFile, String signalBamFile, String bedFile, String chrSizeFile, int window, int step) throws IOException, DrmaaException {
		this(null, bedFile, chrSizeFile, window, step);
		initializeWithSingleSample(expressionBamFile, signalBamFile, chrSizeFile);
		initializeScoreMaps();
		initializeProcessors();
	}
	
	// Methods to initalize MultiSamplePeakCaller
	// Populate signal and control data, initialize score hashes, initialize processors, and initialize filter reject writers
	
	/**
	 * Populate MultiSamplePeakCaller with single sample and expression data
	 * @param expressionBamFile
	 * @param signalBamFile
	 * @param chrSizeFile
	 * @throws IOException
	 */
	private void initializeWithSingleSample(String expressionBamFile, String signalBamFile, String chrSizeFile) throws IOException {
		controlSamples = new ArrayList<SampleData>();
		numControls = controlSamples.size();
		signalSamples = new ArrayList<SampleData>();
		SampleData signalSample = new SampleData(signalBamFile, firstReadTranscriptionStrand, genes, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF, true, false);
		signalSamples.add(signalSample);
		numSignals = signalSamples.size();
		expressionData = new GenomeSpaceSampleData(expressionBamFile, chrSizeFile, genes, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF, true);
		allSamples = new ArrayList<SampleData>();
		allSamples.addAll(controlSamples);
		allSamples.addAll(signalSamples);
		numSamples = allSamples.size();
		ctrl = new SampleData(expressionBamFile, firstReadTranscriptionStrand, genes, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF, true, false);
	}
	
	/**
	 * Populate MultiSamplePeakCaller from sample list file
	 * @throws IOException
	 */
	private void initializeSamplesFromSampleListFile() throws IOException {
		SampleFileParser p = new SampleFileParser(sampleFile);
		controlSamples = p.getControlDatasets();
		numControls = controlSamples.size();
		signalSamples = p.getSignalDatasets();
		expressionData = p.getExpressionData();
		ctrl = p.getCtrlData();
		numSignals = signalSamples.size();
		allSamples = new ArrayList<SampleData>();
		allSamples.addAll(controlSamples);
		allSamples.addAll(signalSamples);
		numSamples = allSamples.size();
	}
	
	/**
	 * Initialize hashes to store: 
		 * window enrichments over gene background for each gene in each sample
		 * Significant peaks for each sample for each score for each gene
		 * t-statistic scores for window of each gene
	 */
	private void initializeScoreMaps() {
		
		tStatisticWindowScores = new TreeMap<Gene, Map<Annotation, Double>>();
		singleSampleWindowEnrichmentOverGene = new HashMap<SampleData, Map<Gene, Map<Annotation, Double>>>();
		singleSampleScanPeaks = new HashMap<SampleData,Map<Gene,Map<String, Collection<Gene>>>>();
		windowScores = new HashMap<SampleData, Map<Gene, Map<Annotation, MultiScore>>>();
		finalWindowScores = new HashMap<SampleData, Map<Gene, Map<Annotation, MultiScore>>>();
		
		for(SampleData sample : allSamples) {
			singleSampleScanPeaks.put(sample,new HashMap<Gene, Map<String, Collection<Gene>>>());
			singleSampleWindowEnrichmentOverGene.put(sample, new TreeMap<Gene, Map<Annotation, Double>>());
			windowScores.put(sample, new HashMap<Gene,Map<Annotation,MultiScore>>());
			finalWindowScores.put(sample, new HashMap<Gene,Map<Annotation,MultiScore>>());
		}
	}
	
	/**
	 * Initialize MultiScore processor for each sample
	 * @throws IOException
	 */
	private void initializeProcessors() throws IOException {
		processors = new HashMap<SampleData,WindowProcessor<MultiScore>>();
		for(SampleData sample : allSamples) {
			if (!processors.containsKey(sample)) {
				try {
					WindowProcessor<MultiScore> p = new MultiScore.Processor(sample.getData(),ctrl.getData(),false);
					processors.put(sample, p);
				} catch(Exception e) {
					try {
						sample.getData();
					} catch(NullPointerException f) {
						logger.debug("sample data is null");
					}
					try {
						ctrl.getData();
					} catch(NullPointerException g) {
						logger.debug("ctrl data is null");
					}
				}
			}
		}
	}
	
	/**
	 * Initialize filter reject writers
	 * @param commonSuffix suffix for output bed
	 * @param outDir
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	protected void initializeFilterRejectWriters(String commonSuffix, String outDir) throws IOException {
		File dir = new File(outDir);
		boolean madeDir = dir.mkdir();
		if(!dir.exists()) {
			throw new IllegalStateException("Could not make directory " + outDir + ".");
		}
		for (String scoreName : scoresToUse) {
			File scoreDir = new File(outDir + "/" + scoreName);
			boolean madeScoreDir = scoreDir.mkdir();
			if(!scoreDir.exists()) {
				throw new IllegalStateException("Could not make directory " + scoreDir + ".");
			}
		}
		logger.info("");
		logger.info("Writing regions rejected by filters to files in directory " + outDir);
		expressionFilterRejectWriters = new HashMap<SampleData, FileWriter>();
		windowCountRejectFileWriters = new HashMap<SampleData, FileWriter>();
		windowScanPvalAllFragmentsRejectFileWriters = new HashMap<SampleData,Map<String,FileWriter>>();
		windowScanPvalWithFragmentLengthFilterRejectFileWriters = new HashMap<SampleData,Map<String,FileWriter>>();
		peakScanPvalRejectWriters = new HashMap<SampleData,Map<String,FileWriter>>();
		duplicateRejectWriters = new HashMap<SampleData,Map<String,FileWriter>>();
		strandRejectWriters = new HashMap<SampleData,Map<String,FileWriter>>();
		allRejectFileWriters = new HashMap<String,ArrayList<FileWriter>>();
		
		for(SampleData sample : allSamples) {
			expressionFilterRejectWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_EXPRESSION_FILTER + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			windowCountRejectFileWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_WINDOW_COUNT + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
		}
		
		for(SampleData sample : allSamples) {
			windowScanPvalAllFragmentsRejectFileWriters.put(sample,new HashMap<String, FileWriter>());
			windowScanPvalWithFragmentLengthFilterRejectFileWriters.put(sample,new HashMap<String, FileWriter>());
			peakScanPvalRejectWriters.put(sample,new HashMap<String, FileWriter>());
			duplicateRejectWriters.put(sample,new HashMap<String, FileWriter>());
			strandRejectWriters.put(sample,new HashMap<String, FileWriter>());
		
			for (String scoreName : scoresToUse) {
				windowScanPvalAllFragmentsRejectFileWriters.get(sample).put(scoreName, new FileWriter(outDir + "/" + scoreName + "/" + REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_ALL_FRAGMENTS + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
				windowScanPvalWithFragmentLengthFilterRejectFileWriters.get(sample).put(scoreName, new FileWriter(outDir + "/" + scoreName + "/" + REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_WITH_FRAGMENT_LENGTH + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
				peakScanPvalRejectWriters.get(sample).put(scoreName, new FileWriter(outDir + "/" + scoreName + "/" + REJECT_FILE_PREFIX_PEAK_SCAN_PVAL + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
				duplicateRejectWriters.get(sample).put(scoreName, new FileWriter(outDir + "/" + scoreName + "/" + REJECT_FILE_PREFIX_DUPLICATE + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
				strandRejectWriters.get(sample).put(scoreName, new FileWriter(outDir + "/" + scoreName + "/" + REJECT_FILE_PREFIX_STRAND + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
				
				if (!allRejectFileWriters.containsKey(scoreName)){
					allRejectFileWriters.put(scoreName,new ArrayList<FileWriter>());
				}
				allRejectFileWriters.get(scoreName).add(expressionFilterRejectWriters.get(sample));
				allRejectFileWriters.get(scoreName).add(windowCountRejectFileWriters.get(sample));
				allRejectFileWriters.get(scoreName).add(windowScanPvalAllFragmentsRejectFileWriters.get(sample).get(scoreName));
				allRejectFileWriters.get(scoreName).add(windowScanPvalWithFragmentLengthFilterRejectFileWriters.get(sample).get(scoreName));
				allRejectFileWriters.get(scoreName).add(peakScanPvalRejectWriters.get(sample).get(scoreName));
				allRejectFileWriters.get(scoreName).add(duplicateRejectWriters.get(sample).get(scoreName));
				allRejectFileWriters.get(scoreName).add(strandRejectWriters.get(sample).get(scoreName));
			}
		}
	}
	
	/**
	 * Close all reject filter writers for each score
	 * @throws IOException
	 */
	protected void closeFilterRejectWriters() throws IOException {
		for(String scoreName : scoresToUse) {
			for(FileWriter w : allRejectFileWriters.get(scoreName)) {
				w.close();
			}
		}
	}
	
	// Methods to set parameters
	
	/**
	 * Copy parameters from other MultiSamplePeakCaller
	 * @param other
	 */
	private void copyParameters(MultiSamplePeakCaller other) {
		setPeakTrimQuantile(other.trimQuantile);
		setPeakWindowCountCutoff(other.peakWindowCountCutoff);
		setPeakWindowScanPvalCutoff(other.peakWindowScanPvalCutoff);
		setPeakCutoffMostCommonReplicate(other.peakCutoffMaxReplicatePct);
		setFirstReadTranscriptionStrand(other.firstReadTranscriptionStrand);
		setExpressionScanPvalueCutoff(other.expressionData.getExpressionScanPvalueCutoff());
		setFilterByStrand(other.filterByStrand);
		setExtraFields(other.extraFields);
		setScoreList(other.scoresToUse);
	}
	
	/**
	 * Set whether to filter by strand info
	 * @param useStrandFilter
	 */
	public void setFilterByStrand(boolean useStrandFilter) { filterByStrand = useStrandFilter; }
	
	/**
	 * Set whether to print extra fields in output beds
	 * @param useExtraFields
	 */
	public void setExtraFields(boolean useExtraFields) { extraFields = useExtraFields; }
	
	/**
	 * Set which scores to compute peaks for
	 * @param scoreList List of keys of scores
	 */
	public void setScoreList(ArrayList<String> scoreList) { 
		scoresToUse = scoreList; 
		String listString = "";
		for(String score : scoresToUse) {
			listString += score + " ";
		}
		logger.info("USING SCORES: " + listString);
	}
	
	/**
	 * Set cutoff for the percentage of fragments overlapping a peak that come from the most common replicate fragment
	 * @param maxPct The max percentage
	 */
	public void setPeakCutoffMostCommonReplicate(double maxPct) { peakCutoffMaxReplicatePct = maxPct; }
	
	/**
	 * Set the scheduler
	 * @param sched Scheduler
	 * @throws DrmaaException 
	 */
	public void setScheduler(Scheduler sched) throws DrmaaException { 
		scheduler = sched;
		if(scheduler.equals(Scheduler.OGS)) {
			drmaaSession = OGSUtils.getDrmaaSession();
		}
	}
	
	/**
	 * Set quantile for trim max contiguous algorithm for peak calling
	 * @param trimPeakQuantile The quantile of data values
	 */
	public void setPeakTrimQuantile(double trimPeakQuantile) { trimQuantile = trimPeakQuantile; }
	
	/**
	 * Set minimum number of fragments overlapping a window to be considered for possible peak
	 * @param peakCountCutoff Minimum fragment count for each window
	 */
	public void setPeakWindowCountCutoff(double peakCountCutoff) { peakWindowCountCutoff = peakCountCutoff; }
	
	/**
	 * Set scan P value cutoff for window within transcript to be considered for possible peak
	 * @param peakScanPvalCutoff Max scan P value
	 */
	public void setPeakWindowScanPvalCutoff(double peakScanPvalCutoff) { peakWindowScanPvalCutoff = peakScanPvalCutoff; }
	
	/**
	 * Set whether read 1 is transcription strand
	 * @param firstReadIsTranscriptionStrand True if read 1 is transcription strand
	 */
	public void setFirstReadTranscriptionStrand(boolean firstReadIsTranscriptionStrand) { firstReadTranscriptionStrand = firstReadIsTranscriptionStrand; }
	
	/**
	 * Set genome wide scan P value cutoff for expression of transcript
	 * @param expressionScanPvalCutoff P value cutoff for transcript expression against genomic background
	 */
	public void setExpressionScanPvalueCutoff(double expressionScanPvalCutoff) {
		for(SampleData sample : allSamples) {
			sample.setExpressionScanPvalueCutoff(expressionScanPvalCutoff);
		}
	}

	// Methods to get expression information
	
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
	
	// Batch write peak files - submit job for each chromosome for each sample
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSamplePeaksAllSamples(String[] commandArgs) throws IOException, InterruptedException, DrmaaException {
		batchWriteSingleSamplePeaksAllSamples(commandArgs, null, DEFAULT_BATCH_MEM_REQUEST);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSamplePeaksAllSamples(String[] commandArgs, String chrListFile) throws IOException, InterruptedException, DrmaaException {
		batchWriteSingleSamplePeaksAllSamples(commandArgs, chrListFile, DEFAULT_BATCH_MEM_REQUEST);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSamplePeaksAllSamples(String[] commandArgs, int memRequestGb) throws IOException, InterruptedException, DrmaaException {
		batchWriteSingleSamplePeaksAllSamples(commandArgs, null, memRequestGb);
	}
	
	private void batchWriteSingleSamplePeaksAllSamples(String[] commandArgs, String chrListFile, int memRequestGb) throws IOException, InterruptedException, DrmaaException {
		
		logger.info("");
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
					r.close();
					b.close();
					throw new IllegalArgumentException("Chromosome name " + line + " not recognized.");
				}
				chrs.add(line);
			}
			r.close();
			b.close();
		}
		
		String jar = commandLineBatchJar(commandArgs);
		ArrayList<Job> jobs = new ArrayList<Job>();
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
				String[] batchedCmmdArgs = BatchedMultiSamplePeakCaller.extendSuperArgsForSampleAndChr(commandArgs, sample.getSampleName(), chr);
				String args = "";
				for(int i=0; i < batchedCmmdArgs.length; i++) {
					args += batchedCmmdArgs[i] + " ";
				}
				String cmmd = "java -jar -Xmx" + xmx + "g -Xms" + xms + "g -Xmn" + xmn + "g " + jar + " " + args;
				logger.info("Running command: " + cmmd);
				String jobID = sample.getSampleName() + "_" + chr + "_" + Long.valueOf(System.currentTimeMillis()).toString();
				switch(scheduler) {
				case LSF:
					LSFJob lsfJob = new LSFJob(Runtime.getRuntime(), jobID, cmmd, outDir + "/" + jobID + ".bsub", "week", memRequestGb);
					jobs.add(lsfJob);
					cmmds.put(jobID, cmmd);
					logger.info("LSF job ID is " + jobID + ".");
					// Submit job
					lsfJob.submit();
					break;
				case OGS:
					OGSJob ogsJob = new OGSJob(drmaaSession, cmmd);
					jobs.add(ogsJob);
					ogsJob.submit();
					logger.info("OGS job ID is " + ogsJob.getID());
					cmmds.put(ogsJob.getID(), cmmd);
					break;
				default:
					throw new IllegalStateException("Case fall through in switch on scheduler value");
				}
			}
		}
		
		logger.info("");
		logger.info("Waiting for jobs to finish...");
		JobUtils.waitForAll(jobs);
		
		logger.info("\nAll jobs finished.\n");
		
	}
	
	// Methods to scan windows in all genes, identify significant peaks, and write to output beds

	/**
	 * Write scan peaks for all samples to separate bed files
	 * @throws IOException
	 */
	private void writeSingleSamplePeaksAllSamples(String outDir) throws IOException {
		logger.info("");
		logger.info("Writing single sample scan peaks for each sample...");
		File o = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = o.mkdir();
		if(!o.exists()) {
			throw new IOException("Could not create directory " + outDir);
		}
		for(SampleData signal : signalSamples) {
			Map<String,String> outfiles = getPeakBedFileName(signal, outDir, null);
			writeSingleSamplePeaks(signal, outfiles);
		}		
		for(SampleData control : controlSamples) {
			Map<String,String> outfiles = getPeakBedFileName(control, outDir, null);
			writeSingleSamplePeaks(control, outfiles);
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
	private void writeSingleSamplePeaks(SampleData sample, Map<String,String> outFiles) throws IOException {
		writeSingleSamplePeaks(sample, outFiles, null);
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
	protected void writeSingleSamplePeaks(SampleData sample, Map<String,String> outFiles, String chrName) throws IOException {
		logger.info("Writing single sample scan peaks for sample " + sample.getSampleName() + "...");
		for (String scoreName : scoresToUse) {
			logger.info(scoreName + " to file " + outFiles.get(scoreName));
		}
		Map<String,FileWriter> w = new HashMap<String,FileWriter>();
		for (String scoreName : scoresToUse) {
			w.put(scoreName,new FileWriter(outFiles.get(scoreName)));
		}
		for(String chr : genes.keySet()) {
			if(chrName != null) {
				if(!chr.equals(chrName)) {
					continue;
				}
			}
			for(Gene gene : genes.get(chr)) {
				Map<String,Collection<Gene>> peaks = getSingleSamplePeaks(sample, gene);
				for (String scoreName : scoresToUse) {
					for(Gene window : peaks.get(scoreName)) {
						GeneWindow geneWindow = new GeneWindow(window); 
						int r = RGB_RED_UNKNOWN;
						int g = RGB_GREEN_UNKNOWN;
						int b = RGB_BLUE_UNKNOWN;
						if(window.getOrientation().equals(Strand.UNKNOWN)) {
							String name = window.getName();
							name += "_STRAND_UNKNOWN";
							geneWindow.setName(name);						
						}
						if(window.getOrientation().equals(gene.getOrientation()) && !window.getOrientation().equals(Strand.UNKNOWN)) {
							r = RGB_RED_WITH_GENE;
							g = RGB_GREEN_WITH_GENE;
							b = RGB_BLUE_WITH_GENE;
						}
						if(!window.getOrientation().equals(gene.getOrientation()) && !window.getOrientation().equals(Strand.UNKNOWN)) {
							String name = window.getName();
							name += "_STRAND_AGAINST_GENE";
							geneWindow.setName(name);
							r = RGB_RED_AGAINST_GENE;
							g = RGB_GREEN_AGAINST_GENE;
							b = RGB_BLUE_AGAINST_GENE;
						}
						geneWindow.setBedScore(window.getScore());
						if (extraFields) {
							w.get(scoreName).write(geneWindow.toBED(true, r, g, b) + "\n");
						} else {
							w.get(scoreName).write(geneWindow.toBED(r, g, b) + "\n");
						}
				}
				}
			}
		}
		for (String scoreName : scoresToUse) {
			w.get(scoreName).close();
		}
		
		logger.info("Done writing scan peaks for sample " + sample.getSampleName() + ".");
	}
	
	/**
	 * Get scan peaks for the sample and the gene
	 * @param sample The sample
	 * @param gene The gene
	 * @return Significant scan peaks
	 * @throws IOException 
	 */
	public Map<String,Collection<Gene>> getSingleSamplePeaks(SampleData sample, Gene gene) throws IOException {
		if(singleSampleScanPeaks.get(sample).containsKey(gene)) {
			return singleSampleScanPeaks.get(sample).get(gene);
		} 
		identifySingleSamplePeaks(sample, gene, expressionFilterRejectWriters.get(sample));
		return singleSampleScanPeaks.get(sample).get(gene);
	}
	
	/**
	 * Identify significant windows for each score within the sample and the gene, 
	 * Merge overlapping windows, trim merged windows, filter again by new scores, filter by strand and duplicates,
	 * Stores final peaks of given sample for each test in the gene in singleSampleScanPeaks
	 * @param sample The sample
	 * @param gene The gene
	 * @throws IOException
	 */
	private void identifySingleSamplePeaks(SampleData sample, Gene gene, FileWriter rejectFileWriterExpression) throws IOException {
		
		Map<String,TreeSet<Annotation>> finalPeaks = new HashMap<String,TreeSet<Annotation>>();
		Map<String,TreeSet<Gene>> rtrnPeaks = new HashMap<String,TreeSet<Gene>>();
		Map<String,TreeSet<Annotation>> significantWindows = new HashMap<String,TreeSet<Annotation>>();
		for (String scoreName : scoresToUse) {
			finalPeaks.put(scoreName,new TreeSet<Annotation>());
			rtrnPeaks.put(scoreName, new TreeSet<Gene>());
			significantWindows.put(scoreName, new TreeSet<Annotation>());
		}
		TranscriptomeSpaceAlignmentModel data = sample.getData();
		
		singleSampleScanPeaks.get(sample).put(gene, new HashMap<String,Collection<Gene>>());
		
		// If gene is not expressed, skip
		if(!isExpressed(gene)) {
			logger.info("Gene " + gene.getName() + " (" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + ") not expressed in expression dataset.");
			for (String scoreName : scoresToUse) {
				singleSampleScanPeaks.get(sample).get(gene).put(scoreName, rtrnPeaks.get(scoreName));
			}
			rejectFileWriterExpression.write(gene.toBED() + "\n");
			return;
		}
				
		logger.info("Finding scan peaks for sample " + sample.getSampleName() + " and gene " + gene.getName() + " (" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + ")");
		
		// Get fixed size windows with sufficient count and significant scan statistic
		significantWindows = findSignificantWindows(sample, gene, windowCountRejectFileWriters.get(sample), windowScanPvalAllFragmentsRejectFileWriters.get(sample), windowScanPvalWithFragmentLengthFilterRejectFileWriters.get(sample));

		for (String scoreName : scoresToUse) {
			// If no significant windows return
			if(significantWindows.get(scoreName).isEmpty()) {
				singleSampleScanPeaks.get(sample).get(gene).put(scoreName, rtrnPeaks.get(scoreName));
				continue;
			}
		
			// Merge overlapping windows
			Collection<Annotation> mergedWindows = mergePeaks(significantWindows.get(scoreName));
			logger.info("Found " + mergedWindows.size() + " merged windows for " + scoreName);
		
			// Trim each window
			TreeSet<Annotation> mergedTree = new TreeSet<Annotation>();
			mergedTree.addAll(mergedWindows);
			TreeSet<Annotation> trimmedMergedWindows = trimWindows(mergedTree, data);
		
			// Filter by score again
			TreeSet<Annotation> scanSigWindows = filterByScore(trimmedMergedWindows, sample, gene, scoreName, peakScanPvalRejectWriters.get(sample).get(scoreName));
			logger.info("Found " + scanSigWindows.size() + " windows after score filter for " + scoreName);
			
			// Filter on strand information
			TreeSet<Annotation> correctStrandWindows = new TreeSet<Annotation>();
			if(filterByStrand) {
				correctStrandWindows.addAll(filterByStrandInformation(scanSigWindows, sample, gene, strandRejectWriters.get(sample).get(scoreName)));
			} else {
				correctStrandWindows.addAll(scanSigWindows);
			}
			logger.info("Found " + correctStrandWindows.size() + " windows after strand filter for " + scoreName);
		
			// Filter on percent duplicates
			TreeSet<Annotation> dupsOk = filterByDuplicates(correctStrandWindows, sample, peakCutoffMaxReplicatePct, duplicateRejectWriters.get(sample).get(scoreName));
			logger.info("Found " + dupsOk.size() + " windows after duplicate filter for " + scoreName);
			
			// Final peaks
			finalPeaks.get(scoreName).addAll(dupsOk);
		
			double geneAvgCoverage = sample.getGeneAverageCoverage(gene);
			double geneCount = sample.getGeneCount(gene);
			int geneSize = coord.getSize(gene);
		
			// Add finishing touches to peaks
			for(Annotation peak : finalPeaks.get(scoreName)) {
			
				Gene window = new Gene(peak);
				// Name peaks
				window.setName(gene.getName() + ":" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd());
			
				MultiScore score = scoreWindow(gene, window, sample);
				double geneEnrichment = sample.getEnrichmentOverGene(gene, window);
				double ctrlEnrichment = score.getEnrichmentOverControl();
				double pval = score.getPvalue(scoreName);
				
				double[] extraFields;
				extraFields = new double[8];
				extraFields[0] = score.getCount();
				extraFields[1] = score.getCtrlCount();
				extraFields[2] = score.getSampleRegionCount();
				extraFields[3] = score.getCtrlRegionCount();
				extraFields[4] = geneEnrichment;
				extraFields[5] = window.size();
				extraFields[6] = (double) geneSize;
				extraFields[7] = pval;

				logger.debug("FINAL_PEAK\t" + gene.getName());
				logger.debug("FINAL_PEAK\t" + window.toBED());
				logger.debug("FINAL_PEAK\tname=" + window.getName());
				logger.debug("FINAL_PEAK\twindow_count=" + score.getCount());
				logger.debug("FINAL_PEAK\twindow_size=" + coord.getSize(window));
				logger.debug("FINAL_PEAK\twindow_avg_coverage=" + score.getAverageCoverage(data));
				logger.debug("FINAL_PEAK\tgene_count=" + geneCount);
				logger.debug("FINAL_PEAK\tgene_size=" + geneSize);
				logger.debug("FINAL_PEAK\tgene_avg_coverage=" + geneAvgCoverage);
				logger.debug("FINAL_PEAK\tenrichment_over_transcript=" + geneEnrichment);
				logger.debug("FINAL_PEAK\torientation=" + window.getOrientation().toString());
			
				window.setScore(ctrlEnrichment);
				window.setExtraFields(extraFields);
			
				rtrnPeaks.get(scoreName).add(window);

			}
		
		singleSampleScanPeaks.get(sample).get(gene).putAll(rtrnPeaks);
		
		}
	}
	
	/**
	 * Identify significant windows for each score within the sample and the gene
	 * Recompute scores for significant windows with fragment length filter and filter again
	 * @param sample The sample
	 * @param gene The gene
	 * @return Hash of score_name : significant_windows
	 * @throws IOException
	 */
	private Map<String,TreeSet<Annotation>> findSignificantWindows(SampleData sample, Gene gene, FileWriter windowCountRejectFileWriter, Map<String,FileWriter> rejectFileWriterAllFragments, Map<String,FileWriter> rejectFileWriterFragmentLengthFilter) throws IOException {
		Map<String,TreeSet<Annotation>> rtrn = new HashMap<String,TreeSet<Annotation>>();
		for (String scoreName : scoresToUse) {
			rtrn.put(scoreName,new TreeSet<Annotation>());
		}
		Map<Annotation, MultiScore> geneWindowScores = getWindowScores(sample,gene);
		for(Annotation window : geneWindowScores.keySet()) {
			MultiScore score = geneWindowScores.get(window);
			double count = score.getCount();
			if(count < peakWindowCountCutoff) {
				windowCountRejectFileWriter.write(window.toBED() + "\n");
				continue;
			}
			score.getDebugInfo("FIXED_SIZE_WINDOW_SCORES: ");
			for (String scoreName : scoresToUse) {
				double pval = score.getPvalue(scoreName);		
				if(pval < peakWindowScanPvalCutoff) {
					logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER");
					MultiScore fragmentLengthFilterScore = scoreWindowWithFragmentLengthFilter(gene, window, sample);
					score.getDebugInfo("FIXED_SIZE_WINDOW_SCORES_AFTER_FRAGMENT_LENGTH_FILTER: ");
					double pval2 = fragmentLengthFilterScore.getPvalue(scoreName);
					if(pval2 < peakWindowScanPvalCutoff) {
						window.setScore(pval2);
						logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER");
						rtrn.get(scoreName).add(window);
					} else {
						logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER");
						rejectFileWriterFragmentLengthFilter.get(scoreName).write(window.toBED() + "\n");
					}
				} else {
					logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER");
					rejectFileWriterAllFragments.get(scoreName).write(window.toBED() + "\n");
				}
			}
		}
		logger.info("Finished computing significant windows for " + gene.getName());
		for (String scoreName : scoresToUse) {
			logger.info("Found " + rtrn.get(scoreName).size() + " significant windows for " + scoreName + " score");	
		}
		return rtrn;
	}
	
	/**
	 * Check for gene in computed window scores, else compute window scores for the gene
	 * @param sample The sample
	 * @param gene The gene
	 * @return Hash of window : MultiScore
	 */
	private Map<Annotation, MultiScore> getWindowScores(SampleData sample,Gene gene) {
		if (windowScores.get(sample).containsKey(gene)) {
			return windowScores.get(sample).get(gene);
		} else {
			computeWindowScores(sample,gene);
			return windowScores.get(sample).get(gene);
		}
	}
	
	/**
	 * Make windows over gene, iterate over windows using MultiScore processor
	 * Store scores in windowScores
	 * @param sample The sample
	 * @param gene The gene
	 */
	private void computeWindowScores(SampleData sample, Gene gene) {
		Map<Annotation, MultiScore> scores = new TreeMap<Annotation, MultiScore>();
		if(gene.getSize() < windowSize) {
			logger.info(gene.getName() + " is smaller than window size. Not computing window binding site scores.");
			windowScores.get(sample).put(gene, scores);
			return;
		}
		WindowScoreIterator<MultiScore> iter = null;
		try {
			iter = sample.getData().scan(gene,windowSize,windowSize-stepSize,processors.get(sample));
		} catch(NullPointerException e) {
			logger.debug("Gene: " + gene.getName());
			logger.debug("Gene: " + gene.toBED());
		}
		
		double sampleGeneCount = sample.getGeneCount(gene);
		double ctrlGeneCount = ctrl.getGeneCount(gene);
		while (iter.hasNext()) {
			MultiScore score = iter.next();
			Annotation window = score.getAnnotation();
			score.setGene(gene);
			score.setSampleRegionCount(sampleGeneCount);
			score.setCtrlRegionCount(ctrlGeneCount);
			score.setRegionLength(gene.size());
			score.updateScores();
			//score.getDebugInfo("SCORE_ALL_WINDOWS_IN_GENE\t");
			//new ScanStatisticScore(score).getDebugInfo("SCORE_ALL_WINDOWS_IN_GENE\t");
			scores.put(window, score);
		}
		windowScores.get(sample).put(gene, scores);
	}
	
	/**
	 * Score a given window
	 * @param sample The sample
	 * @param window The window
	 * @param gene The gene
	 * @return MultiScore
	 */
	private MultiScore scoreWindow(Gene gene, Annotation window, SampleData sample) {
		if (finalWindowScores.get(sample).get(gene).containsKey(window)) {
			return finalWindowScores.get(sample).get(gene).get(window);
		}
		double sampleGeneCount = sample.getGeneCount(gene);
		double ctrlGeneCount = ctrl.getGeneCount(gene);
		double regionLength = gene.getSize();
		MultiScore score = new MultiScore(sample.getData(),ctrl.getData(),window,regionLength);
		score.setSampleRegionCount(sampleGeneCount);
		score.setCtrlRegionCount(ctrlGeneCount);
		score.setGene(gene);
		score.updateScores();
		return score;
	}
	
	/**
	 * Score a given window with fragment length filter
	 * @param sample The sample
	 * @param window The window
	 * @param gene The gene
	 * @return MultiScore
	 */
	private MultiScore scoreWindowWithFragmentLengthFilter(Gene gene, Annotation window, SampleData sample) {
		MultiScore score = new MultiScore(sample.getFragmentLengthFilterData(),ctrl.getFragmentLengthFilterData(),window,gene);
		return score;
	}
	
	/**
	 * Trim windows to max contiguous subregion above certain quantile
	 * @param data The alignment data
	 * @return New windows
	 */
	private TreeSet<Annotation> trimWindows(TreeSet<Annotation> untrimmed, TranscriptomeSpaceAlignmentModel data) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		for(Annotation window : untrimmed) {
			List<Double> coverageData = data.getPositionCountList(new Gene(window));
			Annotation trimmed = SampleData.trimMaxContiguous(window, coverageData, trimQuantile);
			rtrn.add(trimmed);
			logger.debug("MERGED_TRIMMED_WINDOW\t" + window.toBED());
		}
		return rtrn;
	}
	
	/**
	 * Rescores windows, filter by pvalue of given score 
	 * @param data The alignment data
	 * @return New windows
	 */
	private TreeSet<Annotation> filterByScore(TreeSet<Annotation> preFilter, SampleData sample, Gene gene, String scoreName, FileWriter rejectFileWriter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		double p;	
		if (!finalWindowScores.get(sample).containsKey(gene)) {
			finalWindowScores.get(sample).put(gene, new HashMap<Annotation,MultiScore>());
		}
		for(Annotation window : preFilter) {
			MultiScore score = scoreWindow(gene,window,sample);
			p = score.getPvalue(scoreName);
			if(p < peakWindowScanPvalCutoff) {
				rtrn.add(window);
				finalWindowScores.get(sample).get(gene).put(window, score);
			} else {
				rejectFileWriter.write(window.toBED() + "\n");
			}
		}
		return rtrn;
	}
	
	/**
	 * Rescores windows, filter by pvalue of given score 
	 * @param data The alignment data
	 * @return New windows
	 */
	private static TreeSet<Annotation> filterByDuplicates(TreeSet<Annotation> preFilter, SampleData sample, double maxPctMostCommonRead, FileWriter rejectFileWriter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		AlignmentModel data = sample.getData();
		for(Annotation window : preFilter) {
			Map<Alignment, Integer> replicateCounts = data.getOverlappingReadReplicateCounts(window, false);
			int total = 0;
			int largestReplicate = 0;
			String mostCommon = "";
			for(Alignment read : replicateCounts.keySet()) {
				int count = replicateCounts.get(read).intValue();
				total += count;
				if(count > largestReplicate) {
					largestReplicate = count;
					mostCommon = read.getChr() + ":" + read.getStart() + "-" + read.getEnd();
				}
			}
			double mostCommonPct = (double) largestReplicate / (double) total;
			if(mostCommonPct > maxPctMostCommonRead) {
				rejectFileWriter.write(window.toBED() + "\n");
				continue;
			}
			rtrn.add(window);
		}
		return rtrn;
	}
	
	/**
	 * Filters regions with less than 90% of reads on same strand as gene
	 * @return New windows
	 */
	private static TreeSet<Annotation> filterByStrandInformation(TreeSet<Annotation> preFilter, SampleData sample, Gene gene, FileWriter rejectFileWriter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		for(Annotation window : preFilter) {
			Strand orientation = AlignmentUtils.assignOrientationToWindow(sample.getOriginalBamFile(), window, sample.firstReadTranscriptionStrand(), 0.9);
			window.setOrientation(orientation);
			if(orientation.equals(gene.getOrientation())) {
				rtrn.add(window);
			} else if(orientation.equals(Strand.UNKNOWN)) {			
				rejectFileWriter.write(window.toBED() + "\n");
			} else {				
				rejectFileWriter.write(window.toBED() + "\n");
			}
		}
		return rtrn;
	}
	
	// Methods to compute window enrichments over gene background
	
	/**
	 * Compute window enrichments in all genes, store in hash.
	 * Note: only used for t-statistic
	 */
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
		logger.info("Done computing single sample window enrichments.");
	}
	
	
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
			Map<Annotation, MultiScore> scores = getWindowScores(sample,gene);
			
			if(gene.getSize() < windowSize) {
				logger.info(gene.getName() + " is smaller than window size. Not computing single sample window enrichments.");
				singleSampleWindowEnrichmentOverGene.get(sample).put(gene, sampleWindowEnrichments);
				continue;
			}
			for(Annotation window : scores.keySet()) {
				double windowAvgCoverage = scores.get(window).getAverageCoverage(sample.getData());
				double enrichment = windowAvgCoverage / geneAvgCoverage;
				logger.info(sample.getSampleName() + "\t" + gene.getName() + "\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + "\tavg_coverage=" + windowAvgCoverage + "\twindow_enrichment=" + enrichment);
				sampleWindowEnrichments.put(window, Double.valueOf(enrichment));
			}
			singleSampleWindowEnrichmentOverGene.get(sample).put(gene, sampleWindowEnrichments);
		}
	}
	
	// Methods for t-statistic testing
	
	/**
	 * Score all genes using t-statistic
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
				e += control.getSampleName() + ":" + getWindowScores(control,gene).get(window).getCount() + ":" + singleSampleWindowEnrichmentOverGene.get(control).get(gene).get(window).toString() + "\t";
			}
			for(SampleData signal : signalSamples) {
				e += signal.getSampleName() + ":" + getWindowScores(signal,gene).get(window).getCount() + ":" + singleSampleWindowEnrichmentOverGene.get(signal).get(gene).get(window).toString() + "\t";
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
	
	// Methods to parse command line arguments and sample file and set parameters for output

	private static CommandLineParser getCommandLineParser(String[] commandArgs) {
		CommandLineParser p = new CommandLineParser(true);
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
		p.addDoubleArg("-ep", "Scan P value cutoff for gene expression", false, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF);
		p.addBooleanArg("-ft", "First read is transcription strand", false, DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND);
		p.addDoubleArg("-r", "Cutoff for percentage of reads in peak coming from the most common replicate fragment", false, DEFAULT_PEAK_MAX_PCT_DUPLICATES);
		p.addBooleanArg("-sf", "Apply strand filter using read strand info", false, DEFAULT_FILTER_BY_STRAND);
		p.addBooleanArg("-ef", "Print additional info in BED file", false, DEFAULT_EXTRA_FIELDS);
		p.addStringListArg("-score", "Score to use in peak calling. Can be used multiple times. Possible values: binomial, scan_statistic", false, new ArrayList<String>(Arrays.asList(DEFAULT_SCORES)));
		p.addStringArg("-sc", "Scheduler (options: " + Scheduler.getCommaSeparatedList() + ")", false, Scheduler.OGS.toString());
		p.parse(commandArgs,true);
		return p;
	}
	
	protected static MultiSamplePeakCaller createFromCommandArgs(String[] commandArgs) throws IOException, DrmaaException {
		CommandLineParser p = getCommandLineParser(commandArgs);
		String sampleListFile = p.getStringArg("-l");
		String bedFile = p.getStringArg("-b");
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-s");
		double scanPvalCutoff = p.getDoubleArg("-sp");
		double trimQuantile = p.getDoubleArg("-q");
		String chrSizeFile = p.getStringArg("-c");
		double windowCountCutoff = p.getDoubleArg("-cp");
		double expressionScanPvalCutoff = p.getDoubleArg("-ep");
		boolean firstReadIsTranscriptionStrand = p.getBooleanArg("-ft");
		double maxPctMostCommonReplicatePerPeak = p.getDoubleArg("-r");
		boolean useStrandFilter = p.getBooleanArg("-sf");
		boolean extraFields =  p.getBooleanArg("-ef");
		Scheduler scheduler = Scheduler.fromString(p.getStringArg("-sc"));
		ArrayList<String> scoresToUse = p.getStringListArg("-score");
		
		ArrayList<String> validScoreNames = new ArrayList<String>(Arrays.asList(DEFAULT_SCORES));
		if (!validScoreNames.containsAll(scoresToUse)) {
			throw new IllegalArgumentException("Invalid score name. Scores must belong to " + validScoreNames.toString());
		}

		MultiSamplePeakCaller m = new MultiSamplePeakCaller(sampleListFile, bedFile, chrSizeFile, windowSize, stepSize);
		m.setExpressionScanPvalueCutoff(expressionScanPvalCutoff);
		m.setFirstReadTranscriptionStrand(firstReadIsTranscriptionStrand);
		m.setPeakCutoffMostCommonReplicate(maxPctMostCommonReplicatePerPeak);
		m.setPeakTrimQuantile(trimQuantile);
		m.setPeakWindowCountCutoff(windowCountCutoff);
		m.setPeakWindowScanPvalCutoff(scanPvalCutoff);
		m.setFilterByStrand(useStrandFilter);
		m.setExtraFields(extraFields);
		m.setScoreList(scoresToUse);
		m.setScheduler(scheduler);
		
		return m;
		 
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
	
	/**
	 * Get name of bed file to write for peaks
	 * @param sample Sample
	 * @param outDir Output directory name or null if current directory
	 * @param chrName Chromosome name or null if all chromosomes
	 * @return File name
	 */
	protected Map<String,String> getPeakBedFileName(SampleData sample, String outDir, String chrName) {
		Map<String,String> rtrn = new HashMap<String,String>();
		for (String scoreName : scoresToUse) {
			String fileName = "";
			if(outDir != null) {
				fileName += outDir + "/";
			}
			fileName += sample.getSampleName() + "_" + scoreName + "_peaks_" + windowSize + "_" + stepSize + "_" + peakWindowScanPvalCutoff + "_" + trimQuantile;
			if(chrName != null) {
				fileName += "_" + chrName;
			}
			fileName += ".bed";
			rtrn.put(scoreName, fileName);
		}
		return rtrn;
	}
	
	/**
	 * Get chromosome list for batch writing
	 * @param commandArgs
	 * @return
	 */
	private static String commandLineBatchChrList(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getStringArg("-cl");
	}
	
	/**
	 * Get memory requirements
	 * @param commandArgs
	 * @return
	 */
	private static int commandLineBatchMemRequest(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getIntArg("-m");
	}
	
	/**
	 * Get debug level
	 * @param commandArgs
	 * @return
	 */
	protected static boolean commandLineHasDebugFlag(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getBooleanArg("-d");
	}
	
	/**
	 * Get whether to use batch writer
	 * @param commandArgs
	 * @return
	 */
	private static boolean commandLineHasBatchFlag(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getBooleanArg("-batch");
	}
	
	/**
	 * Get outdir
	 * @param commandArgs
	 * @return
	 */
	protected static String commandLineOutDir(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getStringArg("-o");		
	}
	
	/**
	 * Get jarfile for batch writer
	 * @param commandArgs
	 * @return
	 */
	private static String commandLineBatchJar(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		String jar = p.getStringArg("-bj");
		if(jar == null) {
			throw new IllegalArgumentException("Must provide batch peak caller jar file with option -bj.");
		}
		return jar;
	}
	
	// Main method
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		MultiSamplePeakCaller m = createFromCommandArgs(args);
		
		if(commandLineHasDebugFlag(args)) {
			m.setLoggerLevel(Level.DEBUG);
		}
			
		if(commandLineHasBatchFlag(args)) {
			m.batchWriteSingleSamplePeaksAllSamples(args, commandLineBatchChrList(args), commandLineBatchMemRequest(args));
		} else {
			m.initializeFilterRejectWriters("all_chr", commandLineOutDir(args) + "/" + FILTER_REJECT_DIR);
			m.writeSingleSamplePeaksAllSamples(commandLineOutDir(args));
			m.closeFilterRejectWriters();
		}
				
		logger.info("");
		logger.info("All done.");	
		
	}
	
	// Inherited methods from nextgen.core.analysis.PeakCaller
	
	/* (non-Javadoc)
	 * @see nextgen.core.analysis.PeakCaller#scoreWindows(java.util.Collection)
	 */
	@Override
	public void scoreWindows(Collection<Annotation> windows) {
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
		TreeSet<Annotation> peakTree = new TreeSet<Annotation>();
		peakTree.addAll(peaks);
		Collection<Annotation> mergedWindows = AnnotationUtils.mergeOverlappingBlocks(peakTree);

		return mergedWindows;
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
	
	// Extra classes for file parser, sample permutation
	
	/**
	 * Helper class to parse file containing list of sample types and bam file names
	 * @author prussell
	 *
	 */
	private class SampleFileParser {
		
		private String EXPRESSION_LABEL = "Expression";
		private String EXPRESSION_PVAL_CUTOFF_LABEL = "Expression_pval_cutoff";
		private String EXPRESSION_AVG_COVERAGE_CUTOFF_LABEL = "Expression_avg_coverage_cutoff";
		private String CONTROL_FOR_SCORES_LABEL = "Control_for_scores";
		private String CONTROL_LABEL = "Control";
		private String SIGNAL_LABEL = "Signal";
		private String FIRST_READ_TRANSCRIPTION_STRAND_LABEL = "first_read_transcription_strand";
		private String SECOND_READ_TRANSCRIPTION_STRAND_LABEL = "second_read_transcription_strand";
		private GenomeSpaceSampleData expressionSampleData;
		private ArrayList<SampleData> controlData;
		private ArrayList<SampleData> signalData;

		
		public SampleFileParser(String sampleFileName) throws IOException {
			controlData = new ArrayList<SampleData>();
			signalData = new ArrayList<SampleData>();
			expressionSampleData = null;
			ctrl = null;
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
		 * Get the  control
		 * @return  control dataset
		 */
		public SampleData getCtrlData() {
			return ctrl;
		}
		
		/**
		 * Parse the sample file and populate data sets
		 * @throws IOException
		 */
		private void parseFile(String sampleFileName) throws IOException {
			boolean foundExpressionData = false;
			boolean foundScoreCtrl = false;
			String scoringCtrlFileName = null;
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
			
			String ctrlLine = b.readLine();
			s.parse(ctrlLine);
			if(s.getFieldCount() != 2) {
				crashWithHelpMessage();
			}
			if(s.asString(0).equals(CONTROL_FOR_SCORES_LABEL)) {
				scoringCtrlFileName = s.asString(1);
				logger.debug("Using for control: " + scoringCtrlFileName);
			} else {
				crashWithHelpMessage();
			}
			
			if(!b.ready()) {
				crashWithHelpMessage();
			}
			
			while(b.ready()) {
				
				String line = b.readLine();
				s.parse(line);
				
				if(s.getFieldCount() == 0) continue;
				if(line.substring(0, 1).equals(COMMENT_CHAR)) continue;
				
				String label = s.asString(0);
				String bamFile = s.asString(1);

				if(label.equals(EXPRESSION_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					if(foundExpressionData) {
						crashWithHelpMessage();
					}
					boolean firstReadIsTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadIsTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					logger.info("Creating sample data object for gene expression from bam file " + bamFile);
					GenomeSpaceSampleData sample = new GenomeSpaceSampleData(bamFile, sizeFile, genes, cutoff, true);
					expressionSampleData = sample;
					foundExpressionData = true;
					if (bamFile.equals(scoringCtrlFileName)) {
						ctrl = new SampleData(bamFile, firstReadIsTranscriptionStrand, genes, cutoff, expByScanPval, false);	
						foundScoreCtrl = true;
						logger.debug("Found control: " + bamFile);
					}
					continue;
				}
				
				if(label.equals(CONTROL_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					logger.info("Creating sample data object for bam file " + bamFile);
					boolean firstReadIsTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadIsTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					SampleData sample = new SampleData(bamFile, firstReadIsTranscriptionStrand, genes, cutoff, expByScanPval, false);
					controlData.add(sample);
					if (bamFile.equals(scoringCtrlFileName)) {
						ctrl = sample;	
						foundScoreCtrl = true;
						logger.debug("Found control: " + bamFile);
					}
					continue;
				}
				
				if(label.equals(SIGNAL_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					logger.info("Creating sample data object for bam file " + bamFile);
					boolean firstReadIsTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadIsTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					SampleData sample = new SampleData(bamFile, firstReadIsTranscriptionStrand, genes, cutoff, expByScanPval, false);
					signalData.add(sample);
					if (bamFile.equals(scoringCtrlFileName)) {
						ctrl = sample;	
						foundScoreCtrl = true;
						logger.debug("Found control: " + bamFile);
					}
					continue;
				}
				
				crashWithHelpMessage();
				
			}
			
			r.close();
			b.close();
			
			if(!foundExpressionData|!foundScoreCtrl) {
				logger.debug("Found expression data: " + foundExpressionData + "\nFound score ctrl: " + foundScoreCtrl);
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
			logger.error("Second line must be:");
			logger.error(CONTROL_FOR_SCORES_LABEL + "\t<bam_file_name>");
			logger.error("Exactly one line must be of the form:");
			logger.error(EXPRESSION_LABEL + "\t<bam_file_name>\t" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("");
			logger.error("Each additional line must be of the form:");
			logger.error(CONTROL_LABEL + "\t<bam_file_name>\t" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("- or -");
			logger.error(SIGNAL_LABEL + "\t<bam_file_name>" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("- or -");
			logger.error(COMMENT_CHAR + "\tcomment");
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

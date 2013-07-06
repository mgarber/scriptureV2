package broad.pda.seq.protection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.utils.AnnotationUtils;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

/**
 * @author prussell
 *
 */
public class TwoSampleScanSkellamAnalysis extends TwoSampleScanSkellamPeakCaller {

	@SuppressWarnings("unused")
	private static int DEFAULT_MAX_FRAGMENT_LENGTH = 2000;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	private boolean hasBackgroundExpressionScores;
	private Map<Gene, Collection<Annotation>> significantPeaks;
	
	
	/**
	 * Construct object from bam files
	 * @param backgroundAlignmentFile Bam file of background sample alignments
	 * @param signalAlignmentFile Bam file of signal sample alignments
	 * @param bedFile Bed file of annotations
	 * @param window Window size
	 * @param step Step size
	 * @param alphaSkellam Max P value for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @param trimQuantile Coverage quantile within untrimmed peak for trimming by max contiguous subsequence
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public TwoSampleScanSkellamAnalysis(String backgroundAlignmentFile, String signalAlignmentFile, String bedFile, int window, int step, double alphaSkellam, double alphaScan, double trimQuantile) throws IllegalArgumentException, IOException {
		this(backgroundAlignmentFile, signalAlignmentFile, BEDFileParser.loadDataByChr(new File(bedFile)), window, step, alphaSkellam, alphaScan, trimQuantile);
	}
	
	
	/**
	 * Construct object from bam files
	 * @param backgroundAlignmentFile Bam file of background sample alignments
	 * @param signalAlignmentFile Bam file of signal sample alignments
	 * @param genesByChr Genes by chromosome
	 * @param window Window size
	 * @param step Step size
	 * @param alphaSkellam Max P value for skellam statistic
	 * @param alphaScan Max P value for scan statistic
	 * @param trimQuantile Coverage quantile within untrimmed peak for trimming by max contiguous subsequence
	 * @throws IllegalArgumentException
	 * @throws IOException
	 */
	public TwoSampleScanSkellamAnalysis(String backgroundAlignmentFile, String signalAlignmentFile, Map<String, Collection<Gene>> genesByChr, int window, int step, double alphaSkellam, double alphaScan, double trimQuantile) throws IllegalArgumentException, IOException {
		genes = genesByChr;
		transcriptomeSpace = new TranscriptomeSpace(genes);
		backgroundData = new ScanStatisticDataAlignmentModel(backgroundAlignmentFile, transcriptomeSpace);
		signalData = new TranscriptomeSpaceAlignmentModel(signalAlignmentFile, transcriptomeSpace);
		windowSize = window;
		stepSize = step;
		pValueCutoffSkellam = alphaSkellam;
		pValueCutoffScan = alphaScan;
		trimPeakByQuantile = trimQuantile;
		hasBackgroundExpressionScores = false;
		backgroundProcessor = new CountScore.Processor(backgroundData, false);
		signalProcessor = new ScanStatisticScore.Processor(signalData, false);
		backgroundName = sampleNameFromBamName(backgroundAlignmentFile);
		signalName = sampleNameFromBamName(signalAlignmentFile);
		significantPeaks = new TreeMap<Gene, Collection<Annotation>>();
		logger.info("Instantiated " + TwoSampleScanSkellamAnalysis.class.getName() + " object.");
	}

	/**
	 * Get sample name from bam file name
	 * @param bamFileName Bam file name
	 * @return The part before .bam extension if exists or the whole file name
	 */
	public static String sampleNameFromBamName(String bamFileName) {
		if(bamFileName.substring(bamFileName.length() - 4).equals(".bam")) {
			return bamFileName.substring(0, bamFileName.length() - 4);
		}
		return bamFileName;
	}
	
	/**
	 * For an annotation get the Poisson parameter describing the null model for the number of reads falling in a window in the background dataset
	 * @param region The annotation of interest
	 * @return Lambda, the expected number of reads in the window under the null Poisson model
	 * @throws IOException
	 */
	private double calculatePoissonLambdaBackground(Annotation region) throws IOException {
		int totalSize = region.getSize();
		return getPoissonLambda(totalSize, backgroundData.getCount(region), windowSize);
	}

	/**
	 * For an annotation get the Poisson parameter describing the null model for the number of reads falling in a window in the signal dataset
	 * @param region The annotation of interest
	 * @return Lambda, the expected number of reads in the window under the null Poisson model
	 * @throws IOException
	 */
	private double calculatePoissonLambdaSignal(Annotation region) throws IOException {
		int totalSize = region.getSize();
		return getPoissonLambda(totalSize, signalData.getCount(region), windowSize);
	}

	
	/**
	 * Compute expression scan p value for each gene in background dataset
	 * Read from and/or write to a cached file
	 * If the file exists, read as many scores as possible from the file, then compute the rest and add to the file
	 * @throws IOException
	 */
	private void computeBackgroundExpressionScores() throws IOException {
		
		logger.info("Computing background gene expression scores...");
		backgroundScanPvalues = new TreeMap<String, Map<Gene, Double>>();
		
		// First try to read scores from file
		String scoreFileName = backgroundName + ".backgroundExpressionScores";
		File scoreFile = new File(scoreFileName);
		if(scoreFile.exists()) {
			logger.info("Trying to read background expression scores from file " + scoreFileName);
			boolean foundAllInFile = true;
			FileReader r = new FileReader(scoreFile);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			Map<String, Double> geneNameToScore = new TreeMap<String, Double>();
			while(b.ready()) {
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() != 2) continue;
				String geneName = s.asString(0);
				double score = s.asDouble(1);
				geneNameToScore.put(geneName, Double.valueOf(score));
			}
			
			for(String chr : genes.keySet()) {
				if(!foundAllInFile) break;
				Map<Gene, Double> chrGeneScores = new TreeMap<Gene, Double>();
				for(Gene gene : genes.get(chr)) {
					if(!foundAllInFile) break;
					String name = gene.getName();
					if(!geneNameToScore.containsKey(name)) {
						foundAllInFile = false;
					}
					chrGeneScores.put(gene, geneNameToScore.get(name));
				}
				backgroundScanPvalues.put(chr, chrGeneScores);
			}
			r.close();
			b.close();
			if(foundAllInFile) {
				logger.info("Got all background scores from file.");
				hasBackgroundExpressionScores = true;
				return;
			}
			logger.info("Couldn't get all background scores from file. Computing remaining scores and overwriting file " + scoreFileName);
		}
		
		// Didn't read from file so compute the scores and write to file
		logger.info("Computing background gene expression scores and writing to file " + scoreFileName);
		FileWriter w = new FileWriter(scoreFile);
		int finished = 0;
		for(String chr : genes.keySet()) {
			Map<Gene, Double> chrGeneScores = new TreeMap<Gene, Double>();
			for(Gene gene : genes.get(chr)) {
				double score = -99;
				if(finished > 0 && finished % 1000 == 0) {
					logger.info("Finished " + finished + " genes.");
				}
				boolean found = false;
				if(backgroundScanPvalues.containsKey(chr)) {
					if(backgroundScanPvalues.get(chr).containsKey(gene.getName())) {
						score = backgroundScanPvalues.get(chr).get(gene.getName()).doubleValue();
						found = true;
					}
				}
				if(!found) {
					score = backgroundData.scoreWindow(gene).getScanPvalue();
				}
				chrGeneScores.put(gene, Double.valueOf(score));
				w.write(gene.getName() + "\t" + score + "\n");
				finished++;
			}
			backgroundScanPvalues.put(chr, chrGeneScores);
		}
		hasBackgroundExpressionScores = true;
		w.close();
	}
	
	/**
	 * Get the gene set
	 * @return The gene set
	 */
	public Map<String, Collection<Gene>> getGenes() {
		return genes;
	}
	
	@Override
	public String toString() {
		return backgroundName + "_" + signalName;
	}
	
	@Override
	public int hashCode() {
		return toString().hashCode();
	}
	
	/**
	 * Get the window size
	 * @return The window size
	 */
	public int getWindowSize() {return windowSize;}
	
	/**
	 * Get the step size
	 * @return The step size
	 */
	public int getStepSize() {return stepSize;}
	
	/**
	 * Get the p value cutoff for skellam statistic
	 * @return The p value cutoff
	 */
	public double getPvalueCutoffSkellam() {return pValueCutoffSkellam;}
	
	/**
	 * Get the p value cutoff for scan statistic
	 * @return The p value cutoff
	 */
	public double getPvalueCutoffScan() {return pValueCutoffScan;}

	/**
	 * Get background sample name
	 * @return Background sample name
	 */
	public String getBackgroundName() {
		return backgroundName;
	}
	
	/**
	 * Get the background data alignment model
	 * @return Background alignment data
	 */
	public ScanStatisticDataAlignmentModel getBackgroundData() {
		return backgroundData;
	}
	
	/**
	 * Get the signal data alignment model
	 * @return signal alignment data
	 */
	public ScanStatisticDataAlignmentModel getSignalData() {
		return signalData;
	}

	
	/**
	 * Get signal sample name
	 * @return Signal sample name
	 */
	public String getSignalName() {
		return signalName;
	}

	
	/**
	 * Add fragment length filter for paired reads
	 * @param maxFragmentLength Max fragment length
	 */
	@SuppressWarnings("unused")
	private void addFragmentLengthFilter(int maxFragmentLength) {
		backgroundData.addFilter(new FragmentLengthFilter(transcriptomeSpace, maxFragmentLength));
		signalData.addFilter(new FragmentLengthFilter(transcriptomeSpace, maxFragmentLength));
	}
	
	/**
	 * Add genomic span filter for paired reads
	 * @param maxGenomicSpan Max genomic span
	 */
	private void addGenomicSpanFilter(int maxGenomicSpan) {
		backgroundData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
		signalData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
	}
	
	/**
	 * Get all windows passing skellam score cutoff and scan P value cutoff
	 * @param gene The gene
	 * @param checkBackgroundExpression Filter for gene expression in background sample
	 * @return All significant windows, or null if can't make Poisson model due to lack of reads or no significant windows
	 * @throws IOException
	 */
	private TreeSet<Annotation> findSignificantFixedSizeWindows(Gene gene, boolean checkBackgroundExpression) throws IOException {
		
		if(checkBackgroundExpression) {
			if(!hasBackgroundExpressionScores) computeBackgroundExpressionScores();
			double expressionPval = backgroundScanPvalues.get(gene.getReferenceName()).get(gene).doubleValue();
				
			if(expressionPval > EXPRESSION_PVALUE_CUTOFF) {
				//logger.info(gene.getName() + "\t" + "gene_not_expressed\texpression_pval=" + expressionPval);
				return null;
			}
		}
		
		// Don't consider genes smaller than window size
		/*
		 * TODO should be able to score genes that are smaller than the window size
		 * This would mean Gene.getWindows() would have to return one small window for genes that are smaller than the window size
		 */
		if(gene.size() < windowSize) {
			//logger.info(gene.getName() + "\t" + "smaller_than_window_size");
			return null;			
		}
		
		TreeSet<Annotation> sigWindows = new TreeSet<Annotation>();
		
		try {
			
			double backgroundLambda = calculatePoissonLambdaBackground(gene);
			double signalLambda = calculatePoissonLambdaSignal(gene);
			
			if(signalLambda == 0) {
				//logger.info(gene.getName() + "\t" + "signal_lambda_is_zero");
				return null;
			}
			
			
			Collection<Gene> baseGenes = new TreeSet<Gene>();
			baseGenes.add(gene);
			
			WindowScoreIterator<CountScore> backgroundScoreIter = transcriptomeSpace.scan(baseGenes, windowSize, windowSize - stepSize, backgroundProcessor);
			WindowScoreIterator<ScanStatisticScore> signalScoreIter = transcriptomeSpace.scan(baseGenes, windowSize, windowSize - stepSize, signalProcessor);
			
			
			while(backgroundScoreIter.hasNext()) {
				
				//String logString = gene.getName() + "\tbackgroundLambda=" + backgroundLambda + "\tsignalLambda=" + signalLambda + "\t";
				
				CountScore backgroundScore = backgroundScoreIter.next();
				ScanStatisticScore signalScore = signalScoreIter.next();
				
				Annotation region = backgroundScore.getAnnotation();
				//String backgroundCoords = region.getChr() + ":" + region.getStart() + "-" + region.getEnd();
				//logString += backgroundCoords + "\t";
				
				
				// Skip if signal read count is too low to be interesting
				int signalCount = (int) Math.round(signalScore.getCount());
				//logString += "signal_count=" + signalCount + "\t";
				if(signalCount < Math.max(MIN_WINDOW_COUNT, signalLambda)) {
					//logString += "signal_count_too_low";
					//logger.info(logString);
					continue; 
				}
				
				// Get background count
				int backgroundCount = (int) Math.round(backgroundScore.getCount());
				//logString += "background_count=" + backgroundCount + "\t";
				if(backgroundCount < 0) {
					//logString += "background_count<0";
					//logger.info(logString);
					continue;
				}
				
				// filter for skellam P value
				double skellamPval = getSkellamPvalue(backgroundLambda, signalLambda, backgroundCount, signalCount);
				//logString += "skellam_pval=" + skellamPval + "\t";
				if(skellamPval > pValueCutoffSkellam) {
					//logString += "skellam_pval_not_significant";
					//logger.info(logString);
					continue;
				}
				
				// filter for scan P value
				double scanPval = signalScore.getScanPvalue();
				//logString += "scan_pval=" + scanPval + "\t";
				if(scanPval > pValueCutoffScan) {
					//logString += "scan_pval_not_significant";
					//logger.info(logString);
					continue;
				}
				
				//logString += "region_significant";
				//logger.info(logString);
				sigWindows.add(region);
				
				
			}
		} catch(IllegalArgumentException e) {
			logger.info(gene.getName() + "\t" + "cant_make_poisson_model");
			return null;
		}

		if(sigWindows.isEmpty()) return null;
		return sigWindows;
		
	}
		
	/**
	 * Get significant peaks for a gene
	 * @param gene The gene
	 * @return The set of significant peaks
	 * @throws IOException
	 */
	public Collection<Annotation> getSignificantPeaks(Gene gene) throws IOException {
		return getSignificantPeaks(gene, true);
	}

	
	/**
	 * Get significant peaks for a gene providing expression P value
	 * @param gene The gene
	 * @param checkBackgroundExpression Filter for gene expression in background sample
	 * @return The set of significant peaks
	 * @throws IOException
	 */
	public Collection<Annotation> getSignificantPeaks(Gene gene, boolean checkBackgroundExpression) throws IOException {
		Collection<Annotation> rtrn = findSignificantPeaks(gene, checkBackgroundExpression);
		significantPeaks.put(gene, rtrn);
		return rtrn;
	}
	
	/**
	 * Get merged and trimmed significant peaks for a gene or null if there are none
	 * @param gene The gene
	 * @param checkBackgroundExpression Filter for gene expression in background sample
	 * @param fixedSizeWriter FileWriter to write fixed size significant windows
	 * @param mergedWriter FileWriter to write merged significant windows
	 * @param trimmedWriter FileWriter to write merged and trimmed significant windows
	 * @return Map associating each gene with its significant windows or empty map if there are none
	 * @throws IOException
	 */
	private Collection<Annotation> findSignificantPeaks(Gene gene, boolean checkBackgroundExpression, FileWriter fixedSizeWriter, FileWriter mergedWriter, FileWriter trimmedWriter) throws IOException {
				

		// Get fixed size windows that pass cutoff
		TreeSet<Annotation> sigWindowsFixedSize = findSignificantFixedSizeWindows(gene, checkBackgroundExpression);
		if(sigWindowsFixedSize == null) {
			//logger.info(gene.getName() + "\tno_significant_windows");
			return null;
		}
		if(fixedSizeWriter != null) {
			for(Annotation a : sigWindowsFixedSize) {
				fixedSizeWriter.write(a.toBED() + "\n");
			}
		}
		
		// Merge overlapping windows
		Collection<Annotation> mergedWindows = AnnotationUtils.mergeOverlappingBlocks(sigWindowsFixedSize);
		//logger.info(gene.getName() + "\tnum_merged_windows=" + mergedWindows.size());
		if(mergedWriter != null) {
			for(Annotation a : mergedWindows) {
				mergedWriter.write(a.toBED() + "\n");
			}
		}
		
		// Trim merged regions
		Collection<Annotation> trimmedWindows = new TreeSet<Annotation>();
		for(Annotation peak : mergedWindows) {
			Gene peakGene = new Gene(peak);
			Annotation trimmedPeak = trimPeak(peakGene);
			trimmedWindows.add(trimmedPeak);
			//logger.info(gene.getName() + "\tadded_trimmed_peak\t" + trimmedPeak.toBED());
		}
		if(trimmedWriter != null) {
			for(Annotation a : trimmedWindows) {
				trimmedWriter.write(a.toBED() + "\n");
			}
		}
		
		return trimmedWindows;
		
	}
	
	/**
	 * Get merged and trimmed significant peaks for a gene or null if there are none
	 * @param gene The gene
	 * @param checkBackgroundExpression Filter for gene expression in background sample
	 * @return Map associating each gene with its significant windows or empty map if there are none
	 * @throws IOException
	 */
	private Collection<Annotation> findSignificantPeaks(Gene gene, boolean checkBackgroundExpression) throws IOException {
				
		return findSignificantPeaks(gene, checkBackgroundExpression, null, null, null);
		
	}

	
	/**
	 * Get significant peaks and write to file
	 * @param bedFile Bed file for significant peaks
	 * @param suppOutFilePrefix File prefix for fixed size, merged, and trimmed significant windows
	 * @throws IOException
	 */
	private void findAndWritePeaks(String bedFile, String suppOutFilePrefix) throws IOException {
		logger.info("Getting peaks and writing to file " + bedFile);
		FileWriter fixedSizeWriter = null;
		FileWriter mergedWriter = null;
		FileWriter trimmedWriter = null;
		if(suppOutFilePrefix != null) {
			fixedSizeWriter = new FileWriter(suppOutFilePrefix + "_fixed_size.bed");
			mergedWriter = new FileWriter(suppOutFilePrefix + "_merged.bed");
			trimmedWriter = new FileWriter(suppOutFilePrefix + "_merged_trimmed.bed");
		}
		FileWriter writer = new FileWriter(bedFile);
		for(String chr : genes.keySet()) {
			logger.info("Analyzing genes on reference sequence: " + chr);
			for(Gene gene : genes.get(chr)) {
				Collection<Annotation> sigPeaks = findSignificantPeaks(gene, true, fixedSizeWriter, mergedWriter, trimmedWriter);
				if(sigPeaks == null) continue;
				for(Annotation window : sigPeaks)	{
					Gene peak = new Gene(window);
					writer.write(peak.toBED(255,0,0) + "\n");
				}
			}
		}
		writer.close();
		if(fixedSizeWriter != null) fixedSizeWriter.close();
		if(mergedWriter != null) mergedWriter.close();
		if(trimmedWriter != null) trimmedWriter.close();
	}
	
	/**
	 * Get significant peaks and write to file
	 * @param bedFile Bed file for significant peaks
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void findAndWritePeaks(String bedFile) throws IOException {
		findAndWritePeaks(bedFile, null);
	}

	/**
	 * Get skellam P value of difference of read counts in region
	 * @param region The region
	 * @param backgroundLambda Lambda for background sample
	 * @param signalLambda Lambda for signal sample
	 * @return The Skellam P value of the difference
	 * @throws IOException
	 */
	public double skellamPval(Annotation region, double backgroundLambda, double signalLambda) throws IOException {
		CountScore backgroundScore = new CountScore(backgroundData, region, false);
		CountScore signalScore = new CountScore(signalData, region, false);
		int backgroundCount = (int) Math.round(backgroundScore.getCount());
		int signalCount = (int) Math.round(signalScore.getCount());
		return getSkellamPvalue(backgroundLambda, signalLambda, backgroundCount, signalCount);
	}
	
	/**
	 * Get scan P value in signal dataset with respect to the transcript
	 * @param parentGene The parent gene
	 * @param window The window
	 * @return The scan P value of counts in the region
	 */
	public double signalScanPval(Gene parentGene, Annotation window) {
		ScanStatisticScore windowScore = new ScanStatisticScore(signalData, window, false);
		ScanStatisticScore geneScore = new ScanStatisticScore(signalData, parentGene, false);
		double windowCount = windowScore.getCount();
		double globalLambda = geneScore.getLocalLambda();
		double windowLength = window.getSize();
		double globalLength = parentGene.getSize();
		return AlignmentDataModelStats.calculatePVal(new Double(windowCount).intValue(), globalLambda, windowLength, globalLength);
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalArgumentException 
	 */
	public static void main(String[] args) throws IllegalArgumentException, IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Background bam file", true);
		p.addStringArg("-s", "Signal bam file", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addIntArg("-w", "Window size", true);
		p.addIntArg("-t", "Step size", true);
		p.addFloatArg("-psk", "P value cutoff for skellam statistic", false, (float)0.05);
		p.addFloatArg("-psc", "P value cutoff for scan statistic", false, (float)0.05);
		p.addFloatArg("-q", "Coverage quantile within untrimmed peak for trimming by max contiguous subsequence", false, (float)0.5);
		//p.addIntegerArg("-mf", "Max fragment length for paired reads", false, Integer.valueOf(DEFAULT_MAX_FRAGMENT_LENGTH));
		p.addIntArg("-mg", "Max genomic span for paired reads", false, DEFAULT_MAX_GENOMIC_SPAN);
		p.addStringArg("-o", "Output file", true);
		p.addStringArg("-outdebug", "Output file prefix for debug bed files", false);
		p.parse(args);
		String backgroundFile = p.getStringArg("-b");
		String signalFile = p.getStringArg("-s");
		String bedFile = p.getStringArg("-g");
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-t");
		double alphaSkellam = p.getFloatArg("-psk");
		double alphaScan = p.getFloatArg("-psc");
		double trimQuantile = p.getFloatArg("-q");
		//int maxFragmentLength = p.getIntegerArg("-mf").intValue();
		int maxGenomicSpan = p.getIntArg("-mg");
		String outFile = p.getStringArg("-o");
		String outDebug = p.getStringArg("-outdebug");
		
		
		TwoSampleScanSkellamAnalysis svbs = new TwoSampleScanSkellamAnalysis(backgroundFile, signalFile, bedFile, windowSize, stepSize, alphaSkellam, alphaScan, trimQuantile);
		
		svbs.addGenomicSpanFilter(maxGenomicSpan);
		svbs.findAndWritePeaks(outFile,outDebug);
		
	}

}

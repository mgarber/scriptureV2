/**
 * 
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
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.math.FisherCombinedProbabilityTest;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 *
 */
public final class MultiSampleScanSkellamAnalysis {
	
	private Collection<TwoSampleScanSkellamAnalysis> pairedSamples;
	private Map<String, Collection<Gene>> genes;
	private Map<Gene, Collection<Annotation>> jointSignificantPeaks;
	private int windowSize;
	private int stepSize;
	private double alphaSkellam;
	private double alphaScan;
	private double trimQuantile;
	
	private static int DEFAULT_WINDOW_SIZE = 30;
	private static int DEFAULT_STEP_SIZE = 1;
	private static double DEFAULT_ALPHA_SKELLAM = 0.01;
	private static double DEFAULT_ALPHA_SCAN = 0.01;
	private static double DEFAULT_TRIM_QUANTILE = 0.7;
	
	private Map<Gene, Map<TwoSampleScanSkellamAnalysis, Collection<Annotation>>> sampleLevelPeaks;
	private boolean hasSampleLevelPeaks;
	
	private static Logger logger = Logger.getLogger(MultiSampleScanSkellamAnalysis.class.getName());
	
	/**
	 * Instantiate with bed file of genes and default peak calling parameters
	 * @param bedFile
	 * @throws IOException
	 */
	private MultiSampleScanSkellamAnalysis(String backgroundBamFile, String signalSampleBamList, String bedFile) throws IOException {
		this(backgroundBamFile, signalSampleBamList, bedFile, DEFAULT_WINDOW_SIZE, DEFAULT_STEP_SIZE, DEFAULT_ALPHA_SKELLAM, DEFAULT_ALPHA_SCAN, DEFAULT_TRIM_QUANTILE);
	}
	
	/**
	 * Instantiate with bed file of genes and peak calling parameters
	 * @param backgroundBamFile Background sample bam file
	 * @param signalSampleBamList File containing list of signal bam file names
	 * @param bedFile Bed file of genes
	 * @param window Window size for peak calling
	 * @param step Step size for peak calling
	 * @param skellamPvalCutoff Skellam P value cutoff for peak calling within paired sample
	 * @param scanPvalCutoff Scan P value cutoff for peak calling in signal sample
	 * @param trimMaxQuantile Quantile parameter to trim max contiguous subsequence for peak calling in signal sample
	 * @throws IOException
	 */
	private MultiSampleScanSkellamAnalysis(String backgroundBamFile, String signalSampleBamList, String bedFile, int window, int step, double skellamPvalCutoff, double scanPvalCutoff, double trimMaxQuantile) throws IOException {
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		jointSignificantPeaks = new TreeMap<Gene, Collection<Annotation>>();
		windowSize = window;
		stepSize = step;
		alphaSkellam = skellamPvalCutoff;
		alphaScan = scanPvalCutoff;
		trimQuantile = trimMaxQuantile;
		pairedSamples = new ArrayList<TwoSampleScanSkellamAnalysis>();
		hasSampleLevelPeaks = false;
		FileReader r = new FileReader(signalSampleBamList);
		BufferedReader b = new BufferedReader(r);
		while(b.ready()) {
			String signalBamFile = b.readLine();
			TwoSampleScanSkellamAnalysis p = new TwoSampleScanSkellamAnalysis(backgroundBamFile, signalBamFile, genes, windowSize, stepSize, alphaSkellam, alphaScan, trimQuantile);
			pairedSamples.add(p);
		}
	}
	
	/**
	 * Get the genes by chromosome
	 * @return Genes by chromosome
	 */
	Map<String, Collection<Gene>> getGenes() {
		return genes;
	}
	
	/**
	 * Get peaks that are jointly significant across samples in scan and skellam statistics
	 * @param gene The gene
	 * @return The peaks
	 * @throws IOException
	 */
	private Collection<Annotation> getJointSignificantPeaks(Gene gene) throws IOException {
		return getJointSignificantPeaks(gene, null);
	}
	
	/**
	 * Get peaks that are jointly significant across samples in scan and skellam statistics, and write sample level peaks to a bed file
	 * @param gene The gene
	 * @param bedFileSampleLevelPeaks File stream to write all sample level peaks to bed file
	 * @return The peaks
	 * @throws IOException
	 */
	private Collection<Annotation> getJointSignificantPeaks(Gene gene, FileWriter bedFileSampleLevelPeaks) throws IOException {
		if(jointSignificantPeaks.containsKey(gene)) {
			return jointSignificantPeaks.get(gene);
		}
		
		// Filter for background expression
		TwoSampleScanSkellamAnalysis firstSample = pairedSamples.iterator().next();
		double expressionPval = firstSample.getBackgroundData().scoreWindow(gene).getScanPvalue();
		if(expressionPval > TwoSampleScanSkellamPeakCaller.EXPRESSION_PVALUE_CUTOFF) {
			logger.info(gene.getName() + "\tnot_expressed");
			return null;
		}
		
		findJointSignificantPeaks(gene, bedFileSampleLevelPeaks);
		return jointSignificantPeaks.get(gene);
	}
	
	/**
	 * Find and merge peaks that are jointly significant across samples in scan and skellam statistics, and cache set of peaks
	 * @param gene The gene
	 * @throws IOException
	 */
	private void findJointSignificantPeaks(Gene gene) throws IOException {
		findJointSignificantPeaks(gene, null);
	}
	
	/**
	 * Find and merge peaks that are jointly significant across samples in scan and skellam statistics, and cache set of peaks, and write sample level peaks to a bed file
	 * @param gene The gene
	 * @param bedFileSampleLevelPeaks File stream to write all sample level peaks to bed file
	 * @throws IOException
	 */
	private void findJointSignificantPeaks(Gene gene, FileWriter bedFileSampleLevelPeaks) throws IOException {
		
		TreeSet<Annotation> jointSignificantSampleLevelPeaks = new TreeSet<Annotation>();
		int geneSize = gene.getSize();
		
		
		// Get significant peaks from each sample and filter for joint significance
		for(TwoSampleScanSkellamAnalysis sample : pairedSamples) {
			
			
			
			double backgroundGeneCount = sample.getBackgroundData().getCount(gene);
			double signalGeneCount = sample.getSignalData().getCount(gene);
			
			Collection<Annotation> sigPeaks = sample.getSignificantPeaks(gene, false);
			if(sigPeaks == null) {
				logger.info(sample.getSignalName() + "\t" + gene.getName() + "\tno_sig_peaks");
				continue;
			}

			for(Annotation peak : sigPeaks) {
				
				if(bedFileSampleLevelPeaks != null) {
					Annotation peakCopy = peak.copy();
					peakCopy.setName(sample.getSignalName());
					bedFileSampleLevelPeaks.write(peakCopy.toBED(0, 0, 0) + "\n");
				}
				
				String logString = sample.getSignalName() + "\t" + gene.getName() + "\t";
				
				logString += peak.getChr() + ":" + peak.getStart() + "-" + peak.getEnd() + "\t";
				int peakSize = peak.getSize();
					
				// Check joint scan P value
				double jointScanPval = jointScanPval(gene, peak);
				if(jointScanPval > alphaScan) {
					logString += "joint_scan_pval_not_significant\t" + jointScanPval;
					logger.info(logString);
					continue;						
				}
					
				// Check joint skellam P value
				double backgroundLambda = backgroundGeneCount * peakSize / geneSize;
				double signalLambda = signalGeneCount * peakSize / geneSize;
				double jointSkellamPval = jointSkellamPval(peak, backgroundLambda, signalLambda);
				if(jointSkellamPval > alphaSkellam) {
					logString += "joint_skellam_pval_not_significant\t" + jointSkellamPval;
					logger.info(logString);
					continue;
				}
					
				// Add peak
				logString += "joint_p_values_significant";
				logger.info(logString);
				jointSignificantSampleLevelPeaks.add(peak);
			}
		}

		// Merge overlapping peaks
		Collection<Annotation> mergedWindows = AnnotationUtils.mergeOverlappingBlocksAnyOrientation(jointSignificantSampleLevelPeaks);
		
		// Add to cached set
		jointSignificantPeaks.put(gene, mergedWindows);
		
	}

	
	/**
	 * Find and store peaks for each gene and each paired sample
	 * @throws IOException
	 */
	private void findSampleLevelPeaks() throws IOException {
		sampleLevelPeaks = new TreeMap<Gene, Map<TwoSampleScanSkellamAnalysis, Collection<Annotation>>>();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				Map<TwoSampleScanSkellamAnalysis, Collection<Annotation>> peaksThisGene = new HashMap<TwoSampleScanSkellamAnalysis, Collection<Annotation>>();
				for(TwoSampleScanSkellamAnalysis sample : pairedSamples) {
					Collection<Annotation> peaks = sample.getSignificantPeaks(gene);
					peaksThisGene.put(sample, peaks);
				}
				sampleLevelPeaks.put(gene, peaksThisGene);
			}
		}
		hasSampleLevelPeaks = true;
	}
	
	/**
	 * Joint P value of skellam statistics for all paired samples
	 * @param region The region
	 * @param backgroundLambda Background lambda for the parent gene
	 * @param signalLambda Signal lambda for the parent gene
	 * @return Joint P value according to Fisher's combined probability test
	 * @throws IOException
	 */
	private double jointSkellamPval(Annotation region, double backgroundLambda, double signalLambda) throws IOException {
		ArrayList<Double> indPvals = new ArrayList<Double>();
		for(TwoSampleScanSkellamAnalysis sample : pairedSamples) {
			double pval = sample.skellamPval(region, backgroundLambda, signalLambda);
			indPvals.add(Double.valueOf(pval));
		}
		FisherCombinedProbabilityTest fisherTest = new FisherCombinedProbabilityTest(indPvals);
		return fisherTest.getCombinedPvalue();
	}
	
	/**
	 * Joint P value of scan statistics for all paired samples
	 * @param parentGene The parent gene
	 * @param window The region
	 * @return Joint P value according to Fisher's combined probability test
	 */
	private double jointScanPval(Gene parentGene, Annotation window) {
		ArrayList<Double> indPvals = new ArrayList<Double>();
		for(TwoSampleScanSkellamAnalysis sample : pairedSamples) {
			double pval = sample.signalScanPval(parentGene, window);
			indPvals.add(Double.valueOf(pval));
		}
		FisherCombinedProbabilityTest fisherTest = new FisherCombinedProbabilityTest(indPvals);
		return fisherTest.getCombinedPvalue();

	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("--background_bam", "Background bam file", true);
		p.addStringArg("--signal_bam_list", "File containing list of signal bam files", true);
		p.addStringArg("--genes", "Bed file of genes", true);
		p.addIntArg("--window", "Window size", false, DEFAULT_WINDOW_SIZE);
		p.addIntArg("--step", "Step size", false, DEFAULT_STEP_SIZE);
		p.addDoubleArg("--alpha_skellam", "Joint P value cutoff for skellam statistic", false, DEFAULT_ALPHA_SKELLAM);
		p.addDoubleArg("--alpha_scan", "Joint P value cutoff for scan statistic", false, DEFAULT_ALPHA_SCAN);
		p.addDoubleArg("--trim_quantile", "Quantile for individual sample peak trimming", false, DEFAULT_TRIM_QUANTILE);
		p.addStringArg("--output", "Output bed file of joint significant peaks", true);
		p.addStringArg("--out_sample_level_peaks", "Output bed file of sample level significant peaks", false, null);
		
		p.parse(args);
		
		String backgroundBamFile = p.getStringArg("--background_bam");
		String signalSampleBamList = p.getStringArg("--signal_bam_list");
		String bedFile = p.getStringArg("--genes");
		int window = p.getIntArg("--window");
		int step = p.getIntArg("--step");
		double skellamPvalCutoff = p.getDoubleArg("--alpha_skellam");
		double scanPvalCutoff = p.getDoubleArg("--alpha_scan");
		double trimMaxQuantile = p.getDoubleArg("--trim_quantile");
		String outBed = p.getStringArg("--output");
		String outSampleLevelBed = p.getStringArg("--out_sample_level_peaks");
		
		MultiSampleScanSkellamAnalysis ra = new MultiSampleScanSkellamAnalysis(backgroundBamFile, signalSampleBamList, bedFile, window, step, skellamPvalCutoff, scanPvalCutoff, trimMaxQuantile);
		FileWriter jointWriter = new FileWriter(outBed);
		FileWriter sampleLevelWriter = null;
		if(outSampleLevelBed != null) {
			sampleLevelWriter = new FileWriter(outSampleLevelBed);
		}
		
		for(String chr : ra.getGenes().keySet()) {
			for(Gene gene : ra.getGenes().get(chr)) {
				Collection<Annotation> sigPeaks = ra.getJointSignificantPeaks(gene, sampleLevelWriter);
				if(sigPeaks == null) continue;
				for(Annotation a : sigPeaks) {
					jointWriter.write(a.toBED(255, 0, 0) + "\n");
				}
			}
		}
		jointWriter.close();
		if(sampleLevelWriter != null) {
			sampleLevelWriter.close();
		}
		
	}

}

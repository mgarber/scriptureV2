/**
 * CLIP/Protect seq dataset for one sample
 */
package broad.pda.seq.clip;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.math.Statistics;
import broad.core.parser.StringParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.MappingQualityFilter;
import nextgen.core.readFilters.NumHitsFilter;

import broad.pda.seq.protection.*;

/**
 * @author shari
 *
 */
public class SampleData {

	protected String sampleName;
	protected TranscriptomeSpaceAlignmentModel data;
	protected TranscriptomeSpaceAlignmentModel maxFragmentLengthData;
	protected Map<Gene, ScanStatisticScore> maxFragmentLengthGeneScores;
	protected Map<Gene, ScanStatisticScore> geneScores;
	protected Map<Gene, Double> geneAvgCoverage;
	protected static Logger logger = Logger.getLogger(SampleData.class.getName());
	protected Map<String, Collection<Gene>> genesByChr;
	protected Map<String, Gene> genesByName;
	protected double expressionCutoffValue;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 100000;
	private static int DEFAULT_MAX_FRAGMENT_LENGTH = 150;
	protected boolean expressionByScanPval;
	private String originalBamFile;
	private boolean read1TranscriptionStrand;
	protected boolean fullyContainedReads;
	
	/**
	 * @param bamFile Bam file
	 * @param firstReadTranscriptionStrand Whether read 1 is in direction of transcription
	 * @param genes Genes by chromosome
	 * @param window Window size
	 * @param step Step size
	 * @param expressionCutoff P value cutoff for scan test of gene expression
	 * @param expByScanPval Expression is assessed by scan P value. If false, uses average read depth
	 * @param fullyContained Count fully contained reads only
	 * @throws IOException 
	 */
	public SampleData(String bamFile, boolean firstReadTranscriptionStrand, Map<String, Collection<Gene>> genes, double expressionCutoff, boolean expByScanPval, boolean fullyContained) throws IOException {
		fullyContainedReads = fullyContained;
		geneAvgCoverage = new TreeMap<Gene, Double>();
		originalBamFile = bamFile;
		read1TranscriptionStrand = firstReadTranscriptionStrand;
		StringParser p = new StringParser();
		p.parse(bamFile, "\\.");
		sampleName = p.asString(0);
		for(int i = 1 ; i < p.getFieldCount() - 1; i++) {
			sampleName += "." + p.asString(i);
		}
		expressionCutoffValue = expressionCutoff;
		expressionByScanPval = expByScanPval;
		genesByChr = genes;
		data = new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(genes));
		maxFragmentLengthData = new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(genes));
		
		// Read filters
		data.addFilter(new GenomicSpanFilter(DEFAULT_MAX_GENOMIC_SPAN));
		data.addFilter(new MappingQualityFilter(5,10));
		data.addFilter(new NumHitsFilter(1));
		
		maxFragmentLengthData.addFilter(new GenomicSpanFilter(DEFAULT_MAX_GENOMIC_SPAN));
		maxFragmentLengthData.addFilter(new MappingQualityFilter(5,10));
		maxFragmentLengthData.addFilter(new NumHitsFilter(1));
		maxFragmentLengthData.addFilter(new FragmentLengthFilter(maxFragmentLengthData.getCoordinateSpace(), DEFAULT_MAX_FRAGMENT_LENGTH));
		
		genesByName = new TreeMap<String, Gene>();
		for(String chr : genesByChr.keySet()) {
			for(Gene gene : genesByChr.get(chr)) {
				genesByName.put(gene.getName(), gene);
			}
		}
		geneScores = new TreeMap<Gene, ScanStatisticScore>();
		maxFragmentLengthGeneScores = new TreeMap<Gene, ScanStatisticScore>();
		logger.info("Instantiated sample data object. Name = " + sampleName);
	}
	
	@Override
	public int hashCode() {
		return sampleName.hashCode();
	}
	
	/**
	 * Set genome wide scan P value cutoff for expression of transcript
	 * @param expressionScanPvalCutoff P value cutoff for transcript expression against genomic background
	 */
	public void setExpressionScanPvalueCutoff(double expressionScanPvalCutoff) {
		expressionCutoffValue = expressionScanPvalCutoff;
	}

	/**
	 * Get genome wide scan P value cutoff for expression of transcript
	 * @return P value cutoff for transcript expression against genomic background
	 */
	public double getExpressionScanPvalueCutoff() {
		return expressionCutoffValue;
	}
	
	/**
	 * Get number of fragments mapping to gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The number of fragments mapping to the gene
	 */
	public double getGeneCount(Gene gene) {
		if(geneScores.containsKey(gene)) {
			return geneScores.get(gene).getCount();
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
		return score.getCount();
	}
	
	/**
	 * Get number of fragments mapping to gene after fragment length filter
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The number of fragments mapping to the gene
	 */
	public double getFragmentLengthGeneCount(Gene gene) {
		if(maxFragmentLengthGeneScores.containsKey(gene)) {
			return maxFragmentLengthGeneScores.get(gene).getCount();
		}
		ScanStatisticScore score = new ScanStatisticScore(maxFragmentLengthData, gene, fullyContainedReads);
		return score.getCount();
	}
	
	/**
	 * Get coordinate space wide scan P value of number of fragments mapping to the gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The scan P value of the number of fragments mapping to the gene with respect to teh coordinate space
	 */
	public double getGeneScanPval(Gene gene) {
		if(geneScores.containsKey(gene)) {
			return geneScores.get(gene).getScanPvalue();
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
		geneScores.put(gene, score);
		return score.getScanPvalue();
	}
	
	
	/**
	 * Get average coverage of gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The average coverage of the gene
	 */
	public double getGeneAverageCoverage(Gene gene) {
		if(geneAvgCoverage.containsKey(gene)) {
			return geneAvgCoverage.get(gene).doubleValue();
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
		geneScores.put(gene, score);
		double avgCoverage = score.getAverageCoverage(data);
		geneAvgCoverage.put(gene, Double.valueOf(avgCoverage));
		return avgCoverage;
	}
	
	/**
	 * Whether the gene is significantly expressed in the sample
	 * @param gene The gene
	 * @return Whether the gene is expressed at the given significance level
	 */
	public boolean isExpressed(Gene gene) {
		if(!geneScores.containsKey(gene)) {
			ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
			logger.debug("CHECK_GENE_EXPRESSION\t" + gene.getName());
			geneScores.put(gene, score);			
		}
		ScanStatisticScore score = geneScores.get(gene);
		if(expressionByScanPval) {
			return score.getScanPvalue() <= expressionCutoffValue;
		}
		double avgDepth = getGeneAverageCoverage(gene);
		logger.debug("cutoff=" + expressionCutoffValue + "\tcount=" + score.getCount() + "\tsize=" + score.getCoordinateSpace().getSize(score.getAnnotation()) + "\tscore=" + avgDepth);
		return avgDepth >= expressionCutoffValue;
	}
	
	
	/**
	 * Whether read 1 is in direction of transcription
	 * @return True iff first read transcription strand is true
	 */
	public boolean firstReadTranscriptionStrand() {
		return read1TranscriptionStrand;
	}
	
	/**
	 * Get enrichment of a window over a gene
	 * @param gene The gene
	 * @param window Window contained in the gene
	 * @return Enrichment of window over gene background
	 */
	public double getEnrichmentOverGene(Gene gene, Annotation window) {
		if(!gene.contains(window)) {
			throw new IllegalArgumentException("Gene must contain window.");
		}
		double geneAvgCov = getGeneAverageCoverage(gene);
		ScanStatisticScore windowScore = new ScanStatisticScore(data, window, data.getCount(gene), gene.getSize());
		double windowAvgCoverage = windowScore.getAverageCoverage(data);
		double enrichment = windowAvgCoverage / geneAvgCov;
		return enrichment;
	}

	/**
	 * Get the logger
	 * @return The logger
	 */
	@SuppressWarnings("static-method")
	public Logger getLogger() {
		return logger;
	}
	
	/**
	 * Get the sample name
	 * @return Sample name
	 */
	public String getSampleName() {
		return sampleName;
	}
	
	/**
	 * Get name of original bam file
	 * @return Bam file name
	 */
	public String getOriginalBamFile() {
		return originalBamFile;
	}
	
	/**
	 * Set the sample name
	 * @param name Sample name
	 */
	public void setSampleName(String name) {
		sampleName = name;
	}
	
/*	*//**
	 * Write the window scores to a file for future use
	 * @param m Multi sample binding site caller that this belongs to
	 * @throws IOException
	 *//*
	public void writeWindowScoresToFile(MultiSampleBindingSiteCaller m) throws IOException {
		windowScoreFile.writeWindowScoresToFile(m);
	}
*/	
	/**
	 * Get the alignment data
	 * @return Alignment model
	 */
	public TranscriptomeSpaceAlignmentModel getData() {
		return data;
	}
	
	/**
	 * Get the alignment data filtered by max fragment length
	 * @return Alignment model
	 */
	public TranscriptomeSpaceAlignmentModel getFragmentLengthFilterData() {
		return maxFragmentLengthData;
	}
	
	/**
	 * Trim the region to max contiguous subregion above a certain quantile
	 * @param window The region
	 * @param data Position level list of counts within the region
	 * @param quantile Quantile for trim max contiguous
	 * @return Trimmed region
	 */
	public static Annotation trimMaxContiguous(Annotation window, List<Double> data, double quantile) {
	
		String coverageString = "";
		for(Double d : data) {
			coverageString += d.toString() + " ";
		}
		
		logger.debug("WINDOW_TO_TRIM\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + "\t" + coverageString);
		
		if(window.getSize() != data.size()) {
			throw new IllegalArgumentException("Annotation and data must have same size. Name=" + window.getName() + " " + window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + " size=" + window.getSize() + " data_size=" + data.size());
		}
		
		
		double[] array = new double[data.size()];
		for(int i=0; i < data.size(); i++) {
			array[i] = data.get(i).doubleValue();
		}
		Collections.sort(data);
		
		double cutoff = Statistics.quantile(data, quantile);
		for(int j=0; j<array.length; j++){
			double d = array[j] - cutoff;
			array[j] = d;
		}

		double[] maxSum = MaximumContiguousSubsequence.maxSubSum3(array);
		
		logger.debug("TRIMMED_BOUNDARIES\t" + maxSum[1] + "-" + maxSum[2]);
	
		if(maxSum[0] > 0){
			int deltaStart = new Double(maxSum[1]).intValue();
			int deltaEnd =  new Double(data.size() - 1 - maxSum[2]).intValue();
			if(window.getStrand().equals(Strand.NEGATIVE)) {
			    int tmpStart = deltaStart;
			    deltaStart = deltaEnd;
			    deltaEnd = tmpStart;
			}
			window = window.trim(deltaStart, deltaEnd);
		}
		
		logger.debug("TRIMMED_WINDOW\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd());
		
		return window;
	}

	
}

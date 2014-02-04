/**
 * 
 */
package nextgen.core.normalize;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.protection.GenomeSpaceSampleData;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.utils.FileUtil;

/**
 * @author prussell
 *
 */
public class ComparativeCountsMatrix {
	
	/*
	 * Sample names
	 */
	private List<String> sampleNames;
	
	/*
	 * Map of sample name to normalization object (normalize by control counts)
	 */
	private Map<String, NormalizedCount> normalizationBySampleName;
	
	/**
	 * For each gene and each sample, the count ratio of sample to control
	 */
	private Map<Gene, Map<String, Double>> ratiosByGeneAndSampleName;
	
	/**
	 * For each sample, data object that can report whether a gene is expressed
	 */
	private Map<String, GenomeSpaceSampleData> sampleExpressionData;
	
	/**
	 * Data object that can report whether a gene is expressed in the control sample
	 */
	private GenomeSpaceSampleData controlExpressionData;
	
	/**
	 * Logger
	 */
	private static Logger logger = Logger.getLogger(ComparativeCountsMatrix.class.getName());
	
	/**
	 * Default genome wide scan P value cutoff for expression
	 */
	private static double DEFAULT_EXPRESSION_PVAL_CUTOFF = 0.01;
	
	/**
	 * @param controlBamFile Bam file of control alignments
	 * @param sampleBamListFile Bam file of alignments
	 * @param chrSizeFile Chromosome size file
	 * @param transcriptomeSpaceBedFile Gene bed file for transcriptome space
	 * @param fullyContained Only count fully contained alignments in annotations
	 * @param normalizationName Name of normalization to use for numerator and denominator
	 * @throws IOException
	 */
	public ComparativeCountsMatrix(String controlBamFile, String sampleBamListFile, String chrSizeFile, String transcriptomeSpaceBedFile, boolean fullyContained, String normalizationName) throws IOException {
		this(controlBamFile, FileUtil.fileLinesAsList(sampleBamListFile), chrSizeFile, BEDFileParser.loadDataByChr(new File(transcriptomeSpaceBedFile)), fullyContained, normalizationName);
	}
	
	/**
	 * @param controlBamFile Bam file of control alignments
	 * @param sampleBamFiles Bam file of alignments
	 * @param chrSizeFile Chromosome size file
	 * @param genes Genes by chromosome
	 * @param fullyContained Only count fully contained alignments in annotations
	 * @param normalizationName Name of normalization to use for numerator and denominator
	 * @throws IOException
	 */
	public ComparativeCountsMatrix(String controlBamFile, List<String> sampleBamFiles, String chrSizeFile, Map<String, Collection<Gene>> genes, boolean fullyContained, String normalizationName) throws IOException {
		
		TranscriptomeSpace coordSpace = new TranscriptomeSpace(genes);
		normalizationBySampleName = new TreeMap<String, NormalizedCount>();
		ratiosByGeneAndSampleName = new TreeMap<Gene, Map<String, Double>>();
		
		// Read control alignment data
		logger.info("");
		logger.info("Reading control data from bam file " + controlBamFile + "...");
		AlignmentModel controlData = new AlignmentModel(controlBamFile, coordSpace);
		
		// Establish expression object for control data
		controlExpressionData = new GenomeSpaceSampleData(controlBamFile, chrSizeFile, genes, 1, 1, DEFAULT_EXPRESSION_PVAL_CUTOFF, fullyContained);
		
		// Read alignment data for all samples
		sampleExpressionData = new TreeMap<String, GenomeSpaceSampleData>();
		sampleNames = new ArrayList<String>();
		for(String sampleBam : sampleBamFiles) {
			logger.info("Reading sample data from bam file " + sampleBam + "...");
			String sampleName = sampleBam.replaceAll(".bam", "");
			sampleNames.add(sampleName);
			// Expression object
			sampleExpressionData.put(sampleName, new GenomeSpaceSampleData(sampleBam, chrSizeFile, genes, 1, 1, DEFAULT_EXPRESSION_PVAL_CUTOFF, fullyContained));
			AlignmentModel sampleData = new AlignmentModel(sampleBam, coordSpace);
			NormalizedCount normalization = null;
			if(normalizationName.equals(NormalizedCount.CROSS_SAMPLE_BINOMIAL_ENRICHMENT_SCORE_NAME)) {
				normalization = new CrossSampleBinomialEnrichmentScore(sampleData, controlData, fullyContained);
			} else if (normalizationName.equals(NormalizedCount.MAX_DEPTH_NAME)) {
				normalization = new CrossSampleRatio(sampleData, controlData, new MaxDepthNormalization(sampleData, fullyContained), new MaxDepthNormalization(controlData, fullyContained), fullyContained);
			} else if (normalizationName.equals(NormalizedCount.RAW_COUNTS_NAME)){
				normalization = new CrossSampleRatio(sampleData, controlData, new RawCounts(sampleData, fullyContained), new RawCounts(controlData, fullyContained), fullyContained);
			} else if (normalizationName.equals(NormalizedCount.FRAGMENT_SIZE_NORMALIZATION_NAME)) {
				normalization = new CrossSampleRatio(sampleData, controlData, new FragmentSizeNormalization(sampleData, fullyContained), new FragmentSizeNormalization(controlData, fullyContained), fullyContained);
			} else {
				String possibleNormalizations = "";
				possibleNormalizations += NormalizedCount.RAW_COUNTS_NAME + " ";
				possibleNormalizations += NormalizedCount.CROSS_SAMPLE_BINOMIAL_ENRICHMENT_SCORE_NAME + " ";
				possibleNormalizations += NormalizedCount.MAX_DEPTH_NAME + " ";
				possibleNormalizations += NormalizedCount.FRAGMENT_SIZE_NORMALIZATION_NAME + " ";
				throw new IllegalArgumentException("Invalid normalization name. Options are: " + possibleNormalizations);
			}
			normalizationBySampleName.put(sampleName, normalization);
		}
		logger.info("Using normalization: " + normalizationName);
		logger.info("Done constructing object.");
	}
	
	/**
	 * For each sample get ratio of count to control count over the region
	 * @param region Region
	 * @return Map of sample name to count ratio
	 */
	private Map<String, Double> getRatioBySampleName(Gene region) {
		if(!ratiosByGeneAndSampleName.containsKey(region)) {
			calculateRatioBySampleName(region);
		}
		return ratiosByGeneAndSampleName.get(region);
	}
	
	/**
	 * Calculate and store count ratios for each sample over the region
	 * @param region The region
	 */
	private void calculateRatioBySampleName(Gene region) {
		boolean controlExpressed = controlExpressionData.isExpressed(region);
		Map<String, Double> ratios = new TreeMap<String, Double>();
		for(String sample : sampleNames) {
			boolean sampleExpressed = sampleExpressionData.get(sample).isExpressed(region);
			if(!controlExpressed || !sampleExpressed) {
				ratios.put(sample, Double.valueOf(Double.NaN));
			} else {
				double ratio = normalizationBySampleName.get(sample).getNormalizedCount(region);
				ratios.put(sample, Double.valueOf(ratio));
			}
		}
		ratiosByGeneAndSampleName.put(region, ratios);
	}
	
	/**
	 * Get matrix of count ratios
	 * @param regionBedFile Bed file of regions to use
	 * @return Matrix of count ratios. Columns are samples; rows are regions.
	 * @throws IOException
	 */
	public MatrixWithHeaders getMatrix(String regionBedFile) throws IOException {
		return getMatrix(regionBedFile, null);
	}
	
	
	/**
	 * Get matrix of count ratios for one or all chromosomes
	 * @param regionBedFile Bed file of regions to use
	 * @param chr Single chromosome to include or null if all chromosomes
	 * @return Matrix of count ratios. Columns are samples; rows are regions.
	 * @throws IOException
	 */
	public MatrixWithHeaders getMatrix(String regionBedFile, String chr) throws IOException {
		
		logger.info("");
		logger.info("Getting matrix of count ratios...");
		
		// Get regions from bed file
		Map<String, Collection<Gene>> regions = BEDFileParser.loadDataByChr(new File(regionBedFile));
		
		// Establish set of chromosomes to use
		Collection<String> chrs = new TreeSet<String>();
		if(chr != null) {
			chrs.add(chr);
		} else {
			chrs.addAll(regions.keySet());
		}

		// Establish region names for matrix headers
		Collection<String> regionNames = new TreeSet<String>();
		for(String c : chrs) {
			for(Gene region : regions.get(c)) {
				String regionName = region.getName();
				if(regionNames.contains(regionName)) {
					logger.warn("Duplicate region name " + regionName);
				}
				regionNames.add(region.getName());
			}
		}
		List<String> regionNamesList = new ArrayList<String>();
		regionNamesList.addAll(regionNames);
		
		// Instantiate matrix
		MatrixWithHeaders matrix = new MatrixWithHeaders(regionNamesList, sampleNames);
		
		// Calculate values
		for(String c : chrs) {
			logger.info(c);
			for(Gene region : regions.get(c)) {
				logger.debug("Getting entries for gene " + region.getName());
				String regionName = region.getName();
				Map<String, Double> ratios = getRatioBySampleName(region);
				for(String sampleName : ratios.keySet()) {
					matrix.set(regionName, sampleName, ratios.get(sampleName).doubleValue());
				}
			}
		}
		
		return matrix;
		
	}
	
	/**
	 * Write table of ratios to file for one or all chromosomes
	 * @param matrix Matrix
	 * @param outFile Output file
	 * @throws IOException
	 */
	public static void writeMatrix(MatrixWithHeaders matrix, String outFile) throws IOException {
		logger.info("");
		logger.info("Writing table to file " + outFile + "...");
		
		matrix.write(outFile);
		
		logger.info("Done writing table.");
	}
	
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-cb", "Denominator bam file name", true);
		p.addStringArg("-lb", "File containing list of numerator bam files", true);
		p.addStringArg("-bt", "Bed file for transcriptome space", true);
		p.addStringArg("-br", "Bed file of regions for count ratios", true);
		p.addStringArg("-o", "Output table", true);
		p.addStringArg("-c", "Single chromosome to use", false, null);
		p.addStringArg("-cs", "Chromosome size file", true);
		p.addBooleanArg("-nm", "Normalize matrix by column median", false, true);
		p.addBooleanArg("-l10", "Take log10 of each matrix entry", false, true);
		p.addBooleanArg("-fc", "Only count fully contained reads", true);
		p.addStringArg("-n", "Type of normalization", false, NormalizedCount.RAW_COUNTS_NAME);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.parse(args);
		boolean debug = p.getBooleanArg("-d");
		if(debug) logger.setLevel(Level.DEBUG);
		String controlBam = p.getStringArg("-cb");
		String bamListFile = p.getStringArg("-lb");
		String transcriptomeBed = p.getStringArg("-bt");
		String regionBed = p.getStringArg("-br");
		String outTable = p.getStringArg("-o");
		String chr = p.getStringArg("-c");
		String chrSizeFile = p.getStringArg("-cs");
		boolean normalizeByMedian = p.getBooleanArg("-nm");
		boolean takeLog10 = p.getBooleanArg("-l10");
		boolean fullyContained = p.getBooleanArg("-fc");
		String normalization = p.getStringArg("-n");
		
		ComparativeCountsMatrix c = new ComparativeCountsMatrix(controlBam, bamListFile, chrSizeFile, transcriptomeBed, fullyContained, normalization);	
		MatrixWithHeaders matrix = c.getMatrix(regionBed, chr);
		if(normalizeByMedian) {
			matrix.normalizeColumnsByMedian();
		}
		if(takeLog10) {
			matrix.log10();
		}
		writeMatrix(matrix, outTable);

		logger.info("");
		logger.info("All done.");
	}

}

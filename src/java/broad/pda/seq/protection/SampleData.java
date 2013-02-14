/**
 * Protect seq dataset for one sample
 */
package broad.pda.seq.protection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.math.Statistics;
import broad.core.parser.StringParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.readFilters.GenomicSpanFilter;

/**
 * @author prussell
 *
 */
public class SampleData {

	protected String sampleName;
	protected TranscriptomeSpaceAlignmentModel data;
	private Map<Gene, ScanStatisticScore> geneScores;
	protected Map<Gene, Map<Annotation, ScanStatisticScore>> windowScores;
	protected int windowSize;
	protected int stepSize;
	protected static Logger logger = Logger.getLogger(SampleData.class.getName());
	private WindowProcessor<ScanStatisticScore> processor;
	protected Map<String, Collection<Gene>> genesByChr;
	protected Map<String, Gene> genesByName;
	protected double scanPvalAlpha;
	private CachedScoreFile windowScoreFile;
	private boolean gotWindowScoresFromFile;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	
	/**
	 * @param bamFile Bam file
	 * @param genes Genes by chromosome
	 * @param window Window size
	 * @param step Step size
	 * @param alpha P value cutoff for scan test of gene expression
	 * @throws IOException 
	 */
	public SampleData(String bamFile, Map<String, Collection<Gene>> genes, int window, int step, double alpha) throws IOException {
		StringParser p = new StringParser();
		p.parse(bamFile, "\\.");
		sampleName = p.asString(0);
		for(int i = 1 ; i < p.getFieldCount() - 1; i++) {
			sampleName += "." + p.asString(i);
		}
		scanPvalAlpha = alpha;
		genesByChr = genes;
		data = new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(genes));
		data.addFilter(new GenomicSpanFilter(DEFAULT_MAX_GENOMIC_SPAN));
		processor = new ScanStatisticScore.Processor(data);
		genesByName = new TreeMap<String, Gene>();
		for(String chr : genesByChr.keySet()) {
			for(Gene gene : genesByChr.get(chr)) {
				genesByName.put(gene.getName(), gene);
			}
		}
		geneScores = new TreeMap<Gene, ScanStatisticScore>();
		windowScores = new TreeMap<Gene, Map<Annotation, ScanStatisticScore>>();
		windowSize = window;
		stepSize = step;
		logger.info("Instantiated sample data object. Name = " + sampleName + ", window size = " + windowSize + ", step size = " + stepSize);
		windowScoreFile = new CachedScoreFile(getDefaultWindowScoreFileName());
		gotWindowScoresFromFile = windowScoreFile.readWindowScoresFromFile();
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
		ScanStatisticScore score = new ScanStatisticScore(data, gene);
		geneScores.put(gene, score);
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
		ScanStatisticScore score = new ScanStatisticScore(data, gene);
		geneScores.put(gene, score);
		return score.getScanPvalue();
	}
	
	
	/**
	 * Get coordinate space wide scan P value of number of fragments mapping to the gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The scan P value of the number of fragments mapping to the gene with respect to the coordinate space
	 */
	public double getGeneAverageCoverage(Gene gene) {
		if(geneScores.containsKey(gene)) {
			return geneScores.get(gene).getAverageCoverage();
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene);
		geneScores.put(gene, score);
		return score.getAverageCoverage();
	}
	
	/**
	 * Whether the gene is significantly expressed in the sample
	 * @param gene The gene
	 * @return Whether the gene is expressed at the given significance level
	 */
	public boolean isExpressed(Gene gene) {
		if(geneScores.containsKey(gene)) {
			return geneScores.get(gene).getScanPvalue() < scanPvalAlpha;
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene);
		geneScores.put(gene, score);
		return geneScores.get(gene).getScanPvalue() < scanPvalAlpha;
	}
	
	/**
	 * Get number of fragments mapping to gene window
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @param window The window
	 * @return The number of fragments mapping to the gene
	 */
	public double getWindowCount(Gene gene, Annotation window) {
		if(windowScores.containsKey(gene)) {
			return windowScores.get(gene).get(window).getCount();
		}
		computeWindowScores(gene);
		return windowScores.get(gene).get(window).getCount();
	}
	
	
	
	/**
	 * Get scores for each window in the gene
	 * @param gene The gene
	 * @return The scores by window
	 */
	public Map<Annotation, ScanStatisticScore> getWindowScores(Gene gene) {
		if(windowScores.containsKey(gene)) {
			return windowScores.get(gene);
		}
		computeWindowScores(gene);
		return windowScores.get(gene);		
	}
	
	/**
	 * Get the default name of the window score file in the current directory
	 * @return The file name
	 */
	private String getDefaultWindowScoreFileName() {
		return getDefaultWindowScoreFileName(".");
	}
	
	/**
	 * Whether the window scores were read from an existing file
	 * @return Whether the window scores were read from an existing file
	 */
	public boolean gotWindowScoresFromFile() {
		return gotWindowScoresFromFile;
	}
	
	/**
	 * Get the default name of the window score file in a specified directory
	 * @param directory The directory
	 * @return The file path
	 */
	private String getDefaultWindowScoreFileName(String directory) {
		String name = "window_scores_" + sampleName + "_" + windowSize + "_" + stepSize;
		return directory + "/" + name;
	}
	

	
	/**
	 * Compute and store scores for each window of gene
	 * @param gene The gene
	 */
	private void computeWindowScores(Gene gene) {
		logger.info("Computing window scores for sample " + sampleName + " and gene " + gene.getName());
		Map<Annotation, ScanStatisticScore> scores = new TreeMap<Annotation, ScanStatisticScore>();
		if(gene.getSize() < windowSize) {
			logger.info(gene.getName() + " is smaller than window size. Not computing window binding site scores.");
			windowScores.put(gene, scores);
			return;
			// TODO should be able to score genes that are smaller than window size
		}		
		WindowScoreIterator<ScanStatisticScore> iter = data.scan(gene, windowSize, windowSize - stepSize, processor);
		double geneTotal = getGeneCount(gene);
		double geneLength = gene.getSize();
		while(iter.hasNext()) {
			ScanStatisticScore score = iter.next();
			Annotation window = score.getAnnotation();
			double regionLength = window.getSize();
			double regionTotal = data.getCount(window);
			score.setGlobalLength(geneLength);
			score.setRegionLength(regionLength);
			score.setTotal(geneTotal);
			score.setRegionTotal(regionTotal);
			score.refreshScanPvalue(data);
			scores.put(window, score);
		}
		windowScores.put(gene, scores);
	}
	
	
	/**
	 * Get the sample name
	 * @return Sample name
	 */
	public String getSampleName() {
		return sampleName;
	}
	
	/**
	 * Set the sample name
	 * @param name Sample name
	 */
	public void setSampleName(String name) {
		sampleName = name;
	}
	
	/**
	 * Write the window scores to a file for future use
	 * @throws IOException
	 */
	public void writeWindowScoresToFile(MultiSampleBindingSiteCaller m) throws IOException {
		windowScoreFile.writeWindowScoresToFile(m);
	}
	
	/**
	 * Get the alignment data
	 * @return Alignment model
	 */
	public TranscriptomeSpaceAlignmentModel getData() {
		return data;
	}
	
	/**
	 * Trim the region to max contiguous subregion above a certain quantile
	 * @param window The region
	 * @param data Position level list of counts within the region
	 * @param quantile Quantile for trim max contiguous
	 * @return Trimmed region
	 */
	public static Annotation trimMaxContiguous(Annotation window, List<Double> data, double quantile) {
	
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
	
		if(maxSum[0] > 0){
			int deltaStart = new Double(maxSum[1]).intValue();
			int deltaEnd =  new Double(data.size() - 1 - maxSum[2]).intValue();
			if(window.getStrand().equals(Strand.NEGATIVE)) {
			    int tmpStart = deltaStart;
			    deltaStart = deltaEnd;
			    deltaEnd = tmpStart;
			}
			window.trim(deltaStart, deltaEnd);
		}
		return window;
	}

	
	private class CachedScoreFile {
		
		
		private String TOTAL_READS_IDENTIFIER = "total_reads";
		private String WINDOW_SIZE_IDENTIFIER = "window_size";
		private String STEP_SIZE_IDENTIFIER = "step_size";
		private String TOTAL_GENE_LENGTH_IDENTIFIER = "total_gene_length";
		private String NUM_GENES_IDENTIFIER = "num_genes";
		private String GENE_IDENTIFIER = "gene";
		private String NOT_EXPRESSED_IDENTIFIER = "not_expressed";
		private String DELIMITER = "=";
		private String GLOBAL_LENGTH_IDENTIFIER = "global_length";
		private String REGION_LENGTH_IDENTIFIER = "region_length";
		private String GLOBAL_TOTAL_IDENTIFIER = "global_total";
		private String REGION_TOTAL_IDENTIFIER = "region_total";
		private String SCAN_P_VALUE_IDENTIFIER = "scan_pval";
		private String NO_WINDOWS_IDENTIFIER = "no_windows";
		
		private StringParser stringParser;
		private String fileName;
		
		public CachedScoreFile(String file) {
			fileName = file;
			stringParser = new StringParser();
		}
		
		/**
		 * Read and store all window scores from file
		 * First check if file is valid for this dataset and parameters
		 * @return Whether the scores were read and stored from the file
		 * @throws IOException 
		 */
		public boolean readWindowScoresFromFile() throws IOException {
			
			logger.info("Trying to read window scores from file.");
			
			// Check if cached file is valid
			boolean validated = validateFile();
			if(!validated) {
				logger.warn("Window score file " + fileName + " not valid. Not reading scores from file.");
				return false;
			}
			
			windowScores.clear();
			FileReader r = new FileReader(fileName);
			BufferedReader b = new BufferedReader(r);
			
			TreeSet<String> genesDone = new TreeSet<String>();
			
			while(b.ready()) {
				String line = b.readLine();
				stringParser.parse(line);
				if(!stringParser.asString(0).equals(GENE_IDENTIFIER)) continue;
				Gene gene = getGeneFromLine(line);
				genesDone.add(gene.getName());
				if(genesDone.size() % 1000 == 0) {
					logger.info("Finished " + genesDone.size() + " genes.");
				}
				if(!windowScores.containsKey(gene)) {
					Map<Annotation, ScanStatisticScore> map = new TreeMap<Annotation, ScanStatisticScore>();
					windowScores.put(gene, map);
				}
				if(noWindows(line)) {
					continue;
				}
				Annotation window = getWindowFromLine(line);
				ScanStatisticScore score = getScoreFromLine(line);
				windowScores.get(gene).put(window, score);
				
			}
			
			r.close();
			b.close();
			logger.info("Got window scores from file " + fileName + ".");
			return true;
		}
		
		/**
		 * Write scores to file
		 * @throws IOException 
		 */
		public void writeWindowScoresToFile(MultiSampleBindingSiteCaller m) throws IOException {
			logger.info("Writing window scores to file " + fileName);
			FileWriter w = new FileWriter(fileName);
			w.write(TOTAL_READS_IDENTIFIER + "\t" + data.getGlobalNumReads() + "\n");
			w.write(NUM_GENES_IDENTIFIER + "\t" + data.getReferenceNames().size() + "\n");
			w.write(TOTAL_GENE_LENGTH_IDENTIFIER + "\t" + data.getGlobalLength() + "\n");
			w.write(WINDOW_SIZE_IDENTIFIER + "\t" + windowSize + "\n");
			w.write(STEP_SIZE_IDENTIFIER + "\t" + stepSize + "\n");
			for(String geneName : genesByName.keySet()) {
				Gene gene = genesByName.get(geneName);
				if(!m.isExpressed(gene)) {
					w.write(GENE_IDENTIFIER + "\t" + geneName + "\t" + NOT_EXPRESSED_IDENTIFIER + "\n");
					continue;
				}
				if(windowScores.get(gene).isEmpty()) {
					w.write(GENE_IDENTIFIER + "\t" + geneName + "\t" + NO_WINDOWS_IDENTIFIER + "\n");
					continue;
				}
				for(Annotation window : windowScores.get(gene).keySet()) {
					ScanStatisticScore score = windowScores.get(gene).get(window);
					String lineToWrite = GENE_IDENTIFIER + "\t" + geneName + "\t" + window.getChr() + "\t" + window.getStart() + "\t" + window.getEnd() + "\t";
					lineToWrite += GLOBAL_LENGTH_IDENTIFIER + DELIMITER + score.getGlobalLength() + "\t";
					lineToWrite += REGION_LENGTH_IDENTIFIER + DELIMITER + score.getRegionLength() + "\t";
					lineToWrite += GLOBAL_TOTAL_IDENTIFIER + DELIMITER + score.getTotal() + "\t";
					lineToWrite += REGION_TOTAL_IDENTIFIER + DELIMITER + score.getRegionTotal() + "\t";
					lineToWrite += SCAN_P_VALUE_IDENTIFIER + DELIMITER + score.getScanPvalue();
					w.write(lineToWrite + "\n");
				}
			}
			w.close();
		}
		
		/**
		 * Get the gene described on the line of scores
		 * @param line The line in file
		 * @return The gene
		 */
		private Gene getGeneFromLine(String line) {
			stringParser.parse(line);
			if(!stringParser.asString(0).equals(GENE_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid gene score: " + line);
			}
			return genesByName.get(stringParser.asString(1));
		}
		
		/**
		 * Whether the line indicates that the gene has no windows
		 * @param line The line
		 * @return Whether the line contains the no windows identifier
		 */
		private boolean noWindows(String line) {
			stringParser.parse(line);
			if(!stringParser.asString(0).equals(GENE_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid gene score: " + line);
			}
			
			if(stringParser.asString(2).equals(NO_WINDOWS_IDENTIFIER)) {
				return true;
			}
			return false;
			
		}
		
		/**
		 * Get the window described on the line of scores
		 * @param line The line in file
		 * @return The window
		 */
		private Annotation getWindowFromLine(String line) {
			
			stringParser.parse(line);
			if(!stringParser.asString(0).equals(GENE_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid gene score: " + line);
			}
			
			if(stringParser.asString(2).equals(NOT_EXPRESSED_IDENTIFIER)) {
				return null;
			}
			
			String chr = stringParser.asString(2);
			int start = stringParser.asInt(3);
			int end = stringParser.asInt(4);
			return new BasicAnnotation(chr, start, end);
			
		}
		
		/**
		 * Get scan statistic score described on the line of scores
		 * @param line The line in file
		 * @return The score
		 */
		private ScanStatisticScore getScoreFromLine(String line) {

			stringParser.parse(line);
			if(!stringParser.asString(0).equals(GENE_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid gene score: " + line);
			}
			
			if(stringParser.asString(2).equals(NOT_EXPRESSED_IDENTIFIER)) {
				return null;
			}
			
			StringParser p = new StringParser();
			
			p.parse(stringParser.asString(5),DELIMITER);
			if(!p.asString(0).equals(GLOBAL_LENGTH_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid global length: " + line);
			}
			double globalLength = p.asDouble(1);
			
			p.parse(stringParser.asString(6),DELIMITER);
			if(!p.asString(0).equals(REGION_LENGTH_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid region length: " + line);
			}
			double regionLength = p.asDouble(1);
			
			p.parse(stringParser.asString(7),DELIMITER);
			if(!p.asString(0).equals(GLOBAL_TOTAL_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid global total: " + line);
			}
			double globalTotal = p.asDouble(1);
			
			p.parse(stringParser.asString(8),DELIMITER);
			if(!p.asString(0).equals(REGION_TOTAL_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid region total: " + line);
			}
			double regionTotal = p.asDouble(1);
			
			p.parse(stringParser.asString(9),DELIMITER);
			if(!p.asString(0).equals(SCAN_P_VALUE_IDENTIFIER)) {
				throw new IllegalArgumentException("Line does not contain valid scan P value: " + line);
			}
			double scanPval = p.asDouble(1);

			ScanStatisticScore score = new ScanStatisticScore(getWindowFromLine(line));
			score.setGlobalLength(globalLength);
			score.setRegionLength(regionLength);
			score.setRegionTotal(regionTotal);
			score.setCount(regionTotal);
			score.setTotal(globalTotal);
			score.setScanPvalue(scanPval);
			
			return score;
			
		}
		
		/**
		 * Whether the file exists and represents the same data and parameters as this run
		 * @return True iff the file is valid
		 * @throws IOException 
		 */
		private boolean validateFile() throws IOException {
			
			File file = new File(fileName);
			if(!file.exists()) return false;
			
			FileReader r = new FileReader(fileName);
			BufferedReader b = new BufferedReader(r);

			try {
				
				// Compare total number of reads to the data
				String totalReadsLine = b.readLine();
				stringParser.parse(totalReadsLine);
				if(!stringParser.asString(0).equals(TOTAL_READS_IDENTIFIER)) crashIncorrectFileFormat();
				double totalReads = stringParser.asDouble(1);
				if(totalReads != data.getGlobalNumReads()) {
					logger.warn("Cached score file is not valid because total reads do not match data.");
					return false;
				}
				
				// Compare total number of genes to the transcriptome space
				String numGenesLine = b.readLine();
				stringParser.parse(numGenesLine);
				if(!stringParser.asString(0).equals(NUM_GENES_IDENTIFIER)) crashIncorrectFileFormat();
				long numGenes = stringParser.asInt(1);
				if(numGenes != data.getReferenceNames().size()) {
					logger.warn("Cached score file is not valid because number of genes does not match transcriptome space.");
					return false;
				}
				
				// Compare total length of genes to the transcriptome space
				String totalGeneLengthLine = b.readLine();
				stringParser.parse(totalGeneLengthLine);
				if(!stringParser.asString(0).equals(TOTAL_GENE_LENGTH_IDENTIFIER)) crashIncorrectFileFormat();
				long totalLength = Double.valueOf(stringParser.asDouble(1)).longValue();
				if(totalLength != data.getGlobalLength()) {
					logger.warn("Cached score file is not valid because total length of genes does not match transcriptome space.");
					return false;
				}
				
				// Compare window size
				String windowSizeLine = b.readLine();
				stringParser.parse(windowSizeLine);
				if(!stringParser.asString(0).equals(WINDOW_SIZE_IDENTIFIER)) crashIncorrectFileFormat();
				int window = stringParser.asInt(1);
				if(window != windowSize) {
					logger.warn("Cached score file is not valid because window size does not match.");
					return false;
				}
				
				// Compare step size
				String stepSizeLine = b.readLine();
				stringParser.parse(stepSizeLine);
				if(!stringParser.asString(0).equals(STEP_SIZE_IDENTIFIER)) crashIncorrectFileFormat();
				int step = stringParser.asInt(1);
				if(step != stepSize) {
					logger.warn("Cached score file is not valid because step size does not match.");
					return false;
				}
				
				Collection<String> geneNames = data.getReferenceNames();
				Collection<String> alreadyChecked = new TreeSet<String>();
				while(b.ready()) {
					String line = b.readLine();
					stringParser.parse(line);
					if(!stringParser.asString(0).equals(GENE_IDENTIFIER)) crashIncorrectFileFormat();
					String geneName = stringParser.asString(1);
					if(alreadyChecked.contains(geneName)) continue;
					if(!geneNames.contains(geneName)) {
						logger.warn("Cached score file is not valid because gene " + geneName + " is not in transcriptome space.");
						return false;
					}
					geneNames.remove(geneName);
					alreadyChecked.add(geneName);
				}
				if(!geneNames.isEmpty()) {
					logger.warn("Cached score file is not valid because genes are missing.");
					return false;
				}
				
			} catch(Exception e) {
				e.printStackTrace();
				crashIncorrectFileFormat();
			}
			b.close();
			r.close();
			return true;
			
		}
		
		/**
		 * Throw exception and print correct file format
		 */
		private void crashIncorrectFileFormat() {
			logger.error("Wrong format for cached score file. Sample format:");
			logger.error("total_reads	23434888");
			logger.error("num_genes	20000");
			logger.error("total_gene_length	5000000");
			logger.error("window_size	100");
			logger.error("step_size	10");
			logger.error("gene	NM12345	not_expressed");
			logger.error("gene NM123456	not_expressed");
			logger.error("gene NM123457	chr10	12345678	12345778	global_length=1000	region_length=1000	global_total=5555	region_total=6678	scan_pval=0.05");
			throw new IllegalArgumentException("Cached score file is in wrong format.");
		}
		
		
	}
	
	
	
}

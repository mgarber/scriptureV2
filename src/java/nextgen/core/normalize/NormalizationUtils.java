package nextgen.core.normalize;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.job.LSFJob;
import nextgen.core.writers.WigWriter;

/**
 * @author prussell
 *
 */
public class NormalizationUtils {
	
	private static Logger logger = Logger.getLogger(NormalizationUtils.class.getName());
	
	/**
	 * Write a bed file with score set to the normalized count for the feature
	 * @param inputFeatureBed Bed file of features
	 * @param normalization Normalization
	 * @param outFeatureBed Output bed file
	 * @throws IOException
	 */
	public static void writeFeatureEnrichmentBed(String inputFeatureBed, NormalizedCount normalization, String outFeatureBed) throws IOException {
		writeFeatureEnrichmentBed(BEDFileParser.loadDataByChr(inputFeatureBed), normalization, outFeatureBed);
	}
	
	/**
	 * Write a bed file with score set to the normalized count for the feature
	 * @param features Features
	 * @param normalization Normalization
	 * @param outFeatureBed Output bed file
	 * @throws IOException
	 */
	public static void writeFeatureEnrichmentBed(Map<String, Collection<Gene>> features, NormalizedCount normalization, String outFeatureBed) throws IOException {
		logger.info("");
		logger.info("Writing normalized enrichment for each feature to file " + outFeatureBed + "...");
		FileWriter w = new FileWriter(outFeatureBed);
		int numDone = 0;
		int skipped = 0;
		for(String chr : features.keySet()) {
			logger.info(chr);
			if(!features.containsKey(chr)) {
				continue;
			}
			for(Gene feature : features.get(chr)) {
				if(feature == null) {
					continue;
				}
				try {
					double count = normalization.getNormalizedCount(feature);
					Gene featureWithScore = feature.copy();
					featureWithScore.setBedScore(count);
					featureWithScore.setScore(count);
					w.write(featureWithScore.toBED() + "\n");
				} catch (IllegalArgumentException e) {
					skipped++;
					numDone++;
					if(numDone % 100 == 0) {
						logger.info("Finished " + numDone + " features. Skipped " + skipped + " that do not have parent.");
					}
					continue;
				}
				numDone++;
				if(numDone % 100 == 0) {
					logger.info("Finished " + numDone + " features. Skipped " + skipped + " that do not have parent.");
				}
			}
		}
		w.close();
	}
	
	public static void makeBigWig(String wigPrefix, String wigToBigWig, String chrSizeFile) throws IOException, InterruptedException {
		logger.info("");
		logger.info("Making bigwig file for wig file " + wigPrefix + "...");
		String wig = wigPrefix + ".wig";
		String bw = wigPrefix + ".bw";
		String cmmd = wigToBigWig + " " + wig + " " + chrSizeFile + " " + bw;
		String jobID = Long.valueOf(System.currentTimeMillis()).toString();
		String bsubOut = "make_bigwig_" + jobID + ".bsub";
		logger.info("Running command " + cmmd);
		LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOut, "hour", 4);
		job.submit();
		logger.info("Waiting for job to finish...");
		job.waitFor();
		logger.info("Done writing bigwig to file " + bw + ".");
	}
	

	
	public static void writePositionLevelDataAllGenes(NormalizedCount normalizedCount, Map<String, Collection<Gene>> features, String outFilePrefix, String wigToBigWig, String chrSizeFile) throws IOException, InterruptedException {
		String outFile = outFilePrefix + ".wig";
		logger.info("");
		logger.info("Writing position level data for all genes to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		int numDone = 0;
		for(String c : features.keySet()) {
			logger.info("Calculating data for chromosome " + c + "...");
			TreeMap<Integer, Double> toWrite = new TreeMap<Integer, Double>();
			for(Gene gene : features.get(c)) {
				numDone++;
				if(numDone % 100 == 0) {
					logger.info("Finished " + numDone + " genes.");
				}
				toWrite.putAll(normalizedCount.getNormalizedCountsByPosition(gene));
			}
			logger.info("Done calculating chromosome " + c + ". Writing...");
			WigWriter.write(w, c, toWrite, false);
		}
		w.close();
		logger.info("Done writing wig file.");
		if(wigToBigWig != null && chrSizeFile != null) {
			makeBigWig(outFilePrefix, wigToBigWig, chrSizeFile);
		}
	}


	
}

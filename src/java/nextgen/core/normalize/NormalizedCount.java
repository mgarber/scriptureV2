/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;

import nextgen.core.annotation.Annotation;

/**
 * @author prussell
 *
 */
public interface NormalizedCount {
	
	public static String RAW_COUNTS_NAME = "raw_counts";
	public static String TRANSCRIPT_AVERAGE_NORMALIZATION_NAME = "transcript_average";
	public static String CROSS_SAMPLE_RATIO_NAME = "cross_sample_ratio";
	public static String CROSS_SAMPLE_TRANSCRIPT_AVERAGE = "cross_sample_transcript_average";
	public static String MAX_DEPTH_NAME = "max_depth";
	public static String FRAGMENT_SIZE_NORMALIZATION_NAME = "fragment_size";
	public static String CROSS_SAMPLE_BINOMIAL_ENRICHMENT_SCORE_NAME = "cross_sample_binomial_enrichment";
	
	/**
	 * Get normalized count over a region
	 * @param region The region
	 * @return The normalized count of the entire region
	 */
	public double getNormalizedCount(Annotation region);
	
	/**
	 * Get position level normalized counts across a region
	 * @param region The region
	 * @return Normalized position level count by position
	 */
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region);
	
	/**
	 * Get a string identifier of the normalization
	 * @return The normalization name
	 */
	public String getNormalizationName();

}


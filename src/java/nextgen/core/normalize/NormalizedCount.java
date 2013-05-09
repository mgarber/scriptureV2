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
	 * Get raw count for the region
	 * @param region The region
	 * @return The raw count over the region
	 */
	public double getRawCount(Annotation region);
	
	/**
	 * Get position level raw counts across a region
	 * @param region The region
	 * @return Position level raw count by position
	 */
	public Map<Integer, Double> getRawCountsByPosition(Annotation region);
	
	
	
	
}


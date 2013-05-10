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

}


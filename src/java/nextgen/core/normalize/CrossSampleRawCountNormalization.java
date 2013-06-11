/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;
import java.util.TreeMap;

import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;

/**
 * @author prussell
 *
 */
public class CrossSampleRawCountNormalization implements NormalizedCount {
	
	private RawCounts sampleRawCounts;
	private RawCounts otherSampleRawCounts;
	
	/**
	 * @param alignmentData Alignment data to normalize
	 * @param otherAlignmentData Alignment data to normalize against
	 */
	public CrossSampleRawCountNormalization(AlignmentModel alignmentData, AlignmentModel otherAlignmentData) {
		sampleRawCounts = new RawCounts(alignmentData);
		otherSampleRawCounts = new RawCounts(otherAlignmentData);
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		return sampleRawCounts.getNormalizedCount(region) / otherSampleRawCounts.getNormalizedCount(region);
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		Map<Integer, Double> sampleCountsByPos = sampleRawCounts.getNormalizedCountsByPosition(region);
		Map<Integer, Double> otherSampleCountsByPos = otherSampleRawCounts.getNormalizedCountsByPosition(region);
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		for(Integer i : sampleCountsByPos.keySet()) {
			rtrn.put(i, Double.valueOf(sampleCountsByPos.get(i).doubleValue() / otherSampleCountsByPos.get(i).doubleValue()));
		}
		return rtrn;
	}

}

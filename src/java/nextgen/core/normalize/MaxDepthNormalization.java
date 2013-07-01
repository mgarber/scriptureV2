/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;

/**
 * @author prussell
 * TODO document here
 */
public class MaxDepthNormalization implements NormalizedCount {

	private RawCounts rawCounts;
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(MaxDepthNormalization.class.getName());

	/**
	 * @param alignmentData Alignment data
	 * @param fullyContained Only count fully contained alignments in annotations
	 */
	public MaxDepthNormalization(AlignmentModel alignmentData, boolean fullyContained) {
		rawCounts = new RawCounts(alignmentData, fullyContained);
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		Map<Integer, Double> rawPositionCounts = rawCounts.getNormalizedCountsByPosition(region);
		double rtrn = 0;
		for(Double d : rawPositionCounts.values()) {
			double c = d.doubleValue();
			if(c > rtrn) {
				rtrn = c;
			}
		}
		return rtrn;
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		double maxCount = getNormalizedCount(region);
		Map<Integer, Double> rawPositionCounts = rawCounts.getNormalizedCountsByPosition(region);
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		for(Integer i : rawPositionCounts.keySet()) {
			rtrn.put(i, Double.valueOf(rawPositionCounts.get(i).doubleValue() / maxCount));
		}
		return rtrn;
	}

	@Override
	public String getNormalizationName() {
		return NormalizedCount.MAX_DEPTH_NAME;
	}

}

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
 * TODO document here
 */
public class FragmentSizeNormalization implements NormalizedCount {
	
	private AlignmentModel alignmentData;
	private RawCounts rawCountData;

	/**
	 * @param data Alignment data
	 * @param fullyContained Count fully contained reads only
	 */
	public FragmentSizeNormalization(AlignmentModel data, boolean fullyContained) {
		alignmentData = data;
		rawCountData = new RawCounts(alignmentData, fullyContained);
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		int regionSize = region.getSize();
		double rawCount = rawCountData.getNormalizedCount(region);
		double medianFragmentSize;
		double normalizationFactor = 1;
		try {
			medianFragmentSize = alignmentData.getMedianReadSize(region, alignmentData.getCoordinateSpace(), regionSize, 200);
			normalizationFactor = regionSize / medianFragmentSize;
		} catch (IllegalArgumentException e) {
			// No overlapping fragments
		}
		return rawCount / normalizationFactor;
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		double normalizedCount = getNormalizedCount(region);
		Map<Integer, Double> rawPositionCounts = rawCountData.getNormalizedCountsByPosition(region);
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		for(Integer i : rawPositionCounts.keySet()) {
			rtrn.put(i, Double.valueOf(rawPositionCounts.get(i).doubleValue() / normalizedCount));
		}
		return rtrn;
	}

	@Override
	public String getNormalizationName() {
		return NormalizedCount.FRAGMENT_SIZE_NORMALIZATION_NAME;
	}

}

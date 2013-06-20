/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;
import java.util.TreeMap;

import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.RatioScore;

/**
 * @author prussell
 *
 */
public class CrossSampleRawCountNormalization implements NormalizedCount {
	
	private RawCounts numeratorRawCounts;
	private RawCounts denominatorRawCounts;
	private AlignmentModel numeratorData;
	private AlignmentModel denominatorData;
	
	/**
	 * @param numeratorAlignmentData Alignment data to normalize
	 * @param denominatorAlignmentData Alignment data to normalize against
	 * @param fullyContained Only count fully contained alignments in annotations
	 */
	public CrossSampleRawCountNormalization(AlignmentModel numeratorAlignmentData, AlignmentModel denominatorAlignmentData, boolean fullyContained) {
		numeratorData = numeratorAlignmentData;
		denominatorData = denominatorAlignmentData;
		numeratorRawCounts = new RawCounts(numeratorData, fullyContained);
		denominatorRawCounts = new RawCounts(denominatorData, fullyContained);
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		return numeratorRawCounts.getNormalizedCount(region) / denominatorRawCounts.getNormalizedCount(region);
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		Map<Integer, Double> sampleCountsByPos = numeratorRawCounts.getNormalizedCountsByPosition(region);
		Map<Integer, Double> otherSampleCountsByPos = denominatorRawCounts.getNormalizedCountsByPosition(region);
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		for(Integer i : sampleCountsByPos.keySet()) {
			rtrn.put(i, Double.valueOf(sampleCountsByPos.get(i).doubleValue() / otherSampleCountsByPos.get(i).doubleValue()));
		}
		return rtrn;
	}
	
	/**
	 * Get ratio score object for numerator count and denominator count over region
	 * @param region Region of interest
	 * @return Ratio score object
	 */
	public RatioScore getRatioScore(Annotation region) {
		RatioScore score = new RatioScore(region);
		score.setNumeratorTotal(numeratorData.getGlobalNumReads());
		score.setDenominatorTotal(denominatorData.getGlobalNumReads());
		score.setNumeratorCount(numeratorRawCounts.getNormalizedCount(region));
		score.setDenominatorCount(denominatorRawCounts.getNormalizedCount(region));
		score.setNumeratorRegionTotal(numeratorRawCounts.getNormalizedCount(region));
		score.setDenominatorRegionTotal(denominatorRawCounts.getNormalizedCount(region));
		return score;
	}
	
}

/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.RatioScore;

/**
 * @author prussell
 *
 */
public class CrossSampleRatio implements NormalizedCount {
	
	private NormalizedCount numeratorCountData;
	private NormalizedCount denominatorCountData;
	private AlignmentModel numeratorData;
	private AlignmentModel denominatorData;
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(CrossSampleRatio.class.getName());
	
	/**
	 * @param numeratorAlignmentData Alignment data to normalize
	 * @param denominatorAlignmentData Alignment data to normalize against
	 * @param numeratorCounts Normalized counts object for numerator data
	 * @param denominatorCounts Normalized counts object for denominator data
	 * @param fullyContained Only count fully contained alignments in annotations
	 */
	public CrossSampleRatio(AlignmentModel numeratorAlignmentData, AlignmentModel denominatorAlignmentData, NormalizedCount numeratorCounts, NormalizedCount denominatorCounts, boolean fullyContained) {
		numeratorData = numeratorAlignmentData;
		denominatorData = denominatorAlignmentData;
		numeratorCountData = numeratorCounts;
		denominatorCountData = denominatorCounts;
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		double n = numeratorCountData.getNormalizedCount(region);
		double d = denominatorCountData.getNormalizedCount(region);
		//logger.info("Region " + region.getName() + " numerator normalized count " + n + " denominator normalized count " + d + " ratio " + Double.valueOf(n/d).toString());
		return n / d;
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		Map<Integer, Double> sampleCountsByPos = numeratorCountData.getNormalizedCountsByPosition(region);
		Map<Integer, Double> otherSampleCountsByPos = denominatorCountData.getNormalizedCountsByPosition(region);
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
		score.setNumeratorTotal(numeratorData.getGlobalNumReadsReferenceSeqs()); 
		score.setDenominatorTotal(denominatorData.getGlobalNumReadsReferenceSeqs()); 
		score.setNumeratorCount(numeratorCountData.getNormalizedCount(region));
		score.setDenominatorCount(denominatorCountData.getNormalizedCount(region));
		score.setNumeratorRegionTotal(numeratorCountData.getNormalizedCount(region));
		score.setDenominatorRegionTotal(denominatorCountData.getNormalizedCount(region));
		return score;
	}

	@Override
	public String getNormalizationName() {
		return NormalizedCount.CROSS_SAMPLE_RATIO_NAME;
	}
	
}

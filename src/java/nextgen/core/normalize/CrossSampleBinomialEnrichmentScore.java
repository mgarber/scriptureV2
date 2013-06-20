/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;

import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.RatioScore;

/**
 * @author prussell
 *
 */
public class CrossSampleBinomialEnrichmentScore implements NormalizedCount {
	
	private CrossSampleRawCountNormalization comparativeCounts;
		
	/**
	 * @param numeratorAlignmentData Alignment data to normalize
	 * @param denominatorAlignmentData Alignment data to normalize against
	 * @param fullyContained Only count fully contained alignments in annotations
	 */
	public CrossSampleBinomialEnrichmentScore(AlignmentModel numeratorAlignmentData, AlignmentModel denominatorAlignmentData, boolean fullyContained) {
		comparativeCounts = new CrossSampleRawCountNormalization(numeratorAlignmentData, denominatorAlignmentData, fullyContained);
	}


	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		RatioScore score = comparativeCounts.getRatioScore(region);
		return score.getEnrichmentBinomialScore();
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		throw new UnsupportedOperationException("Method not implemented");
	}

}

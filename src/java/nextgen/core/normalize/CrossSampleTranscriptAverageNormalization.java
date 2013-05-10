/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;

/**
 * @author prussell
 *
 */
public class CrossSampleTranscriptAverageNormalization implements NormalizedCount {

	private TranscriptAverageNormalization withinTranscriptNormalizedData;
	private TranscriptAverageNormalization otherWithinTranscriptNormalizedData;
	private Map<Annotation, Double> normalizedCounts;
	private static Logger logger = Logger.getLogger(CrossSampleTranscriptAverageNormalization.class.getName());
	
	/**
	 * @param transcriptNormalizedData Dataset to normalize. Already normalized within transcript.
	 * @param otherTranscriptNormalizedData Dataset to normalize to. Already normalized within transcript.
	 */
	public CrossSampleTranscriptAverageNormalization(TranscriptAverageNormalization transcriptNormalizedData, TranscriptAverageNormalization otherTranscriptNormalizedData) {
		logger.info("");
		logger.info("Instantiating cross sample normalization object. Single sample normalized counts will be further normalized to other sample.");
		withinTranscriptNormalizedData = transcriptNormalizedData;
		otherWithinTranscriptNormalizedData = otherTranscriptNormalizedData;
		normalizedCounts = new TreeMap<Annotation, Double>();
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		if(!normalizedCounts.containsKey(region)) {
			computeNormalizedCount(region);
		}
		return normalizedCounts.get(region).doubleValue();
	}
	
	private void computeNormalizedCount(Annotation region) {
		double count1 = withinTranscriptNormalizedData.getNormalizedCount(region);
		double count2 = otherWithinTranscriptNormalizedData.getNormalizedCount(region);
		normalizedCounts.put(region, Double.valueOf(count1 / count2));
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		Map<Integer, Double> counts1 = withinTranscriptNormalizedData.getNormalizedCountsByPosition(region);
		Map<Integer, Double> counts2 = otherWithinTranscriptNormalizedData.getNormalizedCountsByPosition(region);
		counts1.keySet().retainAll(counts2.keySet());
		for(Integer i : counts1.keySet()) {
			double count1 = counts1.get(i).doubleValue();
			if(count1 == 0) continue;
			double count2 = counts2.get(i).doubleValue();
			if(count2 == 0) continue;
			double count = Math.log10(count1 / count2);
			rtrn.put(i, Double.valueOf(count));
		}
		return rtrn;
	}


}

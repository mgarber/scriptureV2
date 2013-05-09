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
public class TranscriptAverageNormalization implements NormalizedCount {
	
	private AlignmentModel data;
	private Map<Annotation, Double> regionAverages;
	private RawCounts rawCounts;
	
	/**
	 * @param alignmentData Alignment data
	 */
	public TranscriptAverageNormalization(AlignmentModel alignmentData) {
		regionAverages = new TreeMap<Annotation, Double>();
		data = alignmentData;
		rawCounts = new RawCounts(data);
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		throw new UnsupportedOperationException("");
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		double regionAverage = getRegionAverage(region);
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		if(regionAverage == 0) {
			return rtrn;
		}
		Map<Integer, Double> rawCount = getRawCountsByPosition(region);
		for(Integer i : rawCount.keySet()) {
			rtrn.put(i, Double.valueOf(rawCount.get(i).doubleValue() / regionAverage));
		}
		return rtrn;
	}
	
	/**
	 * Get the average position level read coverage over the region
	 * @param region The region
	 * @return The average coverage
	 */
	public double getRegionAverage(Annotation region) {
		if(regionAverages.containsKey(region)) {
			return regionAverages.get(region).doubleValue();
		}
		computeRegionAverage(region);
		return regionAverages.get(region).doubleValue();
	}

	private void computeRegionAverage(Annotation region) {
		Map<Integer, Double> rawCount = getRawCountsByPosition(region);
		int numCounts = rawCount.size();
		double sum = 0;
		for(Double d : rawCount.values()) {
			sum += d.doubleValue();
		}
		regionAverages.put(region, Double.valueOf(sum / numCounts));
	}

	@Override
	public double getRawCount(Annotation region) {
		return rawCounts.getRawCount(region);
	}

	@Override
	public Map<Integer, Double> getRawCountsByPosition(Annotation region) {
		return rawCounts.getNormalizedCountsByPosition(region);
	}
	
}

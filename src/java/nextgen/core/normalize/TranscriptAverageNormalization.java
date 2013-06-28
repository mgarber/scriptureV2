/**
 * 
 */
package nextgen.core.normalize;

import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;


import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 *
 */
public class TranscriptAverageNormalization implements NormalizedCount {
	
	private AlignmentModel data;
	private Map<Annotation, Double> regionAverages;
	private RawCounts rawCounts;
	private Map<Annotation, Annotation> childToParent;
	private static Logger logger = Logger.getLogger(TranscriptAverageNormalization.class.getName());
	
	/**
	 * @param alignmentData Alignment data
	 */
	public TranscriptAverageNormalization(AlignmentModel alignmentData) {
		this(alignmentData, null);
	}
	
	/**
	 * @param alignmentData Alignment data
	 * @param childAnnotationBedFile Bed file of child annotations
	 * @param parentAnnotationBedFile Bed file of parent annotations
	 * @throws IOException
	 */
	public TranscriptAverageNormalization(AlignmentModel alignmentData, String childAnnotationBedFile, String parentAnnotationBedFile) throws IOException {
		this(alignmentData, AnnotationUtils.mapChildToLargestParent(AnnotationUtils.loadAnnotations(childAnnotationBedFile), AnnotationUtils.loadAnnotations(parentAnnotationBedFile)));
	}
	
	/**
	 * @param alignmentData Alignment data
	 * @param childToParentAnnotation Map associating each annotation with a single parent annotation
	 */
	public TranscriptAverageNormalization(AlignmentModel alignmentData, Map<Annotation, Annotation> childToParentAnnotation) {
		logger.info("");
		logger.info("Instantiating normalization object. Each count will be normalized to the average over the transcript.");
		regionAverages = new TreeMap<Annotation, Double>();
		data = alignmentData;
		rawCounts = new RawCounts(data, false);
		childToParent = childToParentAnnotation;
	}
	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		if(childToParent == null) {
			throw new IllegalStateException("Must instantiate with a map of child to parent annotation.");
		}
		if(!childToParent.containsKey(region)) {
			throw new IllegalArgumentException("Child " + region.getName() + " does not have parent.");
		}
		if(childToParent.get(region) == null) {
			throw new IllegalArgumentException("Parent of child " + region.getName() + " is null.");
		}		
		return getNormalizedRegionAverage(region, childToParent.get(region));
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		double regionAverage = getRawRegionAverage(region);
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		if(regionAverage == 0) {
			return rtrn;
		}
		Map<Integer, Double> rawCount = rawCounts.getNormalizedCountsByPosition(region);
		for(Integer i : rawCount.keySet()) {
			rtrn.put(i, Double.valueOf(rawCount.get(i).doubleValue() / regionAverage));
		}
		return rtrn;
	}
	
	/**
	 * Get the average position level read coverage over the region, normalized to the average coverage over other region
	 * @param region The region
	 * @param other The other region to normalize to
	 * @return The ratio of average read coverage over the region to average read coverage over the other region
	 */
	public double getNormalizedRegionAverage(Annotation region, Annotation other) {
		return getRawRegionAverage(region) / getRawRegionAverage(other);
	}
	
	/**
	 * Get the average position level read coverage over the region
	 * @param region The region
	 * @return The average coverage
	 */
	public double getRawRegionAverage(Annotation region) {
		if(regionAverages.containsKey(region)) {
			return regionAverages.get(region).doubleValue();
		}
		computeRegionAverage(region);
		return regionAverages.get(region).doubleValue();
	}

	private void computeRegionAverage(Annotation region) {
		Map<Integer, Double> rawCount = rawCounts.getNormalizedCountsByPosition(region);
		int numCounts = rawCount.size();
		double sum = 0;
		for(Double d : rawCount.values()) {
			sum += d.doubleValue();
		}
		regionAverages.put(region, Double.valueOf(sum / numCounts));
	}

	@Override
	public String getNormalizationName() {
		return NormalizedCount.TRANSCRIPT_AVERAGE_NORMALIZATION_NAME;
	}

	
}

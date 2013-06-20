/**
 * 
 */
package nextgen.core.normalize;

import java.util.Map;
import java.util.TreeMap;

import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;

/**
 * @author prussell
 *
 */
public class RawCounts implements NormalizedCount {

	private AlignmentModel data;
	private boolean fullyContainedOnly;

	/**
	 * @param alignmentData Alignment data
	 * @param fullyContained Only count fully contained alignments in annotations
	 */
	public RawCounts(AlignmentModel alignmentData, boolean fullyContained) {
		data = alignmentData;
		fullyContainedOnly = fullyContained;
	}

	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		return data.getCount(region, fullyContainedOnly);
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		TreeMap<Integer,Double> rtrn = new TreeMap<Integer,Double>();
		WindowProcessor<CountScore> processor = new CountScore.Processor(data);
		WindowScoreIterator<CountScore> scoreIter = data.scan(region, 1, 0, processor);
		while(scoreIter.hasNext()) {
			CountScore score = scoreIter.next();
			int pos = score.getAnnotation().getStart();
			rtrn.put(Integer.valueOf(pos), Double.valueOf(score.getCount()));
		}
		return rtrn;
	}


}

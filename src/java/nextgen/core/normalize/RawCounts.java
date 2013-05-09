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

	/**
	 * @param alignmentData Alignment data
	 */
	public RawCounts(AlignmentModel alignmentData) {
		data = alignmentData;
	}

	
	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getNormalizedCount(Annotation region) {
		return getRawCount(region);
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getNormalizedCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getNormalizedCountsByPosition(Annotation region) {
		return getRawCountsByPosition(region);
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getRawCount(nextgen.core.annotation.Annotation)
	 */
	@Override
	public double getRawCount(Annotation region) {
		return data.getCount(region);
	}

	/* (non-Javadoc)
	 * @see nextgen.core.normalize.NormalizedCount#getRawCountsByPosition(nextgen.core.annotation.Annotation)
	 */
	@Override
	public Map<Integer, Double> getRawCountsByPosition(Annotation region) {
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

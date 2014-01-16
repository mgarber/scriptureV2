package nextgen.editing.crispr;

import org.apache.log4j.Logger;

/**
 * Guide RNA score based on Feng Zhang's algorithm
 * @author prussell
 */
public class GuideOffTargetScore implements GuideRNAScore {

	public static Logger logger = Logger.getLogger(GuideOffTargetScore.class.getName());

	@Override
	public double getScore(GuideRNA guideRNA) {
		// TODO
		throw new UnsupportedOperationException();
	}

}

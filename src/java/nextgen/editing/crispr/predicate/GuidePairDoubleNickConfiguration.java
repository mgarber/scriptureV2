package nextgen.editing.crispr.predicate;

import nextgen.editing.crispr.NickingGuideRNAPair;

import org.apache.log4j.Logger;

/**
 * Check whether a guide RNA pair is correctly laid out for double nick strategy
 * @author prussell
 */
public class GuidePairDoubleNickConfiguration implements GuideRNAPairPredicate {
	
	public static Logger logger = Logger.getLogger(GuidePairDoubleNickConfiguration.class.getName());
	private int maxInnerDistance;
	
	/**
	 * @param maxInnerDist Max inner distance between the two guide RNAs
	 */
	public GuidePairDoubleNickConfiguration(int maxInnerDist) {
		maxInnerDistance = maxInnerDist;
	}
	
	@Override
	public boolean evaluate(NickingGuideRNAPair guideRnaPair) {
		// The two guide RNAs cannot overlap
		if(!guideRnaPair.nonoverlapping()) {
			//logger.debug("INVALID_PAIR\tGuide RNAs overlap:\t" + guideRnaPair.toString());
			return false;
		}
		// The two guide RNAs must be facing outward from each other
		if(!guideRnaPair.facingOutward()) {
			//logger.debug("INVALID_PAIR\tGuide RNAs face inward:\t" + guideRnaPair.toString());
			return false;
		}
		// The distance between the guide RNAs' 5' ends must be at most 40
		int innerDist = guideRnaPair.getInnerDistance();
		if(innerDist < 0) {
			//logger.debug("INVALID_PAIR\tInner distance < 0:\t" + guideRnaPair.toString());
			return false;
		}
		if(innerDist > maxInnerDistance) {
			//logger.debug("INVALID_PAIR\tInner distance > " + maxInnerDistance + ":\t" + guideRnaPair.toString());
			return false;
		}
		// If all the criteria are met return true
		//logger.debug("VALID_PAIR\t" + guideRnaPair.toString());
		return true;
	}

	@Override
	public String getPredicateName() {
		return "valid_double_nick_configuration";
	}

	@Override
	public String getShortFailureMessage(NickingGuideRNAPair g) {
		return("not_valid_double_nick_configuration");
	}

}

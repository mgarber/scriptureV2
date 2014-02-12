package nextgen.editing.crispr.predicate;

import java.io.IOException;

import nextgen.editing.crispr.GuideRNA;
import nextgen.editing.crispr.score.GuideEfficacyScore;

public class GuideSufficientEfficacy implements GuideRNAPredicate {
	
	private GuideEfficacyScore score;
	private double maxScore;
	
	/**
	 * @param guideEfficacyScore Score object
	 * @param maxAllowableScore Max allowable guide efficacy score
	 */
	public GuideSufficientEfficacy(GuideEfficacyScore guideEfficacyScore, double maxAllowableScore) {
		score = guideEfficacyScore;
		maxScore = maxAllowableScore;
	}
	
	@Override
	public boolean evaluate(GuideRNA g) {
		try {
			return score.getScore(g) < maxScore;
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalStateException();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new IllegalStateException();
		}
	}

	@Override
	public String getPredicateName() {
		return "guide_efficacy";
	}

	@Override
	public String getShortFailureMessage(GuideRNA g) throws IOException, InterruptedException {
		return "guide_efficacy_score_" + score.getScore(g);
	}

	

}

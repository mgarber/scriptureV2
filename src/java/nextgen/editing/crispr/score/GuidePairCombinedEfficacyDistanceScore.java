package nextgen.editing.crispr.score;

import java.io.IOException;

import nextgen.editing.crispr.NickingGuideRNAPair;

/**
 * Combined score prioritizing low inner distance first, low (good) combined efficacy scores of both guides second
 * Lower score is better
 * @author prussell
 *
 */
public class GuidePairCombinedEfficacyDistanceScore implements GuideRNAPairScore {
	
	private GuideEfficacyScore efficacyScore;
	
	/**
	 * @param guideEfficacyScore A guide efficacy score object, preferably instantiated with the individual guide RNAs that make up the pairs that will be evaluated
	 */
	public GuidePairCombinedEfficacyDistanceScore(GuideEfficacyScore guideEfficacyScore) {
		efficacyScore = guideEfficacyScore;
	}
	
	@Override
	public double getScore(NickingGuideRNAPair guideRNApair) throws IOException, InterruptedException {
		
		int innerDist = guideRNApair.getInnerDistance();
		double sumEfficacyScores = efficacyScore.getScore(guideRNApair.getLeftGuideRNA()) + efficacyScore.getScore(guideRNApair.getRightGuideRNA());
		
		// Efficacy scores are between 0 and 1 so sum is between 0 and 2
		// Inner distance is more important than efficacy score
		// Lower is better for both
		
		return innerDist + 5 * sumEfficacyScores;
		
	}

}

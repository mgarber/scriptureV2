package nextgen.editing.crispr.score;

import java.io.IOException;

import nextgen.editing.crispr.NickingGuideRNAPair;

public class GuidePairInnerDistanceScore implements GuideRNAPairScore {

	@Override
	public double getScore(NickingGuideRNAPair guideRNApair) throws IOException,InterruptedException {
		return guideRNApair.getInnerDistance();
	}

}

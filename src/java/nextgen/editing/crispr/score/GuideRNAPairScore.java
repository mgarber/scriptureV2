package nextgen.editing.crispr.score;

import java.io.IOException;

import nextgen.editing.crispr.NickingGuideRNAPair;

public interface GuideRNAPairScore {

	/**
	 * Get the value of the score
	 * @param guideRNApair The guide RNA pair
	 * @return The score
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public double getScore(NickingGuideRNAPair guideRNApair) throws IOException, InterruptedException;
	
}

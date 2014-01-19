package nextgen.editing.crispr.score;

import java.io.IOException;

import nextgen.editing.crispr.GuideRNA;

/**
 * A score for a guide RNA object
 * @author prussell
 *
 */
public interface GuideRNAScore {
	
	/**
	 * Get the value of the score
	 * @param guideRNA The guide RNA
	 * @return The score
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public double getScore(GuideRNA guideRNA) throws IOException, InterruptedException;
	
	/**
	 * Get the name of the score
	 * @return Score name
	 */
	public String getScoreName();
	
}

package nextgen.editing.crispr;

import java.io.IOException;

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
	
}

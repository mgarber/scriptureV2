/**
 * 
 */
package nextgen.synbio;

import broad.core.sequence.Sequence;

/**
 * @author prussell
 *
 */
public class TypeIISRestrictionEnzyme {
	
	/**
	 * Recognition sequence on top strand
	 */
	private Sequence topStrandRecognitionSite;
	
	/**
	 * Cleavage position on top strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 */
	private int topStrandCleavageSite;
	
	/**
	 * Cleavage position on bottom strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 */
	private int bottomStrandCleavageSite;
	
	/**
	 * Enzyme name
	 */
	private RestrictionEnzymeFactory.RestrictionEnzymeName name;
	
	/**
	 * @param enzymeName Enzyme name
	 * @param topStrandRecognitionSequence Recognition sequence on top strand
	 * @param topStrandCleavagePosition Cleavage position on top strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 * @param bottomStrandCleavagePosition Cleavage position on bottom strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 */
	public TypeIISRestrictionEnzyme(RestrictionEnzymeFactory.RestrictionEnzymeName enzymeName, String topStrandRecognitionSequence, int topStrandCleavagePosition, int bottomStrandCleavagePosition) {
		if(bottomStrandCleavagePosition > 0) {
			throw new IllegalArgumentException("Enzymes that cleave within or 3' of the recognition site on the bottom strand not supported (you provided bottom strand cleavage postion = " + bottomStrandCleavagePosition + ")");
		}
		Sequence seq = new Sequence("top_strand_recognition_sequence");
		seq.setSequenceBases(topStrandRecognitionSequence.toUpperCase());
		topStrandRecognitionSite = seq;
		topStrandCleavageSite = topStrandCleavagePosition;
		bottomStrandCleavageSite = bottomStrandCleavagePosition;
		name = enzymeName;
	}

	/**
	 * Get enzyme name
	 * @return Enzyme name
	 */
	public String getName() {
		return name.toString();
	}
	
	/**
	 * Get top strand recoginition sequence
	 * @return Top strand recognition sequence
	 */
	public String getTopStrandRecognitionSite() {
		return topStrandRecognitionSite.getSequenceBases();
	}
	
	/**
	 * Get bottom strand recoginition sequence
	 * @return Bottom strand recognition sequence
	 */
	public String getBottomStrandRecognitionSite() {
		return Sequence.reverseSequence(topStrandRecognitionSite).getSequenceBases();
	}
	
	/**
	 * Get top strand cleavage site
	 * @return Cleavage position on top strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' is 1.)
	 */
	public int getTopStrandCleavageSite() {
		return topStrandCleavageSite;
	}
	
	/**
	 * Get bottom strand cleavage site
	 * @return Cleavage position on bottom strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' is 1.)
	 */
	public int getBottomStrandCleavageSite() {
		return bottomStrandCleavageSite;
	}
	
	@Override
	public String toString() {
		return getName();
	}
	

}

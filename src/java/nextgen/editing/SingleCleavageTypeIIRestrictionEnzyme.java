package nextgen.editing;

import broad.core.sequence.Sequence;

/**
 * Restriction enzyme with contiguous recognition site and single cleavage position
 * @author prussell
 */
public class SingleCleavageTypeIIRestrictionEnzyme implements RestrictionEnzyme {
	
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
	public SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeFactory.RestrictionEnzymeName enzymeName, String topStrandRecognitionSequence, int topStrandCleavagePosition, int bottomStrandCleavagePosition) {
		Sequence seq = new Sequence("top_strand_recognition_sequence");
		seq.setSequenceBases(topStrandRecognitionSequence.toUpperCase());
		topStrandRecognitionSite = seq;
		topStrandCleavageSite = topStrandCleavagePosition;
		bottomStrandCleavageSite = bottomStrandCleavagePosition;
		name = enzymeName;
	}

	public String getName() {
		return name.toString();
	}
	
	public String getTopStrandRecognitionSequence() {
		return topStrandRecognitionSite.getSequenceBases();
	}
	
	public String getBottomStrandRecognitionSequence() {
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

package nextgen.editing;

import java.util.ArrayList;
import java.util.Collection;

import broad.core.sequence.Sequence;

/**
 * Restriction enzyme with contiguous recognition site and single cleavage position
 * @author prussell
 */
public class SingleCleavageTypeIIRestrictionEnzyme implements RestrictionEnzyme {
	
	/**
	 * Recognition sequence on top strand
	 */
	private Collection<Sequence> topStrandRecognitionSite;
	
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
		this(enzymeName, asCollection(topStrandRecognitionSequence), topStrandCleavagePosition, bottomStrandCleavagePosition);
	}
	
	protected static Collection<String> asCollection(String s) {
		Collection<String> rtrn = new ArrayList<String>();
		rtrn.add(s);
		return rtrn;
	}
	
	/**
	 * @param enzymeName Enzyme name
	 * @param topStrandRecognitionSequence Recognition sequence on top strand
	 * @param topStrandCleavagePosition Cleavage position on top strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 * @param bottomStrandCleavagePosition Cleavage position on bottom strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 */
	public SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeFactory.RestrictionEnzymeName enzymeName, Collection<String> topStrandRecognitionSequence, int topStrandCleavagePosition, int bottomStrandCleavagePosition) {
		
		topStrandRecognitionSite = new ArrayList<Sequence>();
		for(String s : topStrandRecognitionSequence) {
			Sequence seq = new Sequence("top_strand_recognition_sequence");
			seq.setSequenceBases(s.toUpperCase());
			topStrandRecognitionSite.add(seq);
		}
		topStrandCleavageSite = topStrandCleavagePosition;
		bottomStrandCleavageSite = bottomStrandCleavagePosition;
		name = enzymeName;
	}

	public String getName() {
		return name.toString();
	}
	
	public Collection<String> getTopStrandRecognitionSequence() {
		Collection<String> rtrn = new ArrayList<String>();
		for(Sequence seq : topStrandRecognitionSite) {
			rtrn.add(seq.getSequenceBases());
		}
		return rtrn;
	}
	
	public Collection<String> getBottomStrandRecognitionSequence() {
		Collection<String> top = getTopStrandRecognitionSequence();
		Collection<String> rtrn = new ArrayList<String>();
		for(String s : top) {
			rtrn.add(Sequence.reverseSequence(s));
		}
		return rtrn;
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

	@Override
	public boolean hasTopStrandRecognitionSequence(String seq) {
		for(Sequence s : topStrandRecognitionSite) {
			if(s.getSequenceBases().equalsIgnoreCase(seq)) {
				return true;
			}
		}
		return false;
	}
	

	
}

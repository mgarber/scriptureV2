package nextgen.editing;

import java.util.Collection;

/**
 * A representation of a restriction enzyme
 * The interface contains no information about how the enyzme cuts
 * Implementations should capture the cutting behavior of the enzyme in an appropriate way
 * @author prussell
 */
public interface RestrictionEnzyme {

	/**
	 * Get top strand recoginition sequence
	 * @return Top strand recognition sequence
	 */
	public Collection<String> getTopStrandRecognitionSequence();
	
	/**
	 * Get enzyme name
	 * @return Enzyme name
	 */
	public String getName();
	
	
	/**
	 * Whether the sequence is a top strand recognition sequence of this enzyme
	 * @param seq The sequence
	 * @return True iff this is a recognition sequence of the enzyme
	 */
	public boolean hasTopStrandRecognitionSequence(String seq);
	
	/**
	 * Whether the sequence contains a recognition sequence of this enzyme
	 * @param seq The sequence
	 * @return True iff the sequence contains a recognition sequence of this enzyme
	 */
	public boolean sequenceContainsTopStrandRecognitionSequence(String seq);
	
	
}

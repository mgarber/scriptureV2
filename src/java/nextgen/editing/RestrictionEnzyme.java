package nextgen.editing;

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
	public String getTopStrandRecognitionSequence();
	
	/**
	 * Get enzyme name
	 * @return Enzyme name
	 */
	public String getName();
	
}

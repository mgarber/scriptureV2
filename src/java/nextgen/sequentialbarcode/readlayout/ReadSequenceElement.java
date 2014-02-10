package nextgen.sequentialbarcode.readlayout;

/**
 * An element that can be found in a read sequence
 * @author prussell
 *
 */
public interface ReadSequenceElement {
	
	/**
	 * Get element length
	 * @return Element length
	 */
	public int getLength();
	
	/**
	 * Check whether this element matches a string
	 * @param s The string
	 * @return True if this element matches the string
	 */
	public boolean matchesFullString(String s);
	
	/**
	 * Check whether this element matches a substring of a string
	 * @param s The string
	 * @param startOnString Start position of substring within string
	 * @return True if this element matches the substring starting at the specified position
	 */
	public boolean matchesSubstringOf(String s, int startOnString);
	
	/**
	 * Whether instances of the element can appear an arbitrary number of times in tandem within a read
	 * @return True iff element is repeatable
	 */
	public boolean isRepeatable();
	
	/**
	 * For repeatable read elements, a string whose presence signals the repeat is over
	 * @return A sequence that comes after the repeatable element, or null if not repeatable
	 */
	public String getStopSignalForRepeatable();
	
	/**
	 * Get the "intended" sequence of the match, i.e., the real barcode, not a version with mismatches found in the read
	 * @param s Read subsequence to search
	 * @return The real intended sequence or null if doesn't match
	 */
	public String matchedElementSequence(String s);
	
	/**
	 * Get the read element matching the sequence
	 * @param s Read subsequence to search
	 * @return The element that matches the sequence or null if it doesn't match
	 */
	public ReadSequenceElement matchedElement(String s);
	
	/**
	 * Get the general name of the element
	 * @return General element name
	 */
	public String elementName();
	
	/**
	 * Get the specific ID of this instance
	 * @return Element ID
	 */
	public String getId();
	
}

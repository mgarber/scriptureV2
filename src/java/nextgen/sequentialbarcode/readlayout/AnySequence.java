package nextgen.sequentialbarcode.readlayout;

import org.apache.log4j.Logger;

/**
 * The notion of an endogenous sequence of a particular length
 * Does not store an actual nucleotide sequence
 * @author prussell
 *
 */
public class AnySequence implements ReadSequenceElement {
	
	private int length;
	public static Logger logger = Logger.getLogger(AnySequence.class.getName());

	/**
	 * @param len The sequence length
	 */
	public AnySequence(int len) {
		length = len;
	}
	
	@Override
	public int getLength() {
		return length;
	}

	/**
	 * Only check that the string has this length
	 */
	@Override
	public boolean matchesFullString(String s) {
		return s.length() == length;
	}

	@Override
	public String elementName() {
		return "endogenous_sequence";
	}

	@Override
	public String matchedElementSequence(String s) {
		return s;
	}

	@Override
	public boolean matchesSubstringOf(String s, int startOnString) {
		return startOnString + length <= s.length();
	}

	@Override
	public String getId() {
		return null;
	}

	@Override
	public ReadSequenceElement matchedElement(String s) {
		if(!matchesFullString(s)) {
			return null;
		}
		return this;
	}

	@Override
	public boolean isRepeatable() {
		return false;
	}

	@Override
	public String getStopSignalForRepeatable() {
		return null;
	}

}

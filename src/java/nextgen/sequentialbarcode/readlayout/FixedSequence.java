package nextgen.sequentialbarcode.readlayout;

import org.apache.log4j.Logger;

import nextgen.core.utils.AlignmentUtils;

/**
 * A fixed sequence that should be the same in all reads
 * @author prussell
 *
 */
public class FixedSequence implements ReadSequenceElement {
	
	private String seq;
	private int maxNumMismatches;
	private String name;
	public static Logger logger = Logger.getLogger(FixedSequence.class.getName());

	/**
	 * @param fixedSeqName Name of fixed sequence
	 * @param sequence Nucleotide sequence
	 * @param maxMismatches Max number of mismatches when identifying this sequence in a read
	 */
	public FixedSequence(String fixedSeqName, String sequence, int maxMismatches) {
		seq = sequence;
		name = fixedSeqName;
		maxNumMismatches = maxMismatches;
	}
	
	@Override
	public int getLength() {
		return seq.length();
	}

	/**
	 * Get the sequence
	 * @return The nucleotide sequence
	 */
	public String getSequence() {
		return seq;
	}

	@Override
	public boolean matchesFullString(String s) {
		if(getLength() != s.length()) {
			return false;
		}
		return AlignmentUtils.hammingDistanceAtMost(s, seq, maxNumMismatches, true);
	}

	@Override
	public String elementName() {
		return name;
	}

	@Override
	public String matchedElementSequence(String s) {
		if(!matchesFullString(s)) {
			return null;
		}
		return seq;
	}

	@Override
	public boolean matchesSubstringOf(String s, int startOnString) {
		return matchesFullString(s.substring(startOnString, startOnString + getLength()));
	}

	@Override
	public String getId() {
		return elementName();
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

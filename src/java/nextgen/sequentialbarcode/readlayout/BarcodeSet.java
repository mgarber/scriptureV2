package nextgen.sequentialbarcode.readlayout;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

/**
 * A set of barcodes that are used together in the lab
 * @author prussell
 *
 */
public class BarcodeSet implements ReadSequenceElement {
	
	private Collection<Barcode> barcodes;
	public static Logger logger = Logger.getLogger(BarcodeSet.class.getName());
	private Map<String, Barcode> seqToMatchedElement;
	private Collection<String> noMatch;
	private int length;
	private String id;
	private boolean repeatable;
	private String nextStringForRepeatable;
	
	/**
	 * @param setId Barcode set ID
	 * @param barcodeSet The barcodes
	 * @param maxMismatches Max number of mismatches when identifying each barcode in reads
	 */
	public BarcodeSet(String setId, Collection<Barcode> barcodeSet, int maxMismatches) {
		this(setId, barcodeSet, maxMismatches, false, null);
	}
	
	/**
	 * @param setId Barcode set ID
	 * @param barcodeSet The barcodes
	 * @param maxMismatches Max number of mismatches when identifying each barcode in reads
	 * @param isRepeatable Whether to look for multiple matches in sequence
	 * @param stopSignalForRepeatable String whose presence in a read signals the end of the region that is expected to contain these barcodes
	 */
	public BarcodeSet(String setId, Collection<Barcode> barcodeSet, int maxMismatches, boolean isRepeatable, String stopSignalForRepeatable) {
		id = setId;
		repeatable = isRepeatable;
		nextStringForRepeatable = stopSignalForRepeatable;
		seqToMatchedElement = new TreeMap<String, Barcode>();
		noMatch = new TreeSet<String>();
		int len = barcodeSet.iterator().next().getLength();
		barcodes = new TreeSet<Barcode>();
		for(Barcode b : barcodeSet) {
			if(b.getLength() != len) {
				throw new IllegalArgumentException("All barcode sequences must have the same length");
			}
			barcodes.add(b);
		}
		length = len;
	}
	
	@Override
	public int getLength() {
		return length;
	}

	@Override
	public boolean matchesFullString(String s) {
		if(s.length() != length) {
			return false;
		}
		if(seqToMatchedElement.containsKey(s)) {
			return true;
		}
		if(noMatch.contains(s)) {
			return false;
		}
		for(Barcode barcode : barcodes) {
			if(barcode.matchesFullString(s)) {
				seqToMatchedElement.put(s, barcode);
				return true;
			}
		}
		noMatch.add(s);
		return false;
	}


	@Override
	public String elementName() {
		return "barcode_set";
	}

	@Override
	public String matchedElementSequence(String s) {
		return matchedElement(s).matchedElementSequence(s);
	}
	
	@Override
	public boolean matchesSubstringOf(String s, int startOnString) {
		return matchesFullString(s.substring(startOnString, startOnString + getLength()));
	}

	@Override
	public String getId() {
		return id;
	}

	@Override
	public ReadSequenceElement matchedElement(String s) {
		if(seqToMatchedElement.containsKey(s)) {
			return seqToMatchedElement.get(s);
		}
		if(!matchesFullString(s)) {
			return null;
		}
		return seqToMatchedElement.get(s);
	}

	@Override
	public boolean isRepeatable() {
		return repeatable;
	}

	@Override
	public String getStopSignalForRepeatable() {
		return nextStringForRepeatable;
	}


}

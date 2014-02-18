package nextgen.sequentialbarcode.readlayout;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.parser.StringParser;

import nextgen.core.utils.AlignmentUtils;
import nextgen.core.utils.FileUtil;

/**
 * A single barcode
 * @author prussell
 *
 */
public class Barcode implements ReadSequenceElement, Comparable<Barcode> {
	
	private String sequence;
	private String id;
	private int maxNumMismatches;
	public static Logger logger = Logger.getLogger(Barcode.class.getName());
	private Collection<String> matchedStrings;
	private boolean repeatable;
	private String nextStringForRepeatable;
	
	/**
	 * @param seq The barcode sequence
	 */
	public Barcode(String seq) {
		this(seq, -1, false, null);
	}
	
	/**
	 * @param seq The barcode sequence
	 * @param maxMismatches Max allowable number of mismatches when identifying this barcode in reads
	 */
	public Barcode(String seq, int maxMismatches) {
		this(seq, maxMismatches, false, null);
	}
	
	/**
	 * @param seq The barcode sequence
	 * @param maxMismatches Max allowable number of mismatches when identifying this barcode in reads
	 * @param isRepeatable Whether this barcode can appear multiple times in tandem
	 * @param stopSignalForRepeatable String whose presence in read signals the end of the region where this barcode is expected to be found
	 */
	public Barcode(String seq, int maxMismatches, boolean isRepeatable, String stopSignalForRepeatable) {
		this(seq, null, maxMismatches, isRepeatable, stopSignalForRepeatable);
	}

	/**
	 * @param seq The barcode sequence
	 * @param barcodeId Unique ID for this barcode
	 * @param maxMismatches Max allowable number of mismatches when identifying this barcode in reads
	 */
	public Barcode(String seq, String barcodeId, int maxMismatches) {
		this(seq, barcodeId, maxMismatches, false, null);
	}
	
	/**
	 * @param seq The barcode sequence
	 * @param barcodeId Unique ID for this barcode
	 */
	public Barcode(String seq, String barcodeId) {
		this(seq, barcodeId, -1, false, null);
	}
	
	/**
	 * @param seq The barcode sequence
	 * @param barcodeId Unique ID for this barcode
	 * @param maxMismatches Max allowable number of mismatches when identifying this barcode in reads
	 * @param isRepeatable Whether this barcode can appear multiple times in tandem
	 * @param stopSignalForRepeatable String whose presence in read signals the end of the region where this barcode is expected to be found
	 */
	public Barcode(String seq, String barcodeId, int maxMismatches, boolean isRepeatable, String stopSignalForRepeatable) {;
		sequence = seq;
		id = barcodeId;
		maxNumMismatches = maxMismatches;
		matchedStrings = new TreeSet<String>();
		repeatable = isRepeatable;
		nextStringForRepeatable = stopSignalForRepeatable;
	}
	
	/**
	 * Create a set of barcodes with IDs from a table file
	 * Line format: barcode_id	barcode_sequence
	 * @param tableFile Table file
	 * @param maxMismatches Max allowable mismatches when matching barcodes
	 * @return Collection of barcodes
	 * @throws IOException
	 */
	public static Collection<Barcode> createBarcodesFromTable(String tableFile, int maxMismatches) throws IOException {
		FileReader r = new FileReader(tableFile);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		Collection<Barcode> rtrn = new TreeSet<Barcode>();
		while(b.ready()) {
			s.parse(b.readLine());
			if(s.getFieldCount() == 0) {
				continue;
			}
			if(s.getFieldCount() != 2) {
				r.close();
				b.close();
				throw new IllegalArgumentException("Format: barcode_id  barcode_sequence");
			}
			rtrn.add(new Barcode(s.asString(1), s.asString(0), maxMismatches));
		}
		r.close();
		b.close();
		return rtrn;
	}

	/**
	 * Create a set of barcodes without IDs from a list file
	 * @param listFile File with one barcode sequence per line
	 * @param maxMismatches Max allowable mismatches when matching barcodes
	 * @return Collection of barcodes
	 * @throws IOException
	 */
	public static Collection<Barcode> createBarcodesFromList(String listFile, int maxMismatches) throws IOException {
		return createBarcodes(FileUtil.fileLinesAsList(listFile), maxMismatches);
	}
	
	/**
	 * Create barcode objects from a collection of barcode sequences
	 * @param barcodeSeqs Barcode sequences
	 * @param maxMismatches Max allowable mismatches in each barcode when matching to read sequences
	 * @return Collection of barcode objects
	 */
	public static Collection<Barcode> createBarcodes(Collection<String> barcodeSeqs, int maxMismatches) {
		Collection<Barcode> rtrn = new TreeSet<Barcode>();
		for(String b : barcodeSeqs) {
			rtrn.add(new Barcode(b, maxMismatches));
		}
		return rtrn;
	}
	
	@Override
	public int compareTo(Barcode o) {
		int c = sequence.compareTo(o.getSequence());
		if(c != 0 || id == null) {
			return c;
		}
		return id.compareTo(o.getId());
	}

	@Override
	public int getLength() {
		return sequence.length();
	}

	@Override
	public String getId() {
		return id;
	}
	
	/**
	 * Get the barcode sequence
	 * @return The barcode sequence
	 */
	public String getSequence() {
		return sequence;
	}

	/**
	 * Get the mismatch tolerance for this barcode
	 * @return Max mismatches
	 */
	public int getMismatchTolerance() {
		return maxNumMismatches;
	}
	
	@Override
	public boolean matchesFullString(String s) {
		if(matchedStrings.contains(s)) {
			return true;
		}
		if(getLength() != s.length()) {
			return false;
		}
		boolean rtrn = AlignmentUtils.hammingDistanceAtMost(s, sequence, maxNumMismatches, true);
		if(rtrn) {
			matchedStrings.add(s);
		}
		return rtrn;
	}

	@Override
	public String elementName() {
		return "barcode";
	}

	@Override
	public String matchedElementSequence(String s) {
		if(!matchesFullString(s)) {
			return null;
		}
		return sequence;
	}
	
	
	@Override
	public boolean matchesSubstringOf(String s, int startOnString) {
		return matchesFullString(s.substring(startOnString, startOnString + getLength()));
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
		return repeatable;
	}

	@Override
	public String getStopSignalForRepeatable() {
		return nextStringForRepeatable;
	}


}

package nextgen.sequentialbarcode;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;

import com.sleepycat.persist.model.Persistent;

import broad.core.parser.StringParser;

import net.sf.samtools.SAMRecord;
import nextgen.sequentialbarcode.readlayout.Barcode;

/**
 * A possible set of barcodes identified with a fragment
 * @author prussell
 *
 */
@Persistent
public class BarcodeSequence implements Comparable<BarcodeSequence> {
	
	private List<Barcode> barcodes;
	public static Logger logger = Logger.getLogger(BarcodeSequence.class.getName());
	private String samAttributeString;

	/**
	 * Instantiate with no barcodes
	 */
	public BarcodeSequence() {
		barcodes = new ArrayList<Barcode>();
	}
	
	/**
	 * Instantiate with an initial set of barcodes
	 * @param barcodeList Ordered list of barcodes
	 */
	public BarcodeSequence(List<Barcode> barcodeList) {
		this();
		appendBarcodes(barcodeList);
		resetSamAttributeString();
	}
	
	/**
	 * Add a barcode to the end of the sequence
	 * @param barcode The barcode to add
	 */
	public void appendBarcode(Barcode barcode) {
		barcodes.add(barcode);
		resetSamAttributeString();
	}
	
	/**
	 * Add a list of barcodes to the end of the sequence
	 * @param newBarcodes Ordered list of barcodes to add
	 */
	public void appendBarcodes(List<Barcode> newBarcodes) {
		for(Barcode b : newBarcodes) {
			appendBarcode(b);
		}
		resetSamAttributeString();
	}
	
	/**
	 * Refresh the SAM attribute string representing this sequence of barcodes
	 */
	private void resetSamAttributeString() {
		samAttributeString = getSamAttributeString();
	}
	
	/**
	 * Create by reading barcodes from a sam attribute
	 * @param samRecord SAM record
	 * @return The barcode collection recorded in the sam record
	 */
	public static BarcodeSequence fromSamRecord(SAMRecord samRecord) {
		String attribute = samRecord.getStringAttribute(BarcodedBamWriter.BARCODES_SAM_TAG);
		return fromSamAttributeString(attribute);
	}
	
	/**
	 * Create FragmentBarcodes object from a string representation as produced by toSamAttributeString()
	 * @param s The string representation
	 * @return The barcode collection represented by the string
	 */
	public static BarcodeSequence fromSamAttributeString(String s) {
		return fromSamAttributeString(s, -1);
	}
	
	/**
	 * Create FragmentBarcodes object from a string representation as produced by toSamAttributeString()
	 * @param s The string representation
	 * @param maxMismatches Max mismatches for each barcode
	 * @return The barcode collection represented by the string
	 */
	public static BarcodeSequence fromSamAttributeString(String s, int maxMismatches) {
		List<Barcode> barcodes = new ArrayList<Barcode>();
		StringParser s1 = new StringParser();
		StringParser s2 = new StringParser();
		s1.parse(s, "\\[");
		for(int i = 0; i < s1.getFieldCount(); i++) {
			String b = s1.asString(i);
			if(b.length() == 0) {
				continue;
			}
			s2.parse(b, "\\]");
			if(s2.getFieldCount() != 2) {
				throw new IllegalArgumentException("String " + s + " is not of required form [id1]barcode1[id2]barcode2...");
			}
			barcodes.add(new Barcode(s2.asString(1), s2.asString(0), maxMismatches));
		}
		return new BarcodeSequence(barcodes);
	}
	
	/**
	 * Create and get a SAM attribute string representing this sequence of barcodes
	 * @return The SAM attribute string
	 */
	private String getSamAttributeString() {
		if(barcodes.isEmpty()) {
			return null;
		}
		Iterator<Barcode> iter = barcodes.iterator();
		String rtrn = "";
		while(iter.hasNext()) {
			Barcode bc = iter.next();
			rtrn += "[" + bc.getId() + "]" + bc.getSequence();
		}
		return rtrn;
	}
	
	/**
	 * Get the SAM attribute string representing this sequence of barcodes
	 * @return The SAM attribute string
	 */
	public String toSamAttributeString() {
		return samAttributeString;
	}
	
	public String toString() {
		return toSamAttributeString();
	}
	
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) {
			return false;
		}
		return toString().equals(o.toString());
	}
	
	public int hashCode() {
		return toString().hashCode();
	}
	
	@Override
	public int compareTo(BarcodeSequence o) {
		return toString().compareTo(o.toString());
	}

	/**
	 * Get the number of barcodes in the sequence
	 * @return The number of barcodes
	 */
	public int getNumBarcodes() {
		return barcodes.size();
	}

	/**
	 * Get the sequence of barcodes
	 * @return Ordered list of barcodes
	 */
	public List<Barcode> getBarcodes() {
		return barcodes;
	}
	
	
	
}

package nextgen.sequentialbarcode;


import java.util.List;

import org.apache.log4j.Logger;

import broad.core.parser.StringParser;

import net.sf.samtools.SAMRecord;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.sequentialbarcode.fragmentgroup.FragmentGroup;
import nextgen.sequentialbarcode.fragmentgroup.NamedBarcodedFragmentGroup;
import nextgen.sequentialbarcode.readlayout.Barcode;
import nextgen.sequentialbarcode.readlayout.BarcodeSet;
import nextgen.sequentialbarcode.readlayout.ReadLayout;
import nextgen.sequentialbarcode.readlayout.ReadSequenceElement;

/**
 * A basic implementation of a barcoded fragment
 * @author prussell
 *
 */
public class BarcodedFragmentImpl implements BarcodedFragment {
	
	protected String id;
	protected String read1sequence;
	protected String read2sequence;
	protected Annotation location;
	protected ReadLayout read1layout;
	protected ReadLayout read2layout;
	protected BarcodeSequence barcodes;
	public static Logger logger = Logger.getLogger(BarcodedFragmentImpl.class.getName());
	protected int barcodeMaxMismatches;
	protected FragmentGroup fragmentGroup;
	
	/**
	 * @param fragmentId Fragment ID
	 * @param barcodeSignature Barcodes for the fragment
	 * @param mappedChr Mapped chromosome for the fragment
	 * @param mappedStart Mapped start
	 * @param mappedEnd Mapped end
	 */
	public BarcodedFragmentImpl(String fragmentId, BarcodeSequence barcodeSignature, String mappedChr, int mappedStart, int mappedEnd) {
		this(fragmentId, barcodeSignature, new BasicAnnotation(mappedChr, mappedStart, mappedEnd));
	}
	
	/**
	 * @param fragmentId Fragment ID
	 * @param barcodeSignature Barcodes for the fragment
	 * @param mappedLocation Mapped location of the fragment
	 */
	public BarcodedFragmentImpl(String fragmentId, BarcodeSequence barcodeSignature, Annotation mappedLocation) {
		id = StringParser.firstField(fragmentId);
		barcodes = barcodeSignature;
		location = mappedLocation;
		fragmentGroup = new NamedBarcodedFragmentGroup(barcodes);
	}
	
	/**
	 * Instantiate from a SAM record by reading location and attributes
	 * @param samRecord SAM record
	 */
	public BarcodedFragmentImpl(SAMRecord samRecord) {
		
		String fragmentId = StringParser.firstField(samRecord.getReadName());
		read1sequence = null;
		read2sequence = null;
		try {
			String seq = samRecord.getReadString();
			if(samRecord.getFirstOfPairFlag()) {
				read1sequence = seq;
			} else {
				read2sequence = seq;
			}
		} catch (IllegalStateException e) {}
		String barcodeString = samRecord.getStringAttribute(BarcodedBamWriter.BARCODES_SAM_TAG);
		BarcodeSequence barcodeSignature = BarcodeSequence.fromSamAttributeString(barcodeString);
		
		id = fragmentId;
		barcodes = barcodeSignature;
		location = new BasicAnnotation(samRecord);
		fragmentGroup = NamedBarcodedFragmentGroup.fromSAMRecord(samRecord);
		
	}
	
	/**
	 * @param fragmentId Fragment ID
	 * @param read1seq Read1 sequence
	 * @param read2seq Read2 sequence
	 * @param barcodeSignature Barcodes for the fragment
	 */
	public BarcodedFragmentImpl(String fragmentId, String read1seq, String read2seq, BarcodeSequence barcodeSignature) {
		this(fragmentId, read1seq, read2seq, barcodeSignature, null);
	}

	/**
	 * @param fragmentId Fragment ID
	 * @param read1seq Read1 sequence
	 * @param read2seq Read2 sequence
	 * @param barcodeSignature Barcodes for the fragment
	 * @param mappedLocation Mapped location of the fragment
	 */
	public BarcodedFragmentImpl(String fragmentId, String read1seq, String read2seq, BarcodeSequence barcodeSignature, Annotation mappedLocation) {
		id = StringParser.firstField(fragmentId);
		read1sequence = read1seq;
		read2sequence = read2seq;
		barcodes = barcodeSignature;
		location = mappedLocation;
		fragmentGroup = new NamedBarcodedFragmentGroup(barcodes);
	}

	/**
	 * @param fragmentId Fragment ID
	 * @param read1seq Read1 sequence
	 * @param read2seq Read2 sequence
	 * @param layoutRead1 Read1 layout or null if not specified
	 * @param layoutRead2 Read2 layout or null if not specified
	 * @param maxMismatchesBarcode Max number of mismatches in each barcode when matching to reads
	 */
	public BarcodedFragmentImpl(String fragmentId, String read1seq, String read2seq, ReadLayout layoutRead1, ReadLayout layoutRead2, int maxMismatchesBarcode) {
		id = StringParser.firstField(fragmentId);
		read1sequence = read1seq;
		read2sequence = read2seq;
		read1layout = layoutRead1;
		read2layout = layoutRead2;
		barcodeMaxMismatches = maxMismatchesBarcode;
		fragmentGroup = new NamedBarcodedFragmentGroup(barcodes);
	}
	
	public BarcodeSequence getBarcodes() {
		if(barcodes == null) {
			findBarcodes();
		}
		return barcodes;
	}
	
	public void findBarcodes() {
		barcodes = new BarcodeSequence();
		if(read1layout != null && read1sequence != null) {
			List<List<ReadSequenceElement>> read1elements = read1layout.getMatchedElements(read1sequence);
			if(read1elements != null) {
				for(int i = 0; i < read1elements.size(); i++) {
					ReadSequenceElement parentElement = read1layout.getElements().get(i);
					if(parentElement.getClass().equals(Barcode.class) || parentElement.getClass().equals(BarcodeSet.class)) {
						for(ReadSequenceElement elt : read1elements.get(i)) {
							barcodes.appendBarcode((Barcode)elt);
						}
						continue;
					}
				}
			}
			read1elements = null;
		}
		if(read2layout != null && read2sequence != null) {
			List<List<ReadSequenceElement>> read2elements = read2layout.getMatchedElements(read2sequence);
			if(read2elements != null) {
				for(int i = 0; i < read2elements.size(); i++) {
					ReadSequenceElement parentElement = read2layout.getElements().get(i);
					if(parentElement.getClass().equals(Barcode.class) || parentElement.getClass().equals(BarcodeSet.class)) {
						for(ReadSequenceElement elt : read2elements.get(i)) {
							barcodes.appendBarcode((Barcode)elt);
						}
						continue;
					}
				}
			}
			read2elements = null;
		}
	}
	
	public String getId() {
		return id;
	}
	
	public String getRead1Sequence() {
		return read1sequence;
	}
	
	public String getRead2Sequence() {
		return read2sequence;
	}
	
	public ReadLayout getRead1Layout() {
		return read1layout;
	}
	
	public ReadLayout getRead2Layout() {
		return read2layout;
	}
	
	public Annotation getMappedLocation() {
		return location;
	}
	
	/**
	 * Set mapped location of fragment
	 * @param mappedLocation Mapped location
	 */
	public void setMappedLocation(Annotation mappedLocation) {
		location = mappedLocation;
	}
	
	public int compareTo(BarcodedFragment other) {
		if(location != null && other.getMappedLocation() != null) {
			int l = location.compareTo(other.getMappedLocation());
			if(l != 0) return l;
		}
		return id.compareTo(other.getId());
	}

	@Override
	public FragmentGroup getFragmentGroup() {
		return fragmentGroup;
	}

	@Override
	public void addFragmentWithSameBarcodes(BarcodedFragment fragment) {
		fragmentGroup.addFragment(fragment);
	}

	@Override
	public String getFullInfoString() {
		return location.getFullInfoString();
	}
	
	
}

package nextgen.sequentialbarcode.fragmentgroup;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import broad.core.parser.StringParser;

import net.sf.samtools.SAMRecord;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.sequentialbarcode.BarcodeSequence;
import nextgen.sequentialbarcode.BarcodedBamWriter;
import nextgen.sequentialbarcode.BarcodedFragment;
import nextgen.sequentialbarcode.BarcodedFragmentImpl;

/**
 * A collection of fragments sharing the same barcodes, storing read IDs and mapped locations
 * @author prussell
 *
 */
public class NamedBarcodedFragmentGroup implements FragmentGroup {
	
	private BarcodeSequence barcodes;
	private Collection<BarcodedFragment> fragments;
	
	/**
	 * @param barcodeSignature The shared barcodes
	 */
	public NamedBarcodedFragmentGroup(BarcodeSequence barcodeSignature) {
		this(barcodeSignature, new TreeSet<BarcodedFragment>());
	}
	
	/**
	 * @param barcodeSignature The shared barcodes
	 * @param barcodedFragments Some fragments with the barcodes
	 */
	public NamedBarcodedFragmentGroup(BarcodeSequence barcodeSignature, Collection<BarcodedFragment> barcodedFragments) {
		barcodes = barcodeSignature;
		fragments = new TreeSet<BarcodedFragment>();
		for(BarcodedFragment fragment : barcodedFragments) {
			addFragment(fragment);
		}
	}
	
	public void addFragment(BarcodedFragment fragment) {
		if(!fragment.getBarcodes().equals(barcodes)) {
			throw new IllegalArgumentException("New fragment must have barcodes " + barcodes.toString());
		}
		fragments.add(fragment);
	}
	
	/**
	 * Delimiter separating fragments when this is expressed as a SAM attribute
	 */
	private static String SAM_ATTRIBUTE_FRAGMENT_DELIM = ";";
	
	/**
	 * Read the barcodes and group of fragments from a SAM attributes found in a SAM record
	 * @param record SAM record
	 * @return Fragment group object represented in the SAM attribute
	 */
	public static NamedBarcodedFragmentGroup fromSAMRecord(SAMRecord record) {
		BarcodeSequence barcodeSignature = BarcodeSequence.fromSamRecord(record);
		return new NamedBarcodedFragmentGroup(barcodeSignature, NamedBarcodedFragmentGroup.fromSamAttributeFragmentGroup(record, barcodeSignature));
	}
	
	/**
	 * Read the group of fragments from a SAM attribute string
	 * @param barcodes The barcode set shared by these fragments
	 * @param fragmentGroupAttribute The fragment group attribute string from a SAM record
	 * @return Fragment group object with these barcodes and the set of fragments specified in the SAM attribute string
	 */
	public static NamedBarcodedFragmentGroup fromSamAttributeString(BarcodeSequence barcodes, String fragmentGroupAttribute) {
		return new NamedBarcodedFragmentGroup(barcodes, NamedBarcodedFragmentGroup.fromSamAttributeFragmentGroup(fragmentGroupAttribute, barcodes));
	}
	
	/**
	 * Read the barcode set and fragment group from SAM attribute strings
	 * @param barcodeAttribute Barcode attribute
	 * @param fragmentGroupAttribute Fragment group attribute
	 * @return Fragment group with barcode set and fragments specified in these SAM attributes
	 */
	public static NamedBarcodedFragmentGroup fromSamAttributeStrings(String barcodeAttribute, String fragmentGroupAttribute) {
		BarcodeSequence barcode = BarcodeSequence.fromSamAttributeString(barcodeAttribute);
		return new NamedBarcodedFragmentGroup(barcode, NamedBarcodedFragmentGroup.fromSamAttributeFragmentGroup(fragmentGroupAttribute, barcode));
	}
	
	public BarcodeSequence getBarcodes() {
		return barcodes;
	}
	
	public String toSamAttributeString() {
		if(fragments.isEmpty()) {
			throw new IllegalStateException("Fragment set is empty");
		}
		Iterator<BarcodedFragment> iter = fragments.iterator();
		BarcodedFragment f = iter.next();
		String rtrn = f.getFullInfoString();
		while(iter.hasNext()) {
			BarcodedFragment fr = iter.next();
			rtrn += SAM_ATTRIBUTE_FRAGMENT_DELIM + fr.getFullInfoString();
		}
		return rtrn;
	}

	/**
	 * Read the fragment group attribute from a SAM record and create fragment objects
	 * @param samRecord SAM record
	 * @param barcodeSequence Barcodes for this fragment group
	 * @return The fragment objects
	 */
	private static Collection<BarcodedFragment> fromSamAttributeFragmentGroup(SAMRecord samRecord, BarcodeSequence barcodeSequence) {
		return fromSamAttributeFragmentGroup(samRecord.getStringAttribute(BarcodedBamWriter.FRAGMENT_GROUP_SAM_TAG), barcodeSequence);
	}

	/**
	 * Create fragment objects from a fragment group SAM attribute
	 * @param attributeString The fragment group SAM attribute
	 * @param barcodeSequence The barcodes shared by these fragments
	 * @return The fragment objects
	 */
	private static Collection<BarcodedFragment> fromSamAttributeFragmentGroup(String attributeString, BarcodeSequence barcodeSequence) {
		StringParser s1 = new StringParser();
		Collection<BarcodedFragment> rtrn = new TreeSet<BarcodedFragment>();
		s1.parse(attributeString, SAM_ATTRIBUTE_FRAGMENT_DELIM);
		for(int i = 0; i < s1.getFieldCount(); i++) {
			String coord = s1.asString(i);
			Annotation a = BasicAnnotation.fromFullInfoString(coord);
			BarcodedFragment bf = new BarcodedFragmentImpl(a.getName(), barcodeSequence, a);
			rtrn.add(bf);
		}
		return rtrn;
	}

}

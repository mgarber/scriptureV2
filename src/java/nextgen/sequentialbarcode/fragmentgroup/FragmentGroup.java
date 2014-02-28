package nextgen.sequentialbarcode.fragmentgroup;

import java.util.Collection;

import nextgen.core.annotation.Annotation;
import nextgen.sequentialbarcode.BarcodeSequence;
import nextgen.sequentialbarcode.BarcodedFragment;

/**
 * A group of fragments sharing a barcode sequence
 * @author prussell
 *
 */
public interface FragmentGroup {

	/**
	 * Add another fragment to the group
	 * @param fragment The fragment
	 */
	public void addFragment(BarcodedFragment fragment);
	
	/**
	 * Get the barcodes shared by the fragments in the group
	 * @return The shared barcodes
	 */
	public BarcodeSequence getBarcodes();
	
	/**
	 * Express as a SAM attribute
	 * @return SAM attribute string
	 */
	public String toSamAttributeString();

	/**
	 * Get the regions covered by all the fragments as annotation objects
	 * @return The collection of regions
	 */
	public Collection<Annotation> getRegions();

}

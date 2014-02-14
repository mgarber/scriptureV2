package nextgen.sequentialbarcode;

import nextgen.core.annotation.Annotation;
import nextgen.sequentialbarcode.fragmentgroup.FragmentGroup;
import nextgen.sequentialbarcode.readlayout.ReadLayout;

/**
 * A fragment with a mapped location and collection of barcodes
 * @author prussell
 *
 */
public interface BarcodedFragment extends Comparable<BarcodedFragment> {
	
	/**
	 * Get the fragment ID
	 * @return Fragment ID
	 */
	public String getId();
	
	/**
	 * Get the sequence of read1
	 * @return Sequence of read1 or null if no read1
	 */
	public String getRead1Sequence();
	
	/**
	 * Get the sequence of read2
	 * @return Sequence of read2 or null if no read1
	 */
	public String getRead2Sequence();
	
	/**
	 * Get the expected layout of read1
	 * @return Read1 layout or null if not specified
	 */
	public ReadLayout getRead1Layout();
	
	/**
	 * Get the expected layout of read2
	 * @return Read2 layout or null if not specified
	 */
	public ReadLayout getRead2Layout();
	
	/**
	 * Get the mapped location of the fragment
	 * @return Mapped location as an annotation
	 */
	public Annotation getMappedLocation();
	
	/**
	 * Get a string with no whitespace that specifies all the info about this fragment
	 * @return Full info string
	 */
	public String getFullInfoString();
	
	/**
	 * Find the barcodes in the fragment
	 */
	public void findBarcodes();
		
	/**
	 * Get the sequence of barcodes found in the fragment
	 * @return The sequence of barcodes
	 */
	public BarcodeSequence getBarcodes();
	
	/**
	 * Get the fragment group that this fragment belongs to
	 * @return The fragment group
	 */
	public FragmentGroup getFragmentGroup();
	
	/**
	 * Add a fragment to the fragment group that this fragment belongs to, provided it has the same barcodes as this fragment
	 * @param fragment New fragment to add to group
	 */
	public void addFragmentWithSameBarcodes(BarcodedFragment fragment);

}

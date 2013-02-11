/**
 * 
 */
package nextgen.core.readFilters;

import java.util.Collection;

import nextgen.core.alignment.Alignment;
import nextgen.core.coordinatesystem.CoordinateSpace;

import org.apache.commons.collections15.Predicate;

/**
 * @author prussell
 *
 */
public class FragmentLengthFilter implements Predicate<Alignment> {

	private int minLength;
	private int maxLength;
	private CoordinateSpace coordSpace;
	
	/**
	 * Instantiate the filter with a maximum length and no minimum length
	 * @param coord Coordinate space
	 * @param maxLen Maximum genomic span for a fragment (inclusive)
	 */
	public FragmentLengthFilter(CoordinateSpace coord, int maxLen) {
		this(coord, 0, maxLen);
	}
	
	/**
	 * Instantiate the filter with minimum and maximum lengths
	 * @param coord Coordinate space
	 * @param minLen Minimum genomic span for a fragment (inclusive)
	 * @param maxLen Maximum genomic span for a fragment (inclusive)
	 */
	public FragmentLengthFilter(CoordinateSpace coord, int minLen, int maxLen) {
		minLength = minLen;
		maxLength = maxLen;
		coordSpace = coord;
	}
	
	
	/**
	 * If the region is within a gene, return false if the fragment size is invalid for any isoform, or true otherwise
	 * If the region is not within a gene, filter on genomic distance
	 */
	@Override
	public boolean evaluate(Alignment align) {

		Collection<Integer> sizes = align.getFragmentSize(coordSpace);
		if(sizes.isEmpty()) {
			// The region is not in a gene
			// Filter on genomic distance
			int genomicSize = align.getEnd() - align.getStart() + 1;
			return genomicSize >= minLength && genomicSize <= maxLength;
		}
		for(Integer i : sizes) {
			int iint = i.intValue();
			if(iint > maxLength) return false;
			if(iint < minLength) return false;
		}
		return true;
		
	}

}

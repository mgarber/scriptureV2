/**
 * 
 */
package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;
import org.apache.commons.collections15.Predicate;


/**
 * @author prussell
 *
 */
public class NumHitsFilter implements Predicate<Alignment>  {

	private int maxHits;
	
	/**
	 * @param maxNH Maximum value of NH tag (inclusive)
	 */
	public NumHitsFilter(int maxNH) {
		if(maxNH < 1) {
			throw new IllegalArgumentException("Max NH must be >= 1");
		}
		maxHits = maxNH;
	}
	
	@Override
	public boolean evaluate(Alignment align) {
		
		Object o = align.getAttribute("NH");
		if(o == null) {
			return true;
		}
		int nh = o.getClass().equals(String.class) ? Integer.parseInt((String)o) : ((Integer)o).intValue();
		return nh <= maxHits;
		
	}

}

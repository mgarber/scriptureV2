package xp.core.Utils;
import java.util.Comparator;
import broad.core.annotation.ShortBED;
/**
 *  Created on 2013-3-3  
 */
public class ShortBEDComparator implements Comparator<ShortBED>
{
	@Override
	public int compare(ShortBED a, ShortBED b) {
		return a.compareTo(b);
	}
	
}
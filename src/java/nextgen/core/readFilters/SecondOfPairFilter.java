/**
 * 
 */
package nextgen.core.readFilters;

import net.sf.samtools.SAMRecord;
import nextgen.core.alignment.Alignment;

import org.apache.commons.collections15.Predicate;

/**
 * @author prussell
 *
 */
public class SecondOfPairFilter implements Predicate<Alignment> {

	public SecondOfPairFilter() {}
	
	@Override
	public boolean evaluate(Alignment align) {
		SAMRecord samRecord = align.toSAMRecord();
		return samRecord.getSecondOfPairFlag();
	}

}

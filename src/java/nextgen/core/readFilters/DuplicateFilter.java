package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;

import org.apache.commons.collections15.Predicate;

/**
 * @author shari
 * Remove duplicate reads
 */
public class DuplicateFilter implements Predicate<Alignment> {
	@Override
	public boolean evaluate(Alignment align) {
		return !align.isDuplicate();
	}
}
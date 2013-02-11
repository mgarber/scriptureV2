package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;

import org.apache.commons.collections15.Predicate;

/**
 * @author engreitz
 * Remove chimeric reads
 */
public class ChimeraFilter implements Predicate<Alignment> {
	@Override
	public boolean evaluate(Alignment align) {
		return !align.isChimera();
	}
}

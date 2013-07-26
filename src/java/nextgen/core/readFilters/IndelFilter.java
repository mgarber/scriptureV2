package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;

import org.apache.commons.collections15.Predicate;

public class IndelFilter implements Predicate<Alignment>{

	@Override
	public boolean evaluate(Alignment align) {
		return !align.hasIndel();
	}

}

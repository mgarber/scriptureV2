package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;

import org.apache.commons.collections15.Predicate;

public class ReadEndsInFilter implements Predicate<Alignment> {
	
	private Annotation window;
	
	public ReadEndsInFilter(Annotation w){
		window = w;
	}
	
	@Override
	public boolean evaluate(Alignment read) {
		
		if(window.getStart()<read.getOrientedEnd() && window.getEnd()>read.getOrientedEnd())
			return true;
		else
			return false;
	}
}
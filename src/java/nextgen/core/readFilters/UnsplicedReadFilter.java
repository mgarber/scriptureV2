package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;

import org.apache.commons.collections15.Predicate;

public class UnsplicedReadFilter implements Predicate<Alignment>{

	@Override
	public boolean evaluate(Alignment align) {
		//if has splice junctions --> return false
		if(!align.getSpliceConnections().isEmpty()){
			return false;
		}
		return true;
	}

}

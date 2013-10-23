package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;

import org.apache.commons.collections15.Predicate;

public class NoSpliceFilter implements Predicate<Alignment>{

	@Override
	public boolean evaluate(Alignment align) {
		//if has does not have splice junctions --> return true
		if(align.getSpliceConnections().isEmpty()){
			return true;
		}
		return false;
	}

	public boolean equals(Object other){
		if(other instanceof SplicedReadFilter){
			return true;
		}
		return false;
	}

	@Override
	public int hashCode() {
		throw new UnsupportedOperationException("Class should override hashCode() because it overrides equals()");
	}
	
}

package nextgen.core.readFilters;

import java.util.Iterator;

import org.apache.commons.collections15.Predicate;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;

public class PairedEndFilter implements Predicate<Alignment> {
	@Override
	public boolean evaluate(Alignment read) {
		//check if the read is paired and both mates on same chromosome
		
		if(read.isPaired()){
			//TODO: SK: REMOVE THIS CHECK ONCE IT IS SURE THAT PAIREDENDWRITER IS NOT WRITING THESE READS AT ALL
			Iterator<Annotation> mates = read.getReadAlignments(null).iterator();
			String refName = null;
			while(mates.hasNext()){
				if(refName==null)
					refName = mates.next().getReferenceName();
				else
					if(refName.equals(mates.next().getReferenceName())){
						return true;
					}
			}
		}
		return false;
	}
}

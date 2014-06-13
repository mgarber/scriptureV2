package nextgen.core.model.score;

import org.apache.commons.collections15.Predicate;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.AlignmentModel.AlignmentCount;
import nextgen.core.readFilters.ReadStartsInFilter;

public class StartInWindowCountScore extends WindowScore.AbstractWindowScore implements Comparable<StartInWindowCountScore>{

	private double count;  
	
	public StartInWindowCountScore(Annotation t) {
		super(t);
	}

	public StartInWindowCountScore(AlignmentModel model, Annotation annotation) {
		super(annotation);
		
		Predicate<Alignment> filter=new ReadStartsInFilter(annotation);
		CloseableIterator<AlignmentCount> iter = model.new WrapAlignmentCountIterator(new CloseableFilterIterator<Alignment>(model.getOverlappingReads(annotation,false), filter));
		setCount(model.getCount(iter));
	}
	
	@Override
	public double getScore() { 
		return getCount();
	}
	
	public double getCount() { 
		return count; 
	}
	
	@Override
	public int compareTo(StartInWindowCountScore o) {
		
		// First compare counts
		double otherCount = o.getCount();
		if(count < otherCount) return -1;
		if(count > otherCount) return 1;
		
		// Then compare annotations
		return getAnnotation().compareTo(o.getAnnotation());
	}
	
	public void setCount(double count) {
		this.count = count;
	}

}

package nextgen.core.scripture;

import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;

public class SpliceJunction extends BasicAnnotation {

	/**
	 * The number of spliced reads that support this junction
	 */
	int count;
	/**
	 * The farthest start of the exon at the 5' end
	 */
	int farthestStart;
	/**
	 * The farthest end of the exon at the 3' end
	 */
	int farthestEnd;
	
	public SpliceJunction(Annotation junction) {		
		super(junction);
		count = 1;
	}
	
	public SpliceJunction(Annotation junction,int start,int end) {		
		this(junction);
		setStart(start);
		setEnd(end);
	}

	/**
	 * Returns true if the junctions are the same
	 * TODO: Make sure this criteria is all-inclusive
	 * @param junction
	 * @return
	 */
	public boolean isSameAs(Annotation junction){
	/*	if(this.getStart()==junction.getStart() && this.getEnd()==junction.getEnd()){
			return true;
		}
		return false;
	*/	
		return(this.equals(junction, true));
	}
	
	public void setStart(int start){
		farthestStart=start;
	}
	
	public void setEnd(int end){
		farthestEnd=end;
	}

	public void update(Annotation junction, int start, int end) {
		count++;
		if(start<farthestStart)
			setStart(start);
		if(end>farthestEnd)
			setEnd(end);
	}
	
	public Collection<Annotation> getFlankingExons(){
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		rtrn.add(new BasicAnnotation(this.getChr(),farthestStart,this.getStart(),this.getOrientation()));
		rtrn.add(new BasicAnnotation(this.getChr(),this.getEnd(),farthestEnd,this.getOrientation()));
		return rtrn;
	}
	
	public Assembly toAssembly(){
		Assembly rtrn = new Assembly(getFlankingExons());
		rtrn.setOrientation(getOrientation());
		rtrn.setName(getName());
		return (rtrn);
	}
	
	public Annotation getJunction(){
		return(this);
	}
	
	public int getCount(){
		return count;
	}
}

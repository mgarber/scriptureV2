package nextgen.core.scripture;

import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;

/**
 * This class represents a NON-SPLICED annotation - an exon
 * @author skadri
 *
 */
public class Exon extends Assembly {

	boolean isConsidered;
	
	public Exon(Annotation junction) {		
		super(junction);
		this.setOrientation(junction.getOrientation());
		isConsidered = false;
	}
	
	public Exon(Annotation junction,int counter) {		
		this(junction);
		this.setName("exon_"+counter);
	}
	
	/**
	 * Returns true if the exon was considered to be a part of an assembly.
	 * @return
	 */
	public boolean isConsidered(){
		return isConsidered;
	}
	
	public void setConsideredFlag(boolean value){
		isConsidered = value;
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

	public Assembly toAssembly(){
		Assembly rtrn = new Assembly(this);		
		rtrn.setOrientation(getOrientation());
		rtrn.setName(getName());
		return (rtrn);
	}
	
	public Annotation getJunction(){
		return(this);
	}
	

}

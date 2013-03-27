package nextgen.core.scripture;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;

/**
 * This class represents an assembly
 * @author mguttman
 *
 */
public class Assembly extends BasicAnnotation{

	private boolean isPossiblePremature;
	private boolean isSpurious;
	
	public Assembly(Annotation annotation){
		super(annotation);
		this.isPossiblePremature=false;
		this.isSpurious=false;
	}
	
	public Assembly(Annotation annotation, boolean isPremature){
		super(annotation);
		this.isPossiblePremature=isPremature;
		this.isSpurious=false;
	}
	
	public Assembly(Collection<? extends Annotation> blocks, boolean isPremature){
		super(blocks);
		this.isPossiblePremature=isPremature;
		this.isSpurious=false;
	}
	
	public Assembly(Collection<? extends Annotation> blocks){
		super(blocks);
		this.isPossiblePremature=false;
		this.isSpurious=false;
	}
	
	public void setPossiblePremature(boolean premature){
		this.isPossiblePremature=premature;
	}
	
	public boolean getPossiblePremature(){
		return this.isPossiblePremature;
	}
	
	public void setSpurious(boolean spurious){
		this.isSpurious = spurious;
	}
	
	public boolean isSpurious(){
		return this.isSpurious;
	}

	public Iterator<Assembly> trimNodes() {
		Collection<Assembly> rtrn=new ArrayList<Assembly>();
		//walk back one intron at a time from the assembly and test compatibility
		List<? extends Annotation> blocks=this.getBlocks();
		for(int i=blocks.size()-1; i>=0; i--){
			blocks.remove(i);
			if(!blocks.isEmpty()){
				Assembly current=new Assembly(blocks, false);
				rtrn.add(current);
			}
		}
		return rtrn.iterator();
	}
		
}

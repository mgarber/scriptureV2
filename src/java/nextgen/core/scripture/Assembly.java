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
	private int isConfident;
	
	public Assembly(Annotation annotation){
		super(annotation);
		this.isPossiblePremature=false;
		this.isSpurious=false;
		//Unset = 0
		isConfident = 0;
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
		return isSpurious;
	}

	public boolean isConfidentIsSet(){
		if(isConfident==0)
			return false;
		else
			return true;
	}
	
	public void setConfident(boolean value){
		if(value){
			isConfident = 1;
		}
		else{
			isConfident = 2;
		}
	}
	
	public boolean isConfident(){
		
		if(isConfident == 1)
			return true;
		else
			return false;
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
	
	public Assembly trim(Annotation read) {
		//Collection<Assembly> rtrn=new ArrayList<Assembly>();
		int size = getBlocks().size();
		//walk back one intron at a time from the assembly and test compatibility
		List<? extends Annotation> blocks=getBlocks();
		if(read.overlaps(blocks.get(size-1)) && (blocks.get(size-1).getStart())<(read.getBlocks().get(0).getEnd())){
			blocks.get(size-1).setEnd(read.getBlocks().get(0).getEnd());
			Assembly current=new Assembly(blocks, false);
			return current;
		}
		return null;
	}
	
	public Annotation getLastBlock(){
		return getBlocks().get(getBlocks().size()-1);
	}
}

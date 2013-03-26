package nextgen.core.alignment;

import java.util.ArrayList;
import java.util.Collection;

import net.sf.samtools.SAMRecord;
import broad.core.datastructures.Pair;

/**
 * Represents a paired alignment in all detail, with possbily multiple mappings
 * through the genome. 
 * An alignment pair is meant to be built iteratively as mappings of an insert are 
 * found. Once all pieces are assembled, it can do useful things.
 * @author mgarber
 *
 */
public class AlignmentPair extends Pair<Collection<SAMRecord>> {

	private Pair<Integer> numHits = new Pair<Integer>();
	

	public void add(SAMRecord record) {
		if(record.getFirstOfPairFlag()){
			Collection<SAMRecord> set=new ArrayList<SAMRecord>();
			if(hasValue1()){set=getValue1();}
			set.add(record);
			setValue1(set);
		}
		else{
			Collection<SAMRecord> set=new ArrayList<SAMRecord>();
			if(hasValue2()){set=getValue2();}
			set.add(record);
			setValue2(set);
		}
		updateNumHits(record);
	}
	
	public boolean isComplete() {
		//First check that we have both ends
		if(hasValue1() && hasValue2()){
			//if so, check that the size matches the expected numHits
			int size1=getValue1().size();
			int size2=getValue2().size();
			
			int expectedSize1=numHits.getValue1().intValue();
			int expectedSize2=numHits.getValue2().intValue();
			
			if(size1==expectedSize1 && size2==expectedSize2){
				return true;
			}
		}
		
		return false;
	}
	

	public Collection<SAMRecord> makePairs() {
		Collection<SAMRecord> rtrn=new ArrayList<SAMRecord>();
		
		Collection<SAMRecord> pair1=getValue1();
		Collection<SAMRecord> pair2=getValue2();
		
		//match up all pair1 and pair2
		for(SAMRecord r1: pair1){
			for(SAMRecord r2: pair2){
				if(isCompatiblePair(r1, r2)){
					Pair<SAMRecord> p=new Pair<SAMRecord>(r1, r2);
					SAMRecord fragment=makePair(p);
					if(r1.getNotPrimaryAlignmentFlag() && r2.getNotPrimaryAlignmentFlag()){
						fragment.setNotPrimaryAlignmentFlag(true);
					}
					rtrn.add(fragment);
				}
			}
		}
		
		return rtrn;
	}
	
	
	private boolean isCompatiblePair(SAMRecord r1, SAMRecord r2) {
		if(r1.getReferenceName()==r2.getReferenceName()){
			if((r1.getAlignmentStart()==r2.getMateAlignmentStart())&& (r1.getMateAlignmentStart()==r2.getAlignmentStart())){
				if(r1.getMateNegativeStrandFlag()==r2.getReadNegativeStrandFlag() && r1.getReadNegativeStrandFlag()==r2.getMateNegativeStrandFlag()){
					if(r1.getMateReferenceName().equalsIgnoreCase(r2.getReferenceName()) && r2.getMateReferenceName().equalsIgnoreCase(r1.getReferenceName())){
						return true;
					}
				}
			}
		}
		return false;
	}


	private SAMRecord makePair(Pair<SAMRecord> pair) {
		SingleEndAlignment a1 = new SingleEndAlignment(pair.getValue1());
		SingleEndAlignment a2 = new SingleEndAlignment(pair.getValue2());
		SAMRecord rtrn = new FragmentAlignment(a1, a2).toSAMRecord();

		return rtrn;
	}



	
	private void updateNumHits(SAMRecord record) {
		Object nh=record.getAttribute("NH");
		
		if(nh!=null){
			int num=new Integer(nh.toString());
			if(record.getFirstOfPairFlag()){
				numHits.setValue1(num);
			}
			else{
				numHits.setValue2(num);
			}
		}
		else{
			//The NH flag is not set so we will default to 1,1
			numHits=new Pair<Integer>(1,1);
		}
	}
	
}

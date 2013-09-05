package nextgen.core.alignment;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecord.SAMTagAndValue;
import nextgen.core.writers.PairedEndWriter;
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
	static Logger logger = Logger.getLogger(AlignmentPair.class.getName());

	public void add(SAMRecord record) {
		if(record.getFirstOfPairFlag()){
			Collection<SAMRecord> set=new ArrayList<SAMRecord>();
			if(hasValue1()){set=getValue1();}
			set.add(record);
			setValue1(set);
		}
		else{
			if(record.getSecondOfPairFlag()){
				Collection<SAMRecord> set=new ArrayList<SAMRecord>();
				if(hasValue2()){
					set=getValue2();
				}
				set.add(record);
				setValue2(set);
			}
			else{
				logger.error("Record has a problem with the first/second of pair flag");
			}
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
			else{
//				System.out.println("Size1="+size1+" ExpectedSize1="+expectedSize1+" Size2="+size2+" ExpectedSize2="+expectedSize2);
			}
		}
		
		return false;
	}
	

	public Collection<SAMRecord> makePairs() {
		Collection<SAMRecord> rtrn=new ArrayList<SAMRecord>();
		
		Collection<SAMRecord> pair1=getValue1();
		Collection<SAMRecord> pair2=getValue2();
		
		int adjustedNumHits = pair1.size() * pair2.size();
		//match up all pair1 and pair2
		for(SAMRecord r1: pair1){
			for(SAMRecord r2: pair2){
				if(isCompatiblePair(r1, r2)){
					Pair<SAMRecord> p=new Pair<SAMRecord>(r1, r2);
					SAMRecord fragment= makePair(p);
					fragment.setAttribute("NH", adjustedNumHits); //MG moved this from PairedEndWriter to here, not sure why was there.
					fragment.setMateUnmappedFlag(false);
					fragment.setNotPrimaryAlignmentFlag(r1.getNotPrimaryAlignmentFlag() && r2.getNotPrimaryAlignmentFlag());
					rtrn.add(fragment);
				}
			}
		}
		
		return rtrn;
	}
	
	
	private boolean isCompatiblePair(SAMRecord r1, SAMRecord r2) {

		return r1.getReferenceName()==r2.getReferenceName() && 
				r1.getMateNegativeStrandFlag()==r2.getReadNegativeStrandFlag() && r1.getReadNegativeStrandFlag()==r2.getMateNegativeStrandFlag() &&
						r1.getMateReferenceName().equalsIgnoreCase(r2.getReferenceName()) && r2.getMateReferenceName().equalsIgnoreCase(r1.getReferenceName()) ;
	}


	private SAMRecord makePair(Pair<SAMRecord> pair) {
		long start = System.nanoTime();
		SingleEndAlignment a1 = new SingleEndAlignment(pair.getValue1());
		long a1Time = System.nanoTime() - start;
		start = System.nanoTime();
		SingleEndAlignment a2 = new SingleEndAlignment(pair.getValue2());
		long a2Time = System.nanoTime() - start;
		start = System.nanoTime();
		SAMRecord rtrn = new FragmentAlignment(a1, a2).toSAMRecord();
		long frgTime = System.nanoTime() - start;
		
		long ta1 =  Math.round(a1Time/(double)1000000);
		long ta2 = 	Math.round(a2Time/(double)1000000);
		long tFrg = Math.round(frgTime/(double)1000000);
		
		//if(ta1 + ta2 + tFrg > 50) {
		//	logger.info(String.format("Times: a1: %d\ta2: %d\tfrg: %d", ta1, ta2, tFrg));
		//}
		return rtrn;
	}



	
	private void updateNumHits(SAMRecord record) {
		Object nh=record.getAttribute("NH");
		
		if(nh!=null){
			int num=new Integer(nh.toString());
			//if the next hit is on a different chromosome
			Object cc=record.getAttribute("CC");			
			//IF NH > 1
				// IF NEXT ALIGNMENT IS NOT ON THE SAME CHROMOSOME
				//TODO: CHECK IF DISTANCE TO TOO FAR, THEN WRITE THIS ONE
				// OR 
				// IF NH>1 && cc==null meaning this is the last alignment in a multi-mapped read alignment
			if(num>1 && ((cc!=null && !cc.toString().equals("="))||(cc==null))){
				if(record.getFirstOfPairFlag()){
					numHits.setValue1(getValue1().size());
				}
				else{
					numHits.setValue2(getValue2().size());
				}
			}
			// IF NH ==1
			// OR
			// IF NH>1 && NEXT ALIGNMENT IS ON THE SAME CHROMOSOME
			else{
				if(record.getFirstOfPairFlag()){
					numHits.setValue1(num);
				}
				else{
					numHits.setValue2(num);
				}
			}
		}
		else{
			//The NH flag is not set so we will default to 1,1
			numHits=new Pair<Integer>(1,1);
		}
	}
	
	public boolean valuesAreEmpty(){
		boolean isEmpty=true;
		if(this.hasValue1()){
			if(!this.getValue1().isEmpty()){
				return false;
			}
		}
		if(this.hasValue2()){
			if(!this.getValue2().isEmpty()){
				return false;
			}
		}
		return true;
	}
}

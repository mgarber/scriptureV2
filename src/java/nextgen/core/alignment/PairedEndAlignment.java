package nextgen.core.alignment;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import com.mysql.jdbc.StringUtils;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.feature.GenomeWindow;
import nextgen.core.feature.Window;
import nextgen.core.annotation.*;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.readers.PairedEndWriter;

public class PairedEndAlignment extends BasicAnnotation implements Alignment {

	static Logger logger = Logger.getLogger(PairedEndAlignment.class.getName());

    //First in pair alignment
    SingleEndAlignment firstMate;
    //Second in pair alignment
    SingleEndAlignment secondMate;
    //Whether first of pair or second of pair are in direction of transcription. Unstranded option in case of an unstranded library
    TranscriptionRead txnRead;
    Map<String,String> attributeMap;
    boolean isProperPair;
	
    /**
     * Constructs a paired end alignment object from two alignments.
     * @param firstMate
     * @param secondMate
     * @throws ChromosomeInconsistencyException 
     */
    public PairedEndAlignment(SingleEndAlignment firstMate, SingleEndAlignment secondMate) {
    	super(firstMate);   //TODO How to deal with pairs that have different chromosomes?

    	this.firstMate = firstMate;
        this.secondMate = secondMate;
        attributeMap = new HashMap<String,String>();
    }
    
    /**
     * Constructs a paired end alignment object from two alignments and provides the read that is in the direction of transcription
     * @param firstMate
     * @param secondMate
     * @param strand
     */
    public PairedEndAlignment(SingleEndAlignment firstMate, SingleEndAlignment secondMate,TranscriptionRead strand) {
    	super(firstMate);   //TODO How to deal with pairs that have different chromosomes?

    	this.firstMate = firstMate;
        this.secondMate = secondMate;
        this.txnRead = strand;
        attributeMap = new HashMap<String,String>();
    }
     
   	/**
     * Returns the strand for the first of pair
     * @return
     */
    Strand getFirstOfPairStrand(){
    	//Mate will always exist because object is only created if both mates exist.
    	return this.firstMate.getFragmentStrand();
    }

    /**
     * Returns the strand for the second of pair
     * @return
     */
    Strand getSecondOfPairStrand(){
    	//Mate will always exist because object is only created if both mates exist.
    	return this.secondMate.getFragmentStrand();
    }

	/**
	 * Returns the name of this read pair.
	 */
	public String getReadName() {
		return firstMate.getReadName();
	}

	/**
	 * Returns the chromosome for this alignment
	 */
	public String getChromosome() {
		return this.getReferenceName();
	}

	/**
	 * Returns the sum of the mapping qualities of both mates
	 */
	public int getMappingQuality() {
		
		return this.firstMate.getMappingQuality()+this.secondMate.getMappingQuality();
	}

	/**
	 * Always returns true for objects of this class.
	 */
	public boolean isPaired() {
		return true;
	}

	/**
     * For paired end, returns true if the read in direction of transcription is negative stranded
     * @return
     */
	@Override
	public boolean isNegativeStrand() {
		if ("-".equals(this.getFragmentStrand()))
			return true;
		else
			//Even unstranded would be false
			return false;
	}

	/**
	 * Returns true is the alignment pair is a duplicate
	 */
	public boolean isDuplicate() {
		return firstMate.isDuplicate();
	}	

	/**
	 * Sets the duplicate flag
	 */
	public void setDuplicateFlag(boolean duplicateFlag) {
		firstMate.setDuplicateFlag(duplicateFlag);
		secondMate.setDuplicateFlag(duplicateFlag);
	}
	
	/**
	 * Get the first read
	 * @return The first read in pair
	 */
	public Alignment getFirstMate() {
		return this.firstMate;
	}
	
	/**
	 * Get the second read
	 * @return The second read in pair
	 */
	public Alignment getSecondMate() {
		return this.secondMate;
	}
	
	/**
	 * Returns an object for the fragment between the read pair in the specified coordinate space
	 */
	public Collection<? extends Window> getFragment(CoordinateSpace C) {
		if (C==null) {
			Collection<Window> rtrn=new TreeSet<Window>();
			rtrn.add(new GenomeWindow(this.getChr(), this.getFragmentStart(), this.getFragmentEnd()));
			return rtrn;
		}
		return (C.getFragment(this.getChr(), this.getFragmentStart(), this.getFragmentEnd()));
	}

	/**
	 * Returns the strand for the fragment, depending on the transcription read
	 */
	public Strand getFragmentStrand() {
		if(this.txnRead==TranscriptionRead.FIRST_OF_PAIR)
			return this.firstMate.getFragmentStrand();
		else if(this.txnRead==TranscriptionRead.SECOND_OF_PAIR)
			return this.secondMate.getFragmentStrand();
		else
			return Strand.UNKNOWN;
	}
	
	/**
	 * Sets the strand for the fragment, depending on the transcription read
	 * For FIRST, "first"
	 * For SECOND, "second"
	 * For unstranded "none"
	 */
	public void setFragmentStrand(String strand) {
		if(strand.equalsIgnoreCase("first"))
			txnRead=TranscriptionRead.FIRST_OF_PAIR;
		else if(strand.equalsIgnoreCase("second"))
			txnRead=TranscriptionRead.SECOND_OF_PAIR;
		else if(strand.equalsIgnoreCase("none"))
			txnRead=TranscriptionRead.UNSTRANDED;
		else{
			logger.error("Fragment strand set to unknown");
			//TODO:Return some error
		}
			
	}
	
	/**
	 * Sets the strand for the fragment to the strand passed as argument
	 * @param strand
	 */
	public void setFragmentStrand(TranscriptionRead strand) {
		txnRead = strand;
	}
	
    /**
     * This method checks whether the specified location on the specified 
     * chromosome is contained within the alignment
     * @param location
     * @return
     */
	public boolean contains(String chromosome, double location,CoordinateSpace C) {
		//TODO: Check in Coordinate Space: Will implement after populating coordinatespace
		if(!this.getChr().equals(chromosome))
			return false;
		else{
			if(location>=this.getFragmentStart() && location<=this.getFragmentEnd())
				return true;
			else
				return false;
		}
	}

	/**
	 * Returns the size of the fragment
	 */
	public Collection<Integer> getFragmentSize(CoordinateSpace C) {
		Collection<Integer> rtrn=new ArrayList<Integer>();
		Collection<? extends Window> fragments = this.getFragment(C);
		
		if(fragments == null) {
			rtrn.add(Integer.valueOf(this.getEnd() - this.getStart()));
			return rtrn;
		}
		
		for(Window w: fragments){
			rtrn.add(w.getSize());
		}
		
		return rtrn;
	}

	/**
	 * Returns the start of the fragment
	 */
	public int getFragmentStart() {
		return Math.min(this.firstMate.getFragmentStart(), this.secondMate.getFragmentStart());
	}

	/**
	 * Returns the end of the fragment
	 */
	public int getFragmentEnd() {
		return Math.max(this.firstMate.getFragmentEnd(), this.secondMate.getFragmentEnd());
	}

	/**
	 * Sets the specified attribute value
	 */
	public void setAttribute(String attribute, String value) {
		attributeMap.put(attribute, value);
	}

	/**
	 * Returns the value of the specified attribute
	 */
	public String getAttribute(String attribute) {
		return attributeMap.get(attribute);
	}
	
	/**
	 * This enumerator represents with of the mates in the read pair (or neither for unstranded library) is in the direction of transcription
	 * @author skadri
	 *
	 */
	public static enum TranscriptionRead{
		FIRST_OF_PAIR,SECOND_OF_PAIR,UNSTRANDED
	}

	@Override
	public Annotation getReadAlignmentBlocks(CoordinateSpace C) {
		/*Collection<Annotation> blocks=new TreeSet<Annotation>();
		blocks.addAll(this.firstMate.getReadAlignmentBlocks(C).getBlocks());
		blocks.addAll(this.secondMate.getReadAlignmentBlocks(C).getBlocks());
		Annotation alignmentBlock=new Gene(blocks, this.firstMate.getName());
		return alignmentBlock;*/
		return firstMate.union(secondMate);
	}

	public String toString(){
		/*
		 * add name into the bed string by @zhuxp
		 */
		String bed = this.getReadAlignmentBlocks(null).toBED();
		String[] a=bed.split("\t");
		a[3]=this.getName();
		return StringUtil.join("\t",a);
		
	}

	@Override
	public int getStart() {
		return this.getAlignmentStart();
	}

	@Override
	public int getEnd() {
		return this.getAlignmentEnd();
	}

	@Override
	public String getChr() {
		return this.getReferenceName();
	}

	@Override
	public int getAlignmentStart() {
		return Math.min(this.firstMate.getAlignmentStart(), this.secondMate.getAlignmentStart());
	}

	@Override
	public int getAlignmentEnd() {
		return Math.max(this.firstMate.getAlignmentEnd(), this.secondMate.getAlignmentEnd());
	}

	@Override
	public boolean isMapped() {
		return true;
	}
	
	@Override
	public boolean isChimera() {
		return !(firstMate.getReferenceName().equals(secondMate.getReferenceName()));
	}

	@Override
	public String getReadSequence() {
		return this.firstMate.getReadSequence()+this.secondMate.getReadSequence();
	}

	@Override
	public double getWeight() {
		if(this.attributeMap.containsKey("NH")){
			return 1.0/new Double(this.attributeMap.get("NH").toString());
		}
		return 1.0;
	}

	@Override
	public boolean isProperPair() {
		return this.isProperPair;
	}

	@Override
	public void setProperPairFlag(boolean properPairFlag) {
		this.isProperPair=properPairFlag;
	}
	
	@Override
	public void shift(int delta) {
		super.shift(delta);
		firstMate.shift(delta);
		secondMate.shift(delta);
	}
	
	@Override
	public void moveToCoordinate(int coord) {
		shift(coord - getStart());
	}
	
	@Override
	public List<? extends Annotation> getBlocks() {
		List<Annotation> blocks = new ArrayList<Annotation>();
		blocks.addAll(firstMate.getBlocks());
		blocks.addAll(secondMate.getBlocks());
		return blocks;
	}
	
	@Override
	public final SAMRecord toSAMRecord() {
		SAMRecord record = firstMate.toSAMRecord();
		
		record.setAttribute(PairedEndWriter.readStartFlag, firstMate.getSAMStart());
		record.setAttribute(PairedEndWriter.readCigarFlag, firstMate.getCigarString());
		record.setAttribute(PairedEndWriter.mateSequenceFlag, secondMate.getReadSequence());
		record.setAttribute(PairedEndWriter.mateCigarFlag, secondMate.getCigarString());
		
		record.setMateAlignmentStart(secondMate.getSAMStart());
		// add by @zhuxp
        record.setAlignmentStart(this.getAlignmentStart()+1);
     
		// end of add (test version)
        
        Annotation fragment = getReadAlignmentBlocks(null);
		record.setCigarString(fragment.getLengthOnReference() + "M");  // NOTE: losing information about indels in the SingleEndAlignments		
		return record;
	}
	
	@Override
	public Strand getOrientation(){
		return getFragmentStrand();
	}

	@Override
	public Collection<Annotation> getReadAlignments(CoordinateSpace space) {
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		rtrn.add(this.firstMate.getReadAlignmentBlocks(space));
		rtrn.add(this.secondMate.getReadAlignmentBlocks(space));
		//set orientation to strand orientation
		for(Annotation m:rtrn){
			m.setOrientation(this.getFragmentStrand());
		}
		return rtrn;
	}

	@Override
	public Collection<? extends Annotation> getSpliceConnections() {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		rtrn.addAll(this.firstMate.getSpliceConnections());
		rtrn.addAll(this.secondMate.getSpliceConnections());
		return rtrn;
	}

	@Override
	public boolean hasIndel() {
		return this.firstMate.hasIndel() || this.secondMate.hasIndel();
	}
	
	@Override
	public boolean equals (Annotation other){
		if(other instanceof PairedEndAlignment){
			PairedEndAlignment otherPair=(PairedEndAlignment) other;
			if(otherPair.getFirstMate().equals(this.getFirstMate()) && otherPair.getSecondMate().equals(this.getSecondMate())){
				return true;
			}
		}
		return false;
	}

	/**
	 * Get the beginning position of the fragment considering strand
	 * @return Highest position if negative strand, lowest position otherwise
	 */
	@Override
	public int getFirstFragmentPositionStranded() {
		Strand strand = getFragmentStrand();
		if(strand.equals(Strand.NEGATIVE)) return getFragmentEnd();
		return getFragmentStart();
	}

	/**
	 * Get the ending position of the fragment considering strand
	 * @return Lowest position if negative strand, highest position otherwise
	 */
	@Override
	public int getLastFragmentPositionStranded() {
		Strand strand = getFragmentStrand();
		if(strand.equals(Strand.NEGATIVE)) return getFragmentStart();
		return getFragmentEnd();
	}
}

package broad.core.multiplealignment;

import java.util.List;

import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;
import broad.core.sequence.Sequence;

import junit.framework.TestCase;

public class AlignedSequenceTest extends TestCase {
	String seqNoGaps = "ACTGAAAGGCCTGGAGGTAA";
	String seqOneGap = "ACGGAACCTGGA-AACTG";
	String seqFewGaps = "A--CCGTAAT--CGTAAAATTGCCG---CCAG-CCGA---AACCGA--A";
	
	public void testGetUngappedChunks() {
		String seqBases = "A-A";
		Sequence seq = new Sequence("seq");
		seq.setSequenceBases(seqBases);
		AlignedSequence alnSeq = new AlignedSequence(seq);
		
		List<int []> ungappedRegions = alnSeq.findUngappedChunks();
		
		assertEquals(2, ungappedRegions.size());
		assertEquals(0, ungappedRegions.get(0)[0]);
		assertEquals(1, ungappedRegions.get(0)[1]);
		assertEquals(2, ungappedRegions.get(1)[0]);
		assertEquals(3, ungappedRegions.get(1)[1]);
		
		alnSeq.setSequenceBases(seqOneGap);
		ungappedRegions = alnSeq.findUngappedChunks();
		
		assertEquals(2, ungappedRegions.size());
		assertEquals(0, ungappedRegions.get(0)[0]);
		assertEquals(12, ungappedRegions.get(0)[1]);
		assertEquals(13, ungappedRegions.get(1)[0]);
		assertEquals(18, ungappedRegions.get(1)[1]);
	}
	
	public void testGapAdjustedCoordinate() throws Exception {
		Sequence baseSequence = new Sequence("seq");
		baseSequence.setSequenceBases(seqNoGaps);
		AlignedSequence seq = new AlignedSequence(baseSequence);
		
		assertEquals(6, seq.getGapAdjustedCoordinate(6));
		
		seq.setSequenceBases(seqOneGap);
		assertEquals(11, seq.getGapAdjustedCoordinate(11));
		assertEquals(14, seq.getGapAdjustedCoordinate(13));
		assertEquals(13, seq.getGapAdjustedCoordinate(12));
		
		seq.setSequenceBases(seqFewGaps);
		assertEquals(0, seq.getGapAdjustedCoordinate(0));
		assertEquals(1 + 2, seq.getGapAdjustedCoordinate(1));
		assertEquals(2 + 2, seq.getGapAdjustedCoordinate(2));
		assertEquals(9  + 4,seq.getGapAdjustedCoordinate(9));
		assertEquals(18 + 4,seq.getGapAdjustedCoordinate(18));
		assertEquals(23 + 7, seq.getGapAdjustedCoordinate(23));
		assertEquals(26 + 8, seq.getGapAdjustedCoordinate(26));
		assertEquals(30 + 11, seq.getGapAdjustedCoordinate(30));
		assertEquals(35 + 13, seq.getGapAdjustedCoordinate(35));
	
	}
	
	public void testGetSequencePosition() throws Exception {
		Sequence baseSequence = new Sequence("seq");
		baseSequence.setSequenceBases(seqNoGaps);
		AlignedSequence seq = new AlignedSequence(baseSequence);
		
		assertEquals(6, seq.getSequencePosition(6));
		
		seq.setSequenceBases(seqOneGap);
		assertEquals(11, seq.getSequencePosition(11));
		assertEquals(13, seq.getSequencePosition(14));
		assertEquals(12, seq.getSequencePosition(13));
		
		seq.setSequenceBases(seqFewGaps);
		assertEquals(0, seq.getSequencePosition(0));
		assertEquals(1, seq.getSequencePosition(1 + 2));
		assertEquals(2, seq.getSequencePosition(2 + 2));
		assertEquals(9 ,seq.getSequencePosition(9 + 4));
		assertEquals(18,seq.getSequencePosition(18 + 4));
		assertEquals(23, seq.getSequencePosition(23 + 7));
		assertEquals(26, seq.getSequencePosition(26 + 8));
		assertEquals(30, seq.getSequencePosition(30 + 11));
		assertEquals(35, seq.getSequencePosition(35 + 13));
	}
}

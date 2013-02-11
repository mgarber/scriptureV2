package broad.core.multiplealignment;

import java.util.List;

import broad.core.multiplealignment.MultipleAlignment;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;

import junit.framework.TestCase;

public class MultipleAlignmentTest extends TestCase {
	
	
	
	public void testGetUngappedReferenceIslandsNoEncoding() {
		MultipleAlignment ma = getTestMultipleAlignment();
		
		List<int [] > ungappedIslands = ma.getUngappedReferenceIslands();
		assertEquals(2, ungappedIslands.size());
		assertEquals("First Island start coords are wrong",0,ungappedIslands.get(0)[0]);
		assertEquals("First Island end coords are wrong",8,ungappedIslands.get(0)[1]);
		
		assertEquals("Second Island start coords are wrong",11,ungappedIslands.get(1)[0]);
		assertEquals("Second Island end coords are wrong",19,ungappedIslands.get(1)[1]);

	}
	
	public void testGetUngappedReferenceIslandsWithMatrixEncoding() {
		MultipleAlignment ma = getTestMultipleAlignment();
		
		List<int [] > ungappedIslands = ma.getUngappedReferenceIslands();
		ma.encodeAsMatrix();
		assertEquals(2, ungappedIslands.size());
		assertEquals("First Island start coords are wrong",0,ungappedIslands.get(0)[0]);
		assertEquals("First Island end coords are wrong",8,ungappedIslands.get(0)[1]);
		
		assertEquals("Second Island start coords are wrong",11,ungappedIslands.get(1)[0]);
		assertEquals("Second Island end coords are wrong",19,ungappedIslands.get(1)[1]);

	}
	
	public void testGetUnNedReferenceIslandsWithMatrixEncoding() {
		MultipleAlignment ma = getTestMultipleAlignment();
		AlignedSequence ref = ma.getAlignedSequence(ma.getReferenceId());
		String refBases = ref.getSequenceBases();
		ref.setSequenceBases(refBases.replaceAll("-", "N"));
		List<int [] > ungappedIslands = ma.getUngappedReferenceIslands();


		assertEquals(1, ungappedIslands.size());
		assertEquals("First Island start coords are wrong",0,ungappedIslands.get(0)[0]);
		assertEquals("First Island end coords are wrong",19,ungappedIslands.get(0)[1]);
		
		ma.encodeAsMatrix();
		ungappedIslands = ma.getUngappedReferenceIslands();
		assertEquals(2, ungappedIslands.size());
		assertEquals("First Island start coords are wrong",0,ungappedIslands.get(0)[0]);
		assertEquals("First Island end coords are wrong",8,ungappedIslands.get(0)[1]);
		
		assertEquals("Second Island start coords are wrong",11,ungappedIslands.get(1)[0]);
		assertEquals("Second Island end coords are wrong",19,ungappedIslands.get(1)[1]);
	}
	
	public void testAppend() {
		MultipleAlignment ma = new MultipleAlignment();
		
		MultipleAlignment a1 = new MultipleAlignment();
		AlignedSequence r1 = new AlignedSequence("H");
		AlignedSequence s1 = new AlignedSequence("S1");
		AlignedSequence s2 = new AlignedSequence("S2");
		
		r1.setSequenceBases("A");
		s1.setSequenceBases("C");
		s2.setSequenceBases("C");
		
		a1.addSequence(r1);
		a1.addSequence(s1);
		a1.addSequence(s2);
		assertEquals("A",a1.getAlignedSequence("H").getSequenceBases());
		assertEquals("AA",a1.getAlignedSequence("H").getSequenceBases() + a1.getAlignedSequence("H").getSequenceBases());
		assertEquals("C",a1.getAlignedSequence("S1").getSequenceBases());
		assertEquals("C",a1.getAlignedSequence("S2").getSequenceBases());
		
		ma.append(a1);
		
		assertTrue(!ma.isEmpty());
		assertEquals(a1.getAlignedSequence("H").getSequenceBases(), ma.getAlignedSequence("H").getSequenceBases());
		assertEquals(a1.getAlignedSequence("S1").getSequenceBases(), ma.getAlignedSequence("S1").getSequenceBases());
		assertEquals(a1.getAlignedSequence("S2").getSequenceBases(), ma.getAlignedSequence("S2").getSequenceBases());
		
		ma.append(a1);
		assertEquals("AA",a1.getAlignedSequence("H").getSequenceBases() + a1.getAlignedSequence("H").getSequenceBases());
		assertEquals("AA",ma.getAlignedSequence("H").getSequenceBases());
		assertEquals(a1.getAlignedSequence("H").getSequenceBases() + a1.getAlignedSequence("H").getSequenceBases(), ma.getAlignedSequence("H").getSequenceBases());
		assertEquals(a1.getAlignedSequence("S1").getSequenceBases() + a1.getAlignedSequence("S1").getSequenceBases(), ma.getAlignedSequence("S1").getSequenceBases());
		assertEquals(a1.getAlignedSequence("S2").getSequenceBases() + a1.getAlignedSequence("S2").getSequenceBases()  , ma.getAlignedSequence("S2").getSequenceBases());		
	}
	
	
	private MultipleAlignment getTestMultipleAlignment() {
		MultipleAlignment ma = new MultipleAlignment();
		ma.setReferenceId("H");
		AlignedSequence ref = new AlignedSequence("H");
		AlignedSequence s1 = new AlignedSequence("S1");
		AlignedSequence s2 = new AlignedSequence("S2");
		
		ref.setSequenceBases("ACCTGCCT---AGTGGTAA");
		s1.setSequenceBases( "AGCTCACTA-AACAGGTAA");
		s2.setSequenceBases( "AGC-AACT--AACAGG-AA");
		
		ma.addSequence(ref);
		ma.addSequence(s1);
		ma.addSequence(s2);
		return ma;
	}

}

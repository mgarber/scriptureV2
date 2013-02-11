package broad.core.sequence;

import java.util.Iterator;
import java.util.List;


import broad.core.annotation.GenomicAnnotation;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;

import junit.framework.TestCase;

public class SequenceTest extends TestCase {
	
	public void testSoftmaskingIdentification() throws Exception {
		Sequence seq = new Sequence("test1");
		seq.setSequenceBases("ACTTTAAACccccAAATCTactgACCTTGGGGGAATCCGTATGAGGGTaaatttccAgAaCcGggaga");
		
		List<GenomicAnnotation> maskedRegs = seq.getSoftmaskedRegions();
		
		int i = 0;
		
		assertEquals("Unexpected number of masked regions",7, maskedRegs.size());
		assertEquals("Bad start of region " + (i + 1),10,maskedRegs.get(i).getStart());
		assertEquals("Bad start of region one",14,maskedRegs.get(i++).getEnd());
		
		assertEquals("Bad start of region " + (i + 1),20,maskedRegs.get(i).getStart());
		assertEquals("Bad start of region " + + (i + 1),24,maskedRegs.get(i++).getEnd());

	}
	
	public void testMotifFinding() throws Exception {
		Sequence seq = new Sequence("test");
		seq.setSequenceBases("ATGTATGTGTGAAAAGTAAATGAAAGTG");
		SequenceMotif sm = new SequenceMotif("GT-TG",1);
		
		Iterator<SequenceRegion> matchIt = sm.match(seq).iterator();
		
		while(matchIt.hasNext()) {
			System.out.println(matchIt.next().getSequenceBases());
		}
		
		
		
	}

}

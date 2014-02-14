package nextgen.editing.crispr.predicate;

import nextgen.editing.crispr.*;
import nextgen.core.capture.filter.PolyBaseFilter;

import org.apache.commons.lang3.StringUtils;

/**
 * @author engreitz
 * Guide RNAs run into trouble when they contain too many Us because a polyU signal is the termination
 * signal for Pol III.  This class implements semi ad hoc checks to eliminate guides with high U content:
 * 		- Check whether a guide RNA contains >N U bases in the first M bases of the seed region (by the NGG) ... 
 * 		  these may synergize with the UUUU sequence in the scaffold to halt transcription
 * 		- Allow no more than 3 Us in a row anywhere in the sequence
 * 		- Total U content of guide sequence cannot be more than 40%
 */
public class UContentPredicate implements GuideRNAPredicate {
	
	private int maxU = 1;
	private int seedBases = 4;
	public String name = "UContent";
	
	public UContentPredicate() {}
	
	public UContentPredicate(int maxU, int seedBases) {
		this.maxU = maxU;
		this.seedBases = seedBases;
	}
	
	@Override
	public boolean evaluate(GuideRNA g) {
		String seq = g.getSequenceString();
		boolean seedFilter = StringUtils.countMatches(seq.subSequence(seq.length()-seedBases, seq.length()),"T") <= maxU;
		boolean percentFilter = StringUtils.countMatches(seq, "T") < (0.4 * (double) seq.length());
		boolean repeatFilter = new PolyBaseFilter("T",4,4).evaluate(g); 
		return (seedFilter & percentFilter & repeatFilter);
	}
	
	@Override
	public String getPredicateName() {
		return name;
	}
	
	@Override
	public String getShortFailureMessage(GuideRNA g) {
		return "fails_UContent_filter";
	}

}

package nextgen.editing.crispr.predicate;

import nextgen.editing.crispr.*;

import org.apache.commons.lang3.StringUtils;

/**
 * @author engreitz
 * Check whether a guide RNA contains >N U bases in the first M bases of the seed region (by the NGG)
 * (poor man's heuristic of GuideSufficientEfficacy)
 */
public class GuideRNASeedU implements GuideRNAPredicate {
	
	private int maxU = 1;
	private int seedBases = 4;
	public String name = "SeedUCount";
	
	public GuideRNASeedU() {}
	
	public GuideRNASeedU(int maxU, int seedBases) {
		this.maxU = maxU;
		this.seedBases = seedBases;
	}
	
	@Override
	public boolean evaluate(GuideRNA g) {
		String seq = g.getSequenceString();
		return (StringUtils.countMatches(seq.subSequence(seq.length()-seedBases, seq.length()),"T") <= maxU);
	}
	
	@Override
	public String getPredicateName() {
		return name;
	}
	
	@Override
	public String getShortFailureMessage(GuideRNA g) {
		return name + "_more_than_" + maxU + "_in_seed_" + seedBases;
	}

}

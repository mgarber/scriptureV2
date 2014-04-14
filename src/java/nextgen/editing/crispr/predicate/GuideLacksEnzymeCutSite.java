package nextgen.editing.crispr.predicate;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import nextgen.editing.RestrictionEnzyme;
import nextgen.editing.crispr.GuideRNA;

public class GuideLacksEnzymeCutSite implements GuideRNAPredicate {
	
	private Collection<RestrictionEnzyme> enzymes;
	
	public GuideLacksEnzymeCutSite(RestrictionEnzyme enzyme) {
		enzymes = new ArrayList<RestrictionEnzyme>();
		enzymes.add(enzyme);
	}
	
	public GuideLacksEnzymeCutSite(Collection<RestrictionEnzyme> restrictionEnzymes) {
		enzymes = restrictionEnzymes;
	}
	
	@Override
	public boolean evaluate(GuideRNA guideRNA) {
		String guideSeq = guideRNA.getSequenceString();
		for(RestrictionEnzyme enzyme : enzymes) {
			if(enzyme.sequenceContainsTopStrandRecognitionSequence(guideSeq)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public String getPredicateName() {
		return "guide_lacks_enzyme_cut_site";
	}

	@Override
	public String getShortFailureMessage(GuideRNA g) throws IOException, InterruptedException {
		return "guide_contains_enzyme_cut_site";
	}

}

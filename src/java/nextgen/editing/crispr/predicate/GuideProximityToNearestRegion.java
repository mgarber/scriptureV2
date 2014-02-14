package nextgen.editing.crispr.predicate;

import java.util.Collection;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.editing.crispr.GuideRNA;


/**
 * Check whether a guide RNA is within a certain max distance of the nearest region in a collection
 * Does not count regions that overlap the guide RNA
 * @author prussell
 */
public class GuideProximityToNearestRegion implements GuideRNAPredicate {

	private Collection<Annotation> regionSet;
	private int maxDistance;
	private String name;
	
	/**
	 * @param regions The regions to consider for distance to the guide RNA
	 * @param maxDist Max acceptable distance
	 * @param predicateName Name of this predicate
	 */
	public GuideProximityToNearestRegion(Collection<Annotation> regions, int maxDist, String predicateName) {
		regionSet = regions;
		maxDistance = maxDist;
		name = predicateName;
	}
	
	@Override
	public boolean evaluate(GuideRNA g) {
		Gene asGene = new Gene(g);
		int nearest = asGene.distanceToNearestNonOverlapper(regionSet);
		return nearest <= maxDistance;
	}

	@Override
	public String getPredicateName() {
		return name;
	}

	@Override
	public String getShortFailureMessage(GuideRNA g) {
		return(name + "_not_within_" + maxDistance);
	}
}

package nextgen.editing.crispr;

import java.util.Collection;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.apache.commons.collections15.Predicate;

/**
 * Check whether a guide RNA is within a certain max distance of the nearest region in a collection
 * Does not count regions that overlap the guide RNA
 * @author prussell
 */
public class ProximityToNearestRegion implements Predicate<GuideRNA> {

	private Collection<Annotation> regionSet;
	private int maxDistance;
	
	/**
	 * @param regions The regions to consider for distance to the guide RNA
	 * @param maxDist Max acceptable distance
	 */
	public ProximityToNearestRegion(Collection<Annotation> regions, int maxDist) {
		regionSet = regions;
		maxDistance = maxDist;
	}
	
	@Override
	public boolean evaluate(GuideRNA g) {
		Gene asGene = new Gene(g.getAnnotation());
		int nearest = asGene.distanceToNearestNonOverlapper(regionSet);
		return nearest <= maxDistance;
	}

}

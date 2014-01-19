package nextgen.editing;

import java.util.ArrayList;
import java.util.Collection;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;

import broad.core.sequence.Sequence;

/**
 * A pair of restriction enzyme cut sites, one on left and one on right
 * @author prussell
 */
public class RestrictionEnzymeCutSitePair {
	
	private RestrictionEnzymeCutSite leftSite;
	private RestrictionEnzymeCutSite rightSite;
	
	/**
	 * @param left Left cut site
	 * @param right Right cut site
	 */
	public RestrictionEnzymeCutSitePair(RestrictionEnzymeCutSite left, RestrictionEnzymeCutSite right) {
		if(getInnerDistance(left, right) < 0) {
			throw new IllegalArgumentException("Left cut site must be to the left of right cut site. Left: " + leftSite.toString() + " right: " + rightSite.toString());
		}
		leftSite = left;
		rightSite = right;
	}
	
	/**
	 * Get all pairs of sites in a window with a certain enzyme on the left and another enzyme on the right
	 * @param leftEnzyme Left enzyme
	 * @param rightEnzyme Right enzyme
	 * @param chr Chromosome
	 * @param windowStart Start position of window to search
	 * @param windowEnd End position of window to search
	 * @param maxInnerDist Max inner distance between two cut sites
	 * @return All site pairs in the window meeting the criteria
	 */
	public static Collection<RestrictionEnzymeCutSitePair> getAllSitePairs(RestrictionEnzyme leftEnzyme, RestrictionEnzyme rightEnzyme, Sequence chr, int windowStart, int windowEnd, int maxInnerDist) {
		Collection<RestrictionEnzymeCutSitePair> rtrn = new ArrayList<RestrictionEnzymeCutSitePair>();
		Collection<RestrictionEnzymeCutSite> leftSites = RestrictionEnzymeCutSite.getAllSites(leftEnzyme, chr, windowStart, windowEnd);
		Collection<RestrictionEnzymeCutSite> rightSites = RestrictionEnzymeCutSite.getAllSites(rightEnzyme, chr, windowStart, windowEnd);
		for(RestrictionEnzymeCutSite left : leftSites) {
			for(RestrictionEnzymeCutSite right : rightSites) {
				if(validatePair(left, right, maxInnerDist)) {
					rtrn.add(new RestrictionEnzymeCutSitePair(left, right));
				}
			}
		}
		return rtrn;
	}
	
	/**
	 * Check whether the pair is laid out correctly and within a specified distance of each other
	 * @param left Left cut site
	 * @param right Right cut site
	 * @param maxInnerDist Max inner distance between two cut sites
	 * @return Whether the pair is valid
	 */
	public static boolean validatePair(RestrictionEnzymeCutSite left, RestrictionEnzymeCutSite right, int maxInnerDist) {
		int innerDist = getInnerDistance(left, right);
		if(innerDist < 0) return false;
		if(innerDist > maxInnerDist) return false;
		return true;
	}
	
	/**
	 * @param left Left cut site
	 * @param right Right cut site
	 * @return Inner distance between the two cut sites
	 */
	public static int getInnerDistance(RestrictionEnzymeCutSite left, RestrictionEnzymeCutSite right) {
		return right.getStart() - left.getEnd();
	}
	
	/**
	 * @return Inner distance between the two cut sites
	 */
	public int getInnerDistance() {
		return rightSite.getStart() - leftSite.getEnd();
	}
	
	/**
	 * @return Bed line for the blocked annotation representation of the pair of sites
	 */
	public String toBED() {
		return asAnnotation().toBED();
	}
	
	/**
	 * Get annotation representations of several cut site pairs
	 * @param pairs The pairs to get
	 * @return Blocked annotations
	 */
	public static Collection<Annotation> asAnnotations(Collection<RestrictionEnzymeCutSitePair> pairs) {
		Collection<Annotation> rtrn = new ArrayList<Annotation>();
		for(RestrictionEnzymeCutSitePair pair : pairs) {
			rtrn.add(pair.asAnnotation());
		}
		return rtrn;
	}
	
	/**
	 * @return Annotation representation with the two sites as blocks
	 */
	public Annotation asAnnotation() {
		BasicAnnotation rtrn = new BasicAnnotation(leftSite.getSite());
		rtrn.addBlocks(rightSite.getSite());
		rtrn.setName(toString());
		return rtrn;
	}
	
	public String toString() {
		return "left_" + leftSite.toString() + "_right_" + rightSite.toString();
	}
	
}

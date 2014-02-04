package nextgen.editing;

import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;

import org.apache.log4j.Logger;

import broad.core.sequence.Sequence;

/**
 * A representation of a cut site for a particular restriction enzyme
 * The genomic sequence must contain the recognition sequence at the specified site
 * @author prussell
 */
public class RestrictionEnzymeCutSite {
	
	public static Logger logger = Logger.getLogger(RestrictionEnzymeCutSite.class.getName());
	private Sequence chromosome;
	private Annotation location;
	private RestrictionEnzyme enzyme;
	
	/**
	 * @param restrictionEnzyme Restriction enzyme
	 * @param chr Chromosome
	 * @param recognitionSite Annotation representing the recognition sequence itself on the chromosome
	 */
	public RestrictionEnzymeCutSite(RestrictionEnzyme restrictionEnzyme, Sequence chr, Annotation recognitionSite) {
		String annotSeq = chr.getSubsequence(recognitionSite).getSequenceBases();
		if(!restrictionEnzyme.hasTopStrandRecognitionSequence(annotSeq)) {
			throw new IllegalArgumentException(restrictionEnzyme.getName() + " recognition sequence is not present at " + recognitionSite.toUCSC() + ":" + recognitionSite.getOrientation().toString() + "(" + annotSeq + ")");
		}
		enzyme = restrictionEnzyme;
		chromosome = chr;
		location = recognitionSite;
	}
	
	/**
	 * Get all the cut sites for a restriction enzyme in a specified genomic window
	 * @param restrictionEnzyme The enzyme
	 * @param chr The chromosome sequence
	 * @param windowStart Window start
	 * @param windowEnd Window end
	 * @return All cut sites fully contained in the window
	 */
	public static Collection<RestrictionEnzymeCutSite> getAllSites(RestrictionEnzyme restrictionEnzyme, Sequence chr, int windowStart, int windowEnd) {
		Collection<Annotation> locations = new TreeSet<Annotation>();
		for(String s : restrictionEnzyme.getTopStrandRecognitionSequence()) {
			locations.addAll(chr.getMatches(s, windowStart, windowEnd, true));
		}
		Collection<RestrictionEnzymeCutSite> rtrn = new ArrayList<RestrictionEnzymeCutSite>();
		for(Annotation a : locations) {
			RestrictionEnzymeCutSite cutSite = new RestrictionEnzymeCutSite(restrictionEnzyme, chr, a);
			rtrn.add(cutSite);
		}
		return rtrn;
	}
	
	public Collection<String> getRecognitionSequence() {
		return enzyme.getTopStrandRecognitionSequence();
	}
	
	public String getEnzymeName() {
		return enzyme.getName();
	}
	
	public String getChrName() {
		return chromosome.getId();
	}
	
	public Sequence getChr() {
		return chromosome;
	}
	
	public int getStart() {
		return location.getStart();
	}
	
	public int getEnd() {
		return location.getEnd();
	}
	
	/**
	 * @return The cut site as a genomic annotation
	 */
	public Annotation getSite() {
		return location;
	}
	
	/**
	 * Get several cut sties as annotations
	 * @param sites The sites
	 * @return Annotation objects representing the cut sites
	 */
	public static Collection<Annotation> asAnnotations(Collection<RestrictionEnzymeCutSite> sites) {
		Collection<Annotation> rtrn = new ArrayList<Annotation>();
		for(RestrictionEnzymeCutSite site : sites) {
			rtrn.add(site.getSite());
		}
		return rtrn;
	}
	
	/**
	 * @return A bed line for the location of the cut site
	 */
	public String toBED() {
		Annotation a = location.copy();
		a.setName(toString());
		return a.toBED();
	}
	
	public String toString() {
		return getEnzymeName() + ":" + location.toUCSC();
	}
	
}

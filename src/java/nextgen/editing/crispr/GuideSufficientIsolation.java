package nextgen.editing.crispr;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

/**
 * Check whether a guide RNA is sufficiently isolated from a collection of regions
 * @author prussell
 */
public class GuideSufficientIsolation implements GuideRNAPredicate {
	
	private Map<String, Collection<Gene>> genes;
	private int minDistance;
	public static Logger logger = Logger.getLogger(GuideSufficientIsolation.class.getName());
	private boolean sameStrand;
	private String name;
	
	/**
	 * @param geneAnnotation Full genome annotation
	 * @param minDistNearestGene Min required distance to nearest neighboring gene
	 * @param sameStrandOnly Only consider genes on the same strand as the target gene of the guide pair
	 * @param predicateName Name of this predicate
	 */
	public GuideSufficientIsolation(Map<String, Collection<Gene>> geneAnnotation, int minDistNearestGene, boolean sameStrandOnly, String predicateName) {
		this(geneAnnotation, makeCollection(null), minDistNearestGene, sameStrandOnly, predicateName);
	}
	
	/**
	 * @param geneAnnotation Full genome annotation
	 * @param ignoreOverlappersOfTheseGenes Any gene in the gene annotation that overlaps this target gene in any orientation will be ignored for purposes of nearest neighbor
	 * @param minDistNearestGene Min required distance to nearest neighboring gene
	 * @param sameStrandOnly Only consider genes on the same strand as the target gene of the guide pair
	 * @param predicateName Name of this predicate
	 */
	public GuideSufficientIsolation(Map<String, Collection<Gene>> geneAnnotation, Gene targetGeneToIgnoreOverlappers, int minDistNearestGene, boolean sameStrandOnly, String predicateName) {
		this(geneAnnotation, makeCollection(targetGeneToIgnoreOverlappers), minDistNearestGene, sameStrandOnly, predicateName);
	}
	
	/**
	 * @param geneAnnotation Full genome annotation
	 * @param ignoreOverlappersOfTheseGenes Any gene in the gene annotation that overlaps one of these in any orientation will be ignored for purposes of nearest neighbor
	 * @param minDistNearestGene Min required distance to nearest neighboring gene
	 * @param sameStrandOnly Only consider genes on the same strand as the target gene of the guide pair
	 * @param predicateName Name of this predicate
	 */
	private GuideSufficientIsolation(Map<String, Collection<Gene>> geneAnnotation, Collection<Gene> ignoreOverlappersOfTheseGenes, int minDistNearestGene, boolean sameStrandOnly, String predicateName) {
		genes = geneAnnotation;
		name = predicateName;
		sameStrand = sameStrandOnly;
		Map<String, Collection<Gene>> toRemove = new TreeMap<String, Collection<Gene>>();
		for(Gene ignore : ignoreOverlappersOfTheseGenes) {
			String chr = ignore.getChr();
			for(Gene gene : genes.get(chr)) {
				if(gene.overlaps(ignore, true)) {
					if(!toRemove.containsKey(chr)) {
						toRemove.put(chr, new TreeSet<Gene>());
					}
					toRemove.get(chr).add(gene);
				}
			}
		}
		for(String chr : toRemove.keySet()) {
			for(Gene gene : toRemove.get(chr)) {
				genes.get(chr).remove(gene);
			}
		}
		minDistance = minDistNearestGene;
	}
	
	private static Collection<Gene> makeCollection(Gene gene) {
		Collection<Gene> rtrn = new TreeSet<Gene>();
		if(gene != null) rtrn.add(gene);
		return rtrn;
	}
	
	@Override
	public boolean evaluate(GuideRNA guideRNA) {
		String chr = guideRNA.getChr();
		Collection<Gene> genesThisChr = genes.get(chr);
		int start = guideRNA.getStart() - minDistance;
		int end = guideRNA.getEnd() + minDistance;
		Strand strand = guideRNA.getTargetGene().getOrientation();
		Collection<Gene> overlappers = Gene.getOverlappers(genesThisChr, chr, start, end, strand, !sameStrand);
		if(!overlappers.isEmpty()) {
			for(Gene overlapper : overlappers) {
				logger.debug("NOT_SUFFICIENTLY_ISOLATED\t" + guideRNA.getName() + " is within " + minDistance + " of gene " + overlapper.getName());
			}
			return false;
		}
		logger.debug("SUFFICIENTLY_ISOLATED\t" + guideRNA.getName());
		return true;
	}

	@Override
	public String getPredicateName() {
		return name;
	}

	@Override
	public String getShortFailureMessage() {
		return(name + "_within_" + minDistance + "_of_nearest_gene");
	}

}

package nextgen.editing.crispr;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;

/**
 * Check whether a guide RNA is sufficiently isolated from a collection of regions
 * @author prussell
 */
public class SufficientIsolation implements Predicate<GuideRNA> {
	
	private Map<String, Collection<Gene>> genes;
	private int minDistance;
	public static Logger logger = Logger.getLogger(SufficientIsolation.class.getName());
	
	/**
	 * @param geneAnnotation Full genome annotation
	 * @param minDistNearestGene Min required distance to nearest neighboring gene
	 */
	public SufficientIsolation(Map<String, Collection<Gene>> geneAnnotation, int minDistNearestGene) {
		this(geneAnnotation, makeCollection(null), minDistNearestGene);
	}
	
	/**
	 * @param geneAnnotation Full genome annotation
	 * @param ignoreOverlappersOfTheseGenes Any gene in the gene annotation that overlaps this target gene in any orientation will be ignored for purposes of nearest neighbor
	 * @param minDistNearestGene Min required distance to nearest neighboring gene
	 */
	public SufficientIsolation(Map<String, Collection<Gene>> geneAnnotation, Gene targetGeneToIgnoreOverlappers, int minDistNearestGene) {
		this(geneAnnotation, makeCollection(targetGeneToIgnoreOverlappers), minDistNearestGene);
	}
	
	/**
	 * @param geneAnnotation Full genome annotation
	 * @param ignoreOverlappersOfTheseGenes Any gene in the gene annotation that overlaps one of these in any orientation will be ignored for purposes of nearest neighbor
	 * @param minDistNearestGene Min required distance to nearest neighboring gene
	 */
	private SufficientIsolation(Map<String, Collection<Gene>> geneAnnotation, Collection<Gene> ignoreOverlappersOfTheseGenes, int minDistNearestGene) {
		genes = geneAnnotation;
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
		Collection<Gene> overlappers = Gene.getOverlappers(genesThisChr, chr, start, end);
		if(!overlappers.isEmpty()) {
			for(Gene overlapper : overlappers) {
				logger.debug("NOT_SUFFICIENTLY_ISOLATED\t" + guideRNA.getName() + " is within " + minDistance + " of gene " + overlapper.getName());
			}
			return false;
		}
		logger.debug("SUFFICIENTLY_ISOLATED\t" + guideRNA.getName());
		return true;
	}

}

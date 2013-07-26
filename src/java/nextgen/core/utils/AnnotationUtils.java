package nextgen.core.utils;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;


import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class AnnotationUtils {

	public static Logger logger = Logger.getLogger(AnnotationUtils.class.getName());
	
	public static Map<String, Collection<Gene>> collapseOverlappers(Map<String, Collection<Gene>> genes) {		
		return collapseOverlappers(genes, false);
	}
	
	public static Map<String, Collection<Gene>> collapseOverlappers(Map<String, Collection<Gene>> genes, boolean ignoreOrientation) {		
		Map<String, Collection<Gene>> rtrn = new TreeMap<String, Collection<Gene>>();
		for(String chr : genes.keySet()) {
			Collection<Gene> genesThisChr = genes.get(chr);
			Collection<Gene> collapsedThisChr = new TreeSet<Gene>();
			for(Gene gene : genesThisChr) {
				gene.setOrientation(Strand.UNKNOWN);
				if(collapsedThisChr.isEmpty()) {
					collapsedThisChr.add(gene);
					continue;
				}
				boolean added = false;
				for(Gene collapsed : collapsedThisChr) {
					if(gene.overlaps(collapsed)) {
						collapsed.addBlocks(gene);
						collapsed.setCDSStart(Math.min(collapsed.getCDSStart(), gene.getCDSStart()));
						collapsed.setCDSEnd(Math.max(collapsed.getCDSEnd(), gene.getCDSEnd()));
						collapsed.setName(collapsed.getName() + "_" + gene.getName());
						added = true;
						break;
					} 
				}
				if(!added) {
					collapsedThisChr.add(gene);
				}
			}
			rtrn.put(chr, collapsedThisChr);
		}
		return rtrn;
	}

	
	/**
	 * Merge annotations that overlap in the same orientation as each other, and leave singleton windows the same
	 * @param windows Annotation set
	 * @return New set of merged windows
	 */
	public static Collection<Annotation> mergeOverlappingBlocks(TreeSet<Annotation> windows) {
		TreeSet<Annotation> plusStrand = new TreeSet<Annotation>();
		TreeSet<Annotation> minusStrand = new TreeSet<Annotation>();
		for(Annotation window : windows) {
			Strand orientation = window.getOrientation();
			if(orientation.equals(Strand.POSITIVE)) {
				plusStrand.add(window);
				continue;
			}
			if(orientation.equals(Strand.NEGATIVE)) {
				minusStrand.add(window);
				continue;
			}
			throw new IllegalArgumentException("Annotation must have an orientation: " + window.getName());
		}
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		rtrn.addAll(mergeOverlappingBlocksSameOrientation(plusStrand));
		rtrn.addAll(mergeOverlappingBlocksSameOrientation(minusStrand));
		return rtrn;
	}
	
	/**
	 * Merge annotations that overlap, and leave singleton windows the same
	 * Throws exception if any annotations have different orientations
	 * @param windows Annotation set
	 * @return New set of merged windows
	 */
	private static Collection<Annotation> mergeOverlappingBlocksSameOrientation(TreeSet<Annotation> windows) {
		
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		
		// If provided set of regions contains zero or one element, return a copy of the set itself
		if(windows.size() == 0) return rtrn;
		if(windows.size() == 1) {
			rtrn.addAll(windows);
			return rtrn;
		}
		
		Iterator<? extends Annotation> iter = windows.iterator();
		// Begin with the first window in the sorted set
		Annotation firstWindow = iter.next();
		Strand orientation = firstWindow.getOrientation();
		// A changing exon set that expands when the next window overlaps it
		// Reset when the next window does not overlap
		Collection<Annotation> growingExonSet = new TreeSet<Annotation>();
		growingExonSet.addAll(firstWindow.getBlocks());
		
		while(iter.hasNext()) {
			
			Annotation growingWindow = new Gene(growingExonSet);
			Annotation nextWindow = iter.next();
			if(!nextWindow.getOrientation().equals(orientation)) {
				throw new IllegalArgumentException("All annotations must have same orientation.");
			}
			
			if(nextWindow.overlaps(growingWindow)) {
				// Next window span overlaps growing exon set genomic span
				// Merge the exon sets
				growingExonSet.clear();
				growingExonSet.addAll(growingWindow.union(nextWindow).getBlocks());
				
				if(!iter.hasNext()) {
					// This is the last window in the set
					// Add the last element and leave loop
					Annotation lastAdd = new Gene(growingExonSet);
					rtrn.add(lastAdd);
					continue;					
				}
			} else {
				// Next window does not overlap growing exon set
				// Make window from latest version of growing exon set and add to the return set
				rtrn.add(growingWindow);
				if(!iter.hasNext()) {
					// This is the last window in the set
					// Add it and leave loop
					rtrn.add(nextWindow);
					continue;
				}
				// Reset growing exon set to this new gene
				growingExonSet.clear();
				growingExonSet.addAll(nextWindow.getBlocks());
				continue;
			}
		}
		
		return rtrn;

	}
	
	/**
	 * Get largest parent
	 * @param child Child annotation
	 * @param parents Possible parent annotations by chromosome
	 * @return Largest parent (annotation containing child on same strand) or null if none exists
	 */
	public static <T extends Annotation> T getLargestParent(T child, Map<String, Collection<T>> parents) {
		Collection<T> allParents = getParents(child, parents);
		T rtrn = null;
		int largestParentSize = 0;
		for(T parent : allParents) {
			if(parent.getSize() > largestParentSize) {
				rtrn = parent;
				largestParentSize = parent.getSize();
			}
		}		
		return rtrn;
	}
	
	/**
	 * Map each child annotation to its largest parent defined as annotations overlapping and containing child
	 * @param children Child annotations by chromosome
	 * @param parents Parent annotations by chromosome
	 * @return Map of child to largest parent. Does not include children with no parents.
	 */
	public static <T extends Annotation> Map<T, T> mapChildToLargestParent(Map<String, Collection<T>> children, Map<String, Collection<T>> parents) {
		Map<T, T> rtrn = new TreeMap<T, T>();
		logger.info("Finding largest parent for each child annotation...");
		int numDone = 0;
		for(String chr : children.keySet()) {
			logger.info(chr);
			for(T child : children.get(chr)) {
				numDone++;
				if(numDone % 1000 == 0) {
					logger.info("Finished " + numDone + " regions.");
				}
				T largestParent = getLargestParent(child, parents);
				if(largestParent != null) {
					rtrn.put(child, largestParent);
				}
			}
		}
		logger.info("Done finding largest parent for each annotation.");
		return rtrn;
	}
	
	/**
	 * Get parents of child annotation
	 * @param child Child
	 * @param parents Parents by chromosome
	 * @return Set of parents containing child or empty set if none
	 */
	public static <T extends Annotation> Collection<T> getParents(T child, Map<String, Collection<T>> parents) {
		Collection<T> rtrn = new TreeSet<T>();
		String chr = child.getChr();
		if(!parents.containsKey(chr)) {
			return rtrn;
		}
		for(T other : parents.get(chr)) {
			if(other.getEnd() < child.getStart() || other.getStart() > child.getEnd()) {
				continue;
			}
			if(other.contains(child) && other.overlaps(child) && child.getOrientation().equals(other.getOrientation())) {
				rtrn.add(other);
			}
		}
		return rtrn;
	}
	
	
	/**
	 * Get children of parent annotation
	 * @param parent Parent
	 * @param children Potential children by chromosome
	 * @return Set of children contained in parent or empty set if none
	 */
	public static <T extends Annotation> Collection<T> getChildren(T parent, Map<String, Collection<T>> children) {
		Collection<T> rtrn = new TreeSet<T>();
		String chr = parent.getChr();
		if(!children.containsKey(chr)) {
			logger.debug("Child set does not contain key " + chr);
			return rtrn;
		}
		for(T other : children.get(chr)) {
			if(other.getEnd() < parent.getStart() || other.getStart() > parent.getEnd()) {
				//logger.debug(other.toUCSC() + " does not overlap " + parent.toUCSC());
				continue;
			}
			boolean b1 = parent.contains(other);
			boolean b2 = other.overlaps(parent);
			boolean b3 = parent.getOrientation().equals(other.getOrientation());
			logger.debug(parent.toUCSC() + "\t" + other.toUCSC() + "\t" + b1 + "\t" + b2 + "\t" + b3);
			if(b1 && b2 && b3) {
				rtrn.add(other);
			}
		}
		logger.debug(parent.toUCSC() + " " + rtrn.size() + " overlappers");
		return rtrn;
	}
	
	
	
	
	/**
	 * Map each parent annotation to its set of contained children
	 * @param <T>
	 * @param children Child annotations by chromosome
	 * @param parents Parent annotations by chromosome
	 * @return Map of parent to set of children. Annotations with no children map to empty set.
	 */
	public static <T extends Annotation> Map<T, Collection<T>> mapParentsToChildren(Map<String, Collection<T>> children, Map<String, Collection<T>> parents) {
		logger.info("Mapping parent annotations to children...");
		Map<T, Collection<T>> rtrn = new TreeMap<T, Collection<T>>();
		for(String chr : parents.keySet()) {
			logger.info(chr);
			if(!children.containsKey(chr)) {
				continue;
			}
			for(T parent : parents.get(chr)) {
				Collection<T> childSet = getChildren(parent, children);
				logger.debug(parent.toUCSC() + " " + childSet.size() + " overlappers");
				rtrn.put(parent, childSet);
			}
		}
		logger.info("Done mapping parent annotations to child annotations.");
		return rtrn;
	}
	
	
	/**
	 * Map each child annotation to its set of parents defined as annotations overlapping and containing child
	 * @param <T>
	 * @param children Child annotations by chromosome
	 * @param parents Parent annotations by chromosome
	 * @return Map of child to set of parents. Does not include children with no parents.
	 */
	public static <T extends Annotation> Map<T, Collection<T>> mapChildToParents(Map<String, Collection<T>> children, Map<String, Collection<T>> parents) {
		logger.info("Mapping child annotations to parent annotations...");
		Map<T, Collection<T>> rtrn = new TreeMap<T, Collection<T>>();
		for(String chr : children.keySet()) {
			logger.info(chr);
			if(!parents.containsKey(chr)) {
				continue;
			}
			for(T child : children.get(chr)) {
				Collection<T> parentSet = getParents(child, parents);
				if(!parentSet.isEmpty()) {
					rtrn.put(child, parentSet);
				}
			}
		}
		logger.info("Done mapping child annotations to parent annotations.");
		return rtrn;
	}
	
	/**
	 * Load annotations in the bed file as annotation objects
	 * @param bedFile Bed file
	 * @return Map of chromosome to collection of annotations
	 * @throws IOException
	 */
	public static Map<String, Collection<Annotation>> loadAnnotations(String bedFile) throws IOException {
		Map<String, Collection<Annotation>> rtrn = new TreeMap<String, Collection<Annotation>>();
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(bedFile));
		for(String chr : genes.keySet()) {
			Collection<Annotation> a = new TreeSet<Annotation>();
			a.addAll(genes.get(chr));
			rtrn.put(chr, a);
		}
		return rtrn;
	}
	
	/**
	 * Get midpoint along mature transcript of a spliced annotation
	 * @param annot Annotation
	 * @return Genome coordinate of mature transcript midpoint
	 */
	public static int getMidpoint(Annotation annot) {
		Gene gene = new Gene(annot);
		int size = gene.size();
		int midpoint = size / 2;
		return gene.transcriptToGenomicPosition(midpoint);
	}

	/**
	 * Get midpoint of a sub annotation with respect to a possibly spliced parent annotation
	 * Ignores blocks of sub annotation and only considers blocks of parent annotation
	 * @param annot Parent annotation
	 * @param subAnnot Sub annotation contained in parent annotation
	 * @return Midpoint of sub annotation span considering blocks of parent annotation
	 */
	public static int getSubAnnotationMidpointWithinAnnotation(Annotation annot, Annotation subAnnot) {
		int start = subAnnot.getStart();
		int end = subAnnot.getEnd();
		if(start < annot.getStart() || end > annot.getEnd()) {
			throw new IllegalArgumentException("Sub-annotation must fall within annotation.");
		}
		Gene trimmed = new Gene(annot);
		trimmed.setStart(start);
		trimmed.setEnd(end);
		int trimmedSize = trimmed.size();
		int midpoint = trimmedSize / 2;
		return trimmed.transcriptToGenomicPosition(midpoint);
	}


}

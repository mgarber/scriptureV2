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


import nextgen.core.annotation.AbstractAnnotation;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class AnnotationUtils {

	private static Logger logger = Logger.getLogger(AnnotationUtils.class.getName());
	
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
	 * Map each child annotation to its largest parent defined as annotations overlapping and containing child
	 * @param children Child annotations by chromosome
	 * @param parents Parent annotations by chromosome
	 * @return Map of child to largest parent. Does not include children with no parents.
	 */
	public static <T extends Annotation> Map<T, T> mapChildToLargestParent(Map<String, Collection<T>> children, Map<String, Collection<T>> parents) {
		Map<T, T> rtrn = new TreeMap<T, T>();
		Map<T, Collection<T>> allRelationships = mapChildToParents(children, parents);
		for(T child : allRelationships.keySet()) {
			int largestParentSize = 0;
			for(T parent : allRelationships.get(child)) {
				if(parent.getSize() > largestParentSize) {
					rtrn.put(child, parent);
					largestParentSize = parent.getSize();
				}
			}
		}
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
				Collection<T> parentSet = new TreeSet<T>();
				for(T other : parents.get(chr)) {
					if(other.getEnd() < child.getStart() || other.getStart() > child.getEnd()) {
						continue;
					}
					if(other.contains(child) && other.overlaps(child)) {
						parentSet.add(other);
					}
				}
				if(!parentSet.isEmpty()) {
					rtrn.put(child, parentSet);
				}
			}
		}
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



}

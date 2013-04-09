package nextgen.core.utils;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;


import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.CompoundInterval;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.SingleInterval;
import nextgen.core.feature.GeneWindow;
import nextgen.core.feature.Window;

public class AnnotationUtils {

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

}

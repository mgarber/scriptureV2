package nextgen.core.utils;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;


import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.CompoundInterval;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.SingleInterval;
import nextgen.core.feature.GeneWindow;
import nextgen.core.feature.Window;

public class AnnotationUtils {

	/**
	 * Merge block sets of windows whose genomic spans overlap in any orientation, and leave singleton windows the same
	 * @param windows
	 * @return New set of merged windows without orientation information
	 */
	public static Collection<Annotation> mergeOverlappingBlocksAnyOrientation(TreeSet<Annotation> windows) {
		
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
		// A changing exon set that expands when the next window overlaps it
		// Reset when the next window does not overlap
		Collection<Annotation> growingExonSet = new TreeSet<Annotation>();
		growingExonSet.addAll(firstWindow.getBlocks());
		
		while(iter.hasNext()) {
			
			Annotation growingWindow = new Gene(growingExonSet);
			Annotation nextWindow = iter.next();
			
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

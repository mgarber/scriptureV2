package broad.core.annotation;

import java.net.URL;
import java.util.Iterator;

import junit.framework.TestCase;

import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;

public class AnnotationReaderTest extends TestCase{
	public void testIntersect() throws Exception {
		URL fragmented = getClass().getResource("test1_ar.bed");
		URL consolidated = getClass().getResource("test1_arboot.bed");
		
		BEDReader set1 = new BEDReader(fragmented.getFile());
		BEDReader set2 = new BEDReader(consolidated.getFile());
		
		int originalSize = set1.getAnnotationList().size();
		
		IntervalTree<BED> tree = set1.getChromosomeTree("chr7");
		assertNotNull(tree);
		IntervalTree<BED> redoneTree = new IntervalTree<BED>();
		Iterator<BED> it = tree.valueIterator();
		while(it.hasNext()) {
			BED bed = it.next();
			redoneTree.put(bed.getStart(),bed.getEnd(), bed);
		}
		LightweightGenomicAnnotation firstAnnotation = set2.getAnnotationList().get(0);
		int overlappers = 0;
		Iterator<Node<BED>> overlapperIt = tree.overlappers(firstAnnotation.getStart(), firstAnnotation.getEnd());
		while(overlapperIt.hasNext()) {
			Node<BED> node = overlapperIt.next();
			overlappers++;
		}
		assertEquals(4,overlappers);
		
		set1.intersect(set2.getAnnotationList());
		
		assertEquals(originalSize, set1.getAnnotationList().size());
	}
}

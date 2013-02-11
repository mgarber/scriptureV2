package broad.core.datastructures;

import java.util.Iterator;

import broad.core.datastructures.IntervalTree;
import junit.framework.TestCase;

public class IntervalTreeTest extends TestCase {
	public void testOverlappers() throws Exception {
		IntervalTree<String> tree = new IntervalTree<String>();
		
		tree.put(115588208,115588208, "chr7:115588208-115588208");
		tree.put(115588246,115588288, "chr7:115588246-115588288");
		tree.put(115588290,115588398, "chr7:115588290-115588398");
		tree.put(115588428,115588441, "chr7:115588428-115588441");
		tree.put(115588451,115588483, "chr7:115588451-115588483");
		tree.put(115588489,115588511, "chr7:115588489-115588511");
		tree.put(115588515,115588532, "chr7:115588515-115588532");

		Iterator<String> overlaperIt = 
			new IntervalTree.ValuesIterator<String>(tree.overlappers(115588188,115588445));
		
		int overlappers = 0;
		while(overlaperIt.hasNext()) {
			String annot = overlaperIt.next();
			overlappers++;
			System.out.println("overlapper " + overlappers + " " + annot);
		}
		
		assertEquals(4,overlappers);

		
		overlaperIt = new IntervalTree.ValuesIterator<String>(tree.overlappers(115588240,115588250));
		
		overlappers = 0;
		while(overlaperIt.hasNext()) {
			String annot = overlaperIt.next();
			overlappers++;
			System.out.println("overlapper " + overlappers + " " + annot);
		}
		
		assertEquals(1,overlappers);
		
		
	}
	
	public void testOverlappersOpenClose() throws Exception {
		IntervalTree<String> tree = new IntervalTree<String>();
		
		tree.put(115588208,115588208, "chr7:115588208-115588208");
		tree.put(115588246,115588288, "chr7:115588246-115588288");
		tree.put(115588290,115588398, "chr7:115588290-115588398");
		tree.put(115588428,115588441, "chr7:115588428-115588441");
		tree.put(115588451,115588483, "chr7:115588451-115588483");
		tree.put(115588489,115588511, "chr7:115588489-115588511");
		tree.put(115588515,115588532, "chr7:115588515-115588532");
		
		
		Iterator<String> overlaperIt = new IntervalTree.ValuesIterator<String>(tree.overlappers(115588240,115588246));
		
		int overlappers = 0;
		while(overlaperIt.hasNext()) {
			String annot = overlaperIt.next();
			overlappers++;
			System.out.println("overlapper " + overlappers + " " + annot);
		}
		
		assertEquals(0,overlappers);
		
		overlaperIt = new IntervalTree.ValuesIterator<String>(tree.overlappers(115588288,115588289));
		
		overlappers = 0;
		while(overlaperIt.hasNext()) {
			String annot = overlaperIt.next();
			overlappers++;
			System.out.println("overlapper " + overlappers + " " + annot);
		}
		
		assertEquals(0,overlappers);
		
	}

}

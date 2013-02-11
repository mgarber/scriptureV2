package broad.core.annotation;

import broad.core.annotation.BasicGenomicAnnotation;
import junit.framework.TestCase;

public class TestAnnotationOperations extends TestCase {
	
	public void testGetOverlap() {
		BasicGenomicAnnotation a1 = new BasicGenomicAnnotation("a1","chr1",0,10); 
		BasicGenomicAnnotation a2 = new BasicGenomicAnnotation("a2","chr1",10,20); 
		BasicGenomicAnnotation a3 = new BasicGenomicAnnotation("a3","chr1",9,15); 
		BasicGenomicAnnotation a4 = new BasicGenomicAnnotation("a4","chr2",9,15); 
		
		assertEquals("a1 does not overlap a2",0,a1.getOverlap(a2));
		assertEquals("a1 overlaps 1 base a3",1,a1.getOverlap(a3));
		assertEquals("a2 overlaps 5 bases a3", 5, a2.getOverlap(a3));
		assertEquals("a2 overlaps 5 bases a3", 5, a3.getOverlap(a2));
		assertEquals("a4 is in a different chr", 0, a2.getOverlap(a4));
		assertEquals("a4 is in a different chr", 0, a4.getOverlap(a2));
	}

}

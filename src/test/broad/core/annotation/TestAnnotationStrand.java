package broad.core.annotation;

import junit.framework.TestCase;
import nextgen.core.annotation.BasicAnnotation;

public class TestAnnotationStrand  extends TestCase {
	public void testIntersect() throws Exception {
		BasicAnnotation negAnnot1 = new BasicAnnotation("chr1", 5, 10, "-");
		BasicAnnotation negAnnot2 = new BasicAnnotation("chr1", 5, 10, "-");
		BasicAnnotation posAnnot1 = new BasicAnnotation("chr1", 5, 10, "+");
		
		assertTrue(negAnnot1.isNegativeStrand());
		assertTrue(negAnnot1.equals(negAnnot2));
		assertFalse(negAnnot1.equals(posAnnot1));
	}
}

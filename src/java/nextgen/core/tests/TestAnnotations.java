package nextgen.core.tests;

import nextgen.core.annotation.*;
import junit.framework.TestCase;

/**
 * @author engreitz
 *
 */
public class TestAnnotations extends TestCase {

	// This is a test comment.
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		testCompoundInterval();
		testBasicAnnotation();
	}
	
	// This is a test comment made online
	
	public static void testBasicAnnotation() {
		System.out.println("testBasicAnnotation");
		Annotation a = new BasicAnnotation("chr1:500-600");
		Annotation b = new BasicAnnotation("chr1", 600, 700, "-");
		
		System.out.println(a);
		System.out.println(b);
		
		BasicAnnotation c = new BasicAnnotation(a);
		c.addBlocks(b);
		System.out.println(c);
		assert c.length() == 200;
		assert c.getBlocks().size() == 1;
		
		c.addBlocks(new BasicAnnotation("chr1", 800, 900));
		assert c.numBlocks() == 2;
		assert c.getBlocks(true).get(0).getStart() == 800;
		assert c.getBlocks(false).get(0).getEnd() == 700;
		
		assert !a.hasOrientation();
		assert b.hasOrientation();
		assert b.getStrand() == Annotation.Strand.NEGATIVE;
		
		assert c.getReferenceCoordinateAtPosition(1) == 501;
		c.setOrientation('-');
		assert c.getReferenceCoordinateAtPosition(1, true) == 501;
		assert c.getReferenceCoordinateAtPosition(1) == 898;
		assert c.getPositionAtReferenceCoordinate(898) == 1;
		assert c.getPositionAtReferenceCoordinate(501, true) == 1;
		
		c.setEnd(1000);
		System.out.println(c);
		assert c.length() == 300;
		
		c.setOrientedStart(900);
		assert c.length() == 200;
		c.setOrientedStart(700);
		assert c.numBlocks() == 1;
		assert c.length() == 200;
		c.setOrientedEnd(550);
		assert c.length() == 150;
		c.setReferenceName("chr2");
		
		Annotation d = new BasicAnnotation(c);
		assert c.equals(d);
		assert d.equals(c);
		assert c.overlaps(d);
		assert c.contains(d);
		assert c.equals(c.union(d));
		assert c.equals(c.intersect(d));
		
		
	}
	
	public static void testCompoundInterval() {
		System.out.println("testIntervals()");
		CompoundInterval a = new CompoundInterval();
		a.addInterval(new SingleInterval(101,201));
		a.addInterval(new SingleInterval(401,501));
		a.addInterval(new SingleInterval(601,701));
		
		CompoundInterval b = new CompoundInterval();
		b.addInterval(new SingleInterval(1,101));
		b.addInterval(new SingleInterval(601,701));
		
		assert a.length() == 300 : "length";
		assert a.getSpan() == 600 : "getSpan";
		assert a.numBlocks() == 3 : "numBlocks";
		assert a.getStart() == 101 : "getStart";
		assert a.getEnd() == 701 : "getEnd";
		assert a.getCoordinateAtPosition(200) == 601 : "getCoordinateAtPosition(200)";
		assert a.getPositionAtCoordinate(401) == 100 : "getPositionAtCoordinate";
		
		assert a.containsExactInterval(new SingleInterval(601,701));
		assert a.containsInterval(new SingleInterval(601, 701));
		
		assert a.contains(b) == false;
		
		b.removeInterval(new SingleInterval(1,101));
		assert b.numBlocks() == 1;
		assert b.length() == 100;
		assert a.contains(b);
		assert a.containsExactInterval(b.getBlocks().first());
		assert a.containsExactInterval(b.getBlocks().last());
		
		CompoundInterval c = new CompoundInterval(b);
		c.setStart(501);
		assert c.getSpan() == 200;
		assert c.length() == 200;
		assert b.length() == 100;
		
		CompoundInterval d = new CompoundInterval();
		d.addInterval(new SingleInterval(100,501));
		assert !(c.overlaps(d));
		d.addInterval(new SingleInterval(503,550));
		assert c.overlaps(d);
		
		System.out.println(c);
		System.out.println(d);
		System.out.println(c.intersect(d));
		System.out.println(c.union(d));
		
		d.setEnd(502);
		System.out.println(d);
		System.out.println(d.intersect(c));
		assert d.length() == 402;
		assert d.intersect(c).numBlocks() == 1;
		
		d.setEnd(400);
		System.out.println(d);
		assert d.length() == 300;
		
		a.addInterval(new SingleInterval(100, 650));
		System.out.println(a);
		assert a.overlaps(new SingleInterval(700,700));
		assert !a.overlaps(new SingleInterval(701,701));
		
		CompoundInterval e = new CompoundInterval();
		e.addInterval(502,502);
		e.addInterval(600,700);
		assert !e.overlaps(c);
		assert e.overlaps(c, 1);
		assert e.overlaps(c,100);
		assert e.overlaps(c, true);
	}

}

package nextgen.core.tests;

import java.util.ArrayList;
import java.util.Collection;

import nextgen.core.feature.GeneWindow;
import nextgen.core.feature.GenomeWindow;
import nextgen.core.feature.Window;

import nextgen.core.annotation.*;
import junit.framework.TestCase;

/**
 * This class tests the nextgen.core.feature.Window class 
 * @author skadri
 *
 */
public class TestWindow extends TestCase{

	/**
	 * Testing All Window constructors
	 */
	public void testConstructors() {
		Collection<Annotation> vertices = new ArrayList<Annotation>();
		BasicAnnotation v1 = new BasicAnnotation("chr2", 10,20,"+");
		BasicAnnotation v2 = new BasicAnnotation("chr2", 50,60,"+");
		vertices.add(v1);
		vertices.add(v2);
		/*
		 * Testing constructor
		 * Window(LightweightGenomicAnnotation region) 
		 */
		Window X = new GenomeWindow(v1);
		System.out.println(X.toBED());
		X = new GenomeWindow(v2);
		System.out.println(X.toBED());
		
		/*
		 * Testing constructor
		 * Window(Collection<LightweightGenomicAnnotation> blocks)
		 */
		X = new GeneWindow(new Gene(vertices));
		System.out.println(X.toBED());

		/*
		 * Testing constructor
		 *  Window(String chr,int start,int end)
		 */
		X = new GenomeWindow("chr1",400,600);
		System.out.println(X.toBED());
		
		/*
		 * Testing constructor
		 * Window(Gene gene)
		 */
		X = new GenomeWindow(new BasicAnnotation(vertices));
		System.out.println(X.toBED());
	}
}

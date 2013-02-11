package nextgen.core.tests;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.Window;

import junit.framework.TestCase;

public class TestCoordinateSpace extends TestCase{

	/**
	 * Testing All GenomicSpace constructors and functions
	 * @throws IOException 
	 */
	public void testGenomicSpace() throws IOException {
		
		System.err.println("Testing Genomic space constructors");
		Map<String,Integer> a = new HashMap<String,Integer>();
		a.put("chr1", 10000);
		a.put("chr2", 50000);
		/*
		 * Constructor
		 */
		GenomicSpace G = new GenomicSpace(a);
		/*
		 * Test getFragment()
		 */
		Collection<? extends Window> frags = G.getFragment("chr1", 20, 4000);
		for(Window w: frags){
			System.out.println(w.toBED());
		}
	
		/*
		 * What if the fragment is outside the coordinate space
		 */
		frags = G.getFragment("chr1", 20, 40000);
		for(Window w: frags){
			System.out.println(w.toBED());
		}
		
		/*
		 * TESTING WINDOW ITERATORS
		 * getWindowIterator(int windowSize, String chr,int start, int end,int overlap)
		 */
		Iterator<Window> iter = G.getWindowIterator(100, "chr1",10, 600,50);
		System.err.println("Starting first window iterator:");
		while(iter.hasNext()){
			Window w= iter.next();
			System.out.println(w.toBED());
		}
		/*
		 * getWindowIterator(String chr)
		 */
		iter = G.getWindowIterator("chr1", 100, 99);
		System.err.println("Starting chromosome window iterator:");
		while(iter.hasNext()){
			Window w= iter.next();
			System.out.println(w.toBED());
		}
		
	}
	
	/**
	 * Testing All GenomicSpace constructors
	 * @throws IOException 
	 */
	public void testTranscriptomeSpace() throws IOException {
		
		System.err.println("Testing the transcriptome space now:");
		Map<String,Collection<Gene>> gmap = new HashMap<String,Collection<Gene>>();
		
		//Read the transcriptome space
		BufferedReader br = null;
		br = new BufferedReader(new FileReader("test.bed"));
		String nextLine;
    	while ((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)) {
    		System.err.println(nextLine);
    		Gene X = new Gene(nextLine,false);
    		if(!gmap.containsKey(X.getChr())){
    			gmap.put(X.getChr(), new ArrayList<Gene>());
    		}
    		gmap.get(X.getChr()).add(X);
    	}
    	
		/*
		 * Constructor
		 */
		TranscriptomeSpace G = new TranscriptomeSpace(gmap);
		System.err.println("Space constructed");
		/*
		 * Test getFragment()
		 */
		Collection<? extends Window> frags = G.getFragment("chr4", 86822, 87864);
		for(Window w: frags){
			System.out.println(w.toBED());
		}
	
		/*
		 * What if the fragment is outside the coordinate space
		 */
		frags = G.getFragment("chr1", 20, 40000);
		if(frags!=null)
			for(Window w: frags){
				System.out.println(w.toBED());
			}
		
		/*
		 * TESTING WINDOW ITERATORS
		 * getWindowIterator(int windowSize, String chr,int start, int end,int overlap)
		 */
		Iterator<Window> iter = G.getWindowIterator(100, "chr4",85682, 87822,0);
		System.err.println("Starting first window iterator: getWindowIterator(int windowSize, String chr,int start, int end,int overlap)");
		while(iter.hasNext()){
			Window w= iter.next();
			System.out.println(w.toBED());
		}
		
		/*
		 * TESTING WINDOW ITERATORS
		 * getWindowIterator(int windowSize, String chr,int start, int end,int overlap)
		 */
		iter = G.getWindowIterator(100, "chr4",1, 87822,0);
		System.err.println("Same iterator with defined region outside of the gene");
		while(iter.hasNext()){
			Window w= iter.next();
			System.out.println(w.toBED());
		}
		/*
		 * getWindowIterator(String chr)
		 */
/*		iter = G.getWindowIterator("chr4");
		System.out.println("Starting chromosome window iterator: getWindowIterator(String chr)");
		while(iter.hasNext()){
			Window w= iter.next();
			System.out.println(w.toBED());
		}*/
		
		/*
		 * Transcriptome space with overlapping exons
		 */
		List<LightweightGenomicAnnotation> vertices = new ArrayList<LightweightGenomicAnnotation>();
		BasicLightweightAnnotation v1 = new BasicLightweightAnnotation("chrM", 10,20);
		BasicLightweightAnnotation v2 = new BasicLightweightAnnotation("chrM", 40,50);
		BasicLightweightAnnotation v3 = new BasicLightweightAnnotation("chrM", 70,100);
		vertices.add(v1);
		vertices.add(v2);
		vertices.add(v3);
		
		List<LightweightGenomicAnnotation> vertices2 = new ArrayList<LightweightGenomicAnnotation>();
		BasicLightweightAnnotation v4 = new BasicLightweightAnnotation("chrM", 10,20);
		BasicLightweightAnnotation v5 = new BasicLightweightAnnotation("chrM", 30,60);
		BasicLightweightAnnotation v6 = new BasicLightweightAnnotation("chrM", 80,100);
		vertices2.add(v4);
		vertices2.add(v5);
		vertices2.add(v6);
		Collection<Gene> genes = new ArrayList<Gene>();
		genes.add(new Gene(vertices));
		genes.add(new Gene(vertices2));
		Map<String,Collection<Gene>> xx = new HashMap<String,Collection<Gene>>();
		xx.put("chrM", genes);
		
		System.err.println("Testing the transcriptome space with overlapping genes with hand coded genes:");
		TranscriptomeSpace C2 = new TranscriptomeSpace(xx);
		/*
		 * TESTING WINDOW ITERATORS
		 * getWindowIterator(int windowSize, String chr,int start, int end,int overlap)
		 */
		Iterator<Window> iter2 = C2.getWindowIterator(10,"chrM",0,10000,5);
		System.err.println("Starting window iterator: getWindowIterator(String chr)");
		while(iter2.hasNext()){
			Window w= iter2.next();
			System.out.println(w.toBED());
		}
	}
	public static String whitespaceDelimiter = "\\s++"; //$NON-NLS-1$
}

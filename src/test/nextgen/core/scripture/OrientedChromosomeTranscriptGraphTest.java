qpackage nextgen.core.scripture;


import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.scripture.OrientedChromosomeTranscriptGraph.TranscriptGraphEdge;

import org.jgrapht.GraphPath;
import org.junit.Test;

import junit.framework.TestCase;

public class OrientedChromosomeTranscriptGraphTest extends TestCase {

	@Test
	public void testGraphConstruction() {
		Annotation exon1 = new BasicAnnotation("chr1", 1000, 1200,"+");
		Annotation exon2 = new BasicAnnotation("chr1", 1500, 1700,"+");
		Annotation exon3 = new BasicAnnotation("chr1", 1900, 2100,"+");
		List<Annotation> exonList1 = new ArrayList<Annotation> (3);
		
		exonList1.add(exon1); exonList1.add(exon2); exonList1.add(exon3);
		
		Annotation threeExonAnnotation = new BasicAnnotation(exonList1);
		
		List<Annotation> exonList2 = new ArrayList<Annotation> (2);
		exonList2.add(exon1);  exonList2.add(exon3);
		Annotation twoExonAnnotation = new BasicAnnotation(exonList2);
		
		OrientedChromosomeTranscriptGraph graph = new OrientedChromosomeTranscriptGraph("test","+");
		
		boolean couldAdd = graph.addAnnotationToGraph(threeExonAnnotation);
		assertTrue("Could not add positive annotation to positive graph!!!", couldAdd);
		couldAdd = graph.addAnnotationToGraph(twoExonAnnotation);
		assertTrue("Could not add positive annotation to positive graph!!!", couldAdd);
		//graph.addEdge(exon1, exon2);
		//graph.addEdge(exon2, exon3);
		//graph.addEdge(exon1, exon3);
		
		Collection<GraphPath<Annotation, TranscriptGraphEdge>> paths = graph.getPaths();
		
		for (GraphPath<Annotation, TranscriptGraphEdge> p : paths) {
			System.err.println(OrientedChromosomeTranscriptGraph.pathToGene(p));
		}
		
	}
	
	public void testAddingOrphanVertex() {
		Annotation exon1 = new BasicAnnotation("chr1", 1000, 1200,"+");
		Annotation exon2 = new BasicAnnotation("chr1", 1500, 1700,"+");
		Annotation exon3 = new BasicAnnotation("chr1", 1900, 2100,"+");
		Annotation exon4 = new BasicAnnotation("chr1", 3900, 4500,"+");
		List<Annotation> exonList = new ArrayList<Annotation> (3);
		
		exonList.add(exon1); exonList.add(exon2); exonList.add(exon3);
		
		Annotation threeExonAnnotation = new BasicAnnotation(exonList);
		
		List<Annotation> exonList2 = new ArrayList<Annotation> (2);
		exonList2.add(exon1);  exonList.add(exon3);
		Annotation twoExonAnnotation = new BasicAnnotation(exonList2);
		
		OrientedChromosomeTranscriptGraph graph = new OrientedChromosomeTranscriptGraph("test","+");
		graph.addAnnotationToGraph(threeExonAnnotation);
		graph.addAnnotationToGraph(twoExonAnnotation);
		//graph.addEdge(exon1, exon2);
		//graph.addEdge(exon2, exon3);
		//graph.addEdge(exon1, exon3);
		
		graph.addVertex(exon4);
		
		Collection<GraphPath<Annotation, TranscriptGraphEdge>> paths = graph.getPaths();
		
		for (GraphPath<Annotation, TranscriptGraphEdge> p : paths) {
			System.err.println(OrientedChromosomeTranscriptGraph.pathToGene(p));
		}
		
		Collection<Annotation> orphans = graph.getOrphanVertices();
		
		for (Annotation ov : orphans) {
			System.err.println(ov);
		}
	}
	
	public void testAddingOverlappingExons() {
		Annotation exon1 = new BasicAnnotation("chr1", 1000, 1200,"+");
		Annotation exon2 = new BasicAnnotation("chr1", 1500, 1700,"+");
		Annotation exon3 = new BasicAnnotation("chr1", 1100, 1200,"+");

		List<Annotation> exonList1 = new ArrayList<Annotation> (2);
		exonList1.add(exon1); exonList1.add(exon2); 
		Annotation exonAnnotation1 = new BasicAnnotation(exonList1);
		
		List<Annotation> exonList2 = new ArrayList<Annotation> (2);
		exonList2.add(exon3); exonList2.add(exon2); 
		Annotation exonAnnotation2 = new BasicAnnotation(exonList2);
		
		OrientedChromosomeTranscriptGraph graph = new OrientedChromosomeTranscriptGraph("test","+");
		
		graph.addAnnotationToGraph(exonAnnotation1);
		graph.addAnnotationToGraph(exonAnnotation2);
		Collection<GraphPath<Annotation, TranscriptGraphEdge>> paths = graph.getPaths();
		System.err.println("Adding overlapping exons test");
		for (GraphPath<Annotation, TranscriptGraphEdge> p : paths) {
			System.err.println(OrientedChromosomeTranscriptGraph.pathToGene(p));
		}
	}
	
	public void testAddingReverseOrientedSubgraphs() {
		Annotation exon1 = new BasicAnnotation("chr1", 1000, 1200,"-");
		Annotation exon2 = new BasicAnnotation("chr1", 1500, 1700,"-");
		Annotation exon3 = new BasicAnnotation("chr1", 1100, 1200,"-");

		List<Annotation> exonList1 = new ArrayList<Annotation> (2);
		exonList1.add(exon1); exonList1.add(exon2); 
		Annotation exonAnnotation1 = new BasicAnnotation(exonList1);
		
		List<Annotation> exonList2 = new ArrayList<Annotation> (2);
		exonList2.add(exon3); exonList2.add(exon2); 
		Annotation exonAnnotation2 = new BasicAnnotation(exonList2);
		System.err.println(exonAnnotation1);
		System.err.println(exonAnnotation2);
		OrientedChromosomeTranscriptGraph graph = new OrientedChromosomeTranscriptGraph("test","-");
		
		graph.addAnnotationToGraph(exonAnnotation1);
		graph.addAnnotationToGraph(exonAnnotation2);
		Collection<GraphPath<Annotation, TranscriptGraphEdge>> paths = graph.getPaths();
		
		System.err.println("testAddingReverseOrientedSubgraphs");
		for (GraphPath<Annotation, TranscriptGraphEdge> p : paths) {
			System.err.println(OrientedChromosomeTranscriptGraph.pathToGene(p));
		}
	}
	
}

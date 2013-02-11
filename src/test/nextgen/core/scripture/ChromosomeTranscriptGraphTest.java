package nextgen.core.scripture;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.jgrapht.GraphPath;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.scripture.OrientedChromosomeTranscriptGraph.TranscriptGraphEdge;
import junit.framework.TestCase;

	public class ChromosomeTranscriptGraphTest extends TestCase{
		public void testCycleHandling() {
		Annotation exon1 = new BasicAnnotation("chr1", 1000, 1200,"-");
		Annotation exon2 = new BasicAnnotation("chr1", 1500, 1700,"-");
		Annotation exon3 = new BasicAnnotation("chr1", 2100, 2300,"-");
		Annotation exon4 = new BasicAnnotation("chr1", 1100, 1200,"+");
		Annotation exon5 = new BasicAnnotation("chr1", 1500, 1700,"+");
	
		List<Annotation> exonList1 = new ArrayList<Annotation> (2);
		exonList1.add(exon1); exonList1.add(exon2); exonList1.add(exon3); 
		Annotation exonAnnotation1 = new BasicAnnotation(exonList1);
		
		List<Annotation> exonList2 = new ArrayList<Annotation> (2);
		exonList2.add(exon4); exonList2.add(exon5); 
		Annotation exonAnnotation2 = new BasicAnnotation(exonList2);
		System.err.println(exonAnnotation1);
		System.err.println(exonAnnotation2);
		ChromosomeTranscriptGraph graph = new ChromosomeTranscriptGraph("test");
		
		graph.addAnnotationToGraph(exonAnnotation1);
		graph.addAnnotationToGraph(exonAnnotation2);
		Collection<GraphPath<Annotation, TranscriptGraphEdge>> paths = graph.getPaths();
		
		System.err.println("Cycles test");
		for (GraphPath<Annotation, TranscriptGraphEdge> p : paths) {
			System.err.println(graph.pathToGene(p));
		}
	}


	public void testConnectVertexToGraph() {
		Annotation exon1 = new BasicAnnotation("chr1", 1000, 1200,"-");
		Annotation exon2 = new BasicAnnotation("chr1", 1500, 1700,"-");
		Annotation exon3 = new BasicAnnotation("chr1", 2100, 2300,"-");
		Annotation exon4 = new BasicAnnotation("chr1", 1100, 1200,"+");
		Annotation exon5 = new BasicAnnotation("chr1", 1500, 1700,"+");
	
		List<Annotation> exonList1 = new ArrayList<Annotation> (2);
		exonList1.add(exon1); exonList1.add(exon2); exonList1.add(exon3); 
		Annotation exonAnnotation1 = new BasicAnnotation(exonList1);
		
		List<Annotation> exonList2 = new ArrayList<Annotation> (2);
		exonList2.add(exon4); exonList2.add(exon5); 
		Annotation exonAnnotation2 = new BasicAnnotation(exonList2);
		System.err.println(exonAnnotation1);
		System.err.println(exonAnnotation2);
		ChromosomeTranscriptGraph graph = new ChromosomeTranscriptGraph("test");
		
		graph.addAnnotationToGraph(exonAnnotation1);
		graph.addAnnotationToGraph(exonAnnotation2);
		
		graph.connectVertexToGraph(new BasicAnnotation("chr1", 1150, 1200,"."));
		Collection<GraphPath<Annotation, TranscriptGraphEdge>> paths = graph.getPaths();
		
		System.err.println("testConnectVertexToGraph");
		for (GraphPath<Annotation, TranscriptGraphEdge> p : paths) {
			System.err.println(graph.pathToGene(p));
		}
	}

}

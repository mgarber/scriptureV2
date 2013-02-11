package nextgen.core.scripture;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.jgrapht.GraphPath;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.scripture.OrientedChromosomeTranscriptGraph.TranscriptGraphEdge;

public class ChromosomeTranscriptGraph implements Serializable{

	private static final long serialVersionUID = -6516913539814555024L;
	private String name;
	private OrientedChromosomeTranscriptGraph plusGraph;
	private OrientedChromosomeTranscriptGraph negativeGraph;
	
	public ChromosomeTranscriptGraph(String name) {
		this.name = name;
		this.plusGraph = new OrientedChromosomeTranscriptGraph(name, "+");
		this.negativeGraph = new OrientedChromosomeTranscriptGraph(name, "-");
	}

	public void connectVertexToGraph(Annotation annotation) {
		if(annotation.getStrand() == Strand.POSITIVE || annotation.isUnoriented()) {
			plusGraph.connectVertexToGraph(annotation);
		} 
		if(annotation.getStrand() == Strand.NEGATIVE || annotation.isUnoriented()) {
			negativeGraph.connectVertexToGraph(annotation);
		}
		
	}

	public void addAnnotationToGraph(Annotation annotation) {
		if(annotation.getStrand() == Strand.POSITIVE || annotation.isUnoriented()) {
			plusGraph.addAnnotationToGraph(annotation);
		} 
		if(annotation.getStrand() == Strand.NEGATIVE || annotation.isUnoriented()) {
			negativeGraph.addAnnotationToGraph(annotation);
		}
		
	}

	public List<GraphPath<Annotation, TranscriptGraphEdge>> getPaths() {
		List<GraphPath<Annotation, TranscriptGraphEdge>> paths = new ArrayList<GraphPath<Annotation, TranscriptGraphEdge>>();
		paths.addAll(plusGraph.getPaths());
		paths.addAll(negativeGraph.getPaths());
		return paths;
	}

	public Gene pathToGene(GraphPath<Annotation, TranscriptGraphEdge> path) {
		return OrientedChromosomeTranscriptGraph.pathToGene(path);
	}

}

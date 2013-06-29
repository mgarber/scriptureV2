package nextgen.core.feature;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.pda.datastructures.Alignments;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.CompoundInterval;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.SingleInterval;


public interface Window extends Annotation{

	/**
	 * Add a pointer to the annotation from which this Window originated
	 * @param annotation
	 */
	void addSourceAnnotation(Annotation annotation);

	/**
	 * @return a collection of annotations from which this originated
	 */
	Collection<? extends Annotation> getSourceAnnotations();
	
	/**
	 * Get collection of windows spanning the annotation
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @return The collection of windows
	 */
	//public Collection<Window> getWindows(int windowSize, int stepSize);
	
}

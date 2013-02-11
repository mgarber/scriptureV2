package nextgen.core.annotation;

import java.io.*;
import java.util.Collection;

import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.general.Predicates;
import org.apache.commons.collections15.Predicate;

/**
 * @author engreitz
 * This class is responsible for reading annotation files and generating AnnotationLists from them.
 * Parsing logic is partially copied from broad.core.annotation.AnnotationReader
 */
public class AnnotationFileReader {

	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, AnnotationFactory<? extends T> factory) {
		return load(file, clazz, factory, null);
	}
	
	
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, AnnotationFactory<? extends T> factory, CoordinateSpace cs) {
		return load(file, clazz, factory, cs, Predicates.alwaysTrue());	
	}
	
	
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, AnnotationFactory<? extends T> factory, CoordinateSpace cs, Collection<Predicate<? super T>> filters) {
		return load(file, clazz, factory, cs, Predicates.and(filters));
	}

	
	/**
	 * Reads an Annotation file line by line. 
	 * @param <T>
	 * @param file File containing annotations
	 * @param clazz	The class of the desired AnnotationList
	 * @param cs CoordinateSpace with which to initialize the AnnotationList. Can be null.
	 * @param factory Factory for parsing and creating the annotation type
	 * @param filter Annotation filters to control the subset of annotations that are stored in the AnnotationList.
	 * @return AnnotationList containing all annotations from file that pass the filter
	 */
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, AnnotationFactory<? extends T> factory, CoordinateSpace cs, Predicate<? super T> filter) {
		AnnotationList<T> annotations = new AnnotationList<T>(cs);
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			
			String line;
			while((line = br.readLine()) != null) {
			
				if (line.toLowerCase().startsWith("track")) {
					// TODO
					throw new IllegalArgumentException("AnnotationFileReader does not support files with track headers (TODO)");
				}
				
				line = line.trim();
				if(line.startsWith("#") || line.length() == 0) {
					continue;
				}
				
				String[] lineSplit = line.split("\t");
				
				T annotation = factory.create(lineSplit);
				if (filter.evaluate(annotation)) {
					annotations.add(annotation);
				}
			}
				
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e.getMessage());
		}
		
		return annotations;
	}
}

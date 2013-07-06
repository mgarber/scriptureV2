package nextgen.core.annotation;

import java.io.*;
import java.util.Collection;
import java.util.Iterator;

import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.general.TabbedReader;
import nextgen.core.general.Predicates;
import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.iterators.FilterIterator;
import org.apache.commons.io.LineIterator;

/**
 * @author engreitz
 * This class is responsible for reading annotation files and generating AnnotationLists from them.
 * Parsing logic is partially copied from broad.core.annotation.AnnotationReader
 */
public class AnnotationFileReader {

	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory) {
		return load(file, clazz, factory, null);
	}
	
	
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, CoordinateSpace cs) {
		return load(file, clazz, factory, cs, Predicates.alwaysTrue());	
	}
	
	
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, CoordinateSpace cs, Collection<Predicate<? super T>> filters) {
		return load(file, clazz, factory, cs, Predicates.and(filters));
	}

	
	/**
	 * Reads an Annotation file line by line and return an AnnotationList containing the results
	 * @param <T>
	 * @param file File containing annotations
	 * @param clazz	The class of the desired AnnotationList
	 * @param cs CoordinateSpace with which to initialize the AnnotationList. Can be null.
	 * @param factory Factory for parsing and creating the annotation type
	 * @param filter Annotation filters to control the subset of annotations that are stored in the AnnotationList.
	 * @return AnnotationList containing all annotations from file that pass the filter
	 */
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, CoordinateSpace cs, Predicate<? super T> filter) {
		AnnotationList<T> annotations = new AnnotationList<T>(cs);
		Iterator<T> itr = read(file, clazz, factory, filter);
		while (itr.hasNext()) {
			annotations.add(itr.next());
		}		
		return annotations;
	}
	
	
	/**
	 * Read an Annotation file line by line in a buffered format and return an iterator over the results.
	 * Should be refactored with the "load" function as appropriate ("load" can call this function and then
	 * add it to its Annotation list)
	 * @param file
	 * @param clazz
	 * @param factory
	 * @param cs
	 * @param filter
	 * @return
	 */
	public static <T extends Annotation> Iterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, Predicate<? super T> filter) {
		return new FilterIterator<T>(new AnnotationIterator<T>(file, factory), filter);
	}
	
	
	public static <T extends Annotation> Iterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory) {
		return read(file, clazz, factory, Predicates.alwaysTrue());
	}
	
	
	private static class AnnotationIterator<T> extends nextgen.core.general.TabbedReader.TabbedIterator<T> {
		
		public AnnotationIterator(File file, TabbedReader.Factory<? extends T> factory) {
			super(file, factory);
		}

		@Override
		protected String getNextLine() {
			while (itr.hasNext()) {
				String line = itr.next();
				if (line.toLowerCase().startsWith("track")) {
					throw new IllegalArgumentException("AnnotationFileReader does not support files with track headers (TODO)");
				}
				
				line = line.trim();
				if (line.startsWith("#") || line.length() == 0) {
					continue;
				} else {
					return line;
				}
			}
			return null;
		}
	}
	
}

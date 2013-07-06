package nextgen.core.general;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Iterator;

import broad.core.error.ParseException;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.iterators.FilterIterator;
import org.apache.commons.io.LineIterator;

public class TabbedReader {

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
	public static <T> Iterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, Predicate<? super T> filter) {
		return new FilterIterator<T>(new TabbedIterator<T>(file, factory), filter);
	}
	
	public static <T> Iterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory) {
		return read(file, clazz, factory, Predicates.alwaysTrue());
	}
	
	public static class TabbedIterator<T> implements Iterator<T> {
		protected LineIterator itr;
		private T curr;
		protected Factory<? extends T> factory;
		
		public TabbedIterator(File file, Factory<? extends T> factory) {
			try {
				itr = new LineIterator(new BufferedReader(new FileReader(file)));
			} catch (Exception e) {
				e.printStackTrace();
				throw new RuntimeException(e.getMessage());
			}
			this.factory = factory;
			advance();
		}
		
		@Override
		public boolean hasNext() {
			return (curr != null);
		}

		@Override
		public T next() {
			T result = curr;
			advance();
			return result;
		}

		private void advance() {
			curr = null;
			String nextLine = getNextLine();
			if (nextLine != null) {
				nextLine.trim();
				curr = factory.create(nextLine.split("\t"));
			}
		}
		
		/**
		 * Override this function to do any custom parsing (e.g., skip commented lines)
		 * @return
		 */
		protected String getNextLine() {
			if (itr.hasNext()) return itr.next();
			else return null;
		}
		
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Remove not supported");
		}
		
	}
	
	
	public interface Factory<T> {
		T create(String[] rawFields) throws ParseException;
	}
}

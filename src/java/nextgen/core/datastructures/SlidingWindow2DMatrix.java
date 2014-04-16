package nextgen.core.datastructures;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;

import nextgen.core.annotation.Annotation;
import nextgen.core.feature.Window;

public interface SlidingWindow2DMatrix {
	
	/**
	 * Get an iterator over all windows corresponding to rows and columns of the matrix
	 * @return Window iterator
	 */
	public Iterator<? extends Window> slideWindows();
	
	/**
	 * Slide windows and populate the entire matrix one row at a time
	 */
	public void populateAll();
	
	/**
	 * Populate a row of the matrix corresponding to one window
	 * @param row The window
	 */
	public void populate(Window row);
	
	/**
	 * Get the matrix entry for a specified row region and column region
	 * Best matrix entry for these regions is determined by implementing class
	 * @param row Row region
	 * @param column Column region
	 * @return Matrix entry in best row corresponding to row region and best column corresponding to column region
	 */
	public double getValue(Annotation row, Annotation column);
	
	/**
	 * Set the matrix entry for a specified row region and column region
	 * Best matrix entry for these regions is determined by implementing class
	 * @param row Row region
	 * @param column Column region
	 * @param value Value to set
	 */
	public void setValue(Annotation row, Annotation column, double value);
	
	/**
	 * Set the matrix entry specified by row and column indices
	 * @param row Row
	 * @param column Column
	 * @param value Value to set
	 */
	public void setValue(int row, int column, double value);
	
	/**
	 * Get the index of the window that most closely overlaps a specified region
	 * @param region The region
	 * @return Closest row/column index of the window in the matrix, or -1 if none
	 */
	public int getBestIndex(Annotation region);
	
	/**
	 * Get matrix indices of all windows that overlap a specified region
	 * @param region The region
	 * @return Collection of row/column indices for all windows overlapping the region
	 */
	public Collection<Integer> getAllIndices(Annotation region);
	
	/**
	 * Get the window corresponding to a row/column index
	 * @param index The index
	 * @return The window represented by this index in the matrix
	 */
	public Window getWindow(int index);
	
	/**
	 * Write the data in the matrix to a file
	 * @param outFile Output file
	 * @param maxEntries Max number of entries in matrix to write file; if this number is exceeded do not write the file
	 * @throws IOException 
	 */
	public void writeToFile(String outFile, int maxEntries) throws IOException;
	
	/**
	 * Paint the matrix to an image file
	 * @param outImageFile Output image file
	 * @throws IOException 
	 */
	public void writeToImage(String outImageFile) throws IOException;
	
	/**
	 * Get the number of rows in the matrix
	 * @return Number of rows in the matrix
	 */
	public int numRows();
	
	/**
	 * Get the number of columns in the matrix
	 * @return Number of columns in the matrix
	 */
	public int numColumns();
	
}

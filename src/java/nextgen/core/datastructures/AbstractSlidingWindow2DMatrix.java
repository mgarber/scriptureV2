package nextgen.core.datastructures;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.feature.GenomeWindow;
import nextgen.core.feature.Window;
import nextgen.core.utils.CountLogger;

public abstract class AbstractSlidingWindow2DMatrix implements SlidingWindow2DMatrix {
	
	protected int windowSize;
	protected int stepSize;
	protected CoordinateSpace coordinateSpace;
	protected SparseDoubleMatrix2D matrix;
	protected Map<String, Integer> chrSizes;
	
	protected static Logger logger = Logger.getLogger(AbstractSlidingWindow2DMatrix.class.getName());
	
	/**
	 * Map of chromsome to interval tree of windows represented by rows and columns in the matrix
	 * The value in the tree is the row/column index of the interval in the large matrix
	 */
	private Map<String, IntervalTree<Integer>> windows;
	
	/**
	 * Generate the interval tree of windows for each chromosome by sliding windows across the coordinate space
	 * The value in the tree will be the row/column index of the interval in the large matrix
	 * Call this in every constructor
	 */
	protected void generateWindowTree() {
		logger.info("");
		logger.info("Generating interval tree of windows...");
		windows = new TreeMap<String, IntervalTree<Integer>>();
		Iterator<? extends Window> windowIter = coordinateSpace.getWindowIterator(windowSize, windowSize - stepSize);
		int index = 0;
		while(windowIter.hasNext()) {
			Window window = windowIter.next();
			String chr = window.getChr();
			if(!windows.containsKey(chr)) {
				logger.info("Found " + chr);
				windows.put(chr, new IntervalTree<Integer>());
			}
			@SuppressWarnings("unused")
			Integer dummy = windows.get(chr).placeInTree(window.getStart(), window.getEnd(), Integer.valueOf(index));
			index++;
		}
		logger.info("Done generating window tree.");
	}
	
	/**
	 * Initialize the matrix
	 * Window tree must have been generated already
	 * Call this in every constructor
	 */
	protected void initializeMatrix() {
		logger.info("");
		logger.info("Initializing matrix...");
		int totalWindows = 0;
		for(String chr : windows.keySet()) {
			totalWindows += windows.get(chr).size();
		}
		matrix = new SparseDoubleMatrix2D(totalWindows, totalWindows);
		logger.info("Matrix has " + matrix.rows() + " rows and " + matrix.columns() + " columns.");
	}
	
	public Iterator<? extends Window> slideWindows() {
		return coordinateSpace.getWindowIterator(windowSize, windowSize - stepSize);
	}
	
	@Override
	public int getBestIndex(Annotation region) {
		//logger.debug("Getting index for region " + region.toUCSC());
		// First check if the exact interval is in the tree
		String chr = region.getChr();
		int start = region.getStart();
		int end = region.getEnd();
		IntervalTree<Integer> tree = windows.get(chr);
		if(tree == null) {
			//logger.warn("Can't get index for " + region.toUCSC() + " because " + chr + " has no windows.");
			return -1;
		}
		Node<Integer> node;
		node = tree.find(start, end);
		if(node != null) {
			return node.getValue().intValue();
		}
		int numOverlappers = tree.numOverlappers(start, end);
		if(numOverlappers == 0) {
			throw new IllegalArgumentException("No overlapping windows for " + region.toUCSC());
		}
		if(numOverlappers % 2 == 1) {
			Iterator<Integer> iter = tree.overlappingValueIterator(start, end);
			int best = (numOverlappers - 1) / 2;
			int c = 0;
			while(iter.hasNext()) {
				Integer i = iter.next();
				if(c == best) {
					return i.intValue();
				}
				c++;
			}
		} else {
			int best1 = numOverlappers / 2;
			int best2 = best1 - 1;
			Node<Integer> node1 = null;
			Node<Integer> node2 = null;
			Iterator<Node<Integer>> iter = tree.overlappers(start, end);
			int c = 0;
			while(iter.hasNext()) {
				Node<Integer> n = iter.next();
				if(c == best1) {
					node1 = n;
				}
				if(c == best2) {
					node2 = n;
				}
				if(node1 != null && node2 != null) {
					break;
				}
				c++;
			}
			int start1 = node1.getStart();
			int end1 = node1.getEnd();
			int start2 = node2.getStart();
			int end2 = node2.getEnd();
			int overlap1 = overlap(start, end, start1, end1);
			int overlap2 = overlap(start, end, start2, end2);
			if(overlap1 > overlap2) {
				return node1.getValue().intValue();
			}
			return node2.getValue().intValue();
		}
		throw new IllegalStateException("Something went wrong.");
	}
	
	private static int overlap(int start1, int end1, int start2, int end2) {
		if(start1 <= start2 && end1 <= end2) {
			return Math.max(0, end1 - start2);
		}
		if(start2 <= start1 && end2 <= end1) {
			return Math.max(0, end2 - start1);
		}
		if(start1 <= start2 && end1 >= end2) {
			return end2 - start2;
		}
		return end1 - start1;
	}
	
	@Override
	public Collection<Integer> getAllIndices(Annotation region) {
		Collection<Integer> rtrn = new TreeSet<Integer>();
		String chr = region.getChr();
		int start = region.getStart();
		int end = region.getEnd();
		Iterator<Node<Integer>> iter = windows.get(chr).overlappers(start, end);
		while(iter.hasNext()) {
			rtrn.add(iter.next().getValue());
		}
		return rtrn;
	}
	
	private String getChr(int index) {
		if(index < 0) {
			throw new IllegalArgumentException("Index must be >= 0");
		}
		int n = 0;
		for(String chr : windows.keySet()) {
			n += windows.get(chr).size();
			if(index < n) {
				return chr;
			}
		}
		throw new IllegalArgumentException("Max index is " + n + " (you provided " + index + ")");
	}

	@Override
	public Window getWindow(int index) {
		String chr = getChr(index);
		IntervalTree<Integer> tree = windows.get(chr);
		Iterator<Node<Integer>> iter = tree.iterator();
		while(iter.hasNext()) {
			Node<Integer> node = iter.next();
			if(node.getValue().intValue() == index) {
				return new GenomeWindow(chr, node.getStart(), node.getEnd());
			}
		}
		throw new IllegalArgumentException("Didn't find window with index " + index);
	}

	
	/**
	 * Get an iterator over the nonzero entries in the matrix, in no special order
	 * @return Iterator over the nonzero entries in the matrix
	 */
	public Iterator<Entry> nonzeroEntries() {
		IntArrayList rowNums = new IntArrayList();
		IntArrayList colNums = new IntArrayList();
		DoubleArrayList vals = new DoubleArrayList();
		matrix.getNonZeros(rowNums, colNums, vals);
		Collection<Entry> entries = new ArrayList<Entry>();
		for(int i = 0; i < rowNums.size(); i++) {
			entries.add(new Entry(rowNums.get(i), colNums.get(i), vals.get(i)));
		}
		return entries.iterator();
	}
	
	public double getValue(Annotation row, Annotation column) {
		int rowIndex = getBestIndex(row);
		int colIndex = getBestIndex(column);
		return matrix.get(rowIndex, colIndex);
	}
	
	public void populateAll() {
		logger.info("");
		logger.info("Populating matrix...");
		CountLogger cl = new CountLogger(numRows(), 100);
		Iterator<? extends Window> windowIterator = slideWindows();
		while(windowIterator.hasNext()) {
			cl.advance();
			Window window = windowIterator.next();
			populate(window);
		}
		logger.info("Done populating matrix.");
	}
	
	public int numRows() {
		return matrix.rows();
	}
	
	public int numColumns() {
		return matrix.columns();
	}
	
	public void setValue(Annotation row, Annotation column, double value) {
		int rowIndex = getBestIndex(row);
		int colIndex = getBestIndex(column);
		setValue(rowIndex, colIndex, value);
	}
	
	public void setValue(int row, int column, double value) {
		matrix.set(row, column, value);
	}
		
	public class Entry {
		
		private int rowIndex;
		private int colIndex;
		private double value;
		
		private Entry(int rowNum, int colNum, double val) {
			rowIndex = rowNum;
			colIndex = colNum;
			val = value;
		}
		
		public int rowIndex() {
			return rowIndex;
		}
		
		public int colIndex() {
			return colIndex;
		}
		
		public double value() {
			return value;
		}
		
	}
	

}

package nextgen.sequentialbarcode;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import javax.imageio.ImageIO;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.datastructures.AbstractSlidingWindow2DMatrix;
import nextgen.core.feature.Window;
import nextgen.sequentialbarcode.fragmentgroup.FragmentGroup;
import nextgen.sequentialbarcode.fragmentgroup.NamedBarcodedFragmentGroup;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class DNAInteractionMatrix2D extends AbstractSlidingWindow2DMatrix {
	
	public static Logger logger = Logger.getLogger(DNAInteractionMatrix2D.class.getName());
	private SAMFileReader samReader;
	private int minBarcodes;
	
	public DNAInteractionMatrix2D(String bamFileWithFragmentGroups, String chrSizeFile, int window, int step, int minBarcodesPerRead) throws IOException {
		windowSize = window;
		stepSize = step;
		minBarcodes = minBarcodesPerRead;
		logger.info("Window size is " + windowSize + " and step size is " + stepSize + ".");
		chrSizes = fromFile(chrSizeFile);
		coordinateSpace = new GenomicSpace(chrSizes);
		generateWindowTree();
		initializeMatrix();
		samReader = new SAMFileReader(new File(bamFileWithFragmentGroups));
		populateAll();
	}
	
	private static Map<String, Integer> fromFile(String file) throws IOException {
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		Map<String, Integer> rtrn = new TreeMap<String, Integer>();
		while(b.ready()) {
			s.parse(b.readLine());
			if(s.getFieldCount() != 2) {
				r.close();
				b.close();
				throw new IllegalArgumentException("Line format: chr_name   chr_size");
			}
			rtrn.put(s.asString(0), Integer.valueOf(s.asInt(1)));
		}
		r.close();
		b.close();
		return rtrn;
	}
	
	@Override
	public void populate(Window row) {
		SAMRecordIterator iter = samReader.query(row.getChr(), row.getStart(), row.getEnd(), false);
		if(!iter.hasNext()) {
			iter.close();
			return;
		}
		int rowIndex = getBestIndex(row);
		while(iter.hasNext()) {
			SAMRecord record = iter.next();
			FragmentGroup fragmentGroup = NamedBarcodedFragmentGroup.fromSAMRecord(record);
			if(fragmentGroup.getBarcodes().getNumBarcodes() < minBarcodes) {
				continue;
			}
			for(Annotation fragment : fragmentGroup.getRegions()) {
				try {
					int colIndex = getBestIndex(fragment);
					setValue(rowIndex, colIndex, 1);
				} catch (Exception e) {
					//logger.warn("Skipping fragment " + fragment.toUCSC() + ": " + e.getMessage());
				}
			}
		}
		iter.close();
	}
	
	public void writeToTable(String outFile) throws IOException {
		logger.info("");
		logger.info("Writing matrix to table " + outFile);
		FileWriter w = new FileWriter(outFile);
		Iterator<Entry> nonzeroEntries = nonzeroEntries();
		while(nonzeroEntries.hasNext()) {
			Entry entry = nonzeroEntries.next();
			String window1 = getWindow(entry.rowIndex()).toUCSC();
			String window2 = getWindow(entry.colIndex()).toUCSC();
			w.write(window1 + "\t" + window2 + "\n");
		}
		w.close();
		logger.info("Done writing matrix.");
	}
	
	@Override
	public void writeToFile(String outFile, int maxEntries) throws IOException {
		logger.info("");
		logger.info("Writing matrix to file " + outFile);
		int rows = matrix.rows();
		int cols = matrix.columns();
		int entries = rows * cols;
		if(entries > maxEntries) {
			logger.warn("NOT WRITING FILE because number of entries (" + entries + ") exceeds threshold (" + maxEntries + ").");
			return;
		}
		FileWriter w = new FileWriter(outFile);
		for(int i = 0; i < rows; i++) {
			String line = "";
			for(int j = 0; j < cols; j++) {
				line += matrix.get(i, j) + "\t";
			}
			w.write(line + "\n");
		}
		w.close();
		logger.info("Done writing matrix.");
	}
	
	@Override
	public void writeToImage(String outPrefix) throws IOException {
		String fileName = outPrefix + ".png";
		logger.info("");
		logger.info("Writing image to " + fileName + "...");
		BufferedImage img = new BufferedImage(numColumns(), numRows(), BufferedImage.TYPE_INT_RGB);
		Iterator<Entry> nonzeroEntries = nonzeroEntries();
		int r = 255;
		int g = 0;
		int b = 0;
		int rgb = (r << 16) | (g << 8) | b;
		while(nonzeroEntries.hasNext()) {
			Entry entry = nonzeroEntries.next();
			int row = entry.rowIndex();
			int col = entry.colIndex();
			img.setRGB(col, row, rgb);
		}
		File file = new File(fileName);
		ImageIO.write(img, "PNG", file);
		logger.info("Done writing image.");
	}
	
	
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Input barcoded grouped bam file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-o", "Output prefix", true);
		p.addBooleanArg("-d", "Debug logging on", false, false);
		p.addIntArg("-mb", "Min number of barcodes to count fragment", true);
		p.addIntArg("-w", "Window size", true);
		p.addIntArg("-s", "Step size", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String chrSizeFile = p.getStringArg("-c");
		String outPrefix = p.getStringArg("-o");
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-s");
		int minBarcodes = p.getIntArg("-mb");
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
			AbstractSlidingWindow2DMatrix.logger.setLevel(Level.DEBUG);
		}
		
		DNAInteractionMatrix2D matrix = new DNAInteractionMatrix2D(bamFile, chrSizeFile, windowSize, stepSize, minBarcodes);
		matrix.writeToFile(outPrefix + ".matrix", 500000000);
		matrix.writeToImage(outPrefix);
		matrix.writeToTable(outPrefix + ".table");
		
		logger.info("");
		logger.info("All done.");
		
	}


}

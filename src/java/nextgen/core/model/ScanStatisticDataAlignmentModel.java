package nextgen.core.model;
import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.parser.StringParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationCollection;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.feature.GenomeWindow;
import nextgen.core.feature.Window;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readers.PairedEndReader;
import nextgen.core.readers.PairedEndWriter;
import nextgen.core.utils.AnnotationUtils;
import nextgen.core.exception.RuntimeIOException;
import org.apache.commons.collections15.Predicate;


/**
 * @author engreitz
 * Analysis module for calculating scan statistics on a collection of reads.
 * This is a hold-over class from the move to a generic scoring function for the
 * DataAlignmentModel, and makes two main changes to the DataAlignmentModel:
 * it makes a PEBAM by default, and uses the ScanStatisticScore analysis module
 * instead of the CountScore analysis module.
 * TODO maybe should implement PeakCaller instead of extending DataAlignmentModel
 */
public class ScanStatisticDataAlignmentModel extends AlignmentModel {

	//TODO Should have a permute function
	static Logger logger = Logger.getLogger(ScanStatisticDataAlignmentModel.class.getName());

	/**
	 * Build with a BAM file
	 * Populate the alignment collection
	 * @param bamFile
	 * @param coordinateSpace 
	 * @param readFilters 
	 * @throws IOException 
	 */
	public ScanStatisticDataAlignmentModel(String bamFile, CoordinateSpace coordinateSpace, Collection<Predicate<Alignment>> readFilters) {
		super(bamFile, coordinateSpace, readFilters);
	}
	
	/**
	 * Build with a BAM file
	 * Populate the alignment collection
	 * @param bamFile
	 * @param coordinateSpace 
	 * @param readFilters 
	 * @param readOrCreatePairedEndBam 
	 * @throws IOException 
	 */
	public ScanStatisticDataAlignmentModel(String bamFile, CoordinateSpace coordinateSpace, Collection<Predicate<Alignment>> readFilters, boolean readOrCreatePairedEndBam) {
		super(bamFile, coordinateSpace, readFilters, readOrCreatePairedEndBam);
	}

	
	public ScanStatisticDataAlignmentModel(String bamFile, CoordinateSpace coordinateSpace) {
		this(bamFile, coordinateSpace, new ArrayList<Predicate<Alignment>>());
	}

	public ScanStatisticDataAlignmentModel(String bamFile, CoordinateSpace coordinateSpace, boolean readOrCreatePairedEndBam) {
		this(bamFile, coordinateSpace, new ArrayList<Predicate<Alignment>>(), readOrCreatePairedEndBam);
	}

	
	
	/**
	 * Iterate through the whole coordinate space in windows and score them using the count function
	 * @param windowSize size of the windows to score
	 */
	public Iterator<ScanStatisticScore> scanAll(int windowSize, int overlap) {
		return super.scan(windowSize, overlap, getScanStatisticProcessor());
	}
	
	private ScanStatisticScore.Processor getScanStatisticProcessor() {
		return new ScanStatisticScore.Processor(this);
	}
	
	/**
	 * Scan across windows and return basic level statistics
	 * @param chr the chromosome to iterate over
	 * @param windowSize the window size
	 * @param overlap the overlap between the windows
	 * @return
	 */
	public Iterator<ScanStatisticScore> scan(String chr, int windowSize, int overlap) {
		return super.scan(getCoordinateSpace().getReferenceAnnotation(chr), windowSize, overlap, getScanStatisticProcessor());
	}
	
	public ScanStatisticScore scoreWindow(Annotation w) {
		return new ScanStatisticScore(this, w);
	}
	
	
	
	/**
	 * Whether the region is significantly expressed
	 * @param gene The annotation
	 * @param pValCutoff Scan P value cutoff
	 * @return Whether the scan P value of the annotation is significant
	 * @throws IOException
	 */
	public boolean isExpressed(Annotation gene, double pValCutoff) throws IOException {
		ScanStatisticScore score = new ScanStatisticScore(this, gene);
		return score.getScanPvalue() < pValCutoff;
	}	
}
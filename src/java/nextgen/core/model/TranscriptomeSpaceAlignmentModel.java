/**
 * 
 */
package nextgen.core.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.Window;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScore;
import nextgen.core.model.score.WindowScoreIterator;

import org.apache.commons.collections15.Predicate;

/**
 * @author prussell
 * Alignment model with special functionality for a transcriptome space
 */
public class TranscriptomeSpaceAlignmentModel extends ScanStatisticDataAlignmentModel {
	

	/**
	 * @param bamAlignmentFile
	 * @param transcriptomeSpace
	 * @param readFilterSet
	 * @param readOrCreatePairedEndBam
	 */
	public TranscriptomeSpaceAlignmentModel(String bamAlignmentFile, TranscriptomeSpace transcriptomeSpace, Collection<Predicate<Alignment>> readFilterSet, boolean readOrCreatePairedEndBam) {
		super(bamAlignmentFile, transcriptomeSpace, readFilterSet, readOrCreatePairedEndBam);
	}

	/**
	 * @param bamAlignmentFile
	 * @param transcriptomeSpace
	 * @param readFilterSet
	 */
	public TranscriptomeSpaceAlignmentModel(String bamAlignmentFile, TranscriptomeSpace transcriptomeSpace, Collection<Predicate<Alignment>> readFilterSet) {
		super(bamAlignmentFile, transcriptomeSpace, readFilterSet);
	}

	/**
	 * @param bamAlignmentFile
	 * @param transcriptomeSpace
	 * @param readOrCreatePairedEndBam
	 */
	public TranscriptomeSpaceAlignmentModel(String bamAlignmentFile, TranscriptomeSpace transcriptomeSpace, boolean readOrCreatePairedEndBam) {
		super(bamAlignmentFile, transcriptomeSpace, readOrCreatePairedEndBam);
	}

	/**
	 * @param bamAlignmentFile
	 * @param transcriptomeSpace
	 */
	public TranscriptomeSpaceAlignmentModel(String bamAlignmentFile, TranscriptomeSpace transcriptomeSpace) {
		super(bamAlignmentFile, transcriptomeSpace);
	}
	
	/**
	 * Get list of position level counts in mature transcript
	 * @param gene The region
	 * @return List of position counts
	 * @throws IOException 
	 */
	public List<Double> getPositionCountList(Gene gene) throws IOException {
		List<Double> rtrn = new ArrayList<Double>();
		WindowProcessor<CountScore> processor = new CountScore.Processor(this);
		WindowScoreIterator<CountScore> scoreIter = scan(gene, 1, 0, processor);
		while(scoreIter.hasNext()) {
			rtrn.add(Double.valueOf(scoreIter.next().getCount()));
		}
		return rtrn;
	}
	
	/**
	 * Scan windows over a gene and score
	 * Only scan exons of this gene
	 * @param gene The gene 
	 * @param windowSize Window size
	 * @param overlap Overlap size
	 * @param processor The window processor
	 * @return A score iterator over windows
	 */
	public <T extends WindowScore> WindowScoreIterator<T> scan(Gene gene, int windowSize, int overlap, WindowProcessor<T> processor) {
		Collection<Gene> baseGenes = new TreeSet<Gene>();
		baseGenes.add(gene);
		Iterator<? extends Window> windowIterator = coordinateSpace.getWindowIterator(baseGenes, windowSize, overlap);
		return new WindowScoreIterator<T>(windowIterator, processor, gene);
	}

	/**
	 * Get the names of genes in the transcriptome space
	 * @return The names of reference sequences (genes) in the transcriptome space
	 */
	public Collection<String> getReferenceNames() {
		return coordinateSpace.getReferenceNames();
	}


}

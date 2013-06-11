/**
 * 
 */
package nextgen.core.readFilters;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

import org.apache.commons.collections15.Predicate;

import broad.core.sequence.TranscribedSequenceComposition;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class FirstNucleotideFilter implements Predicate<Alignment> {
	
	private TranscribedSequenceComposition transcribedSequence;
	private char nucleotide;
	private int offset;
	private Map<String, Collection<Gene>> genes;
	
	/**
	 * @param bamFile Bam file of alignments
	 * @param geneBedFile Bed file of genes
	 * @param genomeFasta Genome fasta file
	 * @param transcriptionRead Transcription read
	 * @param firstNucleotide Nucleotide to require at first fragment position
	 * @param offsetAlongParent Offset along the parent transcript from the beginning position
	 * @throws IOException
	 */
	public FirstNucleotideFilter(String bamFile, String geneBedFile, String genomeFasta, TranscriptionRead transcriptionRead, char firstNucleotide, int offsetAlongParent) throws IOException {
		transcribedSequence = new TranscribedSequenceComposition(bamFile, geneBedFile, genomeFasta, transcriptionRead);
		nucleotide = firstNucleotide;
		offset = offsetAlongParent;
		genes = BEDFileParser.loadDataByChr(new File(geneBedFile));
	}
	
	@Override
	public boolean evaluate(Alignment align) {
		Gene alignGene = new Gene(align);
		Gene parentGene = AnnotationUtils.getLargestParent(alignGene, genes);
		char nuc = transcribedSequence.getTranscribedBaseAtBeginningPosition(parentGene, align, offset);
		return nuc == nucleotide;
	}

}

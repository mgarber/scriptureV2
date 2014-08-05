/**
 * 
 */
package nextgen.core.structure;

import general.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.alignment.FeatureSequenceAlignment;
import nextgen.core.alignment.FeatureSequenceAlignment.UnorderedGenePair;
import nextgen.core.annotation.Gene;
import nextgen.core.feature.GeneWindow;

/**
 * @author prussell
 *
 */
public class PairwiseWindowHybridization {

	private Gene gene;
	private int window;
	private int step;
	private FeatureSequenceAlignment alignments;
	private static Logger logger = Logger.getLogger(PairwiseWindowHybridization.class.getName());

	/**
	 * @param transcript The gene
	 * @param chrs Map of chromosome name to sequence
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @throws IOException
	 */
	public PairwiseWindowHybridization(Gene transcript, Map<String, broad.core.sequence.Sequence> chrs, int windowSize, int stepSize) throws IOException {
		this(transcript, chrs, windowSize, stepSize, FeatureSequenceAlignment.DEFAULT_MATCH_SCORE, FeatureSequenceAlignment.DEFAULT_MISMATCH_SCORE, FeatureSequenceAlignment.DEFAULT_GAP_OPEN_PENALTY, FeatureSequenceAlignment.DEFAULT_GAP_EXTEND_PENALTY);
	}
	
	
	/**
	 * @param transcript The gene
	 * @param genomeFasta Genome fasta file
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @throws IOException
	 */
	public PairwiseWindowHybridization(Gene transcript, String genomeFasta, int windowSize, int stepSize) throws IOException {
		this(transcript, FastaSequenceIO.getChrSequencesFromFasta(genomeFasta), windowSize, stepSize);
	}
	
	/**
	 * @param transcript The gene
	 * @param genomeFasta Fasta file of chromosomes
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @param matchScore Match score for Smith Waterman
	 * @param mismatchScore Mismatch score for Smith Waterman
	 * @param gapOpenPenalty Gap open penalty for Smith Waterman
	 * @param gapExtendPenalty Gap extend penalty for Smith Waterman
	 * @throws IOException
	 */
	public PairwiseWindowHybridization(Gene transcript, String genomeFasta, int windowSize, int stepSize, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) throws IOException {
		this(transcript, FastaSequenceIO.getChrSequencesFromFasta(genomeFasta), windowSize, stepSize, matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
	}
	
	/**
	 * @param transcript The gene
	 * @param chrs Map of chromosome name to sequence
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @param matchScore Match score for Smith Waterman
	 * @param mismatchScore Mismatch score for Smith Waterman
	 * @param gapOpenPenalty Gap open penalty for Smith Waterman
	 * @param gapExtendPenalty Gap extend penalty for Smith Waterman
	 * @throws IOException
	 */
	public PairwiseWindowHybridization(Gene transcript, Map<String, broad.core.sequence.Sequence> chrs, int windowSize, int stepSize, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) throws IOException {
		gene = transcript;
		window = windowSize;
		step = stepSize;
		Collection<GeneWindow> windows = gene.getWindows(window, step, 0);
		Collection<Gene> windowsAsGenes = new TreeSet<Gene>();
		for(GeneWindow w : windows) {
			windowsAsGenes.add(new Gene(w));
		}
		logger.info("Gene " + gene.getName() + " has " + windowsAsGenes.size() + " windows with window size " + window + " and step size " + step + ".");
		alignments = new FeatureSequenceAlignment(windowsAsGenes, chrs, matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
	}
	
	/**
	 * Get all pairs of windows that match in antisense direction
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @return Map of window pair to alignment object
	 */
	public Map<UnorderedGenePair, jaligner.Alignment> getAllHybridizingPairs(float minAlignLength, float minPctIdentity) {
		return alignments.getAllPairwiseAntisenseAlignments(minAlignLength, minPctIdentity);
	}
	
	/**
	 * Get the midpoint of the window in transcript coordinates of the gene
	 * @param region The window
	 * @return The window midpoint in terms of the gene
	 */
	private int getMidpointTranscriptCoords(Gene region) {
		int midpointGenomicCoords = region.getMidpointGenomicCoords();
		return gene.genomicToTranscriptPosition(midpointGenomicCoords);
	}
	
	/**
	 * Write all pairs of windows that match in antisense direction to bed file where matching windows are joined by an "intron"
	 * @param outBedFile Output bed file
	 * @param append Write to end of file rather than beginning
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @throws IOException
	 */
	public void writeAllHybridizingPairsToBed(String outBedFile, boolean append, int minAlignLength, float minPctIdentity) throws IOException {
		alignments.writeAllAntisenseAlignmentsToBed(outBedFile, append, minAlignLength, minPctIdentity);
	}

	
	/**
	 * Write all pairs of windows that match in antisense direction to bed file where matching windows are joined by an "intron"
	 * @param outBedFile Output bed file
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @throws IOException
	 */
	public void writeAllHybridizingPairsToBed(String outBedFile, int minAlignLength, float minPctIdentity) throws IOException {
		writeAllHybridizingPairsToBed(outBedFile, false, minAlignLength, minPctIdentity);
	}

	/**
	 * Write a table of pairs of coordinates that hybridize (midpoint of intervals in transcript coordinates)
	 * @param outTable File to write table to
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @param overlapperBed Bed file of overlapper annotations to check. Table will include a column for neither, one (x coord) or both genes overlap a region from bed file
	 * @throws IOException
	 */
	public void writeAllHybridizingPairsToTableWithOverlaps(String outTable, int minAlignLength, float minPctIdentity, String overlapperBed) throws IOException {
		logger.info("Writing midpoint coordinates of hybridizing pairs to file " + outTable);
		logger.info("Checking for overlaps with regions from file " + overlapperBed);
		Collection<Gene> toCheck = BEDFileParser.loadData(new File(overlapperBed));
		Map<UnorderedGenePair, jaligner.Alignment> hybs = getAllHybridizingPairs(minAlignLength, minPctIdentity);
		FileWriter w = new FileWriter(outTable);
		w.write("first_midpoint\tsecond_midpoint\twhich_has_overlap\n");
		for(UnorderedGenePair pair : hybs.keySet()) {
			Gene firstGene = pair.getFirstGene();
			Gene secondGene = pair.getSecondGene();
			boolean firstOverlaps = false;
			boolean secondOverlaps = false;
			for(Gene g : toCheck) {
				if(firstGene.overlaps(g)) {
					firstOverlaps = true;
				}
				if(secondGene.overlaps(g)) {
					secondOverlaps = true;
				}
			}
			boolean bothOverlap = firstOverlaps && secondOverlaps;
			boolean neitherOverlap = !firstOverlaps && !secondOverlaps;
			int firstMidpoint = getMidpointTranscriptCoords(firstGene);
			int secondMidpoint = getMidpointTranscriptCoords(secondGene);
			String firstFirst = firstMidpoint + "\t" + secondMidpoint + "\t";
			String secondFirst = secondMidpoint + "\t" + firstMidpoint + "\t";
			if(neitherOverlap) {
				firstFirst += "neither";
				secondFirst += "neither";
			} else if(bothOverlap) {
				firstFirst += "both";
				secondFirst += "both";
			} else {
				if(firstOverlaps) {
					firstFirst += "first";
					secondFirst += "second";
				}
				if(secondOverlaps) {
					firstFirst += "second";
					secondFirst += "first";
				}
			}
			w.write(firstFirst + "\n");
			w.write(secondFirst + "\n");
		}
		w.close();
		logger.info("Wrote file.");
	}
	
	/**
	 * Write a table of pairs of coordinates that hybridize (midpoint of intervals in transcript coordinates)
	 * @param outTable File to write table to
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @throws IOException
	 */
	public void writeAllHybridizingPairsToTable(String outTable, int minAlignLength, float minPctIdentity) throws IOException {
		logger.info("Writing midpoint coordinates of hybridizing pairs to file " + outTable);
		Map<UnorderedGenePair, jaligner.Alignment> hybs = getAllHybridizingPairs(minAlignLength, minPctIdentity);
		FileWriter w = new FileWriter(outTable);
		for(UnorderedGenePair pair : hybs.keySet()) {
			int firstMidpoint = getMidpointTranscriptCoords(pair.getFirstGene());
			int secondMidpoint = getMidpointTranscriptCoords(pair.getSecondGene());
			w.write(firstMidpoint + "\t" + secondMidpoint + "\n");
			w.write(secondMidpoint + "\t" + firstMidpoint + "\n");
		}
		w.close();
		logger.info("Wrote file.");
	}

	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of genes", true);
		p.addStringArg("-g", "Genome fasta file", true);
		p.addStringArg("-ob", "Output bed file", false, null);
		p.addStringArg("-ot", "Prefix for output table of coordinate pairs", false, null);
		p.addStringArg("-bo", "Bed file of regions to test for overlap. Requires -ot.", false, null);
		p.addIntArg("-w", "Window size", true);
		p.addIntArg("-s", "Step size", true);
		p.addFloatArg("-ma", "Match score for Smith Waterman", false, FeatureSequenceAlignment.DEFAULT_MATCH_SCORE);
		p.addFloatArg("-mm", "Mismatch score for Smith Waterman", false, FeatureSequenceAlignment.DEFAULT_MISMATCH_SCORE);
		p.addFloatArg("-go", "Gap open penalty for Smith Waterman", false, FeatureSequenceAlignment.DEFAULT_GAP_OPEN_PENALTY);
		p.addFloatArg("-ge", "Gap extend penalty for Smith Waterman", false, FeatureSequenceAlignment.DEFAULT_GAP_EXTEND_PENALTY);
		p.addIntArg("-ml", "Min alignment length", true);
		p.addFloatArg("-mp", "Min percent identity", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		Collection<Gene> genes = BEDFileParser.loadData(new File(bedFile));
		String genomeFasta = p.getStringArg("-g");
		String outBed = p.getStringArg("-ob");
		String outTable = p.getStringArg("-ot");
		Map<String, broad.core.sequence.Sequence> chromosomes = FastaSequenceIO.getChrSequencesFromFasta(genomeFasta);
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-s");
		float matchScore = p.getFloatArg("-ma");
		float mismatchScore = p.getFloatArg("-mm");
		float gapOpenPenalty = p.getFloatArg("-go");
		float gapExtendPenalty = p.getFloatArg("-ge");
		int minAlignLength = p.getIntArg("-ml");
		float minPctIdentity = p.getFloatArg("-mp");
		String bedOverlap = p.getStringArg("-bo");

		boolean first = true;
		for(Gene gene : genes) {
			PairwiseWindowHybridization pwh = new PairwiseWindowHybridization(gene, chromosomes, windowSize, stepSize, matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
			logger.info("Writing hybridizing pairs for gene " + gene.getName() + " to file " + outBed);
			if(outBed != null) {
				pwh.writeAllHybridizingPairsToBed(outBed, !first, minAlignLength, minPctIdentity);
			}
			if(outTable != null) {
				String outFile = outTable + "_" + gene.getName();
				if(bedOverlap == null) {
					pwh.writeAllHybridizingPairsToTable(outFile, minAlignLength, minPctIdentity);
				} else {
					pwh.writeAllHybridizingPairsToTableWithOverlaps(outFile, minAlignLength, minPctIdentity, bedOverlap);
				}
			}
			first = false;
		}
		logger.info("Wrote hybridizing pairs for all genes.");
		
	}

}

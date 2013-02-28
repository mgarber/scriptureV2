package nextgen.core.alignment;

import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;
import jaligner.SmithWatermanGotoh;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class FeatureSequenceAlignment {
	
	private Map<String, broad.core.sequence.Sequence> chromosomes;
	private Map<Gene, jaligner.Sequence> sequences;
	private TreeSet<UnorderedGenePair> featurePairs;
	private Matrix scoringMatrix;
	private float gapOpen;
	private float gapExtend;
	private static Logger logger = Logger.getLogger(FeatureSequenceAlignment.class.getName());
	private TreeMap<UnorderedGenePair, jaligner.Alignment> alignments;
	private TreeMap<UnorderedGenePair, jaligner.Alignment> antisenseAlignments;
	
	private static float DEFAULT_MATCH_SCORE = 5;
	private static float DEFAULT_MISMATCH_SCORE = -4;
	private static float DEFAULT_GAP_OPEN_PENALTY = 8;
	private static float DEFAULT_GAP_EXTEND_PENALTY = 2;
	
	/**
	 * @param featureBedFile Bed file of features
	 * @param genomeFasta Genome fasta file
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(String featureBedFile, String genomeFasta) throws IOException {
		this(featureBedFile, genomeFasta, DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE, DEFAULT_GAP_OPEN_PENALTY, DEFAULT_GAP_EXTEND_PENALTY);
	}
	
	/**
	 * @param genes Set of features
	 * @param genomeFasta Genome fasta file
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(Collection<Gene> genes, String genomeFasta) throws IOException {
		this(genes, genomeFasta, DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE, DEFAULT_GAP_OPEN_PENALTY, DEFAULT_GAP_EXTEND_PENALTY);
	}
	
	/**
	 * @param featureBedFile Bed file of features
	 * @param genomeFasta Genome fasta file
	 * @param matchScore Match score for Smith Waterman
	 * @param mismatchScore Mismatch score for Smith Waterman
	 * @param gapOpenPenalty Gap open penalty for Smith Waterman
	 * @param gapExtendPenalty Gap extend penalty for Smith Waterman
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(String featureBedFile, String genomeFasta, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) throws IOException {
		this(BEDFileParser.loadData(new File(featureBedFile)), genomeFasta, matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
	}

	
	/**
	 * @param genes Set of features
	 * @param genomeFasta Genome fasta file
	 * @param matchScore Match score for Smith Waterman
	 * @param mismatchScore Mismatch score for Smith Waterman
	 * @param gapOpenPenalty Gap open penalty for Smith Waterman
	 * @param gapExtendPenalty Gap extend penalty for Smith Waterman
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(Collection<Gene> genes, String genomeFasta, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) throws IOException {
		
		logger.info("Loading chromosomes from file " + genomeFasta + "...");
		FastaSequenceIO fsio = new FastaSequenceIO(genomeFasta);
		List<broad.core.sequence.Sequence> chrs = fsio.loadAll();
		chromosomes = new TreeMap<String, broad.core.sequence.Sequence>();
		for(broad.core.sequence.Sequence chr : chrs) {
			chromosomes.put(chr.getId(), chr);
		}
		logger.info("Loaded " + chromosomes.size() + " chromosomes.");
		
		// Store genes and sequences
		sequences = new TreeMap<Gene, jaligner.Sequence>();
		for(Gene gene : genes) {
			broad.core.sequence.Sequence broadSeq = chromosomes.get(gene.getChr()).getSubsequence(gene);
			jaligner.Sequence jSeq = new jaligner.Sequence(broadSeq.getSequenceBases());
			jSeq.setId(gene.getName());
			sequences.put(gene, jSeq);
		}
		logger.info("Loaded " + sequences.size() + " features with sequences.");
		
		// Store unordered pairs of genes
		featurePairs = new TreeSet<UnorderedGenePair>();
		for(Gene gene1 : genes) {
			for(Gene gene2 : genes) {
				if(gene1.equals(gene2)) continue;
				featurePairs.add(new UnorderedGenePair(gene1, gene2));
			}
		}
		logger.info("There are " + featurePairs.size() + " unordered pairs of features.");
		
		gapOpen = gapOpenPenalty;
		gapExtend = gapExtendPenalty;
		scoringMatrix = MatrixGenerator.generate(matchScore, mismatchScore);
		alignments = new TreeMap<UnorderedGenePair, jaligner.Alignment>();
		antisenseAlignments = new TreeMap<UnorderedGenePair, jaligner.Alignment>();
		
		logger.info("Match score: " + matchScore);
		logger.info("Mismatch score: " + mismatchScore);
		logger.info("Gap open penalty: " + gapOpen);
		logger.info("Gap extend penalty: " + gapExtend);
	
	}
	
	/**
	 * Get the alignment between two genes
	 * @param genes Gene pair
	 * @return The alignment of transcript sequences
	 */
	private jaligner.Alignment getAlignment(UnorderedGenePair genes) {
		if(!alignments.containsKey(genes)) {
			align(genes);
		}
		return alignments.get(genes);
	}
	
	/**
	 * Get the antisense alignment between two genes
	 * @param genes Gene pair
	 * @return The alignment of transcript sequences
	 */
	private jaligner.Alignment getAntisenseAlignment(UnorderedGenePair genes) {
		if(!antisenseAlignments.containsKey(genes)) {
			align(genes);
		}
		return antisenseAlignments.get(genes);
	}
	/**
	 * Align the sequences of two genes
	 * @param gene1 Gene 1
	 * @param gene2 Gene 2
	 * @return Smith Waterman alignment
	 */
	private jaligner.Alignment align(UnorderedGenePair genes) {
		Gene gene1 = genes.getFirstGene();
		Gene gene2 = genes.getSecondGene();
		jaligner.Sequence seq1 = sequences.get(gene1);
		jaligner.Sequence seq2 = sequences.get(gene2);
		String seq2bases = seq2.getSequence();
		String seq2antisenseBases = Sequence.reverseSequence(seq2bases);
		jaligner.Sequence seq2antisense = new jaligner.Sequence(seq2antisenseBases);
		seq2antisense.setId(seq2.getId() + "_antisense");
		jaligner.Alignment align = SmithWatermanGotoh.align(seq1, seq2, scoringMatrix, gapOpen, gapExtend);
		jaligner.Alignment antisenseAlign = SmithWatermanGotoh.align(seq1, seq2antisense, scoringMatrix, gapOpen, gapExtend);
		int correctedStart2 = seq2antisense.length() - antisenseAlign.getStart2() - antisenseAlign.getLength() + antisenseAlign.getGaps2();
		antisenseAlign.setStart2(correctedStart2);
		alignments.put(genes, align);
		antisenseAlignments.put(genes, antisenseAlign);
		return align;
	}
		
	/**
	 * String of information about the alignment
	 * @param align The alignment
	 * @return Formatted string for printing
	 */
	private static String formattedInfoString(jaligner.Alignment align) {
		String seq1 = new String(align.getSequence1());
		String seq2 = new String(align.getSequence2());
		String markup = new String(align.getMarkupLine());
		String rtrn = align.getSummary();
		rtrn += seq1 + "\n";
		rtrn += markup + "\n";
		rtrn += seq2 + "\n";
		return rtrn;
	}
	
	
	/**
	 * Write all pairwise alignments to file
	 * @param outFile Output file
	 * @throws IOException
	 */
	private void writeAllPairwiseAlignments(String outFile) throws IOException {
		writeAllPairwiseAlignments(outFile, 0, 0);
	}
	
	/**
	 * Get all sense direction alignments passing thresholds
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @return All sense alignments passing criteria
	 */
	public Collection<jaligner.Alignment> getAllPairwiseSenseAlignments(float minAlignLength, float minPctIdentity) {
		Collection<jaligner.Alignment> rtrn = new ArrayList<jaligner.Alignment>();
		for(UnorderedGenePair genes : featurePairs) {
			jaligner.Alignment align = getAlignment(genes);
			if(align.getLength() >= minAlignLength && align.getPercentIdentity() >= minPctIdentity) {
				rtrn.add(align);
			}
		}
		return rtrn;
	}
	
	/**
	 * Get all antisense direction alignments passing thresholds
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @return All antisense alignments passing criteria
	 */
	public Collection<jaligner.Alignment> getAllPairwiseAntisenseAlignments(float minAlignLength, float minPctIdentity) {
		Collection<jaligner.Alignment> rtrn = new ArrayList<jaligner.Alignment>();
		for(UnorderedGenePair genes : featurePairs) {
			jaligner.Alignment align = getAntisenseAlignment(genes);
			if(align.getLength() >= minAlignLength && align.getPercentIdentity() >= minPctIdentity) {
				rtrn.add(align);
			}
		}
		return rtrn;
	}
	
	/**
	 * Write all pairwise alignments (sense direction only) to bed file in genome coordinates
	 * @param outBedFile Output bed file
	 * @param minAlignLength Min alignment length to keep
	 * @param minPctIdentity Min percent identity to keep
	 * @throws IOException
	 */
	private void writeAllSenseAlignmentsToBed(String outBedFile, int minAlignLength, float minPctIdentity) throws IOException {
		logger.info("Writing all sense alignments to bed file " + outBedFile);
		logger.info("Min alignment length: " + minAlignLength);
		logger.info("Min percent identity: " + minPctIdentity);
		FileWriter w = new FileWriter(outBedFile);
		for(UnorderedGenePair genes : featurePairs) {
			jaligner.Alignment align = getAlignment(genes);
			if(align.getLength() >= minAlignLength && align.getPercentIdentity() >= minPctIdentity) {
				try {
					String print = alignmentAsBed(genes);
					w.write(print + "\n");
				} catch (Exception e) {
					logger.error("Caught exception when getting alignment for regions " + genes.getFirstGene().getName() + " and " + genes.getSecondGene().getName());
					logger.error("Alignment:");
					logger.error(formattedInfoString(align));
					e.printStackTrace();
				}				
			}
		}
		w.close();
		logger.info("Done writing sense alignments.");
	}
	
	/**
	 * Write all pairwise alignments (antisense direction only) to bed file in genome coordinates
	 * @param outBedFile Output bed file
	 * @param minAlignLength Min alignment length to keep
	 * @param minPctIdentity Min percent identity to keep
	 * @throws IOException
	 */
	private void writeAllAntisenseAlignmentsToBed(String outBedFile, int minAlignLength, float minPctIdentity) throws IOException {
		logger.info("Writing all antisense alignments to bed file " + outBedFile);
		logger.info("Min alignment length: " + minAlignLength);
		logger.info("Min percent identity: " + minPctIdentity);
		FileWriter w = new FileWriter(outBedFile);
		for(UnorderedGenePair genes : featurePairs) {
			jaligner.Alignment antisenseAlign = getAntisenseAlignment(genes);
			if(antisenseAlign.getLength() >= minAlignLength && antisenseAlign.getPercentIdentity() >= minPctIdentity) {
       				try {
					String print = antisenseAlignmentAsBed(genes);
					w.write(print + "\n");				
				} catch (Exception e) {
					logger.error("Caught exception when getting alignment for regions " + genes.getFirstGene().getName() + " and " + genes.getSecondGene().getName());
					logger.error("Alignment:");
					logger.error(formattedInfoString(antisenseAlign));
					e.printStackTrace();
				}
			}
		}
		w.close();
		logger.info("Done writing antisense alignments.");
	}

	/**
	 * Write all pairwise alignments to file
	 * @param outFile Output file
	 * @param minAlignLength Min alignment length to keep
	 * @param minPctIdentity Min percent identity to keep
	 * @throws IOException
	 */
	private void writeAllPairwiseAlignments(String outFile, int minAlignLength, float minPctIdentity) throws IOException {
		logger.info("Writing all alignments to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		for(UnorderedGenePair genes : featurePairs) {
			jaligner.Alignment align = getAlignment(genes);
			if(align.getLength() >= minAlignLength && align.getPercentIdentity() >= minPctIdentity) {
				String print = formattedInfoString(align);
				w.write(print + "\n");				
			}
			jaligner.Alignment antisenseAlign = getAntisenseAlignment(genes);
			if(antisenseAlign.getLength() >= minAlignLength && antisenseAlign.getPercentIdentity() >= minPctIdentity) {
				String print = formattedInfoString(antisenseAlign);
				w.write(print + "\n");				
			}
		}
		w.close();
		logger.info("Done writing alignments.");
	}
	
	/**
	 * Get the aligned region as a feature in genomic coordinates
	 * @param genes The pair of genes
	 * @return The aligned region
	 */
	private Annotation alignmentAsFeature(UnorderedGenePair genes) {
		return asFeature(getAlignment(genes), genes);
	}
	
	/**
	 * Get the aligned region of the antisense alignment as a feature in genomic coordinates
	 * @param genes The pair of genes
	 * @return The aligned region
	 */
	private Annotation antisenseAlignmentAsFeature(UnorderedGenePair genes) {
		return asFeature(getAntisenseAlignment(genes), genes);
	}

	/**
	 * The alignment as a feature in genomic coordinates
	 * @param align The alignment
	 * @param genes The pair of genes that were aligned
	 * @return The aligned region
	 */
	private static Annotation asFeature(jaligner.Alignment align, UnorderedGenePair genes) {
		Gene gene1 = genes.getFirstGene();
		int gene1start = align.getStart1();
		int gene1end = gene1start + align.getLength() - align.getGaps1() - 1;
		Gene alignedRegion1 = gene1.transcriptToGenomicPosition(gene1start, gene1end);
		int newStart1 = alignedRegion1.isNegativeStrand() ? alignedRegion1.getStart() + 1 : alignedRegion1.getStart() - 1;
		alignedRegion1.setStart(newStart1);
		Gene gene2 = genes.getSecondGene();
		int gene2start = align.getStart2();
		int gene2end = gene2start + align.getLength() - align.getGaps2() - 1;
		Gene alignedRegion2 = gene2.transcriptToGenomicPosition(gene2start, gene2end);
		int newStart2 = alignedRegion2.isNegativeStrand() ? alignedRegion2.getStart() + 1 : alignedRegion2.getStart() - 1;
		alignedRegion2.setStart(newStart2);
		BasicAnnotation rtrn = new BasicAnnotation(alignedRegion1);
		rtrn.addBlocks(alignedRegion2);
		rtrn.setName(gene1.getChr() + ":" + gene1.getStart() + "-" + gene1.getEnd() + ":" + gene2.getStart() + "-" + gene2.getEnd());
		return rtrn;
	}
	
	/**
	 * Get the aligned region as a bed line in genomic coordinates
	 * @param genes The pair of genes
	 * @return
	 */
	private String alignmentAsBed(UnorderedGenePair genes) {
		Annotation alignedRegion = alignmentAsFeature(genes);
		return alignedRegion.toBED(41,144,41);
	}
	
	/**
	 * Get the aligned region of the antisense alignment as a bed line in genomic coordinates
	 * @param genes The pair of genes
	 * @return
	 */
	private String antisenseAlignmentAsBed(UnorderedGenePair genes) {
		Annotation alignedRegion = antisenseAlignmentAsFeature(genes);
		return alignedRegion.toBED(144,59,144);
	}
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of features", true);
		p.addStringArg("-g", "Genome fasta file", true);
		p.addStringArg("-o", "Output file for all pairwise alignments", false, null);
		p.addStringArg("-os", "Output bed file for sense alignments", false, null);
		p.addStringArg("-oa", "Output bed file for antisense alignments", false, null);
		p.addFloatArg("-ma", "Match score", false, DEFAULT_MATCH_SCORE);
		p.addFloatArg("-mi", "Mismatch score", false, DEFAULT_MISMATCH_SCORE);
		p.addFloatArg("-go", "Gap open penalty", false, DEFAULT_GAP_OPEN_PENALTY);
		p.addFloatArg("-ge", "Gap extend penalty", false, DEFAULT_GAP_EXTEND_PENALTY);
		p.addIntArg("-l", "Min alignment length", false, 0);
		p.addFloatArg("-id", "Min alignment idenity", false, 0);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String fastaFile = p.getStringArg("-g");
		String outPairwiseFile = p.getStringArg("-o");
		String outBedFileSense = p.getStringArg("-os");
		String outBedFileAntisense = p.getStringArg("-oa");
		float matchScore = p.getFloatArg("-ma");
		float mismatchScore = p.getFloatArg("-mi");
		float gapOpen = p.getFloatArg("-go");
		float gapExtend = p.getFloatArg("-ge");
		int minLength = p.getIntArg("-l");
		float minIdentity = p.getFloatArg("-id");
		
		
		FeatureSequenceAlignment fsa = new FeatureSequenceAlignment(bedFile, fastaFile, matchScore, mismatchScore, gapOpen, gapExtend);
		
		if(outPairwiseFile != null) {
			fsa.writeAllPairwiseAlignments(outPairwiseFile, minLength, minIdentity);
		}
		
		if(outBedFileSense != null) {
			fsa.writeAllSenseAlignmentsToBed(outBedFileSense, minLength, minIdentity);
		}
		
		if(outBedFileAntisense != null) {
			fsa.writeAllAntisenseAlignmentsToBed(outBedFileAntisense, minLength, minIdentity);
		}
		
	}
	
	
	/**
	 * An unordered pair of genes
	 * @author prussell
	 *
	 */
	private class UnorderedGenePair implements Comparable<UnorderedGenePair> {
		
		private Gene firstGene;
		private Gene secondGene;
		
		/**
		 * Determine which gene is first and which is second
		 * @param gene1 Gene 1
		 * @param gene2 Gene 2
		 */
		public UnorderedGenePair(Gene gene1, Gene gene2) {
			if(gene1.compareTo(gene2) < 0) {
				firstGene = gene1;
				secondGene = gene2;
			}
			else {
				firstGene = gene2;
				secondGene = gene1;
			}
		}
		
		public Gene getFirstGene() {
			return firstGene;
		}
		
		public Gene getSecondGene() {
			return secondGene;
		}

		@Override
		public int compareTo(UnorderedGenePair o) {
			if(firstGene.compareTo(o.getFirstGene()) == 0) {
				return secondGene.compareTo(o.getSecondGene());
			}
			return firstGene.compareTo(o.getFirstGene());
		}
		
		
		@Override
		public boolean equals(Object o) {
			UnorderedGenePair p = (UnorderedGenePair) o;
			return firstGene.equals(p.getFirstGene()) && secondGene.equals(p.getSecondGene());
		}

		@Override
		public int hashCode() {
			String str = firstGene.toString() + secondGene.toString();
			return str.hashCode();
		}
		
	}
	
	
}

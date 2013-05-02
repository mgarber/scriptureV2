/**
 * 
 */
package broad.pda.seq.protection;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.readers.WigReader;

/**
 * @author prussell
 *
 */
public class WigSequenceComposition {

	private WigReader wigReader;
	private Map<String, Map<Integer,Double>> wigData;
	private Map<String, Sequence> chromosomes;
	private Map<String, Collection<Gene>> genes;
	private Map<String, Map<HalfOpenInterval, Strand>> geneStrands;
	private static Logger logger = Logger.getLogger(WigSequenceComposition.class.getName());
	private Map<String, Double> baseCounts;
	
	/**
	 * 
	 * @param wigFile Wig file
	 * @param genomeFasta Genome fasta file
	 * @param bedFile Bed file of genes to get orientations from
	 * @throws IOException
	 */
	private WigSequenceComposition(String wigFile, String genomeFasta, String bedFile) throws IOException {
		// Load chromosomes
		FastaSequenceIO fsio = new FastaSequenceIO(genomeFasta);
		List<Sequence> chrs = fsio.loadAll();
		chromosomes = new TreeMap<String, Sequence>();
		for(Sequence chr : chrs) {
			chromosomes.put(chr.getId(), chr);
		}
		// Load genes
		logger.info("Loading genes from file " + bedFile + "...");
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		int numGenes = 0;
		for(String chr : genes.keySet()) {
			numGenes += genes.get(chr).size();
		}
		logger.info("Loaded " + numGenes + " genes.");
		// Collapse overlapping exons and determine strands
		computeStrands();
		// Read wig file
		wigReader = new WigReader(wigFile);
		wigData = wigReader.getAllValues();
		// Compute counts for each nucleotide
		computeBaseCounts();
	}
	
	/**
	 * Collapse overlapping exons and save strand for each collapsed interval
	 */
	private void computeStrands() {
		logger.info("Collapsing overlapping exons and identifying strand...");
		geneStrands = new TreeMap<String, Map<HalfOpenInterval,Strand>>();
		for(String chr : genes.keySet()) {
			logger.info(chr);
			Map<HalfOpenInterval, Strand> chrStrands = new TreeMap<HalfOpenInterval, Strand>();
			for(Gene gene : genes.get(chr)) {
				Strand strand = gene.getOrientation();
				HalfOpenInterval nextInterval = new HalfOpenInterval(gene.getStart(), gene.getEnd());
				boolean merged = false;
				for(HalfOpenInterval existingInterval : chrStrands.keySet()) {
					if(nextInterval.overlaps(existingInterval)) {
						existingInterval.merge(nextInterval);
						if(!chrStrands.get(existingInterval).equals(strand)) {
							chrStrands.put(existingInterval, Strand.UNKNOWN);
						}
						merged = true;
					}
				}
				if(!merged) {
					chrStrands.put(nextInterval, strand);
				}
			}
			geneStrands.put(chr, chrStrands);
		}
		logger.info("Done identifying intervals and strands.");
	}
	
	/**
	 * Write overlapping intervals to bed file for QC
	 * @param bedFile Output bed file
	 * @throws IOException
	 */
	private void writeOverlappingIntervalsToBed(String bedFile) throws IOException {
		FileWriter w = new FileWriter(bedFile);
		for(String chr : geneStrands.keySet()) {
			for(HalfOpenInterval interval : geneStrands.get(chr).keySet()) {
				Strand strand = geneStrands.get(chr).get(interval);
				Annotation a = new BasicAnnotation(chr, (int)interval.getStart(), (int)interval.getEnd());
				a.setOrientation(strand);
				w.write(a.toBED() + "\n");
			}
		}
		w.close();
	}
	
	/**
	 * Get the transcribed orientation of the position based on the gene set
	 * Returns unknown if the position does not overlap a gene or strand is ambiguous
	 * @param chr Chromosome
	 * @param pos Position
	 * @return The orientation
	 */
	private Strand getOrientation(String chr, int pos) {
		for(HalfOpenInterval interval : geneStrands.get(chr).keySet()) {
			if(interval.contains(pos)) {
				return geneStrands.get(chr).get(interval);
			}
		}
		return Strand.UNKNOWN;
	}
	
	/**
	 * Get the oriented transcribed base at the position based on orientation of overlapping exon
	 * @param chr Chromosome
	 * @param pos Position
	 * @return The transcribed base
	 */
	private char getTranscribedBase(String chr, int pos) {
		Strand strand = getOrientation(chr, pos);
		if(strand.equals(Strand.UNKNOWN)) {
			return 'N';
		}
		Sequence subSeq = chromosomes.get(chr).getSubSequence("", pos, pos+1);
		if(strand.equals(Strand.NEGATIVE)) {
			subSeq.reverse();
		}
		return subSeq.getSequenceBases().charAt(0);
	}
	
	/**
	 * Write base counts to table
	 * @param outFile Output file for table
	 * @throws IOException
	 */
	private void writeBaseCountsToTable(String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		double total = 0;
		for(String base : baseCounts.keySet()) {
			total += baseCounts.get(base).doubleValue();
		}
		for(String base : baseCounts.keySet()) {
			double count = baseCounts.get(base).doubleValue();
			w.write(base + "\t" + count + "\t" + count / total + "\n");
		}
		w.close();
	}
	
	/**
	 * Compute counts for each base
	 * For each position in wig file, first identify the transcribed base at that position using gene annotation
	 * Next increment count for that base by adding the value from the wig file at the position
	 * @throws IOException 
	 */
	private void computeBaseCounts() throws IOException {
		logger.info("Computing base counts...");
		double[] counts = new double[10];
		for(int i = 0; i < counts.length; i++) {
			counts[i] = 0;
		}
		for(String chr : wigData.keySet()) {
			logger.info(chr);
			for(Integer pos : wigData.get(chr).keySet()) {
				double value = wigData.get(chr).get(pos).doubleValue();
				char base = getTranscribedBase(chr, pos.intValue());
				switch(base) {
					case 'A':  counts[Sequence.SHORT_ENCODED_A] += value;
					break;
					case 'C':  counts[Sequence.SHORT_ENCODED_C] += value;
					break;
					case 'G':  counts[Sequence.SHORT_ENCODED_G] += value;
					break;
					case 'T':  counts[Sequence.SHORT_ENCODED_T] += value;
					break;
					case 'a':  counts[Sequence.SHORT_ENCODED_a] += value;
					break;
					case 'c':  counts[Sequence.SHORT_ENCODED_c] += value;
					break;
					case 'g':  counts[Sequence.SHORT_ENCODED_g] += value;
					break;
					case 't':  counts[Sequence.SHORT_ENCODED_t] += value;
					break;
					case 'N':  counts[Sequence.SHORT_ENCODED_N] += value;
					break;
					default: throw new IllegalArgumentException("Base " + base + " not valid.");
				}
			}
		}
		baseCounts = new TreeMap<String, Double>();
		baseCounts.put("A", Double.valueOf(counts[Sequence.SHORT_ENCODED_A]));
		baseCounts.put("C", Double.valueOf(counts[Sequence.SHORT_ENCODED_C]));
		baseCounts.put("G", Double.valueOf(counts[Sequence.SHORT_ENCODED_G]));
		baseCounts.put("T", Double.valueOf(counts[Sequence.SHORT_ENCODED_T]));
		baseCounts.put("a", Double.valueOf(counts[Sequence.SHORT_ENCODED_a]));
		baseCounts.put("c", Double.valueOf(counts[Sequence.SHORT_ENCODED_c]));
		baseCounts.put("g", Double.valueOf(counts[Sequence.SHORT_ENCODED_g]));
		baseCounts.put("t", Double.valueOf(counts[Sequence.SHORT_ENCODED_t]));
		baseCounts.put("N", Double.valueOf(counts[Sequence.SHORT_ENCODED_N]));
		logger.info("Done computing base counts.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Genome fasta", true);
		p.addStringArg("-b", "Bed annotation", true);
		p.addStringArg("-w", "Wig file", true);
		p.addStringArg("-qe", "Output qc bed file for exon overlaps and orientations", false, null);
		p.addStringArg("-o", "Output table file", true);
		p.parse(args);
		String genomeFasta = p.getStringArg("-g");
		String bedFile = p.getStringArg("-b");
		String wigFile = p.getStringArg("-w");
		String qcBedFile = p.getStringArg("-qe");
		String outTable = p.getStringArg("-o");
		
		WigSequenceComposition wsc = new WigSequenceComposition(wigFile, genomeFasta, bedFile);
		wsc.writeBaseCountsToTable(outTable);
		if(qcBedFile != null) {
			wsc.writeOverlappingIntervalsToBed(qcBedFile);
		}
		
	}
	
	
	private class HalfOpenInterval implements Comparable<HalfOpenInterval> {

		private double start;
		private double end;
		
		/**
		 * @param startPos Start position (inclusive)
		 * @param endPos End position (exclusive)
		 */
		public HalfOpenInterval(double startPos, double endPos) {
			start = startPos;
			end = endPos;
		}
		
		/**
		 * Get start position
		 * @return Start position (inclusive)
		 */
		public double getStart() {
			return start;
		}
		
		/**
		 * Get end position
		 * @return End position (exclusive)
		 */
		public double getEnd() {
			return end;
		}
		
		/**
		 * Set start
		 * @param s New start
		 */
		public void setStart(double s) {
			start = s;
		}
		
		/**
		 * Set end
		 * @param e New end
		 */
		public void setEnd(double e) {
			end = e;
		}
		
		/**
		 * Whether this interval overlaps other interval
		 * @param other Other interval
		 * @return True iff the intervals share an inclusive position
		 */
		public boolean overlaps(HalfOpenInterval other) {
			if(start < other.getStart()) {
				return end > other.getStart();
			}
			return start < other.getEnd();
		}
		
		/**
		 * Merge with the other interval if they overlap
		 * Reset start and end positions to the union of the intervals
		 * If they do not overlap, do nothing
		 * @param other Other interval to merge
		 */
		public void merge(HalfOpenInterval other) {
			if(!overlaps(other)) {
				return;
			}
			setStart(Math.min(start, other.getStart()));
			setEnd(Math.max(end, other.getEnd()));
		}
		
		/**
		 * Whether the interval contains the number
		 * @param i Number to check
		 * @return True iff i is strictly contained in the interval
		 */
		public boolean contains(int i) {
			return i >= start && i < end;
		}
		
		/**
		 * Whether the interval contains the number
		 * @param d Number to check
		 * @return True iff d is strictly contained in the interval
		 */
		public boolean contains(double d) {
			return d >= start && d < end;
		}
		
		@Override
		public boolean equals(Object o) {
			HalfOpenInterval oi = (HalfOpenInterval)o;
			return start == oi.getStart() && end == oi.getEnd();
		}
		
		@Override
		public String toString() {
			String s = Double.valueOf(start).toString();
			String e = Double.valueOf(end).toString();
			return s + "-" + e;
		}
		
		@Override
		public int hashCode() {
			return toString().hashCode();
		}
		
		@Override
		public int compareTo(HalfOpenInterval other) {
			if(start == other.getStart()) {
				return 0;
			}
			if(start < other.getStart()) {
				return -1;
			}
			return 1;
		}

	}
	

}

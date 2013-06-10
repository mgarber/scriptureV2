/**
 * 
 */
package broad.core.sequence;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class TranscribedRegions {
	
	private Map<String, Sequence> chromosomes;
	private Map<String, Collection<Gene>> genes;
	private Map<String, Map<HalfOpenInterval, Strand>> geneStrands;
	private static Logger logger = Logger.getLogger(TranscribedRegions.class.getName());
	
	/**
	 * @param genomeFasta Genome fasta file
	 * @param bedFile Bed gene annotation
	 * @throws IOException
	 */
	public TranscribedRegions(String genomeFasta, String bedFile) throws IOException {
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
	}

	/**
	 * Get the transcribed orientation of the position based on the gene set
	 * Returns unknown if the position does not overlap a gene or strand is ambiguous
	 * @param chr Chromosome
	 * @param pos Position
	 * @return The orientation
	 */
	public Strand getOrientation(String chr, int pos) {
		for(HalfOpenInterval interval : geneStrands.get(chr).keySet()) {
			if(interval.contains(pos)) {
				return geneStrands.get(chr).get(interval);
			}
		}
		return Strand.UNKNOWN;
	}

	/**
	 * Get the oriented transcribed base at the position based on orientation provided
	 * or based on orientation of overlapping exon if orientation provided is null
	 * @param chr Chromosome
	 * @param pos Position
	 * @param strand Strand or null if get from annotation
	 * @return The transcribed base
	 */
	public char getTranscribedBase(String chr, int pos, Strand strand) {
		Strand s = strand == null ? getOrientation(chr, pos) : strand;
		if(s.equals(Strand.UNKNOWN)) {
			return 'N';
		}
		Sequence subSeq = chromosomes.get(chr).getSubSequence("", pos, pos+1);
		if(s.equals(Strand.NEGATIVE)) {
			subSeq.reverse();
		}
		return subSeq.getSequenceBases().toUpperCase().charAt(0);		
	}
	
	/**
	 * Get the oriented transcribed base at the position based on orientation of overlapping exon
	 * @param chr Chromosome
	 * @param pos Position
	 * @return The transcribed base
	 */
	public char getTranscribedBase(String chr, int pos) {
		return getTranscribedBase(chr, pos, null);
	}
	

	
	/**
	 * Write overlapping intervals to bed file for QC
	 * @param bedFile Output bed file
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
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
		@SuppressWarnings("unused")
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

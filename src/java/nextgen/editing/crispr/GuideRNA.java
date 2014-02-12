package nextgen.editing.crispr;

import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import broad.core.sequence.Sequence;

/**
 * A guide RNA with a genomic position and a target gene for CRISPR editing
 * @author prussell
 */
public class GuideRNA {
	
	private Sequence chr;
	private Annotation annot;
	private String name;
	public static Logger logger = Logger.getLogger(GuideRNA.class.getName());
	private Sequence sequence20;
	private Gene target;

	
	/**
	 * @param targetGene Target gene
	 * @param chromosome Chromosome
	 * @param start Start position of 20mer
	 * @param end Position after last position of 20mer
	 * @param orientation Strand
	 */
	public GuideRNA(Gene targetGene, Sequence chromosome, int start, int end, Strand orientation) {
		validateStartEnd(start, end);
		validateStrand(orientation);
		chr = chromosome;
		name = chr.getId() + ":" + start + "-" + end + ":" + orientation.toString();
		annot = new BasicAnnotation(chromosome.getId(), start, end, orientation, name);
		target = targetGene;
		
		Sequence sequence23;
		int start20;
		int start23;
		int end20;
		int end23;
		Strand strand;
		
		// Validate
		strand = orientation;
		if(strand.equals(Strand.POSITIVE)) {
			start20 = annot.getStart();
			start23 = start20;
			end20 = annot.getEnd();
			end23 = end20 + 3;
		} else {
			start20 = annot.getStart();
			start23 = annot.getStart() - 3;
			end20 = annot.getEnd();
			end23 = annot.getEnd();
		}
		Sequence shortSeq = chromosome.getSubSequence("", start20, end20);
		shortSeq.setSequenceBases(shortSeq.getSequenceBases().toUpperCase());
		Sequence longSeq = chromosome.getSubSequence("", start23, end23);
		longSeq.setSequenceBases(longSeq.getSequenceBases().toUpperCase());
		sequence20 = strand.equals(Strand.POSITIVE) ? shortSeq : Sequence.reverseSequence(shortSeq);
		sequence23 = strand.equals(Strand.POSITIVE) ? longSeq : Sequence.reverseSequence(longSeq);
		//logger.debug(name + "\tsequence20\t" + sequence20.getSequenceBases() + "\tsequence23\t" + sequence23.getSequenceBases());
		validateSequence20(sequence20);
		validateSequence23(sequence23);
	}

	/**
	 * @return The location of the guide RNA as a genomic annotation
	 */
	public Annotation getAnnotation() {
		return annot;
	}
	
	public String getChr() {
		return chr.getId();
	}
	
	public int getStart() {
		return annot.getStart();
	}
	
	public int getEnd() {
		return annot.getEnd();
	}
	
	public Strand getStrand() {
		return annot.getOrientation();
	}
	
	public boolean isPlusStrand() {
		return getStrand().equals(Strand.POSITIVE);
	}
	
	public boolean isMinusStrand() {
		return getStrand().equals(Strand.NEGATIVE);
	}
	
	public String getSequenceString() {
		return sequence20.getSequenceBases();
	}
	
	public Sequence getSequence() {
		return sequence20;
	}
	
	public String getName() {
		return name;
	}
	
	public Gene getTargetGene() {
		return target;
	}
	
	public void setName(String newName) {
		name = newName;
		annot.setName(name);
	}
	
	/**
	 * Find all guide RNAs on either strand within the window
	 * @param chr Chromosome
	 * @param start Start position of window to look in
	 * @param end Position after last position of window
	 * @return All guide RNAs followed by NGG whose 20nt sequence is fully contained in the window
	 */
	public static Collection<GuideRNA> findAll(Sequence chr, int start, int end, Gene targetGene) {
		//logger.debug("");
		Collection<GuideRNA> rtrn = new ArrayList<GuideRNA>();
		Collection<Annotation> nggs = findAllNGGs(chr, start, end);
		//logger.debug("There are " + nggs.size() + " NGGs in " + chr.getId() + ":" + start + "-" + end);
		for(Annotation ngg : nggs) {
			GuideRNA g = adjacentGuideRNA(chr, start, end, ngg, targetGene);
			if(g != null) {
				//logger.debug("Added " + g.toString());
				rtrn.add(g);
			}
		}
		return rtrn;
	}
	
	public String toString() {
		return chr.getId() + ":" + getStart() + "-" + getEnd() + ":" + getStrand().toString() + ":" + sequence20.getSequenceBases();
	}
	
	/**
	 * Get the correctly oriented guide RNA that is fully contained in the window and ends with the NGG specified
	 * @param windowChr Window chromosome
	 * @param windowStart Window start
	 * @param windowEnd Position after last position of window
	 * @param ngg Oriented NGG location
	 * @param targetGene Target gene
	 * @return The 20nt guide RNA or null if not fully contained in window
	 */
	private static GuideRNA adjacentGuideRNA(Sequence windowChr, int windowStart, int windowEnd, Annotation ngg, Gene targetGene) {
		validateNGG(windowChr, ngg);
		//logger.debug("Getting guide RNA adjacent to " + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
		if(ngg.getOrientation().equals(Strand.POSITIVE)) {
			if(ngg.getStart() - windowStart < 20) {
				// Guide RNA cannot be fully contained in window
				logger.debug("Guide RNA neighboring " + ngg.toUCSC() +" not fully contained in " + windowChr.getId() + ":" + windowStart + "-" + windowEnd);
				return null;
			}
			GuideRNA rtrn = new GuideRNA(targetGene, windowChr, ngg.getStart() - 20, ngg.getStart(), Strand.POSITIVE);
			//logger.debug("Guide RNA neighboring " + ngg.toUCSC() + " is " + rtrn.toString());
			return rtrn;
		}
		if(ngg.getOrientation().equals(Strand.NEGATIVE)) {
			if(ngg.getEnd() + 20 > windowEnd) {
				// Guide RNA cannot be fully contained in window
				logger.debug("Guide RNA neighboring " + ngg.toUCSC() +" not fully contained in " + windowChr.getId() + ":" + windowStart + "-" + windowEnd);
				return null;
			}
			GuideRNA rtrn = new GuideRNA(targetGene, windowChr, ngg.getEnd(), ngg.getEnd() + 20, Strand.NEGATIVE);
			//logger.debug("Guide RNA neighboring " + ngg.toUCSC() + " is " + rtrn.toString());
			return rtrn;
		}
		throw new IllegalArgumentException("Strand must be known");
	}
	
	private static void validateNGG(Sequence chr, Annotation ngg) {
		if(ngg.numBlocks() != 1) {
			throw new IllegalArgumentException("Must have one block");
		}
		Sequence seq = chr.getSubsequence(ngg);
		if(!seq.getSequenceBases().substring(1).toUpperCase().equals("GG")) {
			throw new IllegalArgumentException("Sequence must be NGG. Is " + seq.getSequenceBases());
		}
	}
	
	/**
	 * Get all NGG sequences that are fully contained within the window on either strand
	 * @param chr Chromosome
	 * @param start First position of window
	 * @param end Position after last position of window
	 * @return All NGGs as annotations with strand
	 */
	private static Collection<Annotation> findAllNGGs(Sequence chr, int start, int end) {
		String seq = chr.getSubSequence("", start, end).getSequenceBases();
		Collection<Annotation> rtrn = new TreeSet<Annotation>();
		
		// Find all NGGs
		int i1 = 0;
		while(i1 < seq.length()) {
			int pos = seq.indexOf("GG", i1);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i1++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i1 = pos + 1;
		}
		
		int i2 = 0;
		while(i2 < seq.length()) {
			int pos = seq.indexOf("gg", i2);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i2++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i2 = pos + 1;
		}
		
		int i3 = 0;
		while(i3 < seq.length()) {
			int pos = seq.indexOf("Gg", i3);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i3++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i3 = pos + 1;
		}
		
		int i4 = 0;
		while(i4 < seq.length()) {
			int pos = seq.indexOf("gG", i4);
			if(pos == -1) {
				break;
			}
			if(pos == 0) {
				// Won't be fully contained
				i4++;
				continue;
			}
			if(pos + 2 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ngg = new BasicAnnotation(chr.getId(), start + pos - 1, start + pos + 2, Strand.POSITIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ngg.toUCSC() + ":" + ngg.getOrientation().toString());
			rtrn.add(ngg);
			i4 = pos + 1;
		}
		
		// Find all CCNs
		int j1 = 0;
		while(j1 < seq.length()) {
			int pos = seq.indexOf("CC", j1);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j1 = pos + 1;
		}
		
		int j2 = 0;
		while(j2 < seq.length()) {
			int pos = seq.indexOf("cC", j2);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j2 = pos + 1;
		}
		
		int j3 = 0;
		while(j3 < seq.length()) {
			int pos = seq.indexOf("cc", j3);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j3 = pos + 1;
		}
		
		int j4 = 0;
		while(j4 < seq.length()) {
			int pos = seq.indexOf("Cc", j4);
			if(pos == -1) break;
			if(pos + 3 > seq.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation ccn = new BasicAnnotation(chr.getId(), start + pos, start + pos + 3, Strand.NEGATIVE);
			//logger.debug("NGG in " + chr.getId() + ":" + start + "-" + end + "\t" + ccn.toUCSC() + ":" + ccn.getOrientation().toString());
			rtrn.add(ccn);
			j4 = pos + 1;
		}
		
		return rtrn;
		
	}
	
	private static void validateStartEnd(int start, int end) {
		if(end - start != 20) {
			throw new IllegalArgumentException("end - start must equal 20");
		}
	}
	
	private static void validateStrand(Strand orientation) {
		if(orientation.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand cannot be unknown");
		}
	}
	
	private static void validateSequence20(Sequence sequence) {
		if(sequence.getLength() != 20) {
			throw new IllegalArgumentException("Sequence length must be 20 " + sequence.getSequenceBases());
		}
	}
	
	private static void validateSequence23(Sequence sequence) {
		if(sequence.getLength() != 23) {
			throw new IllegalArgumentException("Sequence length must be 23: " + sequence.getSequenceBases());
		}
		if(!sequence.getSequenceBases().substring(21, 23).equals("GG")) {
			throw new IllegalArgumentException("Sequence must end in GG: " + sequence.getSequenceBases());
		}
	}
	
	public String getBedString() {
		return annot.toBED();
	}
	
	public int hashCode() {
		String s = annot.toBED() + name + sequence20.getSequenceBases() + target.toBED();
		return s.hashCode();
	}
	
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) return false;
		GuideRNA g = (GuideRNA)o;
		return hashCode() == g.hashCode();
	}
	
}

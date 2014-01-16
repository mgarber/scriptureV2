package nextgen.editing.crispr;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

import broad.core.sequence.Sequence;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;

/**
 * A pair of guide RNAs, one on each strand, with a target gene for editing
 * @author prussell
 */
public class GuideRNAPair {
	
	private GuideRNA plusStrandGuideRNA;
	private GuideRNA minusStrandGuideRNA;
	private Gene target;
	public static Logger logger = Logger.getLogger(GuideRNAPair.class.getName());
	
	/**
	 * @param plusStrand The guide RNA on the plus strand
	 * @param minusStrand The guide RNA on the minus strand
	 * @param targetGene The target gene
	 */
	public GuideRNAPair(GuideRNA plusStrand, GuideRNA minusStrand, Gene targetGene) {
		if(!plusStrand.getStrand().equals(Strand.POSITIVE) || !minusStrand.getStrand().equals(Strand.NEGATIVE)) {
			throw new IllegalArgumentException("Strands must be positive,negative");
		}
		plusStrandGuideRNA = plusStrand;
		minusStrandGuideRNA = minusStrand;
		target = targetGene;
	}
	
	/**
	 * Get all valid guide RNA pairs fully contained in the window
	 * @param chr Window chromosome
	 * @param start Window start
	 * @param end Position after last position of window
	 * @param targetGene Target gene
	 * @return All fully contained guide RNA pairs (one member of pair on each strand)
	 */
	public static Collection<GuideRNAPair> findAll(Sequence chr, int start, int end, Gene targetGene) {
		Collection<GuideRNA> all = GuideRNA.findAll(chr, start, end, targetGene);
		logger.debug("There are " + all.size() + " total guide RNAs in " + chr.getId() + ":" + start + "-" + end);
		Collection<GuideRNA> plus = new ArrayList<GuideRNA>();
		Collection<GuideRNA> minus = new ArrayList<GuideRNA>();
		Collection<GuideRNAPair> rtrn = new ArrayList<GuideRNAPair>();
		for(GuideRNA g : all) {
			if(g.isPlusStrand()) {
				plus.add(g);
			}
			if(g.isMinusStrand()) {
				minus.add(g);
			}
		}
		for(GuideRNA p : plus) {
			for(GuideRNA m : minus) {
				GuideRNAPair pair = new GuideRNAPair(p, m, targetGene);
				//logger.debug("Adding pair " + pair.toString());
				rtrn.add(pair);
			}
		}
		return rtrn;
	}
	
	public Gene getTargetGene() {
		return target;
	}
	
	public String toString() {
		return target.getName() + ":" + plusStrandGuideRNA.toString() + "," + minusStrandGuideRNA.toString();
	}
	
	public GuideRNA getPlusStrandGuideRNA() {
		return plusStrandGuideRNA;
	}
	
	public GuideRNA getMinusStrandGuideRNA() {
		return minusStrandGuideRNA;
	}
	
	public boolean leftIsPlusStrand() {
		return plusStrandGuideRNA.getStart() < minusStrandGuideRNA.getStart();
	}
	
	public GuideRNA getLeftGuideRNA() {
		return leftIsPlusStrand() ? plusStrandGuideRNA : minusStrandGuideRNA;
	}
	
	public GuideRNA getRightGuideRNA() {
		return leftIsPlusStrand() ? minusStrandGuideRNA : plusStrandGuideRNA;
	}
	
	/**
	 * Get column headers for the string returned by getOligos()
	 * @return Tab delimited string with field names
	 */
	public static String getOligoFieldNames() {
		String rtrn = "oligo_right_facing_top_strand\t";
		rtrn += "oligo_right_facing_bottom_strand\t";
		rtrn += "oligo_left_facing_top_strand\t";
		rtrn += "oligo_left_facing_bottom_strand";
		return rtrn;
	}
	
	/**
	 * Get a string describing the oligos to order, formatted to be written in a table
	 * The oligo sequences include the extra GTTT or CGGT at the end
	 * Column headers can be obtained with getOligoDescriptionFieldNames()
	 * @return Tab delimited string with the two oligos for each guide RNA
	 */
	public String getOligos() {
		String rightFacingTopStrand = plusStrandGuideRNA.getSequenceString() + "GTTT";
		String rightFacingBottomStrand = Sequence.reverseSequence(plusStrandGuideRNA.getSequenceString()) + "CGGT";
		String leftFacingTopStrand = minusStrandGuideRNA.getSequenceString() + "GTTT";
		String leftFacingBottomStrand = Sequence.reverseSequence(minusStrandGuideRNA.getSequenceString()) + "CGGT";
		String rtrn = rightFacingTopStrand + "\t";
		rtrn += rightFacingBottomStrand + "\t";
		rtrn += leftFacingTopStrand + "\t";
		rtrn += leftFacingBottomStrand;
		return rtrn;
	}
	
	/**
	 * Write oligos to a table
	 * The oligos include the GTTT or CGGT at the end
	 * @param pairs Guide RNA pairs to write
	 * @param outFile Output table
	 * @throws IOException
	 */
	public static void writeOligoTable(Collection<GuideRNAPair> pairs, String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		String header = "target_gene\t";
		header += "oligo_ID\t";
		header += getOligoFieldNames() + "\t";
		w.write(header + "\n");
		for(GuideRNAPair pair : pairs) {
			String line = pair.getTargetGene().getName() + "\t";
			line += pair.toString() + "\t";
			line += pair.getOligos() + "\t";
			w.write(line + "\n");
		}
		w.close();
	}
	
	/**
	 * @return Whether the two guide RNAs are non-overlapping
	 */
	public boolean nonoverlapping() {
		int leftEnd = getLeftGuideRNA().getEnd();
		int rightStart = getRightGuideRNA().getStart();
		return leftEnd <= rightStart;
	}
	
	/**
	 * @return Whether the two guide RNAs face away from each other; i.e., the right RNA is on top strand and left is on bottom strand
	 */
	public boolean facingOutward() {
		boolean rtrn = nonoverlapping() && getLeftGuideRNA().getStrand().equals(Strand.NEGATIVE) && getRightGuideRNA().getStrand().equals(Strand.POSITIVE);
		//if(rtrn) logger.debug("Facing outward: " + toString());
		//else logger.debug("Facing inward: " + toString());
		return rtrn;
	}
	
	/**
	 * Write pairs to a bed file where each pair is displayed as an annotation with two blocks
	 * @param pairs The collection of guide RNA pairs to write
	 * @param bedFile The bed file to write
	 * @throws IOException
	 */
	public static void writeBED(Collection<GuideRNAPair> pairs, String bedFile) throws IOException {
		writeBED(pairs, bedFile, false);
	}
		
	/**
	 * Write pairs to a bed file where each pair is displayed as an annotation with two blocks
	 * @param pairs The collection of guide RNA pairs to write
	 * @param bedFile The bed file to write
	 * @param append If true, write onto the end of existing bed file; if false, overwrite file
	 * @throws IOException
	 */
	public static void writeBED(Collection<GuideRNAPair> pairs, String bedFile, boolean append) throws IOException {
		FileWriter w = new FileWriter(bedFile, append);
		for(GuideRNAPair pair : pairs) {
			w.write(pair.toBED() + "\n");
		}
		w.close();
	}
	
	/**
	 * Write pairs to a bed file where each pair is displayed as an annotation with two blocks
	 * @param pairs The collection of guide RNA pairs to write
	 * @param bedFile FileWriter object for bed output
	 * @throws IOException
	 */
	public static void writeBED(Collection<GuideRNAPair> pairs, FileWriter writer) throws IOException {
		for(GuideRNAPair pair : pairs) {
			writer.write(pair.toBED() + "\n");
		}
	}
	
	/**
	 * Get all the individual guide RNAs that make up the pairs
	 * @param pairs A collection of guide RNA pairs
	 * @return A collection of the individual guide RNAs
	 */
	public static Collection<GuideRNA> getIndividualGuideRNAs(Collection<GuideRNAPair> pairs) {
		Collection<GuideRNA> guides = new ArrayList<GuideRNA>();
		for(GuideRNAPair pair : pairs) {
			guides.add(pair.getLeftGuideRNA());
			guides.add(pair.getRightGuideRNA());
		}
		return guides;
	}
	
	/**
	 * Get the guide pair as a blocked annotation
	 * @return Annotation representation with the two guide RNAs as blocks
	 */
	public Annotation asAnnotation() {
		return asAnnotation(toString());
	}
	
	/**
	 * Get the guide pair as a blocked annotation
	 * @param name Name of annotation to set
	 * @return Annotation representation with the two guide RNAs as blocks
	 */
	public Annotation asAnnotation(String name) {
		Annotation plus = new BasicAnnotation(plusStrandGuideRNA.getChr(), plusStrandGuideRNA.getStart(), plusStrandGuideRNA.getEnd());
		BasicAnnotation minus = new BasicAnnotation(minusStrandGuideRNA.getChr(), minusStrandGuideRNA.getStart(), minusStrandGuideRNA.getEnd());
		BasicAnnotation both = new BasicAnnotation(plus);
		both.addBlocks(minus);
		both.setName(name);
		return both;
	}
	
	/**
	 * Get a BED line for the pair
	 * @return BED format representation of an annotation with the two guide RNAs as blocks
	 */
	public String toBED() {
		return toBED(toString());
	}
	
	/**
	 * Get a BED line for the pair
	 * @param name Name of annotation for bed line
	 * @return BED format representation of an annotation with the two guide RNAs as blocks
	 */
	public String toBED(String name) {
		return asAnnotation(name).toBED();
	}
	
	/**
	 * @return The number of positions between the last position of the left RNA and the first position of the right RNA
	 */
	public int getInnerDistance() {
		int rtrn = getRightGuideRNA().getStart() - getLeftGuideRNA().getEnd();
		//logger.debug("Inner distance of " + toString() + ":\t" + rtrn);
		return rtrn;
	}
	
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) return false;
		GuideRNAPair p = (GuideRNAPair)o;
		return minusStrandGuideRNA.equals(p.getMinusStrandGuideRNA()) && plusStrandGuideRNA.equals(p.getPlusStrandGuideRNA()) && target.equals(p.getTargetGene());
	}
	
	public int hashCode() {
		String h = minusStrandGuideRNA.hashCode() + "_" + plusStrandGuideRNA.hashCode();
		return h.hashCode();
	}
	
	
}

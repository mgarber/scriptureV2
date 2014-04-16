package nextgen.editing.crispr;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.primer3.PrimerPair;
import broad.core.primer3.PrimerUtils;
import broad.core.sequence.Sequence;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

/**
 * A pair of guide RNAs, one on each strand, with a target gene for editing
 * @author prussell
 */
public class NickingGuideRNAPair extends BasicAnnotation {
	

	private GuideRNA plusStrandGuideRNA;
	private GuideRNA minusStrandGuideRNA;
	private Gene target;
	public static Logger logger = Logger.getLogger(NickingGuideRNAPair.class.getName());
	private double score;
	public static int PRIMER_LENGTH = 20;
	public static double OPTIMAL_PRIMER_TM = 67;
	
	/**
	 * @param plusStrand The guide RNA on the plus strand
	 * @param minusStrand The guide RNA on the minus strand
	 * @param targetGene The target gene
	 */
	public NickingGuideRNAPair(GuideRNA plusStrand, GuideRNA minusStrand, Gene targetGene) {
		this(null, plusStrand, minusStrand, targetGene);
	}
	
	/**
	 * @param plusStrand The guide RNA on the plus strand
	 * @param minusStrand The guide RNA on the minus strand
	 * @param targetGene The target gene
	 * @param name Pair name
	 */
	public NickingGuideRNAPair(String name, GuideRNA plusStrand, GuideRNA minusStrand, Gene targetGene) {
		this(name, plusStrand, minusStrand, targetGene, 0);
	}

	/**
	 * @param plusStrand The guide RNA on the plus strand
	 * @param minusStrand The guide RNA on the minus strand
	 * @param targetGene The target gene
	 * @param pairScore Score for guide RNA pair
	 */
	public NickingGuideRNAPair(GuideRNA plusStrand, GuideRNA minusStrand, Gene targetGene, double pairScore) {
		this(null, plusStrand, minusStrand, targetGene, pairScore);
	}
	
	/**
	 * @param plusStrand The guide RNA on the plus strand
	 * @param minusStrand The guide RNA on the minus strand
	 * @param targetGene The target gene
	 * @param pairScore Score for guide RNA pair
	 * @param name Pair name
	 */
	public NickingGuideRNAPair(String name, GuideRNA plusStrand, GuideRNA minusStrand, Gene targetGene, double pairScore) {
		super(plusStrand.getChr(), plusStrand.getStart(), plusStrand.getEnd());
		addBlocks(new BasicAnnotation(minusStrand.getChr(), minusStrand.getStart(), minusStrand.getEnd()));
		setName(name);

		if(!plusStrand.getStrand().equals(Strand.POSITIVE) || !minusStrand.getStrand().equals(Strand.NEGATIVE)) {
			throw new IllegalArgumentException("Strands must be positive,negative");
		}
		plusStrandGuideRNA = plusStrand;
		minusStrandGuideRNA = minusStrand;
		target = targetGene;
		score = pairScore;
	}
	
	/**
	 * Set pair score
	 * @param pairScore Value to set
	 */
	public void setScore(double pairScore) {
		score = pairScore;
	}
	
	/**
	 * Get all valid guide RNA pairs fully contained in the window
	 * @param chr Window chromosome
	 * @param start Window start
	 * @param end Position after last position of window
	 * @param targetGene Target gene
	 * @return All fully contained guide RNA pairs (one member of pair on each strand)
	 */
	public static Collection<NickingGuideRNAPair> findAll(Sequence chr, int start, int end, Gene targetGene) {
		Collection<GuideRNA> all = GuideRNA.findAll(chr, start, end, targetGene);
		//logger.debug("There are " + all.size() + " total guide RNAs in " + chr.getId() + ":" + start + "-" + end);
		Collection<GuideRNA> plus = new ArrayList<GuideRNA>();
		Collection<GuideRNA> minus = new ArrayList<GuideRNA>();
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
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
				NickingGuideRNAPair pair = new NickingGuideRNAPair(p, m, targetGene);
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
		return target.getName() + ":" + getStart() + "-" + getEnd();
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
		throw new UnsupportedOperationException("Not implemented");
	}
	
	/**
	 * Get pair score
	 * @return Pair score
	 */
	public double getScore() {
		return score;
	}
	
	/**
	 * Get a string describing the oligos to order, formatted to be written in a table
	 * Column headers can be obtained with getOligoFieldNames()
	 * @return Tab delimited string with the oligos for each guide RNA
	 */
	public String getOligos() {
		throw new UnsupportedOperationException("Not implemented");
	}
	
	public static void writeOligoTable(Collection<Collection<NickingGuideRNAPair>> pools, String leftFlank, String rightFlank, String outFile, String primer3core) throws IOException {
		FileWriter w = new FileWriter(outFile);
		String header = "target_gene\t";
		header += "oligo_ID_minus_strand\t";
		header += "oligo_ID_plus_strand\t";
		header += "pool\t";
		header += "pool_left_primer\t";
		header += "pool_right_primer\t";
		header += "coords\t";
		header += "minus_strand_guide_RNA\t";
		header += "plus_strand_guide_RNA\t";
		header += "minus_strand_oligo\t";
		header += "plus_strand_oligo\t";
		w.write(header + "\n");
		int poolNum = 0;
		for(Collection<NickingGuideRNAPair> pool : pools) {
			PrimerPair primer = PrimerUtils.getOneSyntheticPrimerPair(PRIMER_LENGTH, primer3core, OPTIMAL_PRIMER_TM, null);
			String leftPrimer = primer.getLeftPrimer();
			String rightPrimer = primer.getRightPrimer();
			String rightPrimerRC = Sequence.reverseSequence(rightPrimer);
			String poolName = "pool" + poolNum;
			for(NickingGuideRNAPair pair : pool) {
				String minus = pair.getMinusStrandGuideRNA().getSequenceString();
				String plus = pair.getPlusStrandGuideRNA().getSequenceString();
				String line = pair.getTargetGene().getName() + "\t";
				line += pair.getMinusStrandGuideRNA().toString() + "\t";
				line += pair.getPlusStrandGuideRNA().toString() + "\t";
				line += poolName + "\t";
				line += leftPrimer + "\t";
				line += rightPrimer + "\t";
				line += pair.toUCSC() + "\t";
				line += minus + "\t";
				line += plus + "\t";
				line += leftPrimer + leftFlank + minus + rightFlank + rightPrimerRC + "\t";
				line += leftPrimer+ leftFlank + plus + rightFlank + rightPrimerRC + "\t";
				w.write(line + "\n");
			}
			poolNum++;
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
	
	public static void writeBED(Collection<Collection<NickingGuideRNAPair>> pools, String bedPrefix) throws IOException {
		int poolNum = 0;
		for(Collection<NickingGuideRNAPair> pool : pools) {
			String outBed = bedPrefix + "_pool" + poolNum + ".bed";
			writeBED(pool, outBed, false);
			poolNum++;
		}
	}
		
	/**
	 * Write pairs to a bed file where each pair is displayed as an annotation with two blocks
	 * @param pairs The collection of guide RNA pairs to write
	 * @param bedFile The bed file to write
	 * @param append If true, write onto the end of existing bed file; if false, overwrite file
	 * @throws IOException
	 */
	public static void writeBED(Collection<NickingGuideRNAPair> pairs, String bedFile, boolean append) throws IOException {
		logger.info("Writing " + pairs.size() + " guide RNA pairs to file " + bedFile + ".");
		FileWriter w = new FileWriter(bedFile, append);
		for(NickingGuideRNAPair pair : pairs) {
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
	public static void writeBED(Collection<NickingGuideRNAPair> pairs, FileWriter writer) throws IOException {
		for(NickingGuideRNAPair pair : pairs) {
			writer.write(pair.toBED() + "\n");
		}
	}
	
	/**
	 * Get all the individual guide RNAs that make up the pairs
	 * @param pairs A collection of guide RNA pairs
	 * @return A collection of the individual guide RNAs
	 */
	public static Collection<GuideRNA> getIndividualGuideRNAs(Collection<NickingGuideRNAPair> pairs) {
		Collection<GuideRNA> guides = new ArrayList<GuideRNA>();
		for(NickingGuideRNAPair pair : pairs) {
			guides.add(pair.getLeftGuideRNA());
			guides.add(pair.getRightGuideRNA());
		}
		return guides;
	}

	/**
	 * Get the collection of pairs as a list with nondescending scores
	 * @param guidePairs The pairs to sort by score
	 * @return List sorted by score
	 */
	public static ArrayList<NickingGuideRNAPair> sortByScoreNondescending(Collection<NickingGuideRNAPair> guidePairs) {
		Map<Double, Collection<NickingGuideRNAPair>> byScore = new TreeMap<Double, Collection<NickingGuideRNAPair>>();
		for(NickingGuideRNAPair g : guidePairs) {
			Double d = new Double(g.getScore());
			if(!byScore.containsKey(d))	{
				byScore.put(d, new ArrayList<NickingGuideRNAPair>());
			}
			byScore.get(d).add(g);
		}
		ArrayList<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		for(Double d : byScore.keySet()) {
			for(NickingGuideRNAPair g : byScore.get(d)) {
				rtrn.add(g);
			}
		}
		return rtrn;
	}
	
	private static Collection<Collection<NickingGuideRNAPair>> getEvenPools(Collection<NickingGuideRNAPair> guidePairs, int numPools) {
		TreeSet<NickingGuideRNAPair> sorted = new TreeSet<NickingGuideRNAPair>();
		sorted.addAll(guidePairs);
		ArrayList<Collection<NickingGuideRNAPair>> rtrn = new ArrayList<Collection<NickingGuideRNAPair>>();
		for(int i = 0; i < numPools; i++) {
			rtrn.add(new TreeSet<NickingGuideRNAPair>());
		}
		int n = 0;
		for(NickingGuideRNAPair g : sorted) {
			rtrn.get(n % numPools).add(g);
			n++;
		}
		return rtrn;
	}
	
	/**
	 * Divide into a specified number of pools of mutually non-overlapping guide pairs
	 * Final pools are nondeterministic
	 * @param guidePairs The guide RNA pairs to pool 
	 * @param numPools The number of pools to make
	 * @return Collection of pools of non-overlapping guide pairs
	 */
	public static Collection<Collection<NickingGuideRNAPair>> poolNonOverlapping(Collection<NickingGuideRNAPair> guidePairs, int numPools) {
		logger.info("");
		logger.info("Distributing " + guidePairs.size() + " guide RNA pairs into " + numPools + " pools of pairwise nonoverlapping pairs.");
		ArrayList<Collection<NickingGuideRNAPair>> rtrn = new ArrayList<Collection<NickingGuideRNAPair>>();
		for(int i = 0; i < numPools; i++) {
			rtrn.add(new TreeSet<NickingGuideRNAPair>());
		}
		Random rand = new Random();
		rand.setSeed(System.currentTimeMillis());
		for(NickingGuideRNAPair pair : guidePairs) {
			boolean inPool = false;
			int r = Math.abs(rand.nextInt());
			for(int i = r; i < r + rtrn.size(); i++) {
				int poolNum = i % rtrn.size();
				Collection<NickingGuideRNAPair> pool = rtrn.get(poolNum);
				boolean ok = true;
				for(NickingGuideRNAPair other : pool) {
					if(pair.overlaps(other)) {
						ok = false;
						break;
					}
				} if(ok) {
					pool.add(pair);
					inPool = true;
					break;
				}
			}
			if(!inPool) {
				// Can't be added to a pool; skip this one
				logger.warn("Couldn't add " + pair.toString() + " to a pool. Skipping.");
				continue;
			}
		}
		return rtrn;
	}
	
	/**
	 * Divide into pools of mutually non-overlapping guide pairs
	 * Try evenly sized pools until they are all non-overlapping
	 * @param guidePairs The guide RNA pairs to pool
	 * @return Collection of pools of non-overlapping guide pairs
	 */
	public static Collection<Collection<NickingGuideRNAPair>> poolNonOverlapping(Collection<NickingGuideRNAPair> guidePairs) {
		logger.info("");
		logger.info("Distributing " + guidePairs.size() + " guide RNA pairs into pools of pairwise nonoverlapping pairs.");
		int numPools = 1;
		while(numPools <= guidePairs.size()) {
			logger.info("Trying " + numPools + " pools");
			Collection<Collection<NickingGuideRNAPair>> pools = getEvenPools(guidePairs, numPools);
			boolean overlaps = false;
			for(Collection<NickingGuideRNAPair> pool : pools) {
				Collection<Annotation> annot = new TreeSet<Annotation>();
				annot.addAll(pool);
				if(!AnnotationUtils.pairwiseNonoverlappingIgnoreStrand(annot)) {
					numPools++;
					overlaps = true;
					break;
				}
			}
			if(!overlaps) {
				logger.info("Success with " + pools.size() + " pools");
				return pools;
			}
		}
		throw new IllegalArgumentException("Can't create pools of pairwise nonoverlapping guide pairs");
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
		NickingGuideRNAPair p = (NickingGuideRNAPair)o;
		return minusStrandGuideRNA.equals(p.getMinusStrandGuideRNA()) && plusStrandGuideRNA.equals(p.getPlusStrandGuideRNA()) && target.equals(p.getTargetGene());
	}
	
	public int hashCode() {
		String h = minusStrandGuideRNA.hashCode() + "_" + plusStrandGuideRNA.hashCode();
		return h.hashCode();
	}
	
	
}

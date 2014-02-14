package nextgen.editing.crispr.score;

import nextgen.editing.crispr.GuideRNA;
import java.io.*;
import java.util.*;

import org.apache.log4j.Logger;

/**
 * Guide RNA score based on Feng Zhang's algorithm
 * Implemented by JE based on http://crispr.mit.edu/about
 * @author prussell
 */
public class GuideOffTargetScore implements GuideRNAScore {

	public static Logger logger = Logger.getLogger(GuideOffTargetScore.class.getName());
	
	private List<String> offTargetSeqs = new ArrayList<String>();
	public static double[] penalty = new double[] {0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};	
	public static int maxMismatch = 5;
	
	/**
	 * Bed file containing off target sites, with the 23-base off-target sequence in the name column
	 * @param offTargetSites
	 */
	public GuideOffTargetScore(File offTargetSites) {
		offTargetSeqs = loadOffTargetSeqs(offTargetSites);
	}
	
	@Override
	public double getScore(GuideRNA guideRNA) {
		return scoreAllHits(guideRNA.getSequenceString());
	}

	@Override
	public String getScoreName() {
		return "guide_off_target_score";
	}
	
	
	private List<String> loadOffTargetSeqs(File file) {
		// For now, load whole thing into memory.  
		// Definitely smarter ways to do this in the future...
		List<String> seqs = new ArrayList<String>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line;
			while ((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				seqs.add(fields[3]);
			}
		} catch (IOException e) {
			logger.error(e.getMessage());
			System.exit(1);
		}
		return seqs;
	}
	
	
	private double scoreAllHits(String guideSequence) {
		double sum = 0;
		for (String offTarget : offTargetSeqs) {
			sum += scoreSingleHit(guideSequence, offTarget, 20);
		}
		return (10000.0 / (100.0 + sum));
	}
	
	
	private double scoreSingleHit(String one, String two, int n) {
		// Assumes two strings of equal length
		// Ask user to pass n to avoid overhead
		int mismatches = 0;
		double mismatchWeightProduct = 1;
		List<Integer> mismatchPositions = new ArrayList<Integer>();
		
		for (int i = 0; i < n; i++) {
			if (one.charAt(i) != two.charAt(i)) {
				mismatches++;
				mismatchWeightProduct = mismatchWeightProduct * (1.0 - penalty[i]);
				mismatchPositions.add(i);
			}
			if (mismatches > maxMismatch) return 0; 
		}
		
		double pairwiseDistancePenalty = 1.0 / (1.0 + (19.0-getMeanPairwiseDistance(mismatchPositions))/19.0 * 4.0);
		
		// Note that the extra factor of 100 is not noted on the web site but is included in their guideRNA score calculations
		double result = (mismatchWeightProduct * pairwiseDistancePenalty / ((double) mismatches * mismatches) * 100.0);
		//logger.info("Single hit score of " + one + " and " + two + " = " + result);
		//logger.info("Mismatch WP = " + mismatchWeightProduct + ", pairwise DP = " + pairwiseDistancePenalty);
		return result;
	}
	
	
	private double getMeanPairwiseDistance(List<Integer> mismatchPositions) {
		int sum = 0;
		int total = 0;
		int length = mismatchPositions.size();
		if (length == 0) return 0;

		for (int i = 0; i < length-1; i++) {
			for (int j = i+1; j < length; j++) {
				total++;
				sum += mismatchPositions.get(j) - mismatchPositions.get(i);
			}
		}

		double result = ((double) sum) / ((double) total);
		//logger.info("Mean pairwise distance of " + mismatchPositions.toString() + " = " + result);
		return result;
	}
}

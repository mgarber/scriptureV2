package nextgen.editing.crispr.score;

import nextgen.editing.crispr.*;
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
	
	private byte[] offTargetSeqs;
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
	
	
	private byte[] loadOffTargetSeqs(File file) {
		// For now, load whole thing into memory.  
		// Definitely smarter ways to do this in the future...
		logger.info("Loading target sequences from bitpacked file ... ");
		byte[] bitSeqs = null;
		try {
			BufferedInputStream bs = new BufferedInputStream(new FileInputStream(file));
			bitSeqs = new byte[bs.available()];
			bs.read(bitSeqs);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);;
		}
		logger.info("Complete.");
		return bitSeqs;
	}
	
	
	private double scoreAllHits(String guideSequence) {
		byte[] guideBytes = BitpackGuideSequences.sequenceToBits(guideSequence);
		double sum = 0;
		int nExact = 0;
		for (int i = 0; i < offTargetSeqs.length; i += 5) {
			double currScore = scoreSingleHit(guideBytes, Arrays.copyOfRange(offTargetSeqs, i, i+5), 20);
			if (Double.isInfinite(currScore)) {
				nExact++;
				if (nExact > 1) return 0.0;
			}
			sum = sum + currScore;
		}
		return (10000.0 / (100.0 + sum));
	}
	
	
	private double scoreSingleHit(byte[] one, byte[] two, int n) {
		// Assumes two byte arrays of equal length
		// Ask user to pass n to avoid overhead
		int mismatches = 0;
		double mismatchWeightProduct = 1;
		List<Integer> mismatchPositions = new ArrayList<Integer>();
		
		for (int i = 0; i < n*2; i += 2) {
			byte mask = (byte) (0b00000011 << (i % 8));
			if ((one[i/8] & mask) != (two[i/8] & mask)) {
				mismatches++;
				mismatchWeightProduct = mismatchWeightProduct * (1.0 - penalty[i/2]);
				mismatchPositions.add(i/2);
			}
			if (mismatches > maxMismatch) return 0; 
		}
		
		double pairwiseDistancePenalty = 1.0 / (1.0 + (19.0-getMeanPairwiseDistance(mismatchPositions))/19.0 * 4.0);
		
		// Note that the extra factor of 100 is not noted on the web site but is included in their guideRNA score calculations
		double result = (mismatchWeightProduct * pairwiseDistancePenalty / ((double) mismatches * mismatches) * 100.0);
		return result;
	}
	
	
	/**
		String version
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
		int n = 0;
		for (String offTarget : offTargetSeqs) {
			sum += scoreSingleHit(guideSequence, offTarget, 20);
			n++;
			if (n % 1000000 == 0) logger.info("Scored " + n + " off target sites.");
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
		logger.info("Single hit score of " + one + " and " + two + " = " + result);
		//logger.info("Mismatch WP = " + mismatchWeightProduct + ", pairwise DP = " + pairwiseDistancePenalty);
		return result;
	}
	*/
	
	
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

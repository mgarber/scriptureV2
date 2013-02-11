/*
 *  Jesse Engreitz
 *  March 12, 2012
 *  Utilities for Methylation Analysis
 */

package broad.pda.methylation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import broad.core.sequence.Sequence;

public class MethylationUtils {
	
	
	/* 
	 *  Returns the indices (0-based) of CpG dinucleotides in a sequence
	 */
	public static List<Integer> findCpGs(String seq) {
		List<Integer> indices = new ArrayList<Integer>();
		
		seq = seq.toUpperCase();
		for (int i = -1; (i = seq.indexOf("CG", i+1)) != -1; ) {
			indices.add(i);
		}
		
		return indices;
	}
	
	/*
	 *  CpG bisulfite conversion.  Assumes that no CpGs are methylated.
	 */
	public static String bisulfiteConvert(String seq) {
		return bisulfiteConvert(seq, false);
	}
	
	
	/*
	 *  CpG bisulfite conversion.  
	 *  @param methylated = true assumes that all CpGs are methylated
	 */
	public static String bisulfiteConvert(String seq, Boolean methylated) {
		if (methylated) {
			seq = seq.replaceAll("CG", "_");
			seq = seq.replaceAll("cg", "Q");
			seq = seq.replaceAll("cG", "X");
			seq = seq.replaceAll("Cg", "Y");
		}
		seq = seq.replaceAll("C", "T");
		seq = seq.replaceAll("c", "t");
		if (methylated) {
			seq = seq.replaceAll("_", "CG");
			seq = seq.replaceAll("Q", "cg");
			seq = seq.replaceAll("X", "cG");
			seq = seq.replaceAll("Y", "Cg");
		}
		
		return seq;
	}
	
	
	/**
	 * @param seq
	 * @return base counts in this order:  A C G T N
	 */
	public static int[] getBaseCounts(String seq) {
		int[] counts = new int[5];
		seq = seq.toUpperCase();
		for (int i = 0; i < seq.length(); i++) {
			if (seq.charAt(i) == 'A') counts[0]++;
			else if (seq.charAt(i) == 'C') counts[1]++;
			else if (seq.charAt(i) == 'G') counts[2]++;
			else if (seq.charAt(i) == 'T') counts[3]++;
			else if (seq.charAt(i) == 'N') counts[4]++;
		}
		return counts;
	}
	
	
	public static double[] getBaseFrequencies(String seq) {
		int[] counts = getBaseCounts(seq);
		double[] freqs = new double[counts.length];
		for (int i = 0; i < counts.length; i++)
			freqs[i] = ((double) counts[i]) / ((double) seq.length());
		return freqs;
	}
	
	
	public static double getGcPct(String seq) {
		double[] freqs = getBaseFrequencies(seq);
		return freqs[1] + freqs[2];
	}
	
	
	public static double getGcRatio(String seq) {
		double[] freqs = getBaseFrequencies(seq);
		return freqs[2] / freqs[1];
	}
	
	
	/**
	 * Returns true if it seems likely that the sequence is not bisulfite converted.
	 * Algorithm:  1. Check number of differences between the sequence and its bisulfite converted form.
	 *             2. Check number of differences between the reverse complement and its bisulfite converted form.
	 *             If both of these numbers >= 2, then guess that it's unconverted.
	 * @param seq
	 * @return
	 */
	public static boolean guessUnconverted(String seq) {
		boolean guess = false;
		String converted = bisulfiteConvert(seq, true);
		
		int senseDiff = numCharDifferences(seq, converted);
		
		String rc = Sequence.reverseSequence(seq);
		String rcConverted = bisulfiteConvert(rc, true);
		int rcDiff = numCharDifferences(rc, rcConverted);
		
		if (senseDiff > 1 && rcDiff > 1) guess = true;
		
		return guess;
	}
	
	private static int numCharDifferences(String s1, String s2) {
		if (s1.length() != s2.length()) 
			throw new UnsupportedOperationException("s1 and s2 must be equal length");
		int diff = 0;
		for (int i = 0; i < s1.length(); i++) {
			if (s1.charAt(i) != s2.charAt(i))
				diff++;
		}
		return diff;
	}
	
	public static List<String> fillArray(List<String> senseOligos, List<String> antisenseOligos, int arraySize) {
		List<String> rtrn = new ArrayList<String>();
		if (senseOligos.size() == 0) {
			for (String line : senseOligos)
				rtrn.add(line);
		}
		else if (senseOligos.size() > arraySize) {
			System.err.println("Error: Number of oligos generated is greater than the array size.");
		} else {
			int totalPrinted = 0;
			boolean printSense = true;
			while (totalPrinted < arraySize) {
				int numToPrint = Math.min(arraySize - totalPrinted, senseOligos.size());
				List<String> listToPrint = printSense ? senseOligos : antisenseOligos;
				
				// If we've reached the last list, shuffle the list so we randomly sample
				// probes from different targets
				if (numToPrint < senseOligos.size()) Collections.shuffle(listToPrint);
				
				for (String line : listToPrint.subList(0, numToPrint))
					rtrn.add(line);
				
				totalPrinted += numToPrint;
				printSense = !printSense;
			}
		}
		return rtrn;
	}
	
}

/**
 * 
 */
package nextgen.core.utils;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;

/**
 * @author prussell
 *
 */
public class AlignmentUtils {
	
	public static boolean hammingDistanceAtMost(String s1, String s2, int maxDistance, boolean ignoreCase) {
		if(s1.length() != s2.length()) {
			throw new IllegalArgumentException("Strings must have same length");
		}
		int dist = 0;
		for(int i = 0; i < s1.length(); i++) {
			String c1 = ignoreCase ? s1.substring(i,i+1).toUpperCase() : s1.substring(i, i+1);
			String c2 = ignoreCase ? s2.substring(i,i+1).toUpperCase() : s2.substring(i, i+1);
			if(!c1.equals(c2)) {
				dist++;
			}
			if(dist > maxDistance) {
				return false;
			}
		}
		return true;
	}
	
	public static int hammingDistance(String s1, String s2, boolean ignoreCase) {
		if(s1.length() != s2.length()) {
			throw new IllegalArgumentException("Strings must have same length");
		}
		int rtrn = 0;
		for(int i = 0; i < s1.length(); i++) {
			String c1 = ignoreCase ? s1.substring(i,i+1).toUpperCase() : s1.substring(i, i+1);
			String c2 = ignoreCase ? s2.substring(i,i+1).toUpperCase() : s2.substring(i, i+1);
			if(!c1.equals(c2)) {
				rtrn++;
			}
		}
		return rtrn;
	}
	
	/**
	 * Assign orientation to window based on number of reads mapping in each orientation
	 * @param bamFile Bam file
	 * @param window The annotation
	 * @param firstReadTranscriptionStrand Whether read 1 is transcription strand
	 * @param cutoff Proportion of reads needed to call a strand one way or the other (between 0.5 and 1)
	 * @return Guessed orientation of annotation based on majority vote of reads
	 */
	public static Strand assignOrientationToWindow(String bamFile, Annotation window, boolean firstReadTranscriptionStrand, double cutoff) {
		if(cutoff < 0.5 || cutoff > 1) {
			throw new IllegalArgumentException("Cutoff must be between 0.5 and 1.");
		}
		SAMFileReader samReader = new SAMFileReader(new File(bamFile));
		SAMRecordIterator iter = samReader.queryOverlapping(window.getChr(), window.getStart(), window.getEnd());
		int numPlus = 0;
		int numMinus = 0;
		while(iter.hasNext()) {
			SAMRecord read = iter.next();
			boolean firstRead = !read.getSecondOfPairFlag();
			boolean forwardStrand = !read.getReadNegativeStrandFlag();
			boolean firstFirstForward = firstReadTranscriptionStrand && firstRead && forwardStrand;
			boolean firstFirstReverse = firstReadTranscriptionStrand && firstRead && !forwardStrand;
			boolean firstSecondForward = firstReadTranscriptionStrand && !firstRead && forwardStrand;
			boolean firstSecondReverse = firstReadTranscriptionStrand && !firstRead && !forwardStrand;
			boolean secondFirstForward = !firstReadTranscriptionStrand && firstRead && forwardStrand;
			boolean secondFirstReverse = !firstReadTranscriptionStrand && firstRead && !forwardStrand;
			boolean secondSecondForward = !firstReadTranscriptionStrand && !firstRead && forwardStrand;
			boolean secondSecondReverse = !firstReadTranscriptionStrand && !firstRead && !forwardStrand;
			if(firstFirstForward || firstSecondReverse || secondFirstReverse || secondSecondForward) {
				numPlus++;
			}
			if(firstFirstReverse || firstSecondForward || secondFirstForward || secondSecondReverse) {
				numMinus++;
			}
		}
		iter.close();
		int total = numPlus + numMinus;
		double plusPct = (double) numPlus / (double) total;
		double minusPct = (double) numMinus / (double) total;
		if(plusPct > cutoff) {
			samReader.close();
			return Strand.POSITIVE;
		}
		else if(minusPct > cutoff) {
			samReader.close();
			return Strand.NEGATIVE;
		}
		samReader.close();
		return Strand.UNKNOWN;
	}


	
}

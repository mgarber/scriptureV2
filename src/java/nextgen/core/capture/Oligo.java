package nextgen.core.capture;

import broad.core.primer3.PrimerPair;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 * The full oligo including primers
 */
public class Oligo implements Comparable<Oligo> {
	
	private Probe probe;
	private PrimerPair primer;
	private ProbeSet probeSet;
	
	
	/**
	 * @param pr Probe
	 */
	public Oligo(Probe pr) {
		this(pr, null, null);
	}
	
	/**
	 * @param pr Probe
	 * @param ps Probe set
	 */
	public Oligo(Probe pr, ProbeSet ps) {
		this(pr, ps, null);
	}
	
	/**
	 * @param pr Probe
	 * @param ps Probe set
	 * @param primerPair Primer pair
	 */
	public Oligo(Probe pr, ProbeSet ps, PrimerPair primerPair) {
		probe = pr;
		probeSet = ps;
		primer = primerPair;
	}
	
	/**
	 * Check whether a primer pair is compatible with a probe
	 * @param primer The primer pair
	 * @param probeSequenceExcludingPrimers Probe sequence without primers
	 * @return True if there are no matches of the primers to the probe
	 */
	public static boolean primerPairCompatibleWithProbe(PrimerPair primer, String probeSequenceExcludingPrimers) {
		String leftPrimer = primer.getLeftPrimer();
		String rightPrimer = primer.getRightPrimer();
		String leftPrimer3primeEnd = leftPrimer.substring(leftPrimer.length() - 8);
		String leftPrimer3primeEndRC = Sequence.reverseSequence(leftPrimer3primeEnd);
		String rightPrimer3primeEnd = rightPrimer.substring(rightPrimer.length() - 8);
		String rightPrimer3primeEndRC = Sequence.reverseSequence(rightPrimer3primeEnd);
		// Check that the two primer ends do not appear in the oligo
		if(probeSequenceExcludingPrimers.contains(leftPrimer3primeEnd)) {
			return false;
		}
		if(probeSequenceExcludingPrimers.contains(rightPrimer3primeEndRC)) {
			return false;
		}
		// Check that the RCed primer ends do not appear at all in the oligo
		if(probeSequenceExcludingPrimers.contains(leftPrimer3primeEndRC)) {
			return false;
		}
		if(probeSequenceExcludingPrimers.contains(rightPrimer3primeEnd)) {
			return false;
		}
		return true;
	}

	/**
	 * Check whether a primer pair is compatible with a full oligo sequence
	 * @param primer The primer pair
	 * @param oligoSequenceIncludingPrimers Full oligo sequence including this primer pair at the ends
	 * @return True if there are no nonspecific matches of the primers to the oligo
	 */
	public static boolean primerPairCompatibleWithFullOligo(PrimerPair primer, String oligoSequenceIncludingPrimers) {
		String leftPrimer = primer.getLeftPrimer();
		String rightPrimer = primer.getRightPrimer();
		String leftPrimer3primeEnd = leftPrimer.substring(leftPrimer.length() - 8);
		String leftPrimer3primeEndRC = Sequence.reverseSequence(leftPrimer3primeEnd);
		String rightPrimer3primeEnd = rightPrimer.substring(rightPrimer.length() - 8);
		String rightPrimer3primeEndRC = Sequence.reverseSequence(rightPrimer3primeEnd);
		// Check the each oligo contains the 3' ends of the primers once
		int firstOccurrenceLeftPrimer3primeEnd = oligoSequenceIncludingPrimers.indexOf(leftPrimer3primeEnd);
		int lastOccurrenceLeftPrimer3primeEnd = oligoSequenceIncludingPrimers.lastIndexOf(leftPrimer3primeEnd);
		int firstOccurrenceRightPrimer3primeEndRC = oligoSequenceIncludingPrimers.indexOf(rightPrimer3primeEndRC);
		int lastOccurrenceRightPrimer3primeEndRC = oligoSequenceIncludingPrimers.lastIndexOf(rightPrimer3primeEndRC);
		// Check that the two primer ends appear in the oligo
		if(firstOccurrenceLeftPrimer3primeEnd == -1) {
			return false;
		}
		if(firstOccurrenceRightPrimer3primeEndRC == -1) {
			return false;
		}
		// Check that the two primer ends appear at most once in the oligo
		if(firstOccurrenceLeftPrimer3primeEnd != lastOccurrenceLeftPrimer3primeEnd) {
			return false;
		}
		if(firstOccurrenceRightPrimer3primeEndRC != lastOccurrenceRightPrimer3primeEndRC) {
			return false;
		}
		// Check that the RCed primer ends do not appear at all in the oligo
		if(oligoSequenceIncludingPrimers.contains(leftPrimer3primeEndRC)) {
			return false;
		}
		if(oligoSequenceIncludingPrimers.contains(rightPrimer3primeEnd)) {
			return false;
		}
		return true;
	}
	
	/**
	 * @return Probe
	 */
	public Probe getProbe() {
		return probe;
	}
	
	/**
	 * @return Primer pair
	 */
	public PrimerPair getPrimer() {
		return primer;
	}
	
	/**
	 * @return Probe set the oligo belongs to
	 */
	public ProbeSet getProbeSet() {
		return probeSet;
	}
	
	/**
	 * @param primerPair Primer pair
	 */
	public void setPrimer(PrimerPair primerPair) {
		primer = primerPair;
	}
	
	/**
	 * @return Full oligo sequence
	 */
	public String getOligoBases() {
		return primer.getLeftPrimer().toUpperCase() + probe.getProbeSequence().toUpperCase() + Sequence.reverseSequence(primer.getRightPrimer()).toUpperCase();
	}

	@Override
	public int compareTo(Oligo o) {
		if(!probe.equals(o.getProbe()) ) {
			return probe.compareTo(o.getProbe());
		}
		return primer.compareTo(o.getPrimer());
	}
	
}

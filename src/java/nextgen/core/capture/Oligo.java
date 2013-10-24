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
		this(pr, null);
	}
	
	/**
	 * @param pr Probe
	 * @param primerPair Primer pair
	 */
	public Oligo(Probe pr, PrimerPair primerPair) {
		probe = pr;
		primer = primerPair;
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

package nextgen.core.capture;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

/**
 * @author prussell
 * A collection of probes that will be assigned a single common primer pair
 */
public class ProbeSet {
	
	private Collection<Probe> probes;
	
	/**
	 * Instantiate empty probe set
	 */
	public ProbeSet() {
		probes = new TreeSet<Probe>();
	}
	
	/**
	 * @param probe Single probe to initialize with
	 */
	public ProbeSet(Probe probe) {
		probes = new TreeSet<Probe>();
		probes.add(probe);
	}
	
	/**
	 * @param probeSet Collection of probes to initialize with
	 */
	public ProbeSet(Collection<Probe> probeSet) {
		probes = new TreeSet<Probe>();
		probes.addAll(probeSet);
	}
	
	/**
	 * @param probe Probe to add
	 */
	public void addProbe(Probe probe) {
		probes.add(probe);
	}
	
	/**
	 * @param other Collection of probes to add
	 */
	public void addAll(ProbeSet other) {
		probes.addAll(other.probes);
	}
	
	/**
	 * @param probe Probe to remove
	 */
	public void remove(Probe probe) {
		if(!probes.contains(probe)) {
			throw new IllegalArgumentException("Probe set does not contain probe " + probe.getID());
		}
		probes.remove(probe);
	}
	
	/**
	 * @return A view of the probes
	 */
	public Collection<Probe> getProbes() {
		return probes;
	}
	
	/**
	 * @return Iterator over probes
	 */
	public Iterator<Probe> iter() {
		return probes.iterator();
	}
	
}

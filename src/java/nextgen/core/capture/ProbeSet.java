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
	private String name;
	
	/**
	 * @param probeSetName Name of this probe set
	 */
	public ProbeSet(String probeSetName) {
		probes = new TreeSet<Probe>();
		name = probeSetName;
	}
	
	/**
	 * @param probe Single probe to initialize with
	 * @param probeSetName Name of this probe set
	 */
	public ProbeSet(Probe probe, String probeSetName) {
		probes = new TreeSet<Probe>();
		probes.add(probe);
		name = probeSetName;
	}
	
	/**
	 * @param probeSet Collection of probes to initialize with
	 * @param probeSetName Name of this probe set
	 */
	public ProbeSet(Collection<Probe> probeSet, String probeSetName) {
		probes = new TreeSet<Probe>();
		probes.addAll(probeSet);
		name = probeSetName;
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
	
	public String getName() {
		return name;
	}
	
	public void setName(String newName) {
		name = newName;
	}
	
}

package nextgen.core.capture.filter;

import java.util.Collection;

import nextgen.core.capture.ArrayFeature;
import nextgen.core.capture.Probe;
import nextgen.core.capture.ProbeSet;

/**
 * @author prussell
 *
 */
public interface ProbeFilter extends ArrayFeature {
	
	/**
	 * Do any pre-processing that is required based on the probes
	 * e.g., BLAT
	 */
	public void setup(Collection<ProbeSet> probeSets);
	
	
	/**
	 * @param probe Probe
	 * @return Whether the filter should remove the probe
	 */
	public boolean rejectProbe(Probe probe);
	
}

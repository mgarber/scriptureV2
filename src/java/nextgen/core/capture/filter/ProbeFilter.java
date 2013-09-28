package nextgen.core.capture.filter;

import nextgen.core.capture.ArrayFeature;
import nextgen.core.capture.Probe;

/**
 * @author prussell
 *
 */
public interface ProbeFilter extends ArrayFeature {
	
	/**
	 * @param probe Probe
	 * @return Whether the filter should remove the probe
	 */
	public boolean rejectProbe(Probe probe);
	
}

package nextgen.core.capture.filter;

import nextgen.core.capture.ArrayFeature;
import nextgen.core.capture.ProbeSet;
import broad.core.primer3.PrimerPair;

/**
 * @author prussell
 *
 */
public interface PrimerFilter extends ArrayFeature {
	
	/**
	 * @param primer Primer pair
	 * @param probes The probes that will potentially be connected to the primer
	 * @return Whether the filter should reject the primer for this probe set
	 */
	public boolean rejectPrimer(PrimerPair primer, ProbeSet probes);

}

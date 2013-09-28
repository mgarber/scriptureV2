package nextgen.core.capture.arrayscheme;

import nextgen.core.capture.ArrayFeature;
import nextgen.core.capture.ProbeSet;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 * Specifies a way to create probes from a transcript
 */
public interface ProbeLayout extends ArrayFeature {
	
	/**
	 * @param transcript Parent transcript
	 * @return Set of probes for the transcript
	 */
	public ProbeSet getProbes(Sequence transcript);
	
}

package nextgen.core.capture.filter;

import java.util.Collection;

import nextgen.core.capture.Probe;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFileOptionValue;

/**
 * @author engreitz
 * Catch-all filter for characters that are not bases (might indicate improper formatting of the FASTA file, for example)
 */
public class SynthesisFilter implements ProbeFilter {

	@Override
	public String name() {
		return "standard_filter";
	}

	/* (non-Javadoc)
	 * @see nextgen.core.capture.ArrayFeature#configFileLineDescription()
	 */
	@Override
	public String configFileLineDescription() {
		return "does not need to be specified in config file";
	}

	/* (non-Javadoc)
	 * @see nextgen.core.capture.ArrayFeature#validConfigFileValue(nextgen.core.pipeline.ConfigFileOptionValue)
	 */
	@Override
	public boolean validConfigFileValue(ConfigFileOptionValue value) {
		// TODO Auto-generated method stub
		return false;
	}

	/* (non-Javadoc)
	 * @see nextgen.core.capture.ArrayFeature#setParametersFromConfigFile(nextgen.core.pipeline.ConfigFileOptionValue)
	 */
	@Override
	public void setParametersFromConfigFile(ConfigFileOptionValue value) {	}

	
	@Override
	public boolean rejectProbe(Probe probe) {
		if (probe.getProbeSequence().matches("[\\p{Print}&&[^ACTG]]")) {
			throw new IllegalArgumentException("Found non-ACTG character in probe: " + probe.getID() + " " + probe.getProbeSequence());
		}
		return false;
	}
	
	@Override
	public void setup(Collection<ProbeSet> probeSets) {}

}

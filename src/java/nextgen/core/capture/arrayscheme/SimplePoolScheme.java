package nextgen.core.capture.arrayscheme;

import java.util.ArrayList;
import java.util.Collection;

import nextgen.core.capture.OligoPool;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFile;
import nextgen.core.pipeline.ConfigFileOptionValue;

import broad.core.sequence.Sequence;

/**
 * @author prussell
 * A single primer for all transcripts
 */
public class SimplePoolScheme implements PoolScheme {
	
	private ProbeLayout layout;
	
	/**
	 * 
	 */
	public SimplePoolScheme() {}
	
	/**
	 * @param probeLayout
	 */
	public SimplePoolScheme(ProbeLayout probeLayout) {
		layout = probeLayout;
	}
	
	@Override
	public Collection<ProbeSet> getProbes(Collection<Sequence> transcripts) {
		ProbeSet probeSet = new ProbeSet(name());
		for(Sequence transcript : transcripts) {
			ProbeSet probesThisSequence = layout.getProbes(transcript);
			probeSet.addAll(probesThisSequence);
		}
		Collection<ProbeSet> rtrn = new ArrayList<ProbeSet>();
		rtrn.add(probeSet);
		return rtrn;
	}

	@Override
	public String name() {
		return "simple_pool_scheme";
	}

	@Override
	public String configFileLineDescription() {
		return OligoPool.poolSchemeOptionFlag + "\t" + name();
	}

	@Override
	public boolean validConfigFileValue(ConfigFileOptionValue value) {
		return value.getActualNumValues() == 2 && value.asString(0).equals(OligoPool.poolSchemeOptionFlag) && value.asString(1).equals(name());
	}

	@Override
	public void setParametersFromConfigFile(ConfigFileOptionValue value) {
		throw new UnsupportedOperationException("Method not applicable.");
	}

	@Override
	public void setFromConfigFile(ConfigFile file) {
		ConfigFileOptionValue value = file.getSingleValue(OligoPool.arraySchemeSection, OligoPool.poolSchemeOption);
		if(!validConfigFileValue(value)) {
			throw new IllegalArgumentException("File line not valid:\n" + value.getFullOptionLine() + "\nlineformat:\n" + configFileLineDescription());
		}
		Collection<ProbeLayout> layouts = OligoPool.getProbeLayoutsFromConfigFile(file);
		if(layouts.size() != 1) {
			throw new IllegalArgumentException("To use pool scheme " + name() + " must provide exactly one probe layout in config file.");
		}
		ProbeLayout probeLayout = layouts.iterator().next();
		layout = probeLayout;
	}

	@Override
	public String toString() {
		return name() + "_" + layout.toString();
	}
	
}

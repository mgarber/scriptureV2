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
 * Stack several simple pool schemes (each simple scheme represents a single primer for all transcripts)
 */
public class StackedSimplePoolScheme implements PoolScheme {
	
	private Collection<SimplePoolScheme> schemes;
	
	/**
	 * 
	 */
	public StackedSimplePoolScheme() {}
	
	/**
	 * @param poolSchemes
	 */
	public StackedSimplePoolScheme(Collection<SimplePoolScheme> poolSchemes) {
		schemes = poolSchemes;
	}
	
	@Override
	public Collection<ProbeSet> getProbes(Collection<Sequence> transcripts) {
		Collection<ProbeSet> rtrn = new ArrayList<ProbeSet>();
		for(SimplePoolScheme scheme : schemes) {
			rtrn.addAll(scheme.getProbes(transcripts));
		}
		return rtrn;
	}

	@Override
	public String name() {
		return "stacked_simple_pool_scheme";
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
		Collection<SimplePoolScheme> simpleSchemes = new ArrayList<SimplePoolScheme>();
		for(ProbeLayout layout : layouts) {
			simpleSchemes.add(new SimplePoolScheme(layout));
		}
		schemes = simpleSchemes;
	}

	@Override
	public String toString() {
		String s = "";
		for(SimplePoolScheme scheme : schemes) {
			s += "_" + scheme.toString();
		}
		return name() + s;
	}
	
}

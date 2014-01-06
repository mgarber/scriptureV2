package nextgen.core.capture.arrayscheme;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

import nextgen.core.capture.OligoPool;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFile;
import nextgen.core.pipeline.ConfigFileOptionValue;

import broad.core.sequence.Sequence;


/**
 * @author engreitz
 * One primer per transcript
 */
public class GenePoolScheme implements PoolScheme {
	
	private ProbeLayout layout;

	private static Logger logger = Logger.getLogger(GenePoolScheme.class.getName());
	
	public GenePoolScheme() {}
	
	/**
	 * @param probeLayout
	 */
	public GenePoolScheme(ProbeLayout probeLayout) {
		layout = probeLayout;
	}
	
	@Override
	public Collection<ProbeSet> getProbes(Collection<Sequence> transcripts) {
		logger.info("");
		Collection<ProbeSet> rtrn = new ArrayList<ProbeSet>();
		for(Sequence transcript : transcripts) {
			if (transcript.getSequenceBases().toUpperCase().indexOf("GNIL") != -1) {
				throw new IllegalArgumentException("Found GNIL");
			}
			ProbeSet probesThisSequence = layout.getProbes(transcript);
			probesThisSequence.setName(layout.name() + "_" + transcript.getId());
			rtrn.add(probesThisSequence);
		}
		return rtrn;
	}

	@Override
	public String name() {
		return "gene_pool_scheme";
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
}

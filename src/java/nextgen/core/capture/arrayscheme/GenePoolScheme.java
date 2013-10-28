package nextgen.core.capture.arrayscheme;

import java.util.ArrayList;
import java.util.Collection;

import nextgen.core.capture.OligoPool;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFile;
import nextgen.core.pipeline.ConfigFileOptionValue;

import broad.core.sequence.Sequence;

/**
 * @author engreitz
 * One primer per transcript
 */
public class GenePoolScheme extends SimplePoolScheme {
	
	public GenePoolScheme() {
		super();
	}
	
	/**
	 * @param probeLayout
	 */
	public GenePoolScheme(ProbeLayout probeLayout) {
		super(probeLayout);
	}
	
	@Override
	public Collection<ProbeSet> getProbes(Collection<Sequence> transcripts) {
		Collection<ProbeSet> rtrn = new ArrayList<ProbeSet>();
		ProbeSet probeSet = new ProbeSet();
		for(Sequence transcript : transcripts) {
			if (transcript.getSequenceBases().toUpperCase().indexOf("GNIL") != -1) {
				throw new IllegalArgumentException("Found GNIL");
			}
			ProbeSet probesThisSequence = layout.getProbes(transcript);
			rtrn.add(probesThisSequence);
		}
		return rtrn;
	}

	@Override
	public String name() {
		return "gene_pool_scheme";
	}
}

package nextgen.core.capture.arrayscheme;

import java.util.Collection;

import nextgen.core.capture.ArrayFeature;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFile;

import broad.core.sequence.Sequence;

/**
 * @author prussell
 * Specifies a way to combine probe sets under different primers
 */
public interface PoolScheme extends ArrayFeature {
	
	/**
	 * @param transcripts Parent transcripts
	 * @return Probes for all transcripts
	 */
	public Collection<ProbeSet> getProbes(Collection<Sequence> transcripts);
	
	/**
	 * Configure the pool scheme from a config file
	 * @param file Config file
	 */
	public void setFromConfigFile(ConfigFile file);
		
}

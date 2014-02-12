package nextgen.core.capture.filter;

import java.util.Collection;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import broad.core.primer3.PrimerPair;

import nextgen.core.capture.OligoPool;
import nextgen.core.capture.Probe;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFileOptionValue;
import nextgen.editing.crispr.GuideRNA;
import nextgen.editing.crispr.predicate.GuideRNAPredicate;

/**
 * @author prussell
 * Filter probes with large number of single base in a row (e.g., PolyA)
 */
public class PolyBaseFilter implements ProbeFilter, PrimerFilter, GuideRNAPredicate {
	
	public String name = "PolyBase";
	private String basesToFilter;
	private int cutoff;
	private int repeatLength;
	private static Logger logger = Logger.getLogger(PolyBaseFilter.class.getName());
	
	/**
	 * 
	 */
	public PolyBaseFilter() {}

	public PolyBaseFilter(String basesToFilter, int cutoff, int repeatLength) {
		this.basesToFilter = basesToFilter;
		this.cutoff = cutoff;
		this.repeatLength = repeatLength;
	}
	
	@Override
	public String name() {
		return "poly_base_filter";
	}

	@Override
	public boolean rejectProbe(Probe probe) {
		if (probe.getProbeSequence().toUpperCase().indexOf("GNIL") != -1) {
			throw new IllegalArgumentException("found GNIL");
		}
		return rejectSequence(probe.getProbeSequence());
	}
	

	public boolean rejectSequence(String s) {
		char[] seq = s.toUpperCase().toCharArray();
		for (char c : basesToFilter.toCharArray()) {
			
			int charMatchesInWindow = 0;
			for (int i = 0; i < seq.length; i++) {
				if (seq[i] == c) {
					charMatchesInWindow++;
				}
				
				if (i >= repeatLength) {
					if (seq[i-repeatLength] == c) {
						charMatchesInWindow--;
					}
				}
				
				if (i >= repeatLength - 1) {
					if (charMatchesInWindow >= cutoff) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
	
	/**
	 * @param primer Primer pair
	 * @param probes The probes that will potentially be connected to the primer
	 * @return Whether the filter should reject the primer for this probe set
	 */
	public boolean rejectPrimer(PrimerPair primer, ProbeSet probes) {
		return rejectSequence(primer.getLeftPrimer()) || rejectSequence(primer.getRightPrimer());
	}
	
	
	@Override
	public String configFileLineDescription() {
		return OligoPool.probeFilterOptionFlag + "\t" + name() + "\t<repeat_length [5 - probe length]>\tmax_pct_repeats [0-1]>\t<repeat_bases [e.g., ACTG]>";
	}

	@Override
	public boolean validConfigFileValue(ConfigFileOptionValue value) {
		if(!value.asString(0).equals(OligoPool.probeFilterOptionFlag) && !value.asString(0).equals(OligoPool.primerFilterOptionFlag)) return false;
		if(!value.asString(1).equals(name())) return false;
		if(value.getActualNumValues() != 4 && value.getActualNumValues() != 5) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		try {

			int length = value.asInt(2);
			double pct = value.asDouble(3);
			if (pct < 0 || pct > 1) {
				logger.error("Correct config file line format: " + configFileLineDescription());
				return false;
			}
			
		} catch(NumberFormatException e) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		return true;
	}
	
	@Override
	public void setParametersFromConfigFile(ConfigFileOptionValue value) {
		if(!validConfigFileValue(value)) {
			throw new IllegalArgumentException("Config file line invalid. Line format:\n" + configFileLineDescription());
		}
		
		repeatLength = value.asInt(2);
		cutoff = (int) Math.floor(value.asDouble(3) * repeatLength);
		basesToFilter = value.asString(4);
	}
	
	@Override
	public String toString() {
		return "poly_base_filter";
	}
	
	@Override
	public void setup(Collection<ProbeSet> probeSets) {}
	

	@Override
	public boolean evaluate(GuideRNA g) {
		String seq = g.getSequenceString();
		return !rejectSequence(seq);
	}
	
	@Override
	public String getPredicateName() {
		return name;
	}
	
	@Override
	public String getShortFailureMessage(GuideRNA g) {
		return name;
	}

}

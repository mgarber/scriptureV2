package nextgen.core.capture.filter;

import java.util.Collection;

import org.apache.log4j.Logger;

import nextgen.core.capture.OligoPool;
import nextgen.core.capture.Probe;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFileOptionValue;

/**
 * @author prussell
 * Filter probes with too many repeat masked bases
 */
public class RepeatFilter implements ProbeFilter {
	
	private boolean includesLower;
	private boolean includesN;
	private double maxPct;
	private static Logger logger = Logger.getLogger(RepeatFilter.class.getName());
	
	/**
	 * 
	 */
	public RepeatFilter() {}
	
	/**
	 * @param maxRepeatPct Maximum allowable percentage of repeat masked bases (0 to 1)
	 * @param includeLowerCase Lower case bases count as repeats
	 * @param includeN Ns count as repeats
	 */
	public RepeatFilter(double maxRepeatPct, boolean includeLowerCase, boolean includeN) {
		if(maxRepeatPct < 0 || maxRepeatPct > 1) {
			throw new IllegalArgumentException("Max repeat percentage must be between 0 and 1");
		}
		if(!includeLowerCase && !includeN) {
			throw new IllegalArgumentException("Must include at least one: lowercase or N");
		}
		maxPct = maxRepeatPct;
		includesLower = includeLowerCase;
		includesN = includeN;
	}
	
	@Override
	public String name() {
		return "repeat_filter";
	}

	@Override
	public boolean rejectProbe(Probe probe) {
		String probeSeq = probe.getProbeSequence();
		int size = probeSeq.length();
		int numRepeats = 0;
		for(int i=0; i<size; i++) {
			char c = probeSeq.charAt(i);
			if(includesLower && !Character.isUpperCase(c)) {
				numRepeats++;
				continue;
			}
			if(includesN && (c == 'n' || c == 'N')) {
				numRepeats++;
				continue;
			}
		}
		double pct = (double) numRepeats / (double) size;
		return pct > maxPct;
	}
	
	private static String INCLUDE_N_FLAG = "N";
	private static String INCLUDE_LOWERCASE_FLAG = "lower";

	@Override
	public String configFileLineDescription() {
		return OligoPool.probeFilterOptionFlag + "\t" + name() + "\t<max_pct_repeats (0-1)>\t<" + INCLUDE_N_FLAG + ">\tand/or<" + INCLUDE_LOWERCASE_FLAG + ">";
	}

	@Override
	public boolean validConfigFileValue(ConfigFileOptionValue value) {
		if(!value.asString(0).equals(OligoPool.probeFilterOptionFlag)) return false;
		if(!value.asString(1).equals(name())) return false;
		if(value.getActualNumValues() != 4 && value.getActualNumValues() != 5) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		try {
			double pct = value.asDouble(2);
			if(pct < 0 || pct > 1) {
				logger.error("Correct config file line format: " + configFileLineDescription());
				return false;
			}
		} catch(NumberFormatException e) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		for(int i = 3; i < value.getActualNumValues(); i++) {
			String s = value.asString(i);
			if(!s.equals(INCLUDE_LOWERCASE_FLAG) && !s.equals(INCLUDE_N_FLAG)) {
				logger.error("Correct config file line format: " + configFileLineDescription());
				return false;
			}
		}
		return true;
	}

	@Override
	public void setParametersFromConfigFile(ConfigFileOptionValue value) {
		if(!validConfigFileValue(value)) {
			throw new IllegalArgumentException("Config file line invalid. Line format:\n" + configFileLineDescription());
		}
		maxPct = value.asDouble(2);
		includesLower = false;
		includesN = false;
		for(int i = 3; i < value.getActualNumValues(); i++) {
			String s = value.asString(i);
			if(s.equals(INCLUDE_LOWERCASE_FLAG)) {
				includesLower = true;
			}
			if(s.equals(INCLUDE_N_FLAG)) {
				includesN = true;
			}
		}		
	}
	
	@Override
	public String toString() {
		String r = "";
		if(includesLower) r += "_lower";
		if(includesN) r += "_N";
		return "repeat_filter_" + maxPct + r;
	}
	
	@Override
	public void setup(Collection<ProbeSet> probeSets) {}

}

package nextgen.core.capture.filter;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.primer3.PrimerPair;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import nextgen.core.capture.OligoPool;
import nextgen.core.capture.Probe;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFileOptionValue;
import nextgen.editing.crispr.GuideRNA;
import nextgen.editing.crispr.predicate.GuideRNAPredicate;




/**
 * @author engreitz
 * Filter out low-complexity sequences.  Simple algorithm:  search probe sequences for perfect matches of
 * length n to all possible low-complexity repeats of length 1 to k by creating a hash of sequences of length n
 * from low-complexity repeats and comparing each n-mer in the probe sequence to this hash.
 */
public class LowComplexityFilter implements ProbeFilter, PrimerFilter, GuideRNAPredicate  {

	public String name = "LowComplexity";
	
	// Setting K = 6 will check for sequencing containing repeats of up to k = 6 bases (e.g., (GGCACA GGCACA ...))
	private int maxK = 6;
	
	// Sets the length of the repeat that must be present to call a match
	private int[] minLength = new int[] { 10, 10, 12, 12, 18, 18 };
	
	private Map<Integer, HashSet<String>> hashes = null;
	
	private static Logger logger = Logger.getLogger(PolyBaseFilter.class.getName());
	
	
	/**
	 * Constructor with default settings
	 */
	public LowComplexityFilter() {
		hashes = getLowComplexityHashes(maxK, minLength);
	}
	
	
	/**
	 * @param k
	 * @param minLengths
	 * @return
	 * Set up hash sets containing all strings of minLengths (n) bases from k-mer repeats
	 */
	private Map<Integer, HashSet<String>> getLowComplexityHashes(int k, int[] minLengths) {
		Map<Integer, HashSet<String>> hashes = new TreeMap<Integer, HashSet<String>>();
		for (int repeatK = 0; repeatK < k; repeatK++) {
			HashSet<String> currHash = hashes.containsKey(minLengths[repeatK]) ?
				currHash = hashes.get(minLengths[repeatK]) : new HashSet<String>();
	
			List<String> allRepeats = new ArrayList<String>();
			addAllKmerRepeats(allRepeats, "", repeatK+1);
			for (String repeatElement : allRepeats) {
				addRepeatToHash(currHash, repeatElement, minLengths[repeatK]);
			}
			
			hashes.put(minLengths[repeatK], currHash);
		}
		return hashes;
	}
	
	
	final char chars[] = new char[] {'A','C','T','G'};
	/**
	 * @param allRepeats
	 * @param base
	 * @param k
	 * Recursive function to generate all possible k-mer sequences
	 */
	private void addAllKmerRepeats(List<String> allRepeats, String base, int k) {
		if (base.length() >= k) {
			//logger.info("Adding " + base + " to list.");
			allRepeats.add(base);
		} else {
			for (int i = 0; i < chars.length; i++)
				addAllKmerRepeats(allRepeats, base + chars[i], k);
		}
	}
	
	
	/**
	 * @param hash
	 * @param repeatElement
	 * @param hashLength
	 * Take a repeat element (e.g., GGA) and add all strings of hashLength bases to the hash
	 */
	private void addRepeatToHash(HashSet<String> hash, String repeatElement, int hashLength) {
		repeatElement = StringUtils.repeat(repeatElement, hashLength*2);
		//logger.info("Adding " + repeatElement);
		for (int i = 0; i < hashLength; i++) {
			hash.add(repeatElement.substring(i,hashLength));
		}
	}

	
	@Override
	public String name() {
		return "low_complexity_filter";
	}

	@Override
	public boolean rejectProbe(Probe probe) {
		return rejectSequence(probe.getProbeSequence());
	}
	

	public boolean rejectSequence(String s) {
		for (int period = 0; period < maxK; period++) {
			if (sequenceContainsNmerInHash(s, minLength[period])) return true;
		}
		return false;
	}


	private boolean sequenceContainsNmerInHash(String seq, int n) {
		for (int charIndex = 0; charIndex < seq.length() - n + 1; charIndex++) {
			String substr = seq.substring(charIndex, charIndex + n);
			if (hashes.get(n).contains(substr)) return true;
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
		return OligoPool.probeFilterOptionFlag + "\t" + name();
	}

	@Override
	public boolean validConfigFileValue(ConfigFileOptionValue value) {
		if(!value.asString(0).equals(OligoPool.probeFilterOptionFlag) && !value.asString(0).equals(OligoPool.primerFilterOptionFlag)) return false;
		if(!value.asString(1).equals(name())) return false;
		if(value.getActualNumValues() != 2) {
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
		// config file currently does not control any parameters
	}
	
	@Override
	public String toString() {
		return name();
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
		return name();
	}
	
	@Override
	public String getShortFailureMessage(GuideRNA g) {
		return name();
	}

	
	/**
	 * Main argument for testing purposes:  Given a FASTA file containing sequences, output a subset of the records that do pass the filters.
	 * @param args
	 * @return
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input FASTA file", true);
		p.addStringArg("-o", "Output FASTA file", true);
		p.addStringArg("-f", "Filtered FASTA file", true);
		p.parse(args);
		String inFile = p.getStringArg("-i");
		String outFile = p.getStringArg("-o");
		String filteredFile = p.getStringArg("-f");
		
		LowComplexityFilter filter = new LowComplexityFilter();
		Map<String,Sequence> inSequences = FastaSequenceIO.getChrSequencesFromFasta(inFile);
		List<Sequence> outSequences = new ArrayList<Sequence>();
		List<Sequence> filteredSequences = new ArrayList<Sequence>();
		for (String key : inSequences.keySet()) {
			if (!filter.rejectSequence(inSequences.get(key).getSequenceBases())) {
				outSequences.add(inSequences.get(key));
			} else {
				filteredSequences.add(inSequences.get(key));
			}
		}
		
		new FastaSequenceIO(outFile).write(outSequences);
		new FastaSequenceIO(filteredFile).write(filteredSequences);
		logger.info(outSequences.size() + " out of " + inSequences.size() + " passed the low complexity filter.");
	}
	
}

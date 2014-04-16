package nextgen.core.capture;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.capture.arrayscheme.*;
import nextgen.core.capture.filter.*;
import nextgen.core.general.TabbedReader;
import nextgen.core.pipeline.ConfigFile;
import nextgen.core.pipeline.ConfigFileOption;
import nextgen.core.pipeline.ConfigFileOptionValue;
import nextgen.core.pipeline.ConfigFileSection;
import broad.core.error.ParseException;
import broad.core.parser.CommandLineParser;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.PrimerUtils;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;


/**
 * @author prussell
 * An oligo pool specified by probe design and pooling schemes
 * Incorporates filters for probes and primers
 */
public class OligoPool {
	
	private static Logger logger = Logger.getLogger(OligoPool.class.getName());
	private ConfigFile configFile;
	private PoolScheme poolScheme;
	private Collection<Sequence> transcripts;
	private List<ProbeFilter> probeFilters;
	private Collection<PrimerFilter> primerFilters;
	private Map<String, List<FullDesignEntry>> probeSetDesign;
	private int primerSize;
	private String primer3corePath;
	/**
	 * BufferedReader for file containing list of primers in format produced by PrimerPair.getPrimerFieldsAsStringForConstructor(), or null if not using
	 */
	private BufferedReader primerReader;
	private double optimalTm;
	
	/**
	 * @param config Config file
	 * @throws IOException 
	 */
	private OligoPool(String config) throws IOException {
		logger.info("");
		logger.info("Instantiating oligo pool...");
		logger.info("");
		logger.info("Establishing config file...");
		configFile = createConfigFile(config);
		logger.info("");
		logger.info("Getting transcripts...");
		transcripts = getTranscriptsFromConfigFile();
		logger.info("");
		logger.info("Establishing array scheme...");
		poolScheme = getPoolSchemeFromConfigFile();
		primerSize = getPrimerSizeFromConfigFile();
		primer3corePath = getPrimer3PathFromConfigFile();
		optimalTm = getOptimalTmFromConfigFile();
		primerReader = getPrimerReaderFromConfigFile();
		logger.info("");
		logger.info("Getting probe filters...");
		probeFilters = getProbeFiltersFromConfigFile();
		logger.info("");
		logger.info("Getting primer filters...");
		primerFilters = getPrimerFiltersFromConfigFile();
		logger.info("");
		logger.info("Done instantiating oligo pool.");
	}
	
	/**
	 * Probe layout option flag
	 */
	public static String probeLayoutOptionFlag = "probe_layout";
	private static ConfigFileOption probeLayoutOption = new ConfigFileOption(probeLayoutOptionFlag, 6, true, true, true);
	/**
	 * Pool scheme option flag
	 */
	public static String poolSchemeOptionFlag = "pool_scheme";
	/**
	 * Pool scheme option
	 */
	public static ConfigFileOption poolSchemeOption = new ConfigFileOption(poolSchemeOptionFlag, 6, true, false, true);
	private static String arraySchemeSectionFlag = "Array_scheme";
	/**
	 * Array scheme section
	 */
	public static ConfigFileSection arraySchemeSection = new ConfigFileSection(arraySchemeSectionFlag, true);
	private static String sequenceFastaOptionFlag = "sequence_fasta";
	public static ConfigFileOption sequenceFastaOption = new ConfigFileOption(sequenceFastaOptionFlag, 2, false, false, true);
	private static String sequencesSectionFlag = "Sequences";
	public static ConfigFileSection sequencesSection = new ConfigFileSection(sequencesSectionFlag, true);
	private static String primerSizeOptionFlag = "primer_size";
	private static ConfigFileOption primerSizeOption = new ConfigFileOption(primerSizeOptionFlag, 2, false, false, true);
	private static String primer3corePathOptionFlag = "primer3_core_path";
	private static ConfigFileOption primer3corePathOption = new ConfigFileOption(primer3corePathOptionFlag, 2, false, false, true);
	private static String optimalTmOptionFlag = "optimal_tm_for_primers";
	private static ConfigFileOption optimalTmOption = new ConfigFileOption(optimalTmOptionFlag, 2, false, false, true);
	private static String primerFileOptionFlag = "primer_file";
	private static ConfigFileOption primerFileOption = new ConfigFileOption(primerFileOptionFlag, 2, false, false, false);
	private static String flank5primeOptionFlag = "extra_flank_5prime";
	private static ConfigFileOption flank5primeOption = new ConfigFileOption(flank5primeOptionFlag, 2, false, false, false, "");
	private static String flank3primeOptionFlag = "extra_flank_3prime";
	private static ConfigFileOption flank3primeOption = new ConfigFileOption(flank3primeOptionFlag, 2, false, false, false, "");
	/**
	 * Probe filter option flag
	 */
	public static String probeFilterOptionFlag = "probe_filter";
	private static ConfigFileOption probeFilterOption = new ConfigFileOption(probeFilterOptionFlag, 5, true, true, false);
	private static String probeFiltersSectionFlag = "Probe_filters";
	private static ConfigFileSection probeFiltersSection = new ConfigFileSection(probeFiltersSectionFlag, false);
	public static String primerFilterOptionFlag = "primer_filter";
	private static ConfigFileOption primerFilterOption = new ConfigFileOption(primerFilterOptionFlag, 5, true, true, false);
	private static String primerFiltersSectionFlag = "Primer_filters";
	private static ConfigFileSection primerFiltersSection = new ConfigFileSection(primerFiltersSectionFlag, false);
	
	private static ConfigFile createConfigFile(String fileName) throws IOException {
		logger.info("Getting config file from " + fileName);
		arraySchemeSection.addAllowableOption(probeLayoutOption);
		arraySchemeSection.addAllowableOption(poolSchemeOption);
		arraySchemeSection.addAllowableOption(primerSizeOption);
		arraySchemeSection.addAllowableOption(primer3corePathOption);
		arraySchemeSection.addAllowableOption(optimalTmOption);
		arraySchemeSection.addAllowableOption(primerFileOption);
		arraySchemeSection.addAllowableOption(flank3primeOption);
		arraySchemeSection.addAllowableOption(flank5primeOption);
		sequencesSection.addAllowableOption(sequenceFastaOption);
		probeFiltersSection.addAllowableOption(probeFilterOption);
		primerFiltersSection.addAllowableOption(primerFilterOption);
		Collection<ConfigFileSection> sections = new ArrayList<ConfigFileSection>();
		sections.add(arraySchemeSection);
		sections.add(sequencesSection);
		sections.add(probeFiltersSection);
		sections.add(primerFiltersSection);
		return new ConfigFile(sections, fileName);
	}
	
	private String getFlank5primeFromConfigFile() {
		return configFile.getSingleValueField(arraySchemeSection, flank5primeOption);
	}
	
	private String getFlank3primeFromConfigFile() {
		return configFile.getSingleValueField(arraySchemeSection, flank3primeOption);
	}
	
	private Collection<Sequence> getTranscriptsFromConfigFile() throws IOException {
		ConfigFileOptionValue seqsValue = configFile.getSingleValue(sequencesSection, sequenceFastaOption);
		String seqsFasta = seqsValue.asString(1);
		FastaSequenceIO fsio = new FastaSequenceIO(seqsFasta);
		Collection<Sequence> rtrn = fsio.loadAll();
		logger.info("Got " + rtrn.size() + " transcripts.");
		return rtrn;
	}
	
	private static ProbeFilter getProbeFilterFromConfigFileValue(ConfigFileOptionValue value) {
		
		// Add additional filter classes to this array:
		ProbeFilter[] filters = new ProbeFilter[] { new RepeatFilter(), new PolyBaseFilter(), new LowComplexityFilter() };
		
		for (ProbeFilter filter : filters) {
			if(filter.validConfigFileValue(value)) {
				filter.setParametersFromConfigFile(value);
				logger.info("Got filter " + filter.toString());
				return filter;
			}
		}
		
		throw new IllegalArgumentException("Probe filter not connected to this method: " + value.getFullOptionLine() + ". Need to implement.");
	}

	private static PrimerFilter getPrimerFilterFromConfigFileValue(ConfigFileOptionValue value) {
		
		// Add additional filter classes to this array:
		PrimerFilter[] filters = new PrimerFilter[] { new PolyBaseFilter() };
		
		for (PrimerFilter filter : filters) {
			if(filter.validConfigFileValue(value)) {
				filter.setParametersFromConfigFile(value);
				logger.info("Got filter " + filter.toString());
				return filter;
			}
		}
		
		throw new IllegalStateException("Primer filter not connected to this method: " + value.getFullOptionLine() + ". Need to implement.");
	}
	
	private static ProbeLayout getProbeLayoutFromConfigFileValue(ConfigFileOptionValue value) {
		// Single tiling probe layout
		SingleTilingProbeLayout singleTilingProbeLayout = new SingleTilingProbeLayout();
		if(singleTilingProbeLayout.validConfigFileValue(value)) {
			singleTilingProbeLayout.setParametersFromConfigFile(value);
			logger.info("Got layout " + singleTilingProbeLayout.toString());
			return singleTilingProbeLayout;
		}
		// Probe layout 2
		
		// Probe layout 3
		
		// ...
		throw new IllegalArgumentException("Invalid config file line: " + value.getFullOptionLine());
	}
	
	/**
	 * @param file Config file
	 * @return Collection of all probe layouts specified in the config file
	 */
	public static Collection<ProbeLayout> getProbeLayoutsFromConfigFile(ConfigFile file) {
		Collection<ProbeLayout> rtrn = new ArrayList<ProbeLayout>();
		Collection<ConfigFileOptionValue> layoutVals = file.getOptionValues(arraySchemeSection, probeLayoutOption);
		for(ConfigFileOptionValue val : layoutVals) {
			rtrn.add(getProbeLayoutFromConfigFileValue(val));
		}
		return rtrn;
	}
	
	private PoolScheme getPoolSchemeFromConfigFile() {
		ConfigFileOptionValue poolSchemeVal = configFile.getSingleValue(arraySchemeSection, poolSchemeOption);
		
		// Add additional scheme classes to this array:
		PoolScheme[] schemes = new PoolScheme[] { new SimplePoolScheme(), new StackedSimplePoolScheme(), new GenePoolScheme(), new GroupedStackedPoolScheme(), new MultipleLayoutGenePoolScheme() };
		
		for (PoolScheme scheme : schemes) {
			if (scheme.validConfigFileValue(poolSchemeVal)) {
				scheme.setFromConfigFile(configFile);
				logger.info("Got pool scheme " + scheme.name());
				return scheme;
			}
		}
		
		throw new IllegalArgumentException("Pool scheme not connected to this method: " + poolSchemeVal.getFullOptionLine() + ". Need to implement.");
	}
	
	private int getPrimerSizeFromConfigFile() {
		ConfigFileOptionValue primerSizeVal = configFile.getSingleValue(arraySchemeSection, primerSizeOption);
		int rtrn = primerSizeVal.asInt(1);
		logger.info("Primer size is " + rtrn + ".");
		return rtrn;
	}
	
	private String getPrimer3PathFromConfigFile() {
		ConfigFileOptionValue pathVal = configFile.getSingleValue(arraySchemeSection, primer3corePathOption);
		String rtrn = pathVal.asString(1);
		logger.info("Primer3 executable is " + rtrn + ".");
		return rtrn;
	}
	
	private BufferedReader getPrimerReaderFromConfigFile() throws FileNotFoundException {
		if(!configFile.hasOption(arraySchemeSection, primerFileOption)) {
			return null;
		}
		ConfigFileOptionValue val = configFile.getSingleValue(arraySchemeSection, primerFileOption);
		String name = val.asString(1);
		FileReader r = new FileReader(name);
		return new BufferedReader(r);
	}
	
	private double getOptimalTmFromConfigFile() {
		ConfigFileOptionValue val = configFile.getSingleValue(arraySchemeSection, optimalTmOption);
		double rtrn = val.asDouble(1);
		logger.info("Optimal TM is " + rtrn + ".");
		return rtrn;
	}
	
	private List<ProbeFilter> getProbeFiltersFromConfigFile() {
		List<ProbeFilter> rtrn = new ArrayList<ProbeFilter>();
		Collection<ConfigFileOptionValue> probeFilterVals = configFile.getOptionValues(probeFiltersSection, probeFilterOption);
		if(probeFilterVals == null) {
			return rtrn;
		}
		for(ConfigFileOptionValue val : probeFilterVals) {
			rtrn.add(getProbeFilterFromConfigFileValue(val));
		}
		rtrn.add(new SynthesisFilter());
		return rtrn;
	}
	
	private Collection<PrimerFilter> getPrimerFiltersFromConfigFile() {
		Collection<PrimerFilter> rtrn = new ArrayList<PrimerFilter>();
		Collection<ConfigFileOptionValue> primerFilterVals = configFile.getOptionValues(primerFiltersSection, primerFilterOption);
		if(primerFilterVals == null) {
			return rtrn;
		}
		for(ConfigFileOptionValue val : primerFilterVals) {
			rtrn.add(getPrimerFilterFromConfigFileValue(val));
		}
		return rtrn;
	}
	
	
	/**
	 * @param outFilePrefix
	 * @param changePrimersInputFile
	 * @throws IOException
	 * Main scripting function to create/read oligos and assign primers
	 */
	private void createOligos(String outFilePrefix, String changePrimersInputFile, boolean reassignPrimers) throws IOException {
		
		probeSetDesign = new TreeMap<String, List<FullDesignEntry>>();
		
		if (changePrimersInputFile != null) {	
			// Read in oligos from file
			logger.info("Reading oligos from file: " + changePrimersInputFile);
			CloseableIterator<FullDesignEntry> itr = TabbedReader.read(new File(changePrimersInputFile), FullDesignEntry.class, new Factory(), 1);  // 1 = number of rows to skip at beginning of file
			
			while (itr.hasNext()) {
				FullDesignEntry entry = itr.next();
				if (probeSetDesign.containsKey(entry.probeSetName)) {
					probeSetDesign.get(entry.probeSetName).add(entry);
				} else {
					List<FullDesignEntry> probeSet = new ArrayList<FullDesignEntry>();
					probeSet.add(entry);
					probeSetDesign.put(entry.probeSetName, probeSet);
				}
			}
		} else {
			// Create probes
			Collection<ProbeSet> probeSets = createProbes(outFilePrefix);
			probeSets = filterProbes(probeSets, outFilePrefix);
			// Convert oligos to text probe sets for easy input/output after this point
			for (ProbeSet probeSet : probeSets) {
				if (!probeSet.getProbes().iterator().hasNext()) continue;
				
				List<FullDesignEntry> entries = new ArrayList<FullDesignEntry>();
				for (Probe probe : probeSet.getProbes()) {
					entries.add(new FullDesignEntry(probe, probeSet.getName(), probeSet.getProbes().size()));
				}
				probeSetDesign.put(probeSet.getName(), entries);
			}		
		}
		
		// Add extra flanking sequences if specified in config file
		String flank5p = getFlank5primeFromConfigFile();
		String flank3p = getFlank3primeFromConfigFile();
		
		if(flank5p != "" || flank3p != "") {
			addFlankingSequences(flank5p, flank3p);
		}
		
		if (reassignPrimers || changePrimersInputFile == null) {
			assignPrimers();
		}
	}
	
		
	/**
	 * @param outFilePrefix
	 * @return
	 * @throws IOException
	 * Generate and filter probe sequences
	 */
	private Collection<ProbeSet> createProbes(String outFilePrefix) throws IOException {
		logger.info("");
		logger.info("Creating oligos...");
		Collection<ProbeSet> probeSets = poolScheme.getProbes(transcripts);
		int numProbes = 0;
		for(ProbeSet probeSet : probeSets) {
			numProbes += probeSet.getProbes().size();
		}
		logger.info("There are " + probeSets.size() + " probe sets with a total of " + numProbes + " probes.");
		return probeSets;
	}
	
	
	private Collection<ProbeSet> filterProbes(Collection<ProbeSet> probeSets, String outFilePrefix) throws IOException {
		// Filter probes
		logger.info("Filtering probes...");
		
		// Output statistics about which filters are removing probes
		// TODO:  Move this to the ProbeLayout or ProbeSet code so you can control how to aggregate the stats
		FileWriter w = new FileWriter(outFilePrefix + "_filter_results.out");
		w.write("ProbeSet\tpassed");
		for (ProbeFilter filter : probeFilters) {
			w.write("\t" + filter.toString());
			filter.setup(probeSets);
		}
		w.write("\n");
		
		int removed = 0;
		int remaining = 0;
		for(ProbeSet probeSet : probeSets) {
			if(!probeSet.getProbes().iterator().hasNext()) {
				continue;
			}
			w.write(probeSet.getProbes().iterator().next().getID());
			Iterator<Probe> iter = probeSet.iter();
			
			int[] counts = new int[probeFilters.size()+1];
			
			while(iter.hasNext()) {
				Probe probe = iter.next();
				boolean rejected = false;
				
				for(int i = 0; i < probeFilters.size(); i++) {
					ProbeFilter filter = probeFilters.get(i);
					if(filter.rejectProbe(probe)) {
						rejected = true;
						counts[i+1]++;
					}
				}
				
				if (rejected) {
					iter.remove();
					removed++;
				} else {
					remaining++;
					counts[0]++;
				}
			}
			for (int count : counts) { 
				w.write("\t" + count);
			}
			w.write("\n");
		}
		w.close();
		logger.info("Done filtering probes. Removed " + removed + " probes. " + remaining + " probes remain.");
		
		return probeSets;
	}
		
	
	private void addFlankingSequences(String flank5prime, String flank3prime) {
		logger.info("Adding flanking sequences " + flank5prime + " and " + flank3prime);
		for(List<FullDesignEntry> entries : probeSetDesign.values()) {
			for(FullDesignEntry oligo : entries) {
				oligo.setFlank5prime(flank5prime);
				oligo.setFlank3prime(flank3prime);
			}
		}
	}
	
	/**
	 * @throws IOException
	 * Assign primers to existing probesets
	 */
	private void assignPrimers() throws IOException {
		// Assign and filter primers
		logger.info("Assigning primers to " + probeSetDesign.keySet().size() + " probe sets...");
		int tried = 0;
		int succeeded = 0;
		for(List<FullDesignEntry> entries : probeSetDesign.values()) {
			boolean foundPrimer = false;
			while(!foundPrimer) {
				boolean rejected = false;
				PrimerPair primer = PrimerUtils.getOneSyntheticPrimerPair(primerSize, primer3corePath, optimalTm, primerReader, null);
				tried++;
				for(PrimerFilter filter : primerFilters) {
					if(filter.rejectPrimer(primer, entries)) {
						rejected = true;
					}
				}
				if(!rejected) {
					// Primer is OK
					for(FullDesignEntry oligo : entries) {
						oligo.leftPrimer = primer.getLeftPrimer();
						oligo.rightPrimer = primer.getRightPrimer();
						oligo.oligoSequence = oligo.getFlank5Prime() + primer.getLeftPrimer().toUpperCase() + oligo.probeSequence.toUpperCase() + Sequence.reverseSequence(primer.getRightPrimer()).toUpperCase() + oligo.getFlank3Prime();
					}
					succeeded++;
					foundPrimer = true;
				}
			}
		}
		logger.info("Done assigning primers. Successfully assigned " + succeeded + " primer pairs after trying " + tried + ".");
		logger.info("Done creating oligos.");
	}
	
	
	private void writeFiles(String outFilePrefix) throws IOException {
		logger.info("");
		logger.info("Writing output files...");
		writeProbeFasta(outFilePrefix + "_probes.fa");
		writeProbeList(outFilePrefix + "_probes.out");
		writeOligoFasta(outFilePrefix + "_oligos.fa");
		writeOligoList(outFilePrefix + "_oligos.out");
		writePrimers(outFilePrefix + "_primers.out");
		writeFullTable(outFilePrefix + "_full_design.out");
		logger.info("Done writing files.");
	}
	
	private void writeProbeFasta(String outFile) throws IOException {
		logger.info("Writing fasta file of probe sequences to file " + outFile);
		List<Sequence> probeSeqs = new ArrayList<Sequence>();
		for (List<FullDesignEntry> probeSet : probeSetDesign.values()) {
			for (FullDesignEntry entry : probeSet) {
				Sequence seq = new Sequence(entry.probeId);
				seq.setSequenceBases(entry.probeSequence);
				probeSeqs.add(seq);
			}
		}
		FastaSequenceIO fsio = new FastaSequenceIO(outFile);
		fsio.write(probeSeqs, 200);
	}
	
	private void writeProbeList(String outFile) throws IOException {
		logger.info("Writing list of probe sequences to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		for (List<FullDesignEntry> probeSet : probeSetDesign.values()) {
			for (FullDesignEntry entry : probeSet) {
				w.write(entry.probeSequence + "\n");
			}
		}
		w.close();
	}
	
	private void writeOligoFasta(String outFile) throws IOException {
		logger.info("Writing fasta file of full oligo sequences to file " + outFile);
		List<Sequence> oligoSeqs = new ArrayList<Sequence>();
		for (List<FullDesignEntry> probeSet : probeSetDesign.values()) {
			for (FullDesignEntry entry : probeSet) {
				Sequence seq = new Sequence(entry.probeId);
				seq.setSequenceBases(entry.oligoSequence);
				oligoSeqs.add(seq);
			}
		}
		FastaSequenceIO fsio = new FastaSequenceIO(outFile);
		fsio.write(oligoSeqs, 200);
	}
	
	private void writeOligoList(String outFile) throws IOException {
		logger.info("Writing list of full oligo sequences to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		
		for (List<FullDesignEntry> probeSet : probeSetDesign.values()) {
			for (FullDesignEntry entry : probeSet) {
				if (entry.oligoSequence.toUpperCase().indexOf("GNIL") != -1) {
					w.close();
					throw new IllegalArgumentException("found GNIL");
				}
				w.write(entry.oligoSequence + "\n");
			}	
		}
		w.close();
	}
	
	private void writePrimers(String outFile) throws IOException {
		logger.info("Writing primers to file " + outFile);

		String header = "Probe_set\t";
		header += "Left_primer\t";
		header += "Right_primer\t";
		FileWriter w = new FileWriter(outFile);
		w.write(header + "\n");
		for(String probeSetName : probeSetDesign.keySet()) {
			// assumes that all probes in the probe set correctly have the same primer assigned
			String line = probeSetName + "\t";
			line += probeSetDesign.get(probeSetName).get(0).leftPrimer + "\t";
			line += probeSetDesign.get(probeSetName).get(0).rightPrimer + "\t";
			w.write(line + "\n");
		}
		w.close();
	}
	
	private void writeFullTable(String outFile) throws IOException {
		logger.info("Writing full design table to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		
		w.write(FullDesignEntry.header + "\n");
		for (List<FullDesignEntry> list : probeSetDesign.values()) {
			for (FullDesignEntry entry : list) {
				w.write(entry.toString() + "\n");
			}	
		}
		w.close();
	}
	
	
	/**
	 * @author engreitz
	 * This class contains all of the information written to the full design output file in String/integer format, allowing for
	 * input/output of this design file without creating all of the objects contained by Oligo (e.g., PrimerPair, ProbeLayout, etc.)
	 */
	public class FullDesignEntry {
		public String probeId, parentTranscriptId, probeSetName, senseOrAntisense, probeLayoutString, leftPrimer, rightPrimer, flank5prime, flank3prime, probeSequence, oligoSequence;
		public int probeSetSize, start, end;
		
		// Note:  Do not change the header column names without also changing the R scripts that read them
		public static final String header = 	"Probe_ID\t" +
										"Parent_sequence\t" +
										"Probe_set\t" +
										"Probe_set_size\t" +
										"Start\t" +
										"End\t" +
										"Orientation\t" +
										"Probe_layout\t" +
										"Left_primer\t" +
										"Right_primer\t" +
										"Probe_sequence\t" +
										"Full_oligo";
		
		public FullDesignEntry(Probe probe, String probeSetName, int probeSetSize) {
			probeId = probe.getID();
			parentTranscriptId = probe.getParentTranscriptId();
			this.probeSetName = probeSetName;
			this.probeSetSize = probeSetSize;
			start = probe.getStartPosOnTranscript();
			end = probe.getEndPosOnTranscript() - 1;
			senseOrAntisense = probe.isAntisenseToTranscript() ? "antisense" : "sense";
			probeLayoutString = probe.getProbeLayoutString();
			leftPrimer = "none";
			rightPrimer = "none";
			flank5prime = "";
			flank3prime = "";
			probeSequence = probe.getProbeSequence();
			oligoSequence = "none";
		}
		
		public void setFlank3prime(String flank3p) {
			flank3prime = flank3p;
		}

		public void setFlank5prime(String flank5p) {
			flank5prime = flank5p;
		}

		public String getFlank3Prime() {
			return flank3prime;
		}

		public String getFlank5Prime() {
			return flank5prime;
		}

		public FullDesignEntry(String probeId, String parentTranscriptId, String probeSetName, int probeSetSize, int start, int end, String senseOrAntisense, String probeLayoutString, String leftPrimer, String rightPrimer, String probeSequence, String oligoSequence) {
			this.probeId = probeId;
			this.parentTranscriptId = parentTranscriptId;
			this.probeSetName = probeSetName;
			this.probeSetSize = probeSetSize;
			this.start = start;
			this.end = end;
			this.senseOrAntisense = senseOrAntisense;
			this.probeLayoutString = probeLayoutString;
			this.leftPrimer = leftPrimer;
			this.rightPrimer = rightPrimer;
			this.flank5prime = "";
			this.flank3prime = "";
			this.probeSequence = probeSequence;
			this.oligoSequence = oligoSequence;
		}
		
		public String toTabbedString() {
			return probeId + "\t" + 
					parentTranscriptId + "\t" + 
					probeSetName + "\t" + 
					probeSetSize + "\t" + 
					start + "\t" + 
					end + "\t" + 
					senseOrAntisense + "\t" + 
					probeLayoutString + "\t" + 
					leftPrimer + "\t" + 
					rightPrimer + "\t" +
					probeSequence + "\t" + 
					oligoSequence;
		}
		public String toString() { return toTabbedString(); }
	}
	
	public class Factory implements TabbedReader.Factory<FullDesignEntry> {
		public FullDesignEntry create(String[] rawFields) throws ParseException {
			if (rawFields.length != 12) throw new ParseException("Incorrect number of fields when reading full design file");
			return new FullDesignEntry(rawFields[0], 
						rawFields[1], 
						rawFields[2], 
						Integer.parseInt(rawFields[3]),
						Integer.parseInt(rawFields[4]),
						Integer.parseInt(rawFields[5]),
						rawFields[6],
						rawFields[7],
						rawFields[8],
						rawFields[9],
						rawFields[10],
						rawFields[11]);
		}
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-c", "Config file", true);
		p.addStringArg("-o", "Output file prefix", true);
		p.addStringArg("-l", "Logger level", false);
		p.addStringArg("-i", "Path to existing design; if using this argument, this designer will use these oligos rather than generating new ones from the FASTA file", false);
		p.addBooleanArg("-p", "Flag to re-assign primers to probesets; must specify -i input", false, false);
		
		p.parse(args);
		String configFile = p.getStringArg("-c");
		String outFilePrefix = p.getStringArg("-o");
		String levelString = p.getStringArg("-l");
		String changePrimersInputFile = p.getStringArg("-i");
		boolean reassignPrimers = p.getBooleanArg("-p");
		
		Level level = Level.toLevel(levelString, Level.INFO);
		logger.setLevel(level);

		// Instantiate OligoPool from config file
		OligoPool oligoPool = new OligoPool(configFile);
		
		// Create oligos
		oligoPool.createOligos(outFilePrefix, changePrimersInputFile, reassignPrimers);
		
		// Write files
		oligoPool.writeFiles(outFilePrefix);
		
		// Close primer file reader
		if(oligoPool.primerReader != null) {
			oligoPool.primerReader.close();
		}
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}

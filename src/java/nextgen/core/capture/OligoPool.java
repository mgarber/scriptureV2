package nextgen.core.capture;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import nextgen.core.capture.arrayscheme.*;
import nextgen.core.capture.filter.*;
import nextgen.core.pipeline.ConfigFile;
import nextgen.core.pipeline.ConfigFileOption;
import nextgen.core.pipeline.ConfigFileOptionValue;
import nextgen.core.pipeline.ConfigFileSection;

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
	private Map<PrimerPair,Integer> primerCounts;
	private Collection<Oligo> oligos;
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
		oligos = new TreeSet<Oligo>();
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
	private static ConfigFileOption sequenceFastaOption = new ConfigFileOption(sequenceFastaOptionFlag, 2, false, false, true);
	private static String sequencesSectionFlag = "Sequences";
	private static ConfigFileSection sequencesSection = new ConfigFileSection(sequencesSectionFlag, true);
	private static String primerSizeOptionFlag = "primer_size";
	private static ConfigFileOption primerSizeOption = new ConfigFileOption(primerSizeOptionFlag, 2, false, false, true);
	private static String primer3corePathOptionFlag = "primer3_core_path";
	private static ConfigFileOption primer3corePathOption = new ConfigFileOption(primer3corePathOptionFlag, 2, false, false, true);
	private static String optimalTmOptionFlag = "optimal_tm_for_primers";
	private static ConfigFileOption optimalTmOption = new ConfigFileOption(optimalTmOptionFlag, 2, false, false, true);
	private static String primerFileOptionFlag = "primer_file";
	private static ConfigFileOption primerFileOption = new ConfigFileOption(primerFileOptionFlag, 2, false, false, false);
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
		ProbeFilter[] filters = new ProbeFilter[] { new RepeatFilter(), new PolyBaseFilter() };
		
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
		PoolScheme[] schemes = new PoolScheme[] { new SimplePoolScheme(), new StackedSimplePoolScheme(), new GenePoolScheme() };
		
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
	
	private void createOligos(String outFilePrefix) throws IOException {
		logger.info("");
		logger.info("Creating oligos...");
		Collection<ProbeSet> probeSets = poolScheme.getProbes(transcripts);
		primerCounts = new TreeMap<PrimerPair,Integer>();
		int numProbes = 0;
		for(ProbeSet probeSet : probeSets) {
			numProbes += probeSet.getProbes().size();
		}
		logger.info("There are " + probeSets.size() + " probe sets with a total of " + numProbes + " probes.");
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
		// Assign and filter primers
		logger.info("Assigning primers to " + probeSets.size() + " probe sets...");
		int tried = 0;
		int succeeded = 0;
		for(ProbeSet probeSet : probeSets) {
			boolean foundPrimer = false;
			while(!foundPrimer) {
				boolean rejected = false;
				PrimerPair primer = PrimerUtils.getOneSyntheticPrimerPair(primerSize, primer3corePath, optimalTm, primerReader);
				tried++;
				for(PrimerFilter filter : primerFilters) {
					if(filter.rejectPrimer(primer, probeSet)) {
						rejected = true;
					}
				}
				if(!rejected) {
					// Primer is OK
					for(Probe probe : probeSet.getProbes()) {
						Oligo oligo = new Oligo(probe, probeSet, primer);
						oligos.add(oligo);
					}
					succeeded++;
					foundPrimer = true;
					primerCounts.put(primer, probeSet.getProbes().size());
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
		for(Oligo oligo : oligos) {
			Sequence seq = new Sequence(oligo.getProbe().getID());
			seq.setSequenceBases(oligo.getProbe().getProbeSequence());
			probeSeqs.add(seq);
		}
		FastaSequenceIO fsio = new FastaSequenceIO(outFile);
		fsio.write(probeSeqs, 200);
	}
	
	private void writeProbeList(String outFile) throws IOException {
		logger.info("Writing list of probe sequences to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		for(Oligo oligo : oligos) {
			w.write(oligo.getProbe().getProbeSequence() + "\n");
		}
		w.close();
	}
	
	private void writeOligoFasta(String outFile) throws IOException {
		logger.info("Writing fasta file of full oligo sequences to file " + outFile);
		List<Sequence> oligoSeqs = new ArrayList<Sequence>();
		for(Oligo oligo : oligos) {
			Sequence seq = new Sequence(oligo.getProbe().getID());
			seq.setSequenceBases(oligo.getOligoBases());
			oligoSeqs.add(seq);
		}
		FastaSequenceIO fsio = new FastaSequenceIO(outFile);
		fsio.write(oligoSeqs, 200);
	}
	
	private void writeOligoList(String outFile) throws IOException {
		logger.info("Writing list of full oligo sequences to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		for(Oligo oligo : oligos) {
			if (oligo.getOligoBases().toUpperCase().indexOf("GNIL") != -1) {
				w.close();
				throw new IllegalArgumentException("found GNIL");
			}
			w.write(oligo.getOligoBases() + "\n");
		}
		w.close();
	}
	
	private void writePrimers(String outFile) throws IOException {
		logger.info("Writing primers to file " + outFile);
		Map<PrimerPair, String> layoutsByPrimer = new TreeMap<PrimerPair, String>();
		for(Oligo oligo : oligos) {
			String layout = oligo.getProbeSet().getName();
			layoutsByPrimer.put(oligo.getPrimer(), layout);
		}
		String header = "Probe_layout\t";
		header += "Left_primer\t";
		header += "Right_primer\t";
		FileWriter w = new FileWriter(outFile);
		w.write(header + "\n");
		for(PrimerPair p : layoutsByPrimer.keySet()) {
			String line = layoutsByPrimer.get(p) + "\t";
			line += p.getLeftPrimer() + "\t";
			line += p.getRightPrimer() + "\t";
			w.write(line + "\n");
		}
		w.close();
	}
	
	private void writeFullTable(String outFile) throws IOException {
		logger.info("Writing full design table to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		String header = "Probe_ID\t";
		header += "Parent_sequence\t";
		header += "Probeset_size\t";
		header += "Start\t";
		header += "End\t";
		header += "Orientation\t";
		header += "Probe_layout\t";
		header += "Left_primer\t";
		header += "Right_primer\t";
		header += "Probe_sequence\t";
		header += "Full_oligo\t";
		w.write(header + "\n");
		for(Oligo oligo : oligos) {
			Probe probe = oligo.getProbe();
			PrimerPair primer = oligo.getPrimer();
			String line = probe.getID() + "\t";
			line += probe.getParentTranscript().getId() + "\t";
			line += primerCounts.get(primer) + "\t";
			line += probe.getStartPosOnTranscript() + "\t";
			int end = probe.getEndPosOnTranscript() - 1;
			line += end + "\t";
			line += probe.isAntisenseToTranscript() ? "antisense\t" : "sense\t";
			line += probe.getProbeLayout().toString() + "\t";
			line += primer.getLeftPrimer() + "\t";
			line += primer.getRightPrimer() + "\t";
			line += probe.getProbeSequence() + "\t";
			line += oligo.getOligoBases() + "\t";
			w.write(line + "\n");
		}
		w.close();
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
		p.parse(args);
		String configFile = p.getStringArg("-c");
		String outFilePrefix = p.getStringArg("-o");
		String levelString = p.getStringArg("-l");
		
		Level level = Level.toLevel(levelString, Level.INFO);
		logger.setLevel(level);
		
		// Instantiate OligoPool from config file
		OligoPool oligoPool = new OligoPool(configFile);
		
		// Create oligos
		oligoPool.createOligos(outFilePrefix);
		
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

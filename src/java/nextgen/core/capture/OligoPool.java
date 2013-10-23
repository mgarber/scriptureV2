package nextgen.core.capture;

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

import nextgen.core.capture.arrayscheme.PoolScheme;
import nextgen.core.capture.arrayscheme.ProbeLayout;
import nextgen.core.capture.arrayscheme.SimplePoolScheme;
import nextgen.core.capture.arrayscheme.SingleTilingProbeLayout;
import nextgen.core.capture.arrayscheme.StackedSimplePoolScheme;
import nextgen.core.capture.filter.PrimerFilter;
import nextgen.core.capture.filter.ProbeFilter;
import nextgen.core.capture.filter.RepeatFilter;
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
	private Collection<ProbeFilter> probeFilters;
	private Collection<PrimerFilter> primerFilters;
	private Collection<Oligo> oligos;
	private int primerSize;
	private String primer3corePath;
	
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
	/**
	 * Probe filter option flag
	 */
	public static String probeFilterOptionFlag = "probe_filter";
	private static ConfigFileOption probeFilterOption = new ConfigFileOption(probeFilterOptionFlag, 5, true, true, false);
	private static String probeFiltersSectionFlag = "Probe_filters";
	private static ConfigFileSection probeFiltersSection = new ConfigFileSection(probeFiltersSectionFlag, false);
	private static String primerFilterOptionFlag = "primer_filter";
	private static ConfigFileOption primerFilterOption = new ConfigFileOption(primerFilterOptionFlag, 5, true, true, false);
	private static String primerFiltersSectionFlag = "Primer_filters";
	private static ConfigFileSection primerFiltersSection = new ConfigFileSection(primerFiltersSectionFlag, false);
	
	private static ConfigFile createConfigFile(String fileName) throws IOException {
		logger.info("Getting config file from " + fileName);
		arraySchemeSection.addAllowableOption(probeLayoutOption);
		arraySchemeSection.addAllowableOption(poolSchemeOption);
		arraySchemeSection.addAllowableOption(primerSizeOption);
		arraySchemeSection.addAllowableOption(primer3corePathOption);
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
		// Repeat filter
		RepeatFilter repeatFilter = new RepeatFilter();
		if(repeatFilter.validConfigFileValue(value)) {
			repeatFilter.setParametersFromConfigFile(value);
			logger.info("Got filter " + repeatFilter.toString());
			return repeatFilter;
		}
		// Probe filter 2 (can pass pointer to this as argument in constructor)
		
		// ...
		throw new IllegalArgumentException("Probe filter not connected to this method: " + value.getFullOptionLine() + ". Need to implement.");
	}

	private static PrimerFilter getPrimerFilterFromConfigFileValue(@SuppressWarnings("unused") ConfigFileOptionValue value) {
		// Primer filter 1 (can pass pointer to this as argument in constructor)
		
		// Primer filter 2 (can pass pointer to this as argument in constructor)
		
		// ...
		throw new IllegalStateException("No primer filters are connected to this method. Need to implement.");
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
		// Simple pool scheme
		SimplePoolScheme simplePoolScheme = new SimplePoolScheme();
		if(simplePoolScheme.validConfigFileValue(poolSchemeVal)) {
			simplePoolScheme.setFromConfigFile(configFile);
			logger.info("Got pool scheme " + simplePoolScheme.toString());
			return simplePoolScheme;
		}
		// Stacked simple pool scheme
		StackedSimplePoolScheme stackedSimplePoolScheme = new StackedSimplePoolScheme();
		if(stackedSimplePoolScheme.validConfigFileValue(poolSchemeVal)) {
			stackedSimplePoolScheme.setFromConfigFile(configFile);
			logger.info("Got pool scheme " + stackedSimplePoolScheme.toString());
			return stackedSimplePoolScheme;
		}
		// Pool scheme 3
		
		// Pool scheme 4
		
		// ...
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
	
	private Collection<ProbeFilter> getProbeFiltersFromConfigFile() {
		Collection<ProbeFilter> rtrn = new ArrayList<ProbeFilter>();
		Collection<ConfigFileOptionValue> probeFilterVals = configFile.getOptionValues(probeFiltersSection, probeFilterOption);
		if(probeFilterVals == null) {
			return rtrn;
		}
		for(ConfigFileOptionValue val : probeFilterVals) {
			rtrn.add(getProbeFilterFromConfigFileValue(val));
		}
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
	
	private void createOligos() throws IOException {
		logger.info("");
		logger.info("Creating oligos...");
		Collection<ProbeSet> probeSets = poolScheme.getProbes(transcripts);
		int numProbes = 0;
		for(ProbeSet probeSet : probeSets) {
			numProbes += probeSet.getProbes().size();
		}
		logger.info("There are " + probeSets.size() + " probe sets with a total of " + numProbes + " probes.");
		// Filter probes
		logger.info("Filtering probes...");
		int removed = 0;
		int remaining = 0;
		for(ProbeSet probeSet : probeSets) {
			Iterator<Probe> iter = probeSet.iter();
			while(iter.hasNext()) {
				Probe probe = iter.next();
				boolean rejected = false;
				for(ProbeFilter filter : probeFilters) {
					if(filter.rejectProbe(probe)) {
						iter.remove();
						rejected = true;
						removed++;
						break;
					}
				}
				if(!rejected) remaining++;
			}
		}
		logger.info("Done filtering probes. Removed " + removed + " probes. " + remaining + " probes remain.");
		// Assign and filter primers
		logger.info("Assigning primers to " + probeSets.size() + " probe sets...");
		int tried = 0;
		int succeeded = 0;
		for(ProbeSet probeSet : probeSets) {
			boolean foundPrimer = false;
			while(!foundPrimer) {
				boolean rejected = false;
				PrimerPair primer = PrimerUtils.getOneSyntheticPrimerPair(primerSize, primer3corePath);
				tried++;
				for(PrimerFilter filter : primerFilters) {
					if(filter.rejectPrimer(primer, probeSet)) {
						rejected = true;
					}
				}
				if(!rejected) {
					// Primer is OK
					for(Probe probe : probeSet.getProbes()) {
						Oligo oligo = new Oligo(probe, primer);
						oligos.add(oligo);
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
			w.write(oligo.getOligoBases() + "\n");
		}
		w.close();
	}
	
	private void writePrimers(String outFile) throws IOException {
		logger.info("Writing primers to file " + outFile);
		Map<PrimerPair, String> layoutsByPrimer = new TreeMap<PrimerPair, String>();
		for(Oligo oligo : oligos) {
			String layout = oligo.getProbe().getProbeLayout().toString();
			layoutsByPrimer.put(oligo.getPrimer(), layout);
		}
		String header = "Probe_layout\t";
		header += "Left_primer\t";
		header += "Right_primer\t";
		FileWriter w = new FileWriter(outFile);
		w.write(header + "\n");
		for(PrimerPair p : layoutsByPrimer.keySet()) {
			String line = layoutsByPrimer.get(p).toString() + "\t";
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
		oligoPool.createOligos();
		
		// Write files
		oligoPool.writeFiles(outFilePrefix);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}

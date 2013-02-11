/**
 * Design tiling oligos and primers that do not cross-hybridize with each other or the rest of the transcriptome
 */
package broad.pda.capture.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.motif.SearchException;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.Primer3Configuration;
import broad.core.primer3.Primer3ConfigurationFactory;
import broad.core.primer3.Primer3IO;
import broad.core.primer3.Primer3SequenceInputTags;
import broad.core.primer3.PrimerPair;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;
import broad.core.util.PipelineUtils;

/**
 * For a given set of sequences, design tiling oligos and primers that do not cross-hybridize with each other or the rest of the transcriptome
 * Basic steps:
 * 1. Load input transcripts for pulldown and depletion
 * 2. Classify transcripts by gene and species
 * 3. Create qPCR primers against pulldown genes
 * 4. Create tiling probes
 * 5. Filter redundant probes
 * 6. Design PCR tails
 * 7. Filter cross priming tails
 * 8. Assign tails to probes
 * 
 * Afterward, you need to filter probes that cross hybridize to another transcript, and remove those probes from the final design.
 * 
 * WARNING: At most stages, the code tries to save time by reading back existing output files from a previous run. Clear previous output in the working directory if you do not want it to be reused.
 * IMPORTANT: The code tries to read previously created qPCR primers and previously created probes. The probe indicates whether it contains a qPCR primer site. If you are allowing the code to read existing files, make sure these two were generated on the same run.
 * 
 * IMPORTANT: This code was used to create two arrays. One array consisted of probes for depletion and pulldown, as well as some short probes for pulldown. The other array consisted of mRNAs for depletion. Some array properties are hard coded. To reuse this code you would need to check some hardcoded properties. Affected lines are marked with comments.
 * IMPORTANT: This code launches LSF processes off of the build in Pam's directory. The location of the build is hard coded. To reuse this code, you will need to hard code the location of the code base, or change it to a configurable parameter.
 * 
 * @author prussell
 */
public final class ArrayDesigner {
	
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Private fields
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	
	/**
	 * The probe size
	 */
	protected static final int OLIGO_SIZE = 120;

	//protected static final int OLIGO_SIZE = 60; // for short probes experiment
	
	//private static final int TILING_DENSITY = 120; // for mrna depletion array
	
	/**
	 * The match cutoff for probes to hybridize to another gene
	 */
	protected static final int MATCH_CUTOFF = 30;

	/**
	 * The number of location domains along pulldown probes
	 */
	private static final int NUM_DOMAINS = 4;

	//private static final int TILING_DENSITY = 120; // for mrna depletion array
	
	/**
	 * The maximum number of repetetive bases in a probe
	 */
	protected static final int REPEAT_CUTOFF = 50;

	/**
	 * The spacing of probes along each transcript
	 */
	private static final int TILING_DENSITY = 15; // for main array

	//protected static final int OLIGO_SIZE = 60; // for short probes experiment
	
	//private static final int TILING_DENSITY = 120; // for mrna depletion array
	
	/**
	 * The length of a perfect match between two primers across genes
	 */
	private int betweenGenePrimerDimerMatchLength;

	/**
	 * The inner left primers for each probe for depletion
	 */
	private Map<Species, HashMap< RnaClass , TreeMap< Probe, String > > > depletionInnerLeftPrimers;

	/**
	 * The inner left primers for each probe for depletion
	 */
	private Map<Species, HashMap< RnaClass , TreeMap< Probe, String > > > depletionInnerRightPrimers;

	/**
	 * The outer left primers for each probe for depletion
	 */
	private Map<Species, HashMap< RnaClass , TreeMap< Probe, String > > > depletionOuterLeftPrimers;

	/**
	 * The outer left primers for each probe for depletion
	 */
	private Map<Species, HashMap< RnaClass , TreeMap< Probe, String > > > depletionOuterRightPrimers;

	/**
	 * The primers for depletion before filtering for cross priming
	 */
	private Map< PrimerPair, ArrayList<PrimerPair> > depletionPrimers;

	/**
	 * Name of file containing depletion sequences to design probes against
	 */
	private String depletionseqfile;

	/**
	 * Whether the main array needs to include design for depletion sequences
	 */
	private boolean mainArrayHasDepletionSequences;
	
	/**
	 * Whether to make nested primers or just one gene primer for main array
	 */
	private boolean makeFullNestedPrimerDesign;
	
	/**
	 * Whether to make the second array for mRNA depletion
	 */
	private boolean makeMrnaDepletionArray;
	
	/**
	 * The length of perfect matches between both primers and a probe
	 */
	private int doublePrimerMatchToProbe;

	/**
	 * Name of file containing transcriptome
	 */
	private String genefile;

	/**
	 * Probes designed against each gene collapsed from isoforms
	 */
	private Map< SpeciesGene, TreeSet<Probe> > geneProbeSet;

	/**
	 * The transcriptome to avoid when designing primers
	 */
	private Collection<Sequence> genes;

	/**
	 * Full set of probes for all sequences
	 * Each sequence is associated with its set of probes
	 * Individual probes know their attributes such as start and end position on the sequence, tiling path, etc.
	 */
	private Map< Sequence, TreeSet<Probe> > isoformProbeSet;

	/**
	 * The set of amplification primers that have been assigned
	 */
	private HashSet<String> keptPrimerKmers;

	/**
	 * The minimum number of times a kmer appears in probes before it needs to be filtered from primers
	 */
	private int minKmerOccurrencesInProbes;

	/**
	 * Name of file containing mRNA sequences to put on depletion array
	 */
	private String mrnadepletionfile;

	/**
	 * The number of tiling paths
	 */
	private int numTilingPaths;

	/**
	 * Name of output file for the full design
	 */
	private String outdesignfile;

	/**
	 * Name of output file for primers
	 */
	private String outprimersfile;

	/**
	 * Name of output file for probes
	 */
	private String outprobesfile;

	/**
	 * The inner left primers for pulldown
	 */
	private Map< Probe, String > pulldownInnerLeftPrimers;

	/**
	 * The inner right primers for pulldown
	 */
	private Map< Probe, String > pulldownInnerRightPrimers;

	/**
	 * The middle left primers for pulldown
	 */
	private Map< Probe, String > pulldownMiddleLeftPrimers;

	/**
	 * The middle right primers for pulldown
	 */
	private Map< Probe, String > pulldownMiddleRightPrimers;

	/**
	 * The outer left primers for pulldown
	 */
	private Map< Probe, String > pulldownOuterLeftPrimers;

	/**
	 * The outer right primers for pulldown
	 */
	private Map< Probe, String > pulldownOuterRightPrimers;

	/**
	 * The primers for pulldown before filtering for cross priming
	 */
	private Map< PrimerPair, Map< PrimerPair, ArrayList<PrimerPair> > > pulldownPrimers;

	/**
	 * Name of file containing pulldown sequences to design probes against
	 */
	private String pulldownseqfile;

	/**
	 * Name of file containing qPCR primers
	 */
	private String qpcrfile;

	/**
	 * The set of sequences to design probes against
	 */
	private Collection<Sequence> sequences;

	/**
	 * Map keeping track of the purpose of each sequence on the array
	 */
	private Map< Sequence, ProbePurpose > sequencesByPurpose;

	/**
	 * Map associating gene names with the set of sequences representing their isoforms
	 */
	private Map< SpeciesGene, Collection<Sequence> > sequencesBySpeciesGene;

	/**
	 * The length of a perfect match between one primer and a probe
	 */
	private int singlePrimerMatchToProbe;

	/**
	 * 
	 */
	private int singlePrimerMatchToProbeWithinGene;

	/**
	 * The length of a perfect match between two primers in the same gene
	 */
	private int withinGenePrimerDimerMatchLength;

	/**
	 * ArrayDesigner should only be instantiated within the few classes that it delegates work to
	 * @throws IOException 
	 */
	protected ArrayDesigner() throws IOException {
		this.genefile = null;
		this.genes = new HashSet<Sequence>();
		this.outdesignfile = null;
		this.outprimersfile = null;
		this.outprobesfile = null;
		this.pulldownseqfile = null;
		this.mrnadepletionfile = null;
		this.depletionseqfile = null;
		this.qpcrfile = null;
		this.sequences = new HashSet<Sequence>();
		this.sequencesBySpeciesGene = new HashMap<SpeciesGene,Collection<Sequence>>();
		this.sequencesByPurpose = new HashMap<Sequence, ProbePurpose>();
		this.isoformProbeSet = new HashMap<Sequence, TreeSet<Probe>>();
		this.geneProbeSet = new HashMap<SpeciesGene, TreeSet<Probe>>();
		BigInteger os = BigInteger.valueOf(OLIGO_SIZE);
		BigInteger td = BigInteger.valueOf(TILING_DENSITY);
		BigInteger tp = os.divide(os.gcd(td));
		this.numTilingPaths = tp.intValue();
		this.depletionPrimers = new HashMap<PrimerPair, ArrayList<PrimerPair>>();
		this.pulldownPrimers = new HashMap< PrimerPair, Map< PrimerPair, ArrayList<PrimerPair> > >();
		this.depletionOuterLeftPrimers = new HashMap< Species, HashMap< RnaClass, TreeMap< Probe, String > > >();
		this.depletionInnerLeftPrimers = new HashMap< Species, HashMap< RnaClass, TreeMap< Probe, String > > >();
		this.pulldownOuterLeftPrimers = new TreeMap< Probe, String >();
		this.pulldownMiddleLeftPrimers = new TreeMap< Probe, String >();
		this.pulldownInnerLeftPrimers = new TreeMap< Probe, String >();
		this.depletionOuterRightPrimers = new HashMap< Species, HashMap< RnaClass, TreeMap< Probe, String > > >();
		this.depletionInnerRightPrimers = new HashMap< Species, HashMap< RnaClass, TreeMap< Probe, String > > >();
		this.pulldownOuterRightPrimers = new TreeMap< Probe, String >();
		this.pulldownMiddleRightPrimers = new TreeMap< Probe, String >();
		this.pulldownInnerRightPrimers = new TreeMap< Probe, String >();
		this.keptPrimerKmers = new HashSet<String>();
		
	}
	
	/**
	 * Get an array of all substrings of a string for a specified length
	 * @param s the string
	 * @param length the length of substrings
	 * @return set of all substrings
	 */
	private static String[] getAllSubstringsArray(String s, int length) {
		String[] rtrn = new String[s.length()-length+1];
		for(int i=0; i <= s.length()-length; i++) {
			rtrn[i] = s.substring(i,i + length);
		}
		return rtrn;
	}

	/**
	 * Get a hashset of all substrings of a string for a specified length
	 * @param s the string
	 * @param length the length of substrings
	 * @return set of all substrings
	 */
	private static HashSet<String> getAllSubstringsHashSet(String s, int length) {
		HashSet<String> rtrn = new HashSet<String>();
		for(int i=0; i <= s.length()-length; i++) {
			rtrn.add(s.substring(i,i + length));
		}
		return rtrn;
	}

	/**
	 * Parse command line and set values of arguments
	 * @param args the command line arguments
	 */
	private void loadCommandArgs(String[] args) {
		// Add command arguments
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("Choose probes and design primers.");
		p.addStringArg("-p", "Fasta file of all pulldown sequences", true);
		p.addStringArg("-d", "Fasta file of all depletion sequences for main array", false);
		p.addBooleanArg("-f", "Full nested primer design", true);
		p.addStringArg("-m", "Fasta file of all mRNA sequences for depletion array", false);
		p.addStringArg("-g", "Fasta file of all genes", true);
		p.addStringArg("-q", "Fasta file of qPCR primer sequences", true);
		p.addStringArg("-design", "Output design file without extension, for table and fasta files", true);
		p.addStringArg("-probes", "Output probes file without extension, for table and fasta files", false);
		p.addStringArg("-primers", "Output primers file", false);
		p.addIntegerArg("-spp", "match length for single primer match to any probe on array", true);
		p.addIntegerArg("-bgpd", "match length for between genes primer dimer", true);
		p.addIntegerArg("-wgpd", "match length for within gene primer dimer", true);
		p.addIntegerArg("-kp","max allowable kmer occurrences in probe for single primer match",true);
		p.addIntegerArg("-dpp","match length for double primer match to probe",true);
		p.addIntegerArg("-sppwg", "match length for single primer match to probe within gene",true);
		p.parse(args);
	
		// Get values of command arguments
		// All pulldown sequences to design probes against
		this.pulldownseqfile = p.getStringArg("-p");
		// All depletion sequences for main array
		this.depletionseqfile = p.getStringArg("-d");
		if(this.depletionseqfile == null) this.mainArrayHasDepletionSequences = false;
		else this.mainArrayHasDepletionSequences = true;
		// All mRNA sequences for depletion array
		this.mrnadepletionfile = p.getStringArg("-m");
		if(this.mrnadepletionfile == null) this.makeMrnaDepletionArray = false;
		else this.makeMrnaDepletionArray = true;
		// Whether to make full nested primer design
		this.makeFullNestedPrimerDesign = p.getBooleanArg("-f").booleanValue();
		// Transcripts to avoid when designing primers
		this.genefile = p.getStringArg("-g");
		// Output design file
		this.outdesignfile = p.getStringArg("-design");
		// Output probe file
		this.outprobesfile = p.getStringArg("-probes");
		// Output primers file
		this.outprimersfile = p.getStringArg("-primers");
		// qPCR primers
		this.qpcrfile = p.getStringArg("-q");
		// single primer match to any probe on array
		this.singlePrimerMatchToProbe = p.getIntegerArg("-spp").intValue();
		// between gene primer dimer
		this.betweenGenePrimerDimerMatchLength = p.getIntegerArg("-bgpd").intValue();
		// within gene primer dimer
		this.withinGenePrimerDimerMatchLength = p.getIntegerArg("-wgpd").intValue();
		// kmer occurrences in probe
		this.minKmerOccurrencesInProbes = p.getIntegerArg("-kp").intValue();
		// double primer match to probe
		this.doublePrimerMatchToProbe = p.getIntegerArg("-dpp").intValue();
		// single primer match to probe within gene
		this.singlePrimerMatchToProbeWithinGene = p.getIntegerArg("-sppwg").intValue();
	}
	
	/**
	 * Identify species and gene names in RefSeq fasta header and collect all isoforms of each gene
	 */
	private void loadSequencesBySpeciesGene() {
		if(this.sequences.isEmpty()) {
			throw new IllegalStateException("Set of sequences is empty.");
		}
		this.sequencesBySpeciesGene.clear();
		
		// Fill in map
		for(Sequence seq : this.sequences) {
			Collection<Sequence> tempseqs = new ArrayList<Sequence>();
			SpeciesGene speciesgene = new SpeciesGene(seq.getRefseqSpeciesName(),seq.getRefseqGeneName());
			// Set the purpose of this gene to the purpose of the first isoform encountered
			speciesgene.setPurpose(this.sequencesByPurpose.get(seq));
			speciesgene.setRnaClass(RnaClass.getRnaClass(seq.getRefseqRnaClassName()));
			if(this.sequencesBySpeciesGene.containsKey(speciesgene)) {
				// Temporarily copy the set, add the new sequence and put back
				tempseqs = this.sequencesBySpeciesGene.get(speciesgene);
				tempseqs.add(seq);
				this.sequencesBySpeciesGene.put(speciesgene, tempseqs);
			}
			else {
				tempseqs.add(seq);
				this.sequencesBySpeciesGene.put(speciesgene, tempseqs);
			}
			tempseqs = null;
		}
	
		// Don't need sequences by purpose any more
		this.sequencesByPurpose = null;
		
	}

	/**
	 * Assign primers to probes, filtering for cross primes simultaneously
	 */
	private void assignAndFilterPrimers() throws FileNotFoundException, IOException {
		
		if(this.mainArrayHasDepletionSequences && this.depletionPrimers.isEmpty()) {
			throw new IllegalStateException("Can't assign and filter primers. Depletion primer set is empty.");
		}
	
		if(this.pulldownPrimers.isEmpty()) {
			throw new IllegalStateException("Can't assign and filter primers. Pulldown primer set is empty.");
		}
	
		// Clear set of kept primers
		this.keptPrimerKmers.clear();
		
		// Count how many genes are for depletion
		int numDepletionGenes = 0;
		if(this.mainArrayHasDepletionSequences) {
			for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
				if(speciesgene.getPurpose() == ProbePurpose.DEPLETE_MAIN_ARRAY || speciesgene.getPurpose() == ProbePurpose.DEPLETE_ALL_MRNA) {
					numDepletionGenes ++;
	            }
			}
			System.out.println("There are " + numDepletionGenes + " genes for depletion.");
		}
	    
	    
		// Count how many genes are for pulldown
	    int numPulldownGenes = 0;
	    for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
	            if(speciesgene.getPurpose() == ProbePurpose.PULLDOWN) {
	                    numPulldownGenes ++;
	            }
	    }
	    System.out.println("There are " + numPulldownGenes + " genes for pulldown.");
	    
	    
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    // Make hash tables for single matches and double matches
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	            
	    System.out.println("\nMaking kmer tables of probe sequences...");
	    HashMap<String,Integer> singleKmers = new HashMap<String,Integer>();
	    HashSet<StringPair> doubleKmers = new HashSet<StringPair>();
	    
	    boolean readFromSingleFile = false;
	    boolean readFromDoubleFile = false;
	    
	    File singleKmerFile = new File("SingleKmersForPerfectMatch_" + this.singlePrimerMatchToProbe + ".out");
	    File doubleKmerFile = new File("DoubleKmersForPerfectMatch_" + this.doublePrimerMatchToProbe + ".out");
	    
	    // Read kmers from file if possible
	    if(singleKmerFile.exists()) {
	    	readFromSingleFile = true;
	    	System.out.println("WARNING: reading single kmers for perfect match from file " + singleKmerFile);
	    	System.out.println("IMPORTANT: code is assuming these files were generated off of the probe set that is being used.");
	    	
	    	StringParser stringparse = new StringParser();
	    	FileReader r = new FileReader(singleKmerFile);
	    	BufferedReader b = new BufferedReader(r);
	    	while(b.ready()) {
	    		String line = b.readLine();
	    		stringparse.parse(line);
	    		String kmer = stringparse.asString(0);
	    		int num = stringparse.asInt(1);
	    		singleKmers.put(kmer,Integer.valueOf(num));
	    	}
	     	System.out.println("Done reading single kmer table from file.");
	    }
	
	    // Read kmers from file if possible
	    if(doubleKmerFile.exists()) {
	    	int numRead = 0;
	    	readFromDoubleFile = true;
	    	System.out.println("WARNING: reading double kmers for perfect match from file " + doubleKmerFile);
	    	System.out.println("IMPORTANT: code is assuming these files were generated off of the probe set that is being used.");
	    	
	    	StringParser stringparse = new StringParser();
	    	FileReader r2 = new FileReader(doubleKmerFile);
	    	BufferedReader b2 = new BufferedReader(r2);
	    	while(b2.ready()) {
	    		numRead++;
	    		if(numRead % 10000000 == 0) System.out.println("Read " + Integer.valueOf(numRead/1000000).toString() + " million kmer pairs.");
	    		String line = b2.readLine();
	    		stringparse.parse(line);
	    		StringPair pair = this.new StringPair(stringparse.asString(0), stringparse.asString(1));
	    		doubleKmers.add(pair);
	    		line = null;
	    		pair = null;
	    	}
	    	
	    	System.out.println("Done reading double kmer table from file.");
	    }
	
	    // Generate single kmers
	   	if(!readFromSingleFile) {
	   		System.out.println("Generating single kmers...");
	   		int numProbes = 0;
	       	int numDone = 0;
	       	boolean ten = false;
	       	boolean twenty = false;
	       	boolean thirty = false;
	       	boolean forty = false;
	       	boolean fifty = false;
	       	boolean sixty = false;
	       	boolean seventy = false;
	       	boolean eighty = false;
	       	boolean ninety = false;
	
	       	for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) numProbes += this.geneProbeSet.get(speciesgene).size();
	
	   	    for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
	
	   	    	for(Probe probe : this.geneProbeSet.get(speciesgene)) {
	   	    		numDone++;
	   	    		if(numDone % 1000 == 0) {
	   	    			if((float)numDone/(float)numProbes > 0.1 && !ten) {ten = true; System.out.println("Finished 10% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.2 && !twenty) {twenty = true; System.out.println("Finished 20% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.3 && !thirty) {thirty = true; System.out.println("Finished 30% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.4 && !forty) {forty = true; System.out.println("Finished 40% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.5 && !fifty) {fifty = true; System.out.println("Finished 50% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.6 && !sixty) {sixty = true; System.out.println("Finished 60% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.7 && !seventy) {seventy = true; System.out.println("Finished 70% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.8 && !eighty) {eighty = true; System.out.println("Finished 80% of probes.");}
	   	    			if((float)numDone/(float)numProbes > 0.9 && !ninety) {ninety = true; System.out.println("Finished 90% of probes.");}
	   	    		}
	
	   	    		HashSet<String> substrings = getAllSubstringsHashSet(probe.getSequenceBases(),this.singlePrimerMatchToProbe);
	   	    		for(String kmer : substrings) {
	   	    			if(singleKmers.containsKey(kmer)) {
	   	    				Integer newCount = Integer.valueOf(singleKmers.get(kmer).intValue() + 1);
	   	    				singleKmers.put(kmer, newCount);
	   	    				newCount = null;
	   	    			} else singleKmers.put(kmer,Integer.valueOf(1));
	   	    		}
	   	    		substrings = null;
	   	    	}
	    	}
			
	   	    FileWriter singleWriter = new FileWriter(singleKmerFile);
	    	for(String s : singleKmers.keySet()) singleWriter.write(s + " " + singleKmers.get(s) + "\n");
	    	singleWriter.close();
	    	System.out.println("Wrote single kmers to file " + singleKmerFile);
	
	   	    
	   	}
	    			
	   	// Generate double kmers
	    if(!readFromDoubleFile) {
	    	System.out.println("Generating double kmers...");
	
	    	// process probes in batches of 10% and write to temp files
	    	for(int k = 0; k<10; k++) {
	    		HashSet<StringPair> tempDoubleKmers = new HashSet<StringPair>();
	    		String tempDoubleKmerFile = "temp_double_kmers_" + Integer.valueOf(k).toString();
	    		for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
	    			if(speciesgene.hashCode() % 10 == k) {
	    				for(Probe probe : this.geneProbeSet.get(speciesgene)) {
	    					String[] substringsForDoubleMatch = getAllSubstringsArray(probe.getSequenceBases(),this.doublePrimerMatchToProbe);
	    					int len = substringsForDoubleMatch.length;
	    					for(int i=0; i<len; i++) {
	    						for(int j=i+1; j<len; j++) {
	    							String stringI = substringsForDoubleMatch[i];
	    							String stringJ = substringsForDoubleMatch[j];
	    							StringPair pair = new StringPair(stringI,stringJ);
	    							tempDoubleKmers.add(pair);
	    							stringI = null;
	    							stringJ = null;
	    							pair = null;
	    						}
	    					}
	    					substringsForDoubleMatch = null;
	    				}
	    			}
	    		}
				FileWriter doubleWriter = new FileWriter(tempDoubleKmerFile);
				for(StringPair p : tempDoubleKmers) doubleWriter.write(p.first() + " " + p.second() + "\n");
				doubleWriter.close();
				System.out.println("Wrote double kmer set " + Integer.valueOf(k+1).toString() + "/10 to file " + tempDoubleKmerFile);        		
				tempDoubleKmers.clear();
				tempDoubleKmers = null;
	    	}
	    	
	   	
	    	System.out.println("Reading double kmers back from temp files...");
	    	
	    	// Read back from temp files
	    	for(int k=0; k<10; k++) {
	    		StringParser stringparse = new StringParser();
	    		FileReader r = new FileReader("temp_double_kmers_" + Integer.valueOf(k).toString());
	    		BufferedReader b = new BufferedReader(r);
	    		while(b.ready()) {
	    			String line = b.readLine();
	    			stringparse.parse(line);
	    			StringPair pair = this.new StringPair(stringparse.asString(0), stringparse.asString(1));
	    			doubleKmers.add(pair);
	    			line = null;
	    			pair = null;
	    		}
	    		System.out.println("Finished reading file temp_double_kmers_" + Integer.valueOf(k).toString());
	    	}
	    	
	    	System.out.println("Writing all double kmers to file " + doubleKmerFile);
			FileWriter doubleWriter = new FileWriter(doubleKmerFile);
	    	for(StringPair p : doubleKmers) doubleWriter.write(p.first() + " " + p.second() + "\n");
	    	doubleWriter.close();
	    	System.out.println("Wrote double kmers to file " + doubleKmerFile);
	    }
	    
	    
	    System.out.println("Done making kmer tables.");
	    System.out.println("There are " + singleKmers.size() + " unique " + this.singlePrimerMatchToProbe + "-mers in the probe set.");
	    System.out.println("There are " + doubleKmers.size() + " unique ordered pairs of " + this.doublePrimerMatchToProbe + "-mers appearing together in probes.");
	    
	    
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    // Assign depletion outer primers to species and class
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    
    	int numDone = 0;
	
	    if(this.mainArrayHasDepletionSequences) {
	    	int numDepletionPrimersSkippedInsufficientSubprimers = 0;
	    	int numDepletionPrimersSkippedMatchesPrimerAcrossGenes = 0;
	    	int numDepletionPrimersSkippedMatchesPrimerWithinGene = 0;
	    	int numDepletionPrimersSkippedSingleMatch = 0;
	    	int numDepletionPrimersSkippedDoubleMatch = 0;
	    	int numDepletionPrimersSkippedMatchesProbeWithinGene = 0;
	    
	    	System.out.println("\nAssigning depletion outer primers:");
	    
	    	// First match primers to species and class
	    	System.out.println("Matching outer primers to species and RNA class...");
	    	HashMap<Species, HashMap<RnaClass, PrimerPair>> depletionOuterPrimersToUse = new HashMap<Species, HashMap<RnaClass, PrimerPair>>();
	    
	    	// Iterator over depletion outer primers
	    	Iterator<PrimerPair> depletionOuterIter = this.depletionPrimers.keySet().iterator();
	    	for(Species species : Species.values()) {
	    		HashMap<RnaClass, PrimerPair> depletionOuterPrimersToUseForThisSpecies = new HashMap<RnaClass, PrimerPair>();
	    		for(RnaClass rnaclass : RnaClass.values()) {
	    			boolean assigned = false;
	    			while(!assigned) {
	    				try {
	    					PrimerPair thisDepletionOuter = depletionOuterIter.next();
	    					int numInnersForThisOuter = this.depletionPrimers.get(thisDepletionOuter).size();
	    					String thisLeftPrimer = thisDepletionOuter.getLeftPrimer();
	    					String thisRightPrimer = thisDepletionOuter.getRightPrimer();
	    					if(numInnersForThisOuter < numDepletionGenes * 1.5) {
	    						//System.out.println("Skipping depletion outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because only has " + numInnersForThisOuter + " inner primers.");
	    						numDepletionPrimersSkippedInsufficientSubprimers++;
	    						continue;
	    					}
	    					// check if the primer matches a primer that has already been kept
	    					if(primerCrossPrimesExistingPrimer(thisLeftPrimer) || primerCrossPrimesExistingPrimer(thisRightPrimer)) {
	    						//System.out.println("Skipping depletion outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because it matches an existing primer.");
	    						numDepletionPrimersSkippedMatchesPrimerAcrossGenes++;
	    						continue;
	    					}
	    					// check if the primer matches any probe
	    					// first check for a perfect longer match
	    					if(hasSingleMatch(thisLeftPrimer, thisRightPrimer, singleKmers)) {
	    						//System.out.println("Skipping depletion outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because it has a perfect " + MIN_SINGLE_PRIMER_MATCH + "-mer match to more than " + MIN_KMER_OCCURENCES + " probes.");
	    						numDepletionPrimersSkippedSingleMatch++;
	    						continue;
	    					}
	    					// now check for both primers to have a perfect shorter match
	    					if(hasDoubleMatch(this, thisLeftPrimer, thisRightPrimer, doubleKmers)) {
	    						//System.out.println("Skipping depletion outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because there is a probe with perfect " + MIN_DOUBLE_PRIMER_MATCH + "-mer matches to both primers.");
	    						numDepletionPrimersSkippedDoubleMatch++;
	    						continue;
	    					}
					
	    					// if still alive at this point, assign this primer to the species and class
	    					// add the primers to the list of kept primers
	    					depletionOuterPrimersToUseForThisSpecies.put(rnaclass, thisDepletionOuter);
	    					//System.out.println("Using primer pair " + thisLeftPrimer + " " + thisRightPrimer + " for species " + species + " and class " + rnaclass);
	    					this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisLeftPrimer,this.betweenGenePrimerDimerMatchLength));
	    					this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisRightPrimer,this.betweenGenePrimerDimerMatchLength));
	    					assigned = true;
							} catch (NoSuchElementException e) {
								System.err.println("Skipped " + numDepletionPrimersSkippedInsufficientSubprimers + " depletion primers due to insufficient inner primers.");
								System.err.println("Skipped " + numDepletionPrimersSkippedMatchesPrimerAcrossGenes + " depletion primers due to matches with existing primers.");
								System.err.println("Skipped " + numDepletionPrimersSkippedSingleMatch + " depletion primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
								System.err.println("Skipped " + numDepletionPrimersSkippedDoubleMatch + " depletion primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");
								throw new NoSuchElementException("Ran out of depletion outer primers during species " + species + " and class " + rnaclass);
							}
	    			}
				
	    		}
			
	    		depletionOuterPrimersToUse.put(species, depletionOuterPrimersToUseForThisSpecies);
	    		
	    	}
		
	    	System.out.println("Done assigning depletion outer primers to species and class.\n");
	    	System.out.println("So far, skipped " + numDepletionPrimersSkippedInsufficientSubprimers + " depletion primers due to insufficient inner primers.");
	    	System.out.println("So far, skipped " + numDepletionPrimersSkippedMatchesPrimerAcrossGenes + " depletion primers due to matches with existing primers.");
	    	System.out.println("So far, skipped " + numDepletionPrimersSkippedSingleMatch + " depletion primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	    	System.out.println("So far, skipped " + numDepletionPrimersSkippedDoubleMatch + " depletion primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");
	
	    	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    	// Assign depletion inner primers to genes
	    	// Assign outer and inner primers to probes
	    	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	    	System.out.println("\nAssigning inner primers to depletion genes and assigning outer and inner primers to probes.");
		
	    	// Make iterators for each set of inner primers
	    	HashMap<Species, HashMap<RnaClass, Iterator<PrimerPair>>> depletionInnerIterators = new HashMap<Species, HashMap<RnaClass, Iterator<PrimerPair>>>();
	    	for(Species species : Species.values()) {
	    		HashMap<RnaClass, Iterator<PrimerPair>> thisSpecies = new HashMap<RnaClass, Iterator<PrimerPair>>();
	    		for(RnaClass rnaclass : RnaClass.values()) {
	    			Iterator<PrimerPair> thisClass = this.depletionPrimers.get(depletionOuterPrimersToUse.get(species).get(rnaclass)).iterator();
	    			thisSpecies.put(rnaclass, thisClass);
	    		}
	    		depletionInnerIterators.put(species, thisSpecies);
	    	}
		
	    	// Assign and filter depletion inner primers
	    	for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
			
	    		// Skip pulldown genes
	    		if(speciesgene.getPurpose() == ProbePurpose.PULLDOWN) continue;
	    	
	    		numDone++;
	    		if(numDone % 100 == 0) {
	    			System.out.println("\nAssigned inner primers to " + numDone + " depletion genes.");
	    			System.out.println("So far, skipped " + numDepletionPrimersSkippedInsufficientSubprimers + " depletion primers due to insufficient inner primers.");
	    			System.out.println("So far, skipped " + numDepletionPrimersSkippedMatchesPrimerWithinGene + " depletion primers due to matches with existing primers within same gene.");
	    			System.out.println("So far, skipped " + numDepletionPrimersSkippedMatchesPrimerAcrossGenes + " depletion primers due to matches with existing primers from other genes.");
	    			System.out.println("So far, skipped " + numDepletionPrimersSkippedSingleMatch + " depletion primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	    			System.out.println("So far, skipped " + numDepletionPrimersSkippedDoubleMatch + " depletion primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
	    			System.out.println("So far, skipped " + numDepletionPrimersSkippedMatchesProbeWithinGene + " depletion primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
	    		}
	    	      	
	    		Species species = speciesgene.getSpecies();
	    		RnaClass rnaclass = speciesgene.getRnaClass();
	    		PrimerPair outerPrimer = depletionOuterPrimersToUse.get(species).get(rnaclass);
	    		String leftOuterPrimer = outerPrimer.getLeftPrimer();
	    		String rightOuterPrimer = outerPrimer.getRightPrimer();
	
	    		ArrayList<String> leftPrimersKeptForThisGene = new ArrayList<String>();
	    		ArrayList<String> rightPrimersKeptForThisGene = new ArrayList<String>();
	    		leftPrimersKeptForThisGene.add(leftOuterPrimer);
	    		rightPrimersKeptForThisGene.add(rightOuterPrimer);
	    	
	    		boolean assigned = false;
	    		while(!assigned) {
	    			try {
	    				PrimerPair thisDepletionInner = depletionInnerIterators.get(species).get(rnaclass).next();
	    		
	    				String thisLeftPrimer = thisDepletionInner.getLeftPrimer();
	    				String thisRightPrimer = thisDepletionInner.getRightPrimer();
	 			
	    				// check if the primer matches a primer from the same gene
	    				if(primerCrossPrimesInSet(thisLeftPrimer,rightPrimersKeptForThisGene,this.withinGenePrimerDimerMatchLength)) {
	    					numDepletionPrimersSkippedMatchesPrimerWithinGene++;
	    					continue;
	    				}
	    				if(primerCrossPrimesInSet(thisRightPrimer,leftPrimersKeptForThisGene,this.withinGenePrimerDimerMatchLength)) {
	    					numDepletionPrimersSkippedMatchesPrimerWithinGene++;
	    					continue;
	    				}        			
	    				// check if the primer matches a primer that has already been kept
	    				if(primerCrossPrimesExistingPrimer(thisLeftPrimer) || primerCrossPrimesExistingPrimer(thisRightPrimer)) {
	    					//System.out.println("Skipping depletion inner primer " + thisLeftPrimer + " " + thisRightPrimer + " because it matches an existing primer.");
	    					numDepletionPrimersSkippedMatchesPrimerAcrossGenes++;
	    					continue;
	    				}
	    				// check if the primer matches any probe
	    				// first check for a perfect longer match
	    				if(hasSingleMatch(thisLeftPrimer, thisRightPrimer, singleKmers)) {
	    					//System.out.println("Skipping depletion inner primer " + thisLeftPrimer + " " + thisRightPrimer + " because it has a perfect " + MIN_SINGLE_PRIMER_MATCH + "-mer match to more than " + MIN_KMER_OCCURENCES + " probes.");
	    					numDepletionPrimersSkippedSingleMatch++;
	    					continue;
	    				}
	    				// now check for both primers to have a perfect shorter match
	    				if(hasDoubleMatch(this, thisLeftPrimer, thisRightPrimer, doubleKmers)) {
	    					//System.out.println("Skipping depletion inner primer " + thisLeftPrimer + " " + thisRightPrimer + " because there is a probe with perfect " + MIN_DOUBLE_PRIMER_MATCH + "-mer matches to both primers.");
	    					numDepletionPrimersSkippedDoubleMatch++;
	    					continue;
	    				}
	    				// check if the primer matches a probe in the gene
	    				if(leftPrimerHasMatchInProbeSet(thisLeftPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	    					numDepletionPrimersSkippedMatchesProbeWithinGene++;
	    					continue;
	    				}
	    				if(rightPrimerHasMatchInProbeSet(thisRightPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	    					numDepletionPrimersSkippedMatchesProbeWithinGene++;
	    					continue;
	    				}
	    			
	    				// if still alive at this point, assign this primer to the gene
	    				// add the primers to the list of kept primers
	    				//System.out.println("Using primer pair " + thisLeftPrimer + " " + thisRightPrimer + " for gene " + species + " " + speciesgene.getGeneName());
	    		
	    				TreeMap<Probe,String> newLeftInners = new TreeMap<Probe,String>();
	    				try {newLeftInners.putAll(this.depletionInnerLeftPrimers.get(species).get(rnaclass));} catch (NullPointerException e) {}
	    		
	    				TreeMap<Probe,String> newLeftOuters = new TreeMap<Probe,String>();
	    				try {newLeftOuters.putAll(this.depletionOuterLeftPrimers.get(species).get(rnaclass));} catch (NullPointerException e) {}
	    		
	    				TreeMap<Probe,String> newRightInners = new TreeMap<Probe,String>();
	    				try {newRightInners.putAll(this.depletionInnerRightPrimers.get(species).get(rnaclass));} catch (NullPointerException e) {}
	    		
	    				TreeMap<Probe,String> newRightOuters = new TreeMap<Probe,String>();
	    				try {newRightOuters.putAll(this.depletionOuterRightPrimers.get(species).get(rnaclass));} catch (NullPointerException e) {}
	
	    				for(Probe probe : this.geneProbeSet.get(speciesgene)) {
	    					newLeftOuters.put(probe, leftOuterPrimer);
	    					newRightOuters.put(probe, rightOuterPrimer);
	    					newLeftInners.put(probe, thisLeftPrimer);
	    					newRightInners.put(probe, thisRightPrimer);
	    				}
	    		
	    				HashMap<RnaClass,TreeMap<Probe,String>> newLeftInnersForThisSpecies = new HashMap<RnaClass,TreeMap<Probe,String>>();
	    				try {newLeftInnersForThisSpecies.putAll(this.depletionInnerLeftPrimers.get(species));} catch (NullPointerException e) {}
	    				newLeftInnersForThisSpecies.put(rnaclass, newLeftInners);
	    				this.depletionInnerLeftPrimers.put(species, newLeftInnersForThisSpecies);
	    		
	    				HashMap<RnaClass,TreeMap<Probe,String>> newRightInnersForThisSpecies = new HashMap<RnaClass,TreeMap<Probe,String>>();
	    				try {newRightInnersForThisSpecies.putAll(this.depletionInnerRightPrimers.get(species));} catch (NullPointerException e) {}
	    				newRightInnersForThisSpecies.put(rnaclass, newRightInners);
	    				this.depletionInnerRightPrimers.put(species, newRightInnersForThisSpecies);
	    		
	    				HashMap<RnaClass,TreeMap<Probe,String>> newLeftOutersForThisSpecies = new HashMap<RnaClass,TreeMap<Probe,String>>();
	    				try {newLeftOutersForThisSpecies.putAll(this.depletionOuterLeftPrimers.get(species));} catch (NullPointerException e) {}
	    				newLeftOutersForThisSpecies.put(rnaclass, newLeftOuters);
	    				this.depletionOuterLeftPrimers.put(species, newLeftOutersForThisSpecies);
	    		
	    				HashMap<RnaClass,TreeMap<Probe,String>> newRightOutersForThisSpecies = new HashMap<RnaClass,TreeMap<Probe,String>>();
	    				try {newRightOutersForThisSpecies.putAll(this.depletionOuterRightPrimers.get(species));} catch (NullPointerException e) {}
	    				newRightOutersForThisSpecies.put(rnaclass, newRightOuters);
	    				this.depletionOuterRightPrimers.put(species, newRightOutersForThisSpecies);
	    		
	    				this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisLeftPrimer,this.betweenGenePrimerDimerMatchLength));
	    				this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisRightPrimer,this.betweenGenePrimerDimerMatchLength));
	    				leftPrimersKeptForThisGene.add(thisLeftPrimer);
	    				rightPrimersKeptForThisGene.add(thisRightPrimer);
	    				assigned = true;
	    			} catch (NoSuchElementException e) {
	    				System.err.println("Skipped " + numDepletionPrimersSkippedInsufficientSubprimers + " depletion primers due to insufficient inner primers.");
	    				System.err.println("Skipped " + numDepletionPrimersSkippedMatchesPrimerWithinGene + " depletion primers due to matches with existing primers within same gene.");
	    				System.err.println("Skipped " + numDepletionPrimersSkippedMatchesPrimerAcrossGenes + " depletion primers due to matches with existing primers from other genes.");
	    				System.err.println("Skipped " + numDepletionPrimersSkippedSingleMatch + " depletion primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	    				System.err.println("Skipped " + numDepletionPrimersSkippedDoubleMatch + " depletion primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");
	    				throw new NoSuchElementException("Ran out of depletion inner primers during species " + speciesgene.getSpecies() + " and gene " + speciesgene.getGeneName());
	    			}
	    		}
	    	}
	    
	    	System.out.println("\nDone assigning depletion inner primers to genes.");
	    	System.out.println("Skipped " + numDepletionPrimersSkippedInsufficientSubprimers + " depletion primers due to insufficient inner primers.");
	    	System.out.println("Skipped " + numDepletionPrimersSkippedMatchesPrimerWithinGene + " depletion primers due to matches with existing primers within same gene.");
	    	System.out.println("Skipped " + numDepletionPrimersSkippedMatchesPrimerAcrossGenes + " depletion primers due to matches with existing primers from other genes.");
	    	System.out.println("Skipped " + numDepletionPrimersSkippedSingleMatch + " depletion primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	    	System.out.println("Skipped " + numDepletionPrimersSkippedDoubleMatch + " depletion primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");
	    	System.out.println("Skipped " + numDepletionPrimersSkippedMatchesProbeWithinGene + " depletion primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
	    	System.out.println("\nDone assigning outer and inner primers to depletion probes.\n");
	    }
	    
		
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    // Assign pulldown outer primers to genes
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	    System.out.println("\nAssigning pulldown outer primers:");
	
	    // First match primers to genes
	    System.out.println("Matching outer primers to genes...");
	    HashMap<SpeciesGene, PrimerPair> pulldownOuterPrimersToUse = new HashMap<SpeciesGene, PrimerPair>();
	    
	    int numMiddlesNeededPerOuter = this.numTilingPaths + 1;
	    
	    numDone = 0;
	    int numPulldownPrimersSkippedInsufficientSubprimers=0;
	    int numPulldownPrimersSkippedMatchesPrimerBetweenGenes=0;
	    int numPulldownPrimersSkippedMatchesPrimerWithinGene=0;
	    int numPulldownPrimersSkippedSingleMatch=0;
	    int numPulldownPrimersSkippedDoubleMatch=0;
	    int numPulldownPrimersSkippedMatchesProbeWithinGene = 0;
	    
	    // Iterator over pulldown outer primers
	    Iterator<PrimerPair> pulldownOuterIter = this.pulldownPrimers.keySet().iterator();
		for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
			
	    		if(speciesgene.getPurpose() != ProbePurpose.PULLDOWN) continue;
			
	        	numDone++;
	        	if(numDone % 100 == 0) {
	            	System.out.println("\nAssigned outer primers to " + numDone + " pulldown genes.");
	                System.out.println("So far, skipped " + numPulldownPrimersSkippedInsufficientSubprimers + " pulldown primers due to insufficient middle primers.");
	                System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesPrimerBetweenGenes + " pulldown primers due to matches with existing primers from other genes.");
	                System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesPrimerWithinGene + " pulldown primers due to matches with existing primer within same gene.");
	                System.out.println("So far, skipped " + numPulldownPrimersSkippedSingleMatch + " pulldown primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	                System.out.println("So far, skipped " + numPulldownPrimersSkippedDoubleMatch + " pulldown primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
	                System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesProbeWithinGene + " pulldown primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
	        	}
	    		
				boolean assigned = false;
				while(!assigned) {
					try {
						PrimerPair thisPulldownOuter = pulldownOuterIter.next();
					
						int numMiddlesForThisOuter = this.pulldownPrimers.get(thisPulldownOuter).size();
						String thisLeftPrimer = thisPulldownOuter.getLeftPrimer();
						String thisRightPrimer = thisPulldownOuter.getRightPrimer();
						if(numMiddlesForThisOuter < numMiddlesNeededPerOuter * 2) {
							//System.out.println("Skipping pulldown outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because only has " + numMiddlesForThisOuter + " middle primers.");
							numPulldownPrimersSkippedInsufficientSubprimers++;
							continue;
						}
						// check if the primer matches a primer that has already been kept
						if(primerCrossPrimesExistingPrimer(thisLeftPrimer) || primerCrossPrimesExistingPrimer(thisRightPrimer)) {
							//System.out.println("Skipping pulldown outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because it matches an existing primer.");
							numPulldownPrimersSkippedMatchesPrimerBetweenGenes++;
							continue;
						}
						// check if the primer matches any probe
						// first check for a perfect longer match
						if(hasSingleMatch(thisLeftPrimer, thisRightPrimer, singleKmers)) {
							//System.out.println("Skipping pulldown outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because it has a perfect " + MIN_SINGLE_PRIMER_MATCH + "-mer match to more than " + MIN_KMER_OCCURENCES + " probes.");
							numPulldownPrimersSkippedSingleMatch++;
							continue;
						}
						// now check for both primers to have a perfect shorter match
						if(hasDoubleMatch(this, thisLeftPrimer, thisRightPrimer, doubleKmers)) {
							//System.out.println("Skipping pulldown outer primer " + thisLeftPrimer + " " + thisRightPrimer + " because there is a probe with perfect " + MIN_DOUBLE_PRIMER_MATCH + "-mer matches to both primers.");
							numPulldownPrimersSkippedDoubleMatch++;
							continue;
						}
	        			// check if the primer matches a probe in the gene
	        			if(leftPrimerHasMatchInProbeSet(thisLeftPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	        				numPulldownPrimersSkippedMatchesProbeWithinGene++;
	        				continue;
	        			}
	        			if(rightPrimerHasMatchInProbeSet(thisRightPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	        				numPulldownPrimersSkippedMatchesProbeWithinGene++;
	        				continue;
	        			}
	
						// if still alive at this point, assign this primer to the gene
						// add the primers to the list of kept primers
						pulldownOuterPrimersToUse.put(speciesgene, thisPulldownOuter);
						this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisLeftPrimer,this.betweenGenePrimerDimerMatchLength));
						this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisRightPrimer,this.betweenGenePrimerDimerMatchLength));
						//System.out.println("For species " + speciesgene.getSpecies() + " and gene " + speciesgene.getGeneName() + " using primer " + thisLeftPrimer + " " + thisRightPrimer);
						assigned = true;
					} catch (NoSuchElementException e) {
		                System.err.println("Skipped " + numPulldownPrimersSkippedInsufficientSubprimers + " pulldown primers due to insufficient middle primers.");
	                    System.err.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerBetweenGenes + " pulldown primers due to matches with existing primers from other genes.");
	                    System.err.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerWithinGene + " pulldown primers due to matches with existing primer within same gene.");
		                System.err.println("Skipped " + numPulldownPrimersSkippedSingleMatch + " pulldown primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
		                System.err.println("Skipped " + numPulldownPrimersSkippedDoubleMatch + " pulldown primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
	                    System.err.println("Skipped " + numPulldownPrimersSkippedMatchesProbeWithinGene + " pulldown primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
		                throw new NoSuchElementException("Ran out of pulldown outer primers during species " + speciesgene.getSpecies() + " and gene " + speciesgene.getGeneName());
					}
				}
				
			}
			
		System.out.println("\nDone assigning pulldown outer primers to genes.");
	    System.out.println("So far, skipped " + numPulldownPrimersSkippedInsufficientSubprimers + " pulldown primers due to insufficient middle primers.");
	    System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesPrimerBetweenGenes + " pulldown primers due to matches with existing primers from other genes.");
	    System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesPrimerWithinGene + " pulldown primers due to matches with existing primer within same gene.");
	    System.out.println("So far, skipped " + numPulldownPrimersSkippedSingleMatch + " pulldown primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	    System.out.println("So far, skipped " + numPulldownPrimersSkippedDoubleMatch + " pulldown primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
	    System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesProbeWithinGene + " pulldown primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
		
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    // Assign pulldown middle primers to tiling paths and qpcr
		// Assign pulldown inner primers to even/odd and domain
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
		System.out.println("\nAssigning pulldown middle and inner primers to tiling paths, qPCR, even/odd and domain.");
		
	    HashMap<SpeciesGene, HashMap<Integer, PrimerPair> > pulldownTilingPathLeftPrimersToUse = new HashMap<SpeciesGene, HashMap<Integer, PrimerPair> >();
	    HashMap<SpeciesGene, HashMap<Integer, PrimerPair> > pulldownQpcrRightPrimersToUse = new HashMap<SpeciesGene, HashMap<Integer, PrimerPair> >();
	    HashMap<SpeciesGene, HashMap<Integer, HashMap<Integer,PrimerPair>>> pulldownEvenOddLeftPrimersToUse = new HashMap<SpeciesGene, HashMap<Integer, HashMap<Integer,PrimerPair>>>();
	    HashMap<SpeciesGene, HashMap<Integer, HashMap<Integer,PrimerPair>>> pulldownDomainRightPrimersToUse = new HashMap<SpeciesGene, HashMap<Integer, HashMap<Integer,PrimerPair>>>();
	
	    // Iterators over the middle primers corresponding to each outer primer
	    HashMap<SpeciesGene, Iterator<PrimerPair>> pulldownTilingPathLeftIterators = new HashMap<SpeciesGene, Iterator<PrimerPair>>();
	    HashMap<SpeciesGene, Iterator<PrimerPair>> pulldownQpcrRightIterators = new HashMap<SpeciesGene, Iterator<PrimerPair>>();
	
	    for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
	    	if(speciesgene.getPurpose() != ProbePurpose.PULLDOWN) continue;
	    	Iterator<PrimerPair> leftIter = this.pulldownPrimers.get(pulldownOuterPrimersToUse.get(speciesgene)).keySet().iterator();
	    	Iterator<PrimerPair> rightIter = this.pulldownPrimers.get(pulldownOuterPrimersToUse.get(speciesgene)).keySet().iterator();
	    	pulldownTilingPathLeftIterators.put(speciesgene, leftIter);
	    	pulldownQpcrRightIterators.put(speciesgene, rightIter);
	    }
	    
	    // Record the tiling path IDs
	    Integer[] theTilingPaths = new Integer[this.numTilingPaths + 1];
	    theTilingPaths[0] = Integer.valueOf(-1);
	    for(int i = 1; i <= this.numTilingPaths; i++) theTilingPaths[i] = Integer.valueOf((i-1) * TILING_DENSITY);
	    
	    numDone = 0;
	    int numMiddlesSkippedBecauseInnersRanOut = 0;
	    
	    // Assign and filter pulldown middle primers
	    for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
	    	
	    	// Only pulldown genes
	    	if(speciesgene.getPurpose() != ProbePurpose.PULLDOWN) continue;
	    	
	    	numDone++;
	    	if(numDone % 100 == 0) {
	        	System.out.println("\nAssigned middle and inner primers to " + numDone + " pulldown genes.");
	            System.out.println("So far, skipped " + numPulldownPrimersSkippedInsufficientSubprimers + " pulldown primers due to insufficient subprimers.");
	            System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesPrimerBetweenGenes + " pulldown primers due to matches with existing primers from other genes.");
	            System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesPrimerWithinGene + " pulldown primers due to matches with existing primer within same gene.");
	            System.out.println("So far, skipped " + numPulldownPrimersSkippedSingleMatch + " pulldown primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	            System.out.println("So far, skipped " + numPulldownPrimersSkippedDoubleMatch + " pulldown primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
				System.out.println("So far, skipped " + numMiddlesSkippedBecauseInnersRanOut + " pulldown middle primers because inner primers ran out.");
	            System.out.println("So far, skipped " + numPulldownPrimersSkippedMatchesProbeWithinGene + " pulldown primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
	    	}
	    	
	    	HashMap<Integer, PrimerPair> tilingPathLeftPrimersForThisGene = new HashMap<Integer, PrimerPair>();
	    	HashMap<Integer, PrimerPair> qpcrRightPrimersForThisGene = new HashMap<Integer, PrimerPair>();
	    	HashMap<Integer, HashMap<Integer, PrimerPair>> evenOddLeftPrimersForThisGene = new HashMap<Integer, HashMap<Integer, PrimerPair>>();
	    	HashMap<Integer, HashMap<Integer, PrimerPair>> domainRightPrimersForThisGene = new HashMap<Integer, HashMap<Integer, PrimerPair>>();
	   	
	    	
	    	ArrayList<String> keptLeftPrimersForThisGene = new ArrayList<String>();
	    	ArrayList<String> keptRightPrimersForThisGene = new ArrayList<String>();
	    	
	    	PrimerPair outerPrimer = pulldownOuterPrimersToUse.get(speciesgene);
	    	String leftOuterPrimer = outerPrimer.getLeftPrimer();
	    	String rightOuterPrimer = outerPrimer.getRightPrimer();
	    	
	    	keptLeftPrimersForThisGene.add(leftOuterPrimer);
	    	keptRightPrimersForThisGene.add(rightOuterPrimer);
	
	    	// assign tiling path left primers
	    	for(int i=0; i < theTilingPaths.length; i++) {
	    		boolean middleAssigned = false;
	    		while(!middleAssigned) {
	    			try {
	    				PrimerPair thisPulldownLeftMiddle = pulldownTilingPathLeftIterators.get(speciesgene).next();
	    			
	    				String thisLeftPrimer = thisPulldownLeftMiddle.getLeftPrimer();
	    				// check if the primer matches a primer that has already been kept
	    				if(primerCrossPrimesInSet(thisLeftPrimer,keptRightPrimersForThisGene,this.withinGenePrimerDimerMatchLength)) {
	    					numPulldownPrimersSkippedMatchesPrimerWithinGene++;
	    					continue;
	    				}
	    				if(primerCrossPrimesExistingPrimer(thisLeftPrimer)) {
	    					//System.out.println("Skipping pulldown left middle primer " + thisLeftPrimer + " because it matches an existing primer.");
	    					numPulldownPrimersSkippedMatchesPrimerBetweenGenes++;
	    					continue;
	    				}
	    				// check if the primer matches any probe
	    				// first check for a perfect longer match
	    				if(hasSingleMatch(thisLeftPrimer, thisLeftPrimer, singleKmers)) {
	    					//System.out.println("Skipping pulldown left middle primer " + thisLeftPrimer + " because it has a perfect " + MIN_SINGLE_PRIMER_MATCH + "-mer match to more than " + MIN_KMER_OCCURENCES + " probes.");
	    					numPulldownPrimersSkippedSingleMatch++;
	    					continue;
	    				}
	    				// check if the primer has a double match with any right primer already kept for this gene
	    				if(leftPrimerHasDoubleMatchWithRightSet(this, thisLeftPrimer, keptRightPrimersForThisGene, doubleKmers)) {
	    					//System.out.println("Skipping pulldown left middle primer " + thisLeftPrimer + " because there is already a right primer for the gene and a probe with perfect " + MIN_DOUBLE_PRIMER_MATCH + "-mer matches to both primers.");
	    					numPulldownPrimersSkippedDoubleMatch++;
	    					continue;
	    				}
	        			// check if the primer matches a probe in the gene
	        			if(leftPrimerHasMatchInProbeSet(thisLeftPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	        				numPulldownPrimersSkippedMatchesProbeWithinGene++;
	        				continue;
	        			}
	
	    				// if still alive at this point, assign this primer to the tiling path
	    				// add the primer to the list of kept primers
	    				// add the primer to the map of left middle primers for the gene
	    				this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisLeftPrimer,this.betweenGenePrimerDimerMatchLength));
	    				tilingPathLeftPrimersForThisGene.put(theTilingPaths[i], thisPulldownLeftMiddle);
	    				keptLeftPrimersForThisGene.add(thisLeftPrimer);
	    				//System.out.println("Assigned left primer " + thisLeftPrimer + " to tiling path " + theTilingPaths[i]);
	    				middleAssigned = true;
	    			
	    			
	    				// assign the even/odd inner primers for this middle primer
	    				Iterator<PrimerPair> innerIter = this.pulldownPrimers.get(outerPrimer).get(thisPulldownLeftMiddle).iterator();
	    				HashMap<Integer,PrimerPair> evenOddPrimersForThisTilingPath = new HashMap<Integer,PrimerPair>();
	    				for(int j=-1; j <= 1; j++) {
	    					boolean innerAssigned = false;
	    					while(!innerAssigned) {
	    						try {
	    							PrimerPair thisPulldownLeftInner = innerIter.next();
	    						
	    							String thisLeftInnerPrimer = thisPulldownLeftInner.getLeftPrimer();
	    							// check if the primer matches a primer that has already been kept
	    	        				if(primerCrossPrimesInSet(thisLeftInnerPrimer,keptRightPrimersForThisGene,this.withinGenePrimerDimerMatchLength)) {
	    	        					numPulldownPrimersSkippedMatchesPrimerWithinGene++;
	    	        					continue;
	    	        				}
	    	        				if(primerCrossPrimesExistingPrimer(thisLeftInnerPrimer)) {
	    								//System.out.println("Skipping pulldown left inner primer " + thisLeftInnerPrimer + " because it matches an existing primer.");
	    								numPulldownPrimersSkippedMatchesPrimerBetweenGenes++;
	    								continue;
	    							}
	    							// check if the primer matches any probe
	    							// first check for a perfect longer match
	    							if(hasSingleMatch(thisLeftInnerPrimer, thisLeftInnerPrimer, singleKmers)) {
	    								//System.out.println("Skipping pulldown left inner primer " + thisLeftInnerPrimer + " because it has a perfect " + MIN_SINGLE_PRIMER_MATCH + "-mer match to more than " + MIN_KMER_OCCURENCES + " probes.");
	    								numPulldownPrimersSkippedSingleMatch++;
	    								continue;
	    							}
	    							// check if the primer has a double match with any right primer already kept for this gene
	    							if(leftPrimerHasDoubleMatchWithRightSet(this, thisLeftInnerPrimer, keptRightPrimersForThisGene, doubleKmers)) {
	    								//System.out.println("Skipping pulldown left inner primer " + thisLeftInnerPrimer + " because there is already a right primer for the gene and a probe with perfect " + MIN_DOUBLE_PRIMER_MATCH + "-mer matches to both primers.");
	    								numPulldownPrimersSkippedDoubleMatch++;
	    								continue;
	    							}
	    	            			// check if the primer matches a probe in the gene
	    	            			if(leftPrimerHasMatchInProbeSet(thisLeftPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	    	            				numPulldownPrimersSkippedMatchesProbeWithinGene++;
	    	            				continue;
	    	            			}
	
	        			
	    							// if still alive at this point, assign this primer to the even/odd
	    							// add the primer to the list of kept primers
	    							// add the primer to the map of left middle primers for the gene
	    							this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisLeftInnerPrimer,this.betweenGenePrimerDimerMatchLength));
	    							evenOddPrimersForThisTilingPath.put(Integer.valueOf(j), thisPulldownLeftInner);
	    							keptLeftPrimersForThisGene.add(thisLeftInnerPrimer);
	    							//System.out.println("Assigned left primer " + thisLeftInnerPrimer + " to even/odd " + j);
	    							innerAssigned = true;
	    						} catch (NoSuchElementException e) {
	    							// if ran out of inner primers, kill this middle primer and start over
	    							middleAssigned = false;
	    							numMiddlesSkippedBecauseInnersRanOut++;
	    							break;
	    						}
	    					}
	    				}
	    			
	    				if(middleAssigned) evenOddLeftPrimersForThisGene.put(theTilingPaths[i], evenOddPrimersForThisTilingPath);
	    			} catch (NoSuchElementException e) {
		                System.err.println("Skipped " + numPulldownPrimersSkippedInsufficientSubprimers + " pulldown primers due to insufficient subprimers.");
	                    System.out.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerBetweenGenes + " pulldown primers due to matches with existing primers from other genes.");
	                    System.out.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerWithinGene + " pulldown primers due to matches with existing primer within same gene.");
		                System.err.println("Skipped " + numPulldownPrimersSkippedSingleMatch + " pulldown primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
		                System.err.println("Skipped " + numPulldownPrimersSkippedDoubleMatch + " pulldown primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
	    				System.err.println("Skipped " + numMiddlesSkippedBecauseInnersRanOut + " pulldown middle primers because inner primers ran out.");
	                    System.err.println("Skipped " + numPulldownPrimersSkippedMatchesProbeWithinGene + " pulldown primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
	    				throw new NoSuchElementException("Ran out of pulldown middle primers during gene " + speciesgene.getSpeciesName() + " " + speciesgene.getGeneName() + " and tiling path " + theTilingPaths[i]);
	    			}
	    			
	    		}
	    	}
	    	
	    	// assign qpcr right primers
	    	for(int i= -1; i <= 1; i++) {
	    		boolean middleAssigned = false;
	    		while(!middleAssigned) {
	    			try {
	    				PrimerPair thisPulldownRightMiddle = pulldownQpcrRightIterators.get(speciesgene).next();
	    				String thisRightPrimer = thisPulldownRightMiddle.getRightPrimer();
	    				// check if the primer matches a primer that has already been kept
	    				if(primerCrossPrimesInSet(thisRightPrimer,keptLeftPrimersForThisGene,this.withinGenePrimerDimerMatchLength)) {
	    					numPulldownPrimersSkippedMatchesPrimerWithinGene++;
	    					continue;
	    				}
	    				if(primerCrossPrimesExistingPrimer(thisRightPrimer)) {
	    					//System.out.println("Skipping pulldown right middle primer " + thisRightPrimer + " because it matches an existing primer.");
	    					numPulldownPrimersSkippedMatchesPrimerBetweenGenes++;
	    					continue;
	    				}
	    				// check if the primer matches any probe
	    				// first check for a perfect longer match
	    				if(hasSingleMatch(thisRightPrimer, thisRightPrimer, singleKmers)) {
	    					//System.out.println("Skipping pulldown right middle primer " + thisRightPrimer + " because it has a perfect " + MIN_SINGLE_PRIMER_MATCH + "-mer match to more than " + MIN_KMER_OCCURENCES + " probes.");
	    					numPulldownPrimersSkippedSingleMatch++;
	    					continue;
	    				}
	    				// check if the primer has a double match with any left primer already kept for this gene
	    				if(rightPrimerHasDoubleMatchWithLeftSet(this, keptLeftPrimersForThisGene, thisRightPrimer, doubleKmers)) {
	    					//System.out.println("Skipping pulldown right middle primer " + thisRightPrimer + " because there is already a left primer for the gene and a probe with perfect " + MIN_DOUBLE_PRIMER_MATCH + "-mer matches to both primers.");
	    					numPulldownPrimersSkippedDoubleMatch++;
	    					continue;
	    				}
	        			// check if the primer matches a probe in the gene
	        			if(rightPrimerHasMatchInProbeSet(thisRightPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	        				numPulldownPrimersSkippedMatchesProbeWithinGene++;
	        				continue;
	        			}
	
	    				// if still alive at this point, assign this primer to the qpcr
	    				// add the primer to the list of kept primers
	    				// add the primer to the map of left middle primers for the gene
	    				this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisRightPrimer,this.betweenGenePrimerDimerMatchLength));
	    				qpcrRightPrimersForThisGene.put(Integer.valueOf(i), thisPulldownRightMiddle);
	    				keptRightPrimersForThisGene.add(thisRightPrimer);
	    				middleAssigned = true;
	    			
	    				// assign the domain inner primers for this middle primer
	    				Iterator<PrimerPair> innerRightIter = this.pulldownPrimers.get(outerPrimer).get(thisPulldownRightMiddle).iterator();
	    				HashMap<Integer,PrimerPair> domainPrimersForThisQpcr = new HashMap<Integer,PrimerPair>();
	    				for(int j=-1; j < NUM_DOMAINS; j++) {
	    					boolean innerAssigned = false;
	    					while(!innerAssigned) {
	    						try {
	    								PrimerPair thisPulldownRightInner = innerRightIter.next();
	    								String thisRightInnerPrimer = thisPulldownRightInner.getRightPrimer();
	    								// check if the primer matches a primer that has already been kept
	    		        				if(primerCrossPrimesInSet(thisRightInnerPrimer,keptLeftPrimersForThisGene,this.withinGenePrimerDimerMatchLength)) {
	    		        					numPulldownPrimersSkippedMatchesPrimerWithinGene++;
	    		        					continue;
	    		        				}
	    		        				if(primerCrossPrimesExistingPrimer(thisRightInnerPrimer)) {
	    									//System.out.println("Skipping pulldown right inner primer " + thisRightInnerPrimer + " because it matches an existing primer.");
	    									numPulldownPrimersSkippedMatchesPrimerBetweenGenes++;
	    									continue;
	    								}
	    								// check if the primer matches any probe
	    								// first check for a perfect longer match
	    								if(hasSingleMatch(thisRightInnerPrimer, thisRightInnerPrimer, singleKmers)) {
	    									//System.out.println("Skipping pulldown right inner primer " + thisRightInnerPrimer + " because it has a perfect " + MIN_SINGLE_PRIMER_MATCH + "-mer match to more than " + MIN_KMER_OCCURENCES + " probes.");
	    									numPulldownPrimersSkippedSingleMatch++;
	    									continue;
	    								}
	    								// check if the primer has a double match with any left primer already kept for this gene
	    								if(rightPrimerHasDoubleMatchWithLeftSet(this, keptLeftPrimersForThisGene, thisRightInnerPrimer, doubleKmers)) {
	    									//System.out.println("Skipping pulldown right inner primer " + thisRightInnerPrimer + " because there is already a left primer for the gene and a probe with perfect " + MIN_DOUBLE_PRIMER_MATCH + "-mer matches to both primers.");
	    									numPulldownPrimersSkippedDoubleMatch++;
	    									continue;
	    								}
	    		            			// check if the primer matches a probe in the gene
	    		            			if(rightPrimerHasMatchInProbeSet(thisRightPrimer,this.geneProbeSet.get(speciesgene),this.singlePrimerMatchToProbeWithinGene)) {
	    		            				numPulldownPrimersSkippedMatchesProbeWithinGene++;
	    		            				continue;
	    		            			}
	
	    								// if still alive at this point, assign this primer to the domain
	    								// add the primer to the list of kept primers
	    								// add the primer to the map of right middle primers for the gene
	    								this.keptPrimerKmers.addAll(getAllSubstringsHashSet(thisRightInnerPrimer,this.betweenGenePrimerDimerMatchLength));
	    								domainPrimersForThisQpcr.put(Integer.valueOf(j), thisPulldownRightInner);
	    								keptRightPrimersForThisGene.add(thisRightInnerPrimer);
	    								innerAssigned = true;
	    						} catch (NoSuchElementException e) {
	    							// if ran out of inner primers, kill this middle primer and start over
	    							middleAssigned = false;
	    							numMiddlesSkippedBecauseInnersRanOut++;
	    							break;
	    						}
	    					}
	    				}
	    			
	    				if(middleAssigned) domainRightPrimersForThisGene.put(Integer.valueOf(i), domainPrimersForThisQpcr);
	    			} catch (NoSuchElementException e) {
		                System.err.println("Skipped " + numPulldownPrimersSkippedInsufficientSubprimers + " pulldown primers due to insufficient subprimers.");
	                    System.err.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerBetweenGenes + " pulldown primers due to matches with existing primers from other genes.");
	                    System.err.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerWithinGene + " pulldown primers due to matches with existing primer within same gene.");
		                System.err.println("Skipped " + numPulldownPrimersSkippedSingleMatch + " pulldown primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
		                System.err.println("Skipped " + numPulldownPrimersSkippedDoubleMatch + " pulldown primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
	    				System.err.println("Skipped " + numMiddlesSkippedBecauseInnersRanOut + " pulldown middle primers because inner primers ran out.");
	                    System.err.println("Skipped " + numPulldownPrimersSkippedMatchesProbeWithinGene + " pulldown primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
	    				throw new NoSuchElementException("Ran out of pulldown middle primers during gene " + speciesgene.getSpeciesName() + " " + speciesgene.getGeneName() + " and qpcr " + i);
	    			}
	    		}
	    	}
	
	    	pulldownTilingPathLeftPrimersToUse.put(speciesgene, tilingPathLeftPrimersForThisGene);
	    	pulldownQpcrRightPrimersToUse.put(speciesgene, qpcrRightPrimersForThisGene);
	    	pulldownEvenOddLeftPrimersToUse.put(speciesgene, evenOddLeftPrimersForThisGene);
	    	pulldownDomainRightPrimersToUse.put(speciesgene, domainRightPrimersForThisGene);
	    	
	    }
	    
		System.out.println("Done assigning pulldown middle and inner primers.");
	    System.out.println("Skipped " + numPulldownPrimersSkippedInsufficientSubprimers + " pulldown primers due to insufficient subprimers.");
	    System.out.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerBetweenGenes + " pulldown primers due to matches with existing primers from other genes.");
	    System.out.println("Skipped " + numPulldownPrimersSkippedMatchesPrimerWithinGene + " pulldown primers due to matches with existing primer within same gene.");
	    System.out.println("Skipped " + numPulldownPrimersSkippedSingleMatch + " pulldown primers due to perfect " + this.singlePrimerMatchToProbe + "-mer matches with a probe.");
	    System.out.println("Skipped " + numPulldownPrimersSkippedDoubleMatch + " pulldown primers due to perfect double " + this.doublePrimerMatchToProbe + "-mer matches with a probe.");        		
		System.out.println("Skipped " + numMiddlesSkippedBecauseInnersRanOut + " pulldown middle primers because inner primers ran out.");
	    System.out.println("Skipped " + numPulldownPrimersSkippedMatchesProbeWithinGene + " pulldown primers due to perfect single " + this.singlePrimerMatchToProbeWithinGene + "-mer matches to a probe within the gene.");        		
	
		System.out.println("\nDone assigning outer and inner primers to depletion probes.");
	
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    // Assign pulldown primers to probes
	    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
		System.out.println("\nAssigning pulldown primers to probes.");
		
		this.pulldownOuterLeftPrimers.clear();
		this.pulldownMiddleLeftPrimers.clear();
		this.pulldownInnerLeftPrimers.clear();
		this.pulldownOuterRightPrimers.clear();
		this.pulldownMiddleRightPrimers.clear();
		this.pulldownInnerRightPrimers.clear();
		
		for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
			if(speciesgene.getPurpose() != ProbePurpose.PULLDOWN) continue;
			for(Probe probe : this.geneProbeSet.get(speciesgene)) {
				int tilingPath = probe.getTilingPath();
				int evenOdd = probe.getEvenOdd();
				int qpcr = probe.getContainsQpcrPrimer();
				int domain = probe.getDomain();
				this.pulldownOuterLeftPrimers.put(probe,pulldownOuterPrimersToUse.get(speciesgene).getLeftPrimer());
				this.pulldownOuterRightPrimers.put(probe, pulldownOuterPrimersToUse.get(speciesgene).getRightPrimer());
				this.pulldownMiddleLeftPrimers.put(probe, pulldownTilingPathLeftPrimersToUse.get(speciesgene).get(Integer.valueOf(tilingPath)).getLeftPrimer());
				this.pulldownMiddleRightPrimers.put(probe, pulldownQpcrRightPrimersToUse.get(speciesgene).get(Integer.valueOf(qpcr)).getRightPrimer());
				this.pulldownInnerLeftPrimers.put(probe,pulldownEvenOddLeftPrimersToUse.get(speciesgene).get(Integer.valueOf(tilingPath)).get(Integer.valueOf(evenOdd)).getLeftPrimer());
				this.pulldownInnerRightPrimers.put(probe,pulldownDomainRightPrimersToUse.get(speciesgene).get(Integer.valueOf(qpcr)).get(Integer.valueOf(domain)).getRightPrimer());
			}
		}
	    
		System.out.println("Done assigning pulldown primers to probes.");
	    
	}

	/**
	 * Design primers for depletion probes
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void designDepletionPrimers() throws IOException, InterruptedException {
		
		// Check if primers have already been created and if so, read from the file
		File depletionPrimersFile = new File("SecondPrimers_depletion.out");
		if(depletionPrimersFile.exists()) {
			System.out.println("\nNot designing new depletion primers. Reading existing primers from file SecondPrimers_depletion.out.");
			FileReader reader = new FileReader(depletionPrimersFile);
			BufferedReader b = new BufferedReader(reader);
			boolean start = false;
			PrimerPair outerPrimer = null;
			ArrayList<PrimerPair> innerPrimers = new ArrayList<PrimerPair>();
			this.depletionPrimers.clear();
			int outerID = -1;
			int innerID = 0;
			int outerPrimers = 0;
			int numInnerPrimers = 0;
			while(b.ready()) {
				String line = b.readLine();
				StringParser p = new StringParser();
				p.parse(line);
				if(p.asString(0).equals("OUTER_PRIMER")) {
					outerPrimers++;
					if(start) {
						// Copy latest primers into new containers and add
						ArrayList<PrimerPair> toAdd = new ArrayList<PrimerPair>();
						toAdd.addAll(innerPrimers);
						String[] fields = {"0", "0", outerPrimer.getLeftPrimer(), "0", outerPrimer.getRightPrimer(), "0", "0", "0", "0"};
						PrimerPair newOuter = new PrimerPair(fields);
						this.depletionPrimers.put(newOuter, toAdd);
					}
					numInnerPrimers += innerPrimers.size();
					start = true;
					innerPrimers.clear();
					innerID = 0;
					outerID++;
					String newOuterID = "Primer_" + Integer.valueOf(outerID).toString();
					String[] fields = p.getStringArray();
					fields[0] = newOuterID;
					PrimerPair tempPrimer = new PrimerPair(fields);
					outerPrimer = tempPrimer;
				} else {
					String newInnerID = "Primer_" + Integer.valueOf(outerID).toString() + "_" + Integer.valueOf(innerID).toString();
					innerID++;
					String[] fields = p.getStringArray();
					fields[0] = newInnerID;
					PrimerPair tempPrimer = new PrimerPair(fields);
					innerPrimers.add(tempPrimer);
				}
			}
			// Copy the last primers into new containers and add
			String[] fields = {"0", "0", outerPrimer.getLeftPrimer(), "0", outerPrimer.getRightPrimer(), "0", "0", "0", "0"};
			PrimerPair newOuter = new PrimerPair(fields);
			ArrayList<PrimerPair> toAdd = new ArrayList<PrimerPair>();
			toAdd.addAll(innerPrimers);
			this.depletionPrimers.put(newOuter, toAdd);
			numInnerPrimers += innerPrimers.size();
			System.out.println("Read " + outerPrimers + " outer primers and " + numInnerPrimers + " inner primers.");
			b.close();
		} else {
			// File does not exist so need to design new primers
			System.out.println("\nDesigning PCR tails for depletion probes.");
			System.out.println("There are " + Species.numberOfSpecies() + " species and " + RnaClass.numberOfClasses() + " RNA classes.");
			int numDepletionGenes = 0;
			int numDepletionProbes = 0;
			// Count how many genes are for depletion
			for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
				if(speciesgene.getPurpose() == ProbePurpose.DEPLETE_MAIN_ARRAY || speciesgene.getPurpose() == ProbePurpose.DEPLETE_ALL_MRNA) {
					numDepletionGenes ++;
					numDepletionProbes += this.geneProbeSet.get(speciesgene).size();
				}
			}
			System.out.println("There are " + numDepletionGenes + " genes for depletion with a total of " + numDepletionProbes + " probes.");
			// Use the constructor for two overlapping primers
			// Create twice as many primers as necessary in each dimension
			PcrTailDesigner design = new PcrTailDesigner(30, 20, 20, 20, 20, 2 * Species.numberOfSpecies()*RnaClass.numberOfClasses(), 2 * numDepletionGenes, "depletion");
			design.createPrimers();
			this.depletionPrimers = design.getSecondPrimers();
		}
	}

	/**
	 * Design pulldown and depletion primers, then combine all together, then assign to probes
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws SearchException
	 */
	private void designPrimers() throws IOException, InterruptedException, SearchException {
		System.out.println("Designing primers...");
		this.designPulldownPrimers();
		if(this.mainArrayHasDepletionSequences) this.designDepletionPrimers();
		System.out.println("\nAssigning and filtering primers...\n");
		this.assignAndFilterPrimers();
		System.out.println("Done designing primers.\n");
	}

	/** 
	 * Design a set of primers for pulldown with extras in case they cross-hyb with pulldown primers
	 */
	private void designPulldownPrimers() throws IOException, InterruptedException{
		
		// Check if primers have already been created and if so, read from the file
		File pulldownPrimersFile = new File("ThirdPrimers_pulldown.out");
		if(pulldownPrimersFile.exists()) {
			System.out.println("\nNot designing new pulldown primers. Reading existing primers from file ThirdPrimers_pulldown.out.");
			FileReader reader = new FileReader(pulldownPrimersFile);
			BufferedReader b = new BufferedReader(reader);
			boolean start = false;
			boolean middlestart = false;
			PrimerPair outerPrimer = null;
			PrimerPair middlePrimer = null;
			ArrayList<PrimerPair> innerPrimers = new ArrayList<PrimerPair>();
			HashMap< PrimerPair, ArrayList<PrimerPair>> middlePrimers = new HashMap< PrimerPair, ArrayList<PrimerPair>>();
			this.pulldownPrimers.clear();
			int outerID = -1;
			int middleID = -1;
			int innerID = 0;
			int numOuterPrimers = 0;
			int numMiddlePrimers = 0;
			int numInnerPrimers = 0;
			while(b.ready()) {
				String line = b.readLine();
				StringParser p = new StringParser();
				p.parse(line);
				if(p.asString(0).equals("OUTER_PRIMER")) {
					if(start) {
						ArrayList<PrimerPair> toAdd = new ArrayList<PrimerPair>();
						toAdd.addAll(innerPrimers);
						String[] fields = {"0", "0", middlePrimer.getLeftPrimer(), "0", middlePrimer.getRightPrimer(), "0", "0", "0", "0"};
						PrimerPair newMiddle = new PrimerPair(fields);
						middlePrimers.put(newMiddle, toAdd);
						HashMap<PrimerPair, ArrayList<PrimerPair>> toAdd2 = new HashMap<PrimerPair, ArrayList<PrimerPair>>();
						toAdd2.putAll(middlePrimers);
						String[] fields2 = {"0", "0", outerPrimer.getLeftPrimer(), "0", outerPrimer.getRightPrimer(), "0", "0", "0", "0"};
						PrimerPair newOuter = new PrimerPair(fields2);
						this.pulldownPrimers.put(newOuter, toAdd2);
						numMiddlePrimers += toAdd2.size();
					}
					numOuterPrimers++;
					start = true;
					innerPrimers.clear();
					middlePrimers.clear();
					innerID = 0;
					middleID = -1;
					outerID++;
					String[] fields = {"0", "0", p.asString(2), "0", p.asString(4), "0", "0", "0", "0"};
					PrimerPair tempPrimer = new PrimerPair(fields);
					outerPrimer = tempPrimer;
				} else if(p.asString(0).equals("SECOND_PRIMER")) {
					if(middlestart) {
						String[] fields = {"0", "0", middlePrimer.getLeftPrimer(), "0", middlePrimer.getRightPrimer(), "0", "0", "0", "0"};
						PrimerPair newMiddle = new PrimerPair(fields);
						ArrayList<PrimerPair> toAdd = new ArrayList<PrimerPair>();
						toAdd.addAll(innerPrimers);					
						middlePrimers.put(newMiddle, toAdd);
					}
					middlestart = true;
					innerPrimers.clear();
					innerID = 0;
					middleID++;
					String[] fields = {"0", "0", p.asString(2), "0", p.asString(4), "0", "0", "0", "0"};
					PrimerPair tempPrimer = new PrimerPair(fields);
					middlePrimer = tempPrimer;
				} else {
					String newInnerID = "Primer_" + Integer.valueOf(outerID).toString() + "_" + Integer.valueOf(middleID).toString() + "_" + Integer.valueOf(innerID).toString();
					innerID++;
					String[] fields = p.getStringArray();
					fields[0] = newInnerID;
					PrimerPair tempPrimer = new PrimerPair(fields);
					innerPrimers.add(tempPrimer);					
					numInnerPrimers++;
				}
			}
			String[] fields = {"0", "0", middlePrimer.getLeftPrimer(), "0", middlePrimer.getRightPrimer(), "0", "0", "0", "0"};
			PrimerPair newMiddle = new PrimerPair(fields);
			ArrayList<PrimerPair> toAdd = new ArrayList<PrimerPair>();
			toAdd.addAll(innerPrimers);					
			middlePrimers.put(newMiddle, toAdd);			
			HashMap<PrimerPair, ArrayList<PrimerPair>> toAdd2 = new HashMap<PrimerPair, ArrayList<PrimerPair>>();
			toAdd2.putAll(middlePrimers);
			String[] fields2 = {"0", "0", outerPrimer.getLeftPrimer(), "0", outerPrimer.getRightPrimer(), "0", "0", "0", "0"};
			PrimerPair newOuter = new PrimerPair(fields2);
			numMiddlePrimers += toAdd2.size();
			this.pulldownPrimers.put(newOuter, toAdd2);
			b.close();
			System.out.println("Read " + numOuterPrimers + " outer primers, " + numMiddlePrimers + " middle primers, and " + numInnerPrimers + " inner primers.");
		} else {
			// File does not exist so need to design new primers
			System.out.println("\nDesigning PCR tails for pulldown probes.");
			System.out.println("There are " + this.geneProbeSet.keySet().size() + " genes, " + this.numTilingPaths + " tiling paths and " + NUM_DOMAINS + " domains.");
			int numPulldownGenes = 0;
			int numPulldownProbes = 0;
			// Count how many genes are for pulldown
			for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
				if(speciesgene.getPurpose() == ProbePurpose.PULLDOWN) {
					numPulldownGenes ++;
					numPulldownProbes += this.geneProbeSet.get(speciesgene).size();
				}
			}
			System.out.println("There are " + numPulldownGenes + " genes for pulldown with a total of " + numPulldownProbes + " probes.");
			// Number of middle primers, for example, is number of tiling paths times 2 for qpcr primer, so each combination can have a primer pair, then times 2 again for a cushion
			PcrTailDesigner design = new PcrTailDesigner(25, 15, 15, 5, 15, 15, 15, 5, 15, 2 * numPulldownGenes, 5 * this.numTilingPaths * 3, 10 * NUM_DOMAINS * 3, "pulldown");
			design.createPrimers();
			this.pulldownPrimers = design.getThirdPrimers();
		}
	}

	/**
	 *  Get a qPCR primer pair for a sequence
	 * @param seq the sequence
	 * @return the primer pair
	 * @throws IOException
	 */
	public static PrimerPair designQpcrPrimer(String seq) throws IOException {
		// Jesse's qPCR primer parameters
		Primer3Configuration best = Primer3ConfigurationFactory.getRAPqPCRConfiguration();
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
		Sequence mRNA = new Sequence("Gene");
		mRNA.setSequenceBases(seq);
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getId());
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		if(pp == null || pp.isEmpty()) {
			p3io.endPrimer3Communications();
			return null;
		}
		p3io.endPrimer3Communications();
		TreeSet<PrimerPair> result = new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){result.add(pair);}
		}
		Iterator<PrimerPair> iter = result.iterator();
		if(iter.hasNext()) return iter.next();
		return null;
	}

	/**
	 * Load genes from file
	 * @throws IOException
	 */
	private void loadGenes() throws IOException {
		this.genes.clear();
		if(this.genefile == null) {
			throw new IllegalStateException("Name of gene file has not been grabbed from command line.");
		}
		System.out.println("Loading genes from file " + this.genefile + "...");
		FastaSequenceIO fsio = new FastaSequenceIO(new File(this.genefile));
		this.genes= fsio.loadAll();
		System.out.println("Loaded " + this.genes.size() + " genes.");
		System.out.println();
	}

	/**
	 * Load sequences from file
	 * @throws IOException thrown by FastaSequenceIO instance
	 */
	private void loadSequences() throws IOException {
		
		this.sequences.clear();
		this.sequencesByPurpose.clear();
		
		// Load pulldown sequences
		if(this.pulldownseqfile == null) {
			throw new IllegalStateException("Name of pulldown sequence file has not been grabbed from command line.");
		}
		System.out.println("Loading pulldown sequences from file " + this.pulldownseqfile + "...");
		FastaSequenceIO fsio = new FastaSequenceIO(new File(this.pulldownseqfile));
		Collection<Sequence> pulldownsequences = null;
		pulldownsequences = fsio.loadAll();
		System.out.println("Loaded " + pulldownsequences.size() + " pulldown sequences.");
		// Keep track of purpose of these sequences
		for(Sequence seq : pulldownsequences) {
			this.sequencesByPurpose.put(seq, ProbePurpose.PULLDOWN);
		}
		// Add to master collection of sequences
		this.sequences.addAll(pulldownsequences);
		pulldownsequences = null;
		
		// Load depletion sequences for main array
		if(this.mainArrayHasDepletionSequences) {
			if(this.depletionseqfile == null) {
				throw new IllegalStateException("Name of depletion sequence file has not been grabbed from command line.");
			}
			System.out.println("Loading depletion sequences from file " + this.depletionseqfile + "...");
			FastaSequenceIO fsio2 = new FastaSequenceIO(new File(this.depletionseqfile));
			Collection<Sequence> depletionsequences = null;
			depletionsequences = fsio2.loadAll();
			System.out.println("Loaded " + depletionsequences.size() + " depletion sequences.");
			// Keep track of purpose of these sequences
			for(Sequence seq : depletionsequences) {
				this.sequencesByPurpose.put(seq, ProbePurpose.DEPLETE_MAIN_ARRAY);
			}
			// Add to master collection of sequences
			this.sequences.addAll(depletionsequences);
			depletionsequences = null;
		}
			
		// Load mrna depletion sequences
		if(this.makeMrnaDepletionArray) {
			if(this.mrnadepletionfile == null) {
				throw new IllegalStateException("Name of mRNA depletion sequence file has not been grabbed from command line.");
			}
			System.out.println("Loading mRNA depletion sequences from file " + this.mrnadepletionfile + "...");
			FastaSequenceIO fsio3 = new FastaSequenceIO(new File(this.mrnadepletionfile));
			Collection<Sequence> mrnadepletionseqs = null;
			mrnadepletionseqs = fsio3.loadAll();
			System.out.println("Loaded " + mrnadepletionseqs.size() + " mRNA depletion sequences.");
			// Keep track of purpose of these sequences
			for(Sequence seq : mrnadepletionseqs) {
				this.sequencesByPurpose.put(seq, ProbePurpose.DEPLETE_ALL_MRNA);
			}
			// Add to master collection of sequences
			this.sequences.addAll(mrnadepletionseqs);
			mrnadepletionseqs = null;
		}
		
		System.out.println();
	}
	
	/**
	 * Make a qPCR primer in each domain of each pulldown sequence
	 * Add to provided input set and write to a file
	 */
	private void makeAndCombineQpcrPrimers() throws IOException {
		if(this.sequencesBySpeciesGene.isEmpty()) {
			throw new IllegalStateException("Can't make qPCR primers because sequence by gene set is empty.");
		}
		
		System.out.println("Designing qPCR primers for pulldown sequences.");
		
		// Save existing qPCR primers from input file
		FastaSequenceIO fsio = new FastaSequenceIO(new File(this.qpcrfile));
		Collection<Sequence> qPcrPrimers = null;
		qPcrPrimers = fsio.loadAll();
		System.out.println("Read " + qPcrPrimers.size() + " from file " + this.qpcrfile + ".");
		
		File qpcrfasta = new File("AllQpcrPrimers.fa");
		File qpcrtable = new File("DomainQpcrPrimers.out");
		System.out.println("Writing all qPCR primers to AllQpcrPrimers.fa");
		System.out.println("Writing new qPCR primers to DomainQpcrPrimers.out");
		
		// Check if qPCR primers exist from a previous run and if so read them in
		HashMap<String, HashMap<Integer, PrimerPair>> existingQpcrPrimers = new HashMap<String, HashMap<Integer, PrimerPair>>();
		if(qpcrtable.exists()) {
			System.out.println("Warning: reading existing qPCR primers from file "+ qpcrtable + ". If also reading gene probes, will assume both were generated on the same run.");
			FileReader reader = new FileReader(qpcrtable);
			BufferedReader b = new BufferedReader(reader);
			StringParser p = new StringParser();
			int numRead = 0;
			while(b.ready()) {
				String line = b.readLine();
				p.parse(line);
				if(p.asString(0).equals("Gene") || p.getFieldCount() != 5) continue;
				String gene = p.asString(0);
				Integer domain = Integer.valueOf(p.asInt(1));
				String leftPrimer = p.asString(2);
				String rightPrimer = p.asString(3);
				
				// Mostly make up data to create a new primer pair
				String[] primerPairData = new String[9];
				primerPairData[0] = gene + "_" + domain.toString();
				primerPairData[1] = "0";
				primerPairData[2] = leftPrimer;
				primerPairData[3] = "100";
				primerPairData[4] = rightPrimer;
				primerPairData[5] = "60";
				primerPairData[6] = "60";
				primerPairData[7] = "100";
				primerPairData[8] = "no_comment";
	
				PrimerPair primer = new PrimerPair(primerPairData);
				
				HashMap<Integer, PrimerPair> tempMap = new HashMap< Integer,PrimerPair>();
				if(existingQpcrPrimers.containsKey(gene)) tempMap.putAll(existingQpcrPrimers.get(gene));
				tempMap.put(domain, primer);
				numRead++;
				existingQpcrPrimers.put(gene, tempMap);
			}
			System.out.println("Read " + numRead + " previously generated qPCR primers from file.");
		}
		
		
		FileWriter outfasta = new FileWriter(qpcrfasta);
		FileWriter outtable = new FileWriter(qpcrtable);
		outtable.write("Gene\tDomain\tLeftPrimer\tRightPrimer\tRightPrimerRC\n");
		
		// Start fasta file by writing the existing primers
		for(Sequence seq : qPcrPrimers) {
			outfasta.write(">" + seq.getId() + "\n");
			outfasta.write(seq.getSequenceBases() + "\n");
		}
		
		
		// Make new primers for each gene
		for(SpeciesGene speciesgene : this.sequencesBySpeciesGene.keySet()) {
			// Only for pulldown lincRNA genes
			if(speciesgene.getPurpose() != ProbePurpose.PULLDOWN) continue;
			//if(speciesgene.getRnaClass() != RnaClass.LINCRNA) continue;
			//System.out.println("Designing qPCR primers for gene " + speciesgene.getGeneName());
			Sequence longestIsoform = new Sequence("");
			int longestIsoformLength = 0;
			// Identify the longest isoform
			for(Sequence isoform : this.sequencesBySpeciesGene.get(speciesgene)) {
				if(isoform.getLength() > longestIsoformLength) {
					longestIsoformLength = isoform.getLength();
					longestIsoform.setSequenceBases(isoform.getSequenceBases());
				}
			}
			//System.out.println("Longest isoform has length " + longestIsoformLength);
			int domainSize = Double.valueOf(Math.floor((double)longestIsoformLength / (double)NUM_DOMAINS)).intValue();
			//System.out.println("Domain size is " + domainSize);
			for(int i=0; i < NUM_DOMAINS; i++) {
				String domain = longestIsoform.getSequenceBases().substring(i * domainSize, (i+1) * domainSize - 1);
				//Check if we already have a primer for this domain
				if(existingQpcrPrimers.containsKey(speciesgene.getGeneName())) {
					if(existingQpcrPrimers.get(speciesgene.getGeneName()).containsKey(Integer.valueOf(i))) {
						PrimerPair primer = existingQpcrPrimers.get(speciesgene.getGeneName()).get(Integer.valueOf(i));
						String leftPrimer = primer.getLeftPrimer();
						String rightPrimer = primer.getRightPrimer();
						String rightPrimerRC = Sequence.reverseSequence(rightPrimer);
						outtable.write(speciesgene.getGeneName() + "\t" + i + "\t" + leftPrimer + "\t" + rightPrimer + "\t" + rightPrimerRC + "\n");
						String fastaheader = ">Gene:";
						fastaheader += speciesgene.getGeneName();
						fastaheader += "_Domain:";
						fastaheader += i;
						String leftfastaheader = fastaheader + "_LeftPrimer\n";
						String rightfastaheader = fastaheader + "_RightPrimer\n";
						outfasta.write(leftfastaheader);
						outfasta.write(primer.getLeftPrimer() + "\n");
						outfasta.write(rightfastaheader);
						outfasta.write(primer.getRightPrimer() + "\n");						
					}
				} else {
					// Otherwise make a new primer
					PrimerPair primer = designQpcrPrimer(domain);
					if(primer == null) {
						System.out.println("Warning: could not design primer for gene " + speciesgene.getGeneName() + " domain " + i);
					} else {
						//System.out.println("Designed primer for gene " + speciesgene.getGeneName() + " and domain " + i + ":\n" + primer.getPrimerFieldsAsStringForConstructor());
						String leftPrimer = primer.getLeftPrimer();
						String rightPrimer = primer.getRightPrimer();
						String rightPrimerRC = Sequence.reverseSequence(rightPrimer);
						outtable.write(speciesgene.getGeneName() + "\t" + i + "\t" + leftPrimer + "\t" + rightPrimer + "\t" + rightPrimerRC + "\n");
						String fastaheader = ">Gene:";
						fastaheader += speciesgene.getGeneName();
						fastaheader += "_Domain:";
						fastaheader += i;
						String leftfastaheader = fastaheader + "_LeftPrimer\n";
						String rightfastaheader = fastaheader + "_RightPrimer\n";
						outfasta.write(leftfastaheader);
						outfasta.write(primer.getLeftPrimer() + "\n");
						outfasta.write(rightfastaheader);
						outfasta.write(primer.getRightPrimer() + "\n");
					}
				}
			}
		}
		outfasta.close();
		outtable.close();
		System.out.println();
	}

	/**
	 * Use parallel LSF jobs to collapse genes and filter probes, then read back in to geneProbeSet
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void makeFinalGeneProbes() throws IOException, InterruptedException {
	
		
		// Get name for temp directory
		File curr = new File(".");
		File dir = new File(curr.getAbsolutePath() + "/TmpGeneAndIsoformProbes");
		Runtime run = Runtime.getRuntime();
		// Make directory
		System.out.println("Creating directory " + dir);
		boolean madedir = dir.mkdirs();
		if(!madedir) {
			System.err.println("Could not create directory " + dir);
			//System.exit(-1);
		}
		
		// Write isoform probes to one file per gene
		System.out.println("Writing isoform probes to directory " + dir + "...");
		for(SpeciesGene speciesgene : this.sequencesBySpeciesGene.keySet()) {
			// Check if there is already a file of collapsed probes for this gene and if so, continue
			File geneDone = new File("TmpGeneAndIsoformProbes/geneProbes_" + speciesgene.toString());
			if(geneDone.exists()) continue;
			// Check if there is already a file of isoforms for this gene and if so, continue
			File isoformsDone = new File("TmpGeneAndIsoformProbes/isoformProbes_" + speciesgene.toString());
			if(isoformsDone.exists()) continue;
			String isoformFile = "TmpGeneAndIsoformProbes/isoformProbes_" + speciesgene.toString();
			this.writeIsoformProbes(speciesgene, isoformFile);
		}
		
		// Make gene probes by submitting parallel jobs to LSF
		
		System.out.println("Building gene probes in directory " + dir + "...");
		
		// Keep track of job IDs to wait for them to finish
		ArrayList<String> jobIDs = new ArrayList<String>();
		String jobID = null;
		int numSubmitted = 0;
		for(SpeciesGene speciesgene : this.sequencesBySpeciesGene.keySet()) {
			
			// Check if there is already a file of collapsed probes for this gene and if so, continue
			File geneDone = new File("TmpGeneAndIsoformProbes/geneProbes_" + speciesgene.toString());
			if(geneDone.exists()) continue;
			
			// Submit 100 jobs at a time
			if(numSubmitted > 0 && numSubmitted % 100 == 0) {
				Thread.sleep(1000);
			}
			
			String isoformFile = "TmpGeneAndIsoformProbes/isoformProbes_" + speciesgene.toString();
			String geneFile = "TmpGeneAndIsoformProbes/geneProbes_" + speciesgene.toString();
			
			// Parameters to submit job to LSF
			// IMPORTANT: BUILD LOCATION IS HARD CODED BELOW AND VARIES ACCORDING TO ARRAY TYPE
			String argstring = speciesgene.toString() + " " + isoformFile + " " + this.genefile + " " + geneFile + " AllQpcrPrimers.fa";
			jobID = speciesgene.toString() + "_" + Long.valueOf(System.currentTimeMillis()).toString();
			// Use the main method which only does this task
			// For main array
			String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools/build:/seq/lincRNA/Pam/Software/gbtools/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/junit-4.4.jar broad.pda.capture.designer.FilteredGeneProbes " + argstring;
			// For mrna depletion array
			//String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/junit-4.4.jar broad.pda.capture.designer.FilteredGeneProbes " + argstring;
			// For short probes
			//String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/junit-4.4.jar broad.pda.capture.designer.FilteredGeneProbes " + argstring;
			
			// Submit the job
			try {
				// Get more memory if there is more than 1 isoform or more than 4 tiling paths
				if(this.sequencesBySpeciesGene.get(speciesgene).size() > 1 || this.numTilingPaths > 4) {
					int exitCode = PipelineUtils.bsubMediumProcess(Runtime.getRuntime(), jobID , cmmd , "TmpGeneAndIsoformProbes/bsub_output_" + speciesgene.toString());
				} else {
					int exitCode = PipelineUtils.bsubSmallProcess(Runtime.getRuntime(), jobID , cmmd , "TmpGeneAndIsoformProbes/bsub_output_" + speciesgene.toString());
				}
				jobIDs.add(jobID);
				numSubmitted++;
			} catch (InterruptedException e) {
				System.err.println("Caught InterruptedException when trying to submit job " + jobID);
				e.printStackTrace();
			} catch (IOException e) {
				System.err.println("Caught IOException when trying to submit job " + jobID);
				e.printStackTrace();					
			}
			
		}
		
		// Wait for last job to finish
		// DOESN'T SEEM TO BE WORKING FOR SOME REASON
		// FOR NOW, JUST KILL AND RESTART. EXISTING OUTPUT WILL BE READ BACK IN.
		System.out.println("Waiting for last job to finish...");
		PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
		
		System.out.println("Done creating gene probe temp files.\n");
		
		// Read all the primers back in from the files
		System.out.println("Reading gene probes back in from temp files...");
	
		this.geneProbeSet.clear();
		// Read one file per SpeciesGene
		for(SpeciesGene speciesgene : this.sequencesBySpeciesGene.keySet())	 {
			String geneFile = "TmpGeneAndIsoformProbes/geneProbes_" + speciesgene.toString();
			TreeSet<Probe> probes = this.readFromFile(geneFile);
			// Add to geneProbeSet
			this.geneProbeSet.put(speciesgene, probes);
		}
		
		int totalGeneProbes = 0;
		// Count the new probes
		for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) totalGeneProbes += this.geneProbeSet.get(speciesgene).size();
		System.out.println("Read in " + totalGeneProbes + " total gene probes.");
	
	}

	/**
	 * Clear isoform probe set and make initial set of probes against each isoform separately
	 * Tiling scheme:
	 * 120-mers tiled every 15bp across transcript
	 * One extra probe at the end
	 */
	private void makeIsoformProbes() {
	
		System.out.println("Making initial probe set...");
	
		if(this.sequences.isEmpty()) {
			throw new IllegalStateException("Set of sequences is empty.");
		}
	
		System.out.println("Oligo size: " + OLIGO_SIZE);
		System.out.println("Tiling density: " + TILING_DENSITY);
		System.out.println("Number of tiling paths not including extra end probe: " + this.numTilingPaths + "\n");
	
		this.isoformProbeSet.clear();
		int numProbes = 0;
		
		BigInteger td = BigInteger.valueOf(TILING_DENSITY);
		BigInteger ntp = BigInteger.valueOf(this.numTilingPaths);
		
		// Create probes for one sequence at a time
		for(Sequence sequence : this.sequences) {
	
			TreeSet<Probe> seqProbes = new TreeSet<Probe>();
				
			// This will be the set of probes for the sequence, sorted by start position
			seqProbes.clear();
			int sliderPos = 0;
			WindowSlider slider = WindowSlider.getSlider(sequence, OLIGO_SIZE, OLIGO_SIZE - TILING_DENSITY);
			
			while(slider.hasNext()){
				Probe probe = new Probe();
				SequenceRegion region = slider.next();
				Sequence seq = region.getSequence();
				// Reverse complement the sequence
				seq.reverse();
				probe.setSequence(seq);
				probe.setStartPos(sliderPos);
				probe.setEndPos(sliderPos + OLIGO_SIZE - 1);
				// Calculate the tiling path of this probe
				BigInteger tp = BigInteger.valueOf(sliderPos).mod(td.multiply(ntp));
				probe.setTilingPath(tp.intValue());
				// EvenOdd is the number of the probe within the tiling path mod 2
				probe.setEvenOdd(((sliderPos - tp.intValue()) / OLIGO_SIZE) % 2);
				float domainSize = sequence.getLength() / NUM_DOMAINS;
				// Domain is the floor of the start position divided by domain size
				probe.setDomain((Double.valueOf(Math.floor(sliderPos / domainSize))).intValue());
				boolean added = seqProbes.add(probe);
				if(added) {
						sliderPos += TILING_DENSITY;
				}
			}
	
			// If end of sequence was not captured, add a final tiling path with one probe
			if(sliderPos != sequence.getLength()) {
				Probe probe = new Probe();
				// If sequence is shorter than oligo size, use entire sequence
				SequenceRegion region = sequence.getRegion(BigInteger.ZERO.max(BigInteger.valueOf(sequence.getEnd() - OLIGO_SIZE + 1)).intValue(), sequence.getEnd());
				Sequence seq = region.getSequence();
				// Reverse complement the sequence
				seq.reverse();
				probe.setSequence(seq);
				probe.setStartPos(region.getRegionStart());
				probe.setEndPos(region.getRegionEnd());
				// This is the extra tiling path for miscellaneous probes
				probe.setTilingPath(-1);
				boolean added = seqProbes.add(probe);
			}
			
			// Add sequence and its probes to the full probe set
			this.isoformProbeSet.put(sequence, seqProbes);
			// Tell garbage collector to free up memory right away
			seqProbes = null;
			//System.out.println("Designed " + this.probeSet.get(sequence).size() + " probes for sequence " + sequence.getId());
			numProbes += this.isoformProbeSet.get(sequence).size();
	
			
		}
	
		System.out.println("Created " + numProbes + " initial probes.\n");
	
	}

	/**
	 * Read probes from a file into a TreeSet
	 * @param file
	 * @return the TreeSet of probes
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	protected TreeSet<Probe> readFromFile(String file) throws FileNotFoundException, IOException {
		TreeSet<Probe> probes = new TreeSet<Probe>();
		FileReader reader = new FileReader(file);
		BufferedReader b = new BufferedReader(reader);
		while(b.ready()) {
			String line = b.readLine();
			Probe probe = new Probe();
			// Parse the data on the line and read into the probe
			probe.readFromString(line);
			probes.add(probe);
		}
		return probes;
	}

	/**
	 * Write all primers, assigned or not, to a file
	 * @throws IOException
	 */
	private void writeAllPrimers() throws IOException {
		System.out.println("Writing all primers to file " + this.outprimersfile);
		if(this.depletionPrimers.isEmpty()) {
			throw new IllegalStateException("Can't write primers to file because depletion primer set is empty.");
		}
		if(this.pulldownPrimers.isEmpty()) {
			throw new IllegalStateException("Can't write primers to file because pulldown primer set is empty.");
		}
		File outfile = new File(this.outprimersfile);
		FileWriter writer = new FileWriter(outfile);
		int numPrimers = 0;
		for(PrimerPair p : this.depletionPrimers.keySet()) {
			writer.write("Depletion_outer\t" + p.getPrimerFieldsAsStringForConstructor() + "\n");
			numPrimers++;
			for(PrimerPair q : this.depletionPrimers.get(p)) {
				writer.write("Depletion_inner\t" + q.getPrimerFieldsAsStringForConstructor() + "\n");
				numPrimers++;
			}
		}
		for(PrimerPair p : this.pulldownPrimers.keySet()) {
			writer.write("Pulldown_outer\t" + p.getPrimerFieldsAsStringForConstructor() + "\n");
			numPrimers++;
			for(PrimerPair q : this.pulldownPrimers.get(p).keySet()) {
				writer.write("Pulldown_middle\t" + q.getPrimerFieldsAsStringForConstructor() + "\n");
				numPrimers++;
				for(PrimerPair r : this.pulldownPrimers.get(p).get(q)) {
					writer.write("Pulldown_inner\t" + r.getPrimerFieldsAsStringForConstructor() + "\n");
					numPrimers++;
				}
			}
		}
		writer.close();
		System.out.println("Wrote " + numPrimers + " to file.\n");
	}

	/**
	 * Write all probe sequences to file
	 * @throws IOException
	 */
	private void writeAllProbes() throws IOException {
		
		// Write flat file of probe data
		String outfile = this.outprobesfile + ".out";
		System.out.println("Writing probe information to file " + outfile);
		FileWriter writer = new FileWriter(outfile);
		writer.write("Species\tGene\tClass\tPurpose\tReferenceIsoformLength\tStartPosition\tEndPosition\tTilingPath\tEvenOdd\tDomain\tContainsQpcrPrimer\tProbeSequence\n");
		for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
			// Get length of longest isoform to print
			int longestIsoform = 0;
			for(Sequence seq : this.sequencesBySpeciesGene.get(speciesgene)) {
				if(seq.getLength() > longestIsoform) longestIsoform = seq.getLength();
			}
			for(Probe probe : this.geneProbeSet.get(speciesgene)) {
				writer.write(speciesgene.getSpeciesName() + "\t" + speciesgene.getGeneName() + "\t" + speciesgene.getRnaClass() + "\t" + speciesgene.getPurpose() + "\t"+ longestIsoform + "\t" + probe + "\n");
			}
		}
		writer.close();
	
		// Write fasta file of probe sequences
		String fastafile = this.outprobesfile + ".fa";
		System.out.println("Writing probe sequences to file " + fastafile +"\n");
		FileWriter fastawriter = new FileWriter(fastafile);
		// Put all the probe information in the fasta header
		for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
			for(Probe probe : this.geneProbeSet.get(speciesgene)) {
				String header = ">";
				header += speciesgene.getSpeciesName();
				header += "_";
				header += speciesgene.getGeneName();
				header += "_";
				header += speciesgene.getPurpose().toString();
				header += "_";
				header += speciesgene.getRnaClass().toString();
				header += "_";
				header += "Start:";
				header += Integer.valueOf(probe.getStartPos()).toString();
				header += "_";
				header += "End:";
				header += Integer.valueOf(probe.getEndPos()).toString();		
				header += "_";
				header += "TilingPath:";
				header += Integer.valueOf(probe.getTilingPath()).toString();
				header += "_";
				header += "EvenOdd:";
				header += Integer.valueOf(probe.getEvenOdd()).toString();
				header += "_";
				header += "Domain:";
				header += Integer.valueOf(probe.getDomain()).toString();
				header += "_";
				header += "ContainsQpcrPrimer:";
				header += Integer.valueOf(probe.getContainsQpcrPrimer()).toString();
				header += "\n";
				fastawriter.write(header);
				fastawriter.write(probe.getSequenceBases() + "\n");
			}
		}
		fastawriter.close();
		
	}

	/**
	 * Write full design to file
	 * @throws IOException
	 */
	private void writeFullDesign() throws IOException {
		if(this.outdesignfile == null) {
			throw new IllegalStateException("Output design file has not been specified.");
		}
		
		String tablefile = this.outdesignfile + ".out";
		String fastafile = this.outdesignfile + ".fa";
		
		// Write design in table format
		System.out.println("Writing full design to file " + tablefile);
		FileWriter writer = new FileWriter(tablefile);
		writer.write("Species\tGene\tRnaClass\tPurpose\tProbeStart\tProbeEnd\tTilingPath\tEvenOdd\tDomain\tContainsQpcrPrimer\tSpeciesClassLeftPrimer\tSpeciesClassRightPrimer\tGeneLeftPrimer\tGeneRightPrimer\tTilingPathLeftPrimer\tEvenOddLeftPrimer\tQpcrRightPrimer\tDomainRightPrimer\tProbeOligo\tFullProbeSequence\n");
	
		// Write design in fasta format
		System.out.println("Writing full probe sequences to file " + fastafile);
		FileWriter fastawriter = new FileWriter(fastafile);
		
		for(SpeciesGene speciesgene : this.geneProbeSet.keySet()) {
			// Get gene attributes
			Species species = speciesgene.getSpecies();
			String gene = speciesgene.getGeneName();
			ProbePurpose purpose = speciesgene.getPurpose();
			RnaClass rnaclass = speciesgene.getRnaClass();
			
			// Depletion probes
			if(purpose == ProbePurpose.DEPLETE_MAIN_ARRAY || purpose == ProbePurpose.DEPLETE_ALL_MRNA) {
								
				for(Probe probe : this.geneProbeSet.get(speciesgene)) {
					
					// Get probe attributes
					int startpos = probe.getStartPos();
					int endpos = probe.getEndPos();
					int tilingpath = -1;
					int evenodd = -1;
					int domain = -1;
					int qpcr = -1;
					String probeseq = probe.getSequenceBases();
					
					// Get primer sequences and full primer sequence
					String leftOuterPrimer = this.depletionOuterLeftPrimers.get(speciesgene.getSpecies()).get(speciesgene.getRnaClass()).get(probe);
					String rightOuterPrimer = this.depletionOuterRightPrimers.get(speciesgene.getSpecies()).get(speciesgene.getRnaClass()).get(probe);
					String leftInnerPrimer = this.depletionInnerLeftPrimers.get(speciesgene.getSpecies()).get(speciesgene.getRnaClass()).get(probe);
					String rightInnerPrimer = this.depletionInnerRightPrimers.get(speciesgene.getSpecies()).get(speciesgene.getRnaClass()).get(probe);
					
					String leftUniqueTail = leftInnerPrimer.substring(10);
					String leftFullPrimer = leftOuterPrimer + leftUniqueTail;
					String rightUniqueTail = rightInnerPrimer.substring(10);
					String rightFullPrimer = rightOuterPrimer + rightUniqueTail;
					
					// RC the right primers
					String rightfullRC = Sequence.reverseSequence(rightFullPrimer);
			
					String line = species.toString() + "\t"; // species
					line += gene + "\t"; // gene
					line += rnaclass.toString() + "\t"; // RNA class
					line += purpose.toString() + "\t"; // purpose
					line += Integer.valueOf(startpos).toString() + "\t"; // probe start
					line += Integer.valueOf(endpos).toString() + "\t"; // probe end
					line += Integer.valueOf(tilingpath).toString() + "\t"; // tiling path
					line += Integer.valueOf(evenodd).toString() + "\t"; // even odd
					line += Integer.valueOf(domain).toString() + "\t"; // domain
					line += Integer.valueOf(qpcr).toString() + "\t"; // qpcr primer site
					line += leftOuterPrimer + "\t"; // species class left primer
					line += rightOuterPrimer + "\t"; // species class right primer
					line += leftInnerPrimer + "\t"; // gene left primer
					line += rightInnerPrimer + "\t"; // gene right primer
					line += "NA\t"; // tiling path left primer
					line += "NA\t"; // even odd left primer
					line += "NA\t"; // qpcr right primer
					line += "NA\t"; // domain right primer
					line += probeseq + "\t"; // probe oligo
					line += leftFullPrimer + probeseq + rightfullRC + "\n"; // full probe sequence
					
					// write to the table
					writer.write(line);
	
					// make the fasta information
					String header = ">";
					header += speciesgene.getSpeciesName();
					header += "_";
					header += speciesgene.getGeneName();
					header += "_";
					header += speciesgene.getPurpose().toString();
					header += "_";
					header += speciesgene.getRnaClass().toString();
					header += "_";
					header += "ProbeStart:";
					header += Integer.valueOf(probe.getStartPos()).toString();
					header += "_";
					header += "ProbeEnd:";
					header += Integer.valueOf(probe.getEndPos()).toString();		
					header += "_";
					header += "TilingPath:";
					header += Integer.valueOf(probe.getTilingPath()).toString();
					header += "_";
					header += "EvenOdd:";
					header += Integer.valueOf(probe.getEvenOdd()).toString();
					header += "_";
					header += "Domain:";
					header += Integer.valueOf(probe.getDomain()).toString();
					header += "_";
					header += "ContainsQpcrPrimer:";
					header += Integer.valueOf(probe.getContainsQpcrPrimer()).toString();
					header += "\n";
					
					// write to fasta file
					fastawriter.write(header);
					fastawriter.write(leftFullPrimer + probe.getSequenceBases() + rightfullRC + "\n");
	
					
				}
			}
			
			// Pulldown probes
			if(purpose == ProbePurpose.PULLDOWN) {
								
				for(Probe probe : this.geneProbeSet.get(speciesgene)) {
					
					// Get probe attributes
					int startpos = probe.getStartPos();
					int endpos = probe.getEndPos();
					int tilingpath = probe.getTilingPath();
					int evenodd = probe.getEvenOdd();
					int domain = probe.getDomain();
					int qpcr = probe.getContainsQpcrPrimer();
					
					// Get primer sequences and trim to specific primer length
					String leftOuterPrimer = this.pulldownOuterLeftPrimers.get(probe);
					String rightOuterPrimer = this.pulldownOuterRightPrimers.get(probe);
					String leftMiddlePrimer = this.pulldownMiddleLeftPrimers.get(probe);
					String rightMiddlePrimer = this.pulldownMiddleRightPrimers.get(probe);
					String leftInnerPrimer = this.pulldownInnerLeftPrimers.get(probe);
					String rightInnerPrimer = this.pulldownInnerRightPrimers.get(probe);
					
					String leftUniqueOuterTail = leftOuterPrimer.substring(0,5);
					String leftUniqueInnerTail = leftInnerPrimer.substring(10);
					
					String rightUniqueOuterTail = rightOuterPrimer.substring(0,5);
					String rightUniqueInnerTail = rightInnerPrimer.substring(10);
					
					String leftFullPrimer = leftUniqueOuterTail + leftMiddlePrimer + leftUniqueInnerTail;
					String rightFullPrimer = rightUniqueOuterTail + rightMiddlePrimer + rightUniqueInnerTail;
					
					// RC the right primers
					String rightfullRC = Sequence.reverseSequence(rightFullPrimer);
			
					String line = species.toString() + "\t"; // species
					line += gene + "\t"; // gene
					line += rnaclass.toString() + "\t"; // RNA class
					line += purpose.toString() + "\t"; // purpose
					line += Integer.valueOf(startpos).toString() + "\t"; // probe start
					line += Integer.valueOf(endpos).toString() + "\t"; // probe end
					line += Integer.valueOf(tilingpath).toString() + "\t"; // tiling path
					line += Integer.valueOf(evenodd).toString() + "\t"; // even odd
					line += Integer.valueOf(domain).toString() + "\t"; // domain
					line += Integer.valueOf(qpcr).toString() + "\t"; // qpcr primer sit
					line += "NA\t"; // species class left primer
					line += "NA\t"; // species class right primer
					line += leftOuterPrimer + "\t"; // gene left primer
					line += rightOuterPrimer + "\t"; // gene right primer
					line += leftMiddlePrimer + "\t"; // tiling path left primer
					line += leftInnerPrimer + "\t"; // even odd left primer
					line += rightMiddlePrimer + "\t"; // qpcr right primer
					line += rightInnerPrimer + "\t"; // domain right primer
					line += probe.getSequenceBases() + "\t"; // probe oligo
					line += leftFullPrimer + probe.getSequenceBases() + rightfullRC + "\n"; // full probe sequence
					
					// write to table file
					writer.write(line);
					
					// make the fasta information
					String header = ">";
					header += speciesgene.getSpeciesName();
					header += "_";
					header += speciesgene.getGeneName();
					header += "_";
					header += speciesgene.getPurpose().toString();
					header += "_";
					header += speciesgene.getRnaClass().toString();
					header += "_";
					header += "ProbeStart:";
					header += Integer.valueOf(probe.getStartPos()).toString();
					header += "_";
					header += "ProbeEnd:";
					header += Integer.valueOf(probe.getEndPos()).toString();		
					header += "_";
					header += "TilingPath:";
					header += Integer.valueOf(probe.getTilingPath()).toString();
					header += "_";
					header += "EvenOdd:";
					header += Integer.valueOf(probe.getEvenOdd()).toString();
					header += "_";
					header += "Domain:";
					header += Integer.valueOf(probe.getDomain()).toString();
					header += "_";
					header += "ContainsQpcrPrimer:";
					header += Integer.valueOf(probe.getContainsQpcrPrimer()).toString();
					header += "\n";
					
					// write to fasta file
					fastawriter.write(header);
					fastawriter.write(leftFullPrimer + probe.getSequenceBases() + rightfullRC + "\n");
	
					
				}
			}
			
		}
		fastawriter.close();
		writer.close();		
	}

	/**
	 * Write all the initial isoform probes for one gene to a file
	 * This is just for reading in by a FilteredGeneProbes object
	 * First field is the entire isoform sequence, then probe data
	 * @param speciesgene
	 * @param outFile
	 * @throws IOException
	 */
	private void writeIsoformProbes(SpeciesGene speciesgene, String outFile) throws IOException {
		FileWriter writer = new FileWriter(outFile);
		for(Sequence isoform : this.sequencesBySpeciesGene.get(speciesgene)) {
			for(Probe probe : this.isoformProbeSet.get(isoform)) {
				writer.write(isoform.getSequenceBases() + " " + probe + "\n");
			}
		}
		writer.close();
	}

	/**
	 * Write a set of probes to a file
	 * @param probes
	 * @param outFile
	 * @throws IOException
	 */
	protected static void writeProbes(TreeSet<Probe> probes, String outFile) throws IOException {
		FileWriter writer = new FileWriter(outFile);
		for(Probe probe : probes) {
			// Write the string representation of the probe that can be read back in by Probe object
			writer.write(probe + "\n");
		}
		writer.close();
	}

	/**
	 * Check if both primers in a pair have perfect matches to a probe
	 * @param leftPrimer left primer sequence
	 * @param rightPrimer right primer sequence
	 * @param kmers the hash set of probe kmers
	 * @param length length of match
	 * @return
	 */
	private boolean hasDoubleMatch(ArrayDesigner a, String leftPrimer, String rightPrimer, HashSet<StringPair> kmers) {
		String leftthreeprime = leftPrimer.substring(leftPrimer.length() - this.doublePrimerMatchToProbe);
		String rightthreeprime = rightPrimer.substring(rightPrimer.length() - this.doublePrimerMatchToProbe);
		String leftthreeprimeRC = Sequence.reverseSequence(leftthreeprime);
		String rightthreeprimeRC = Sequence.reverseSequence(rightthreeprime);
		
		StringPair lrr = a.new StringPair(leftthreeprime,rightthreeprimeRC);
		if(kmers.contains(lrr)) return true;
		
		StringPair rrl = a.new StringPair(rightthreeprime,leftthreeprimeRC);
		if(kmers.contains(rrl)) return true;
	
		return false;
	}

	/**
	 * Check if a primer pair has a single perfect match in hash table
	 * @param primer the primer pair
	 * @param kmers the hash set of kmers
	 * @param length length of match
	 * @return
	 */
	private boolean hasSingleMatch(String leftPrimer, String rightPrimer, HashMap<String,Integer> kmers) {
		String leftthreeprime = leftPrimer.substring(leftPrimer.length() - this.singlePrimerMatchToProbe);
		if(kmers.containsKey(leftthreeprime)) {
			if(kmers.get(leftthreeprime).intValue() >= this.minKmerOccurrencesInProbes) return true;
		}
		String rightthreeprime = rightPrimer.substring(rightPrimer.length() - this.singlePrimerMatchToProbe);
		if(kmers.containsKey(rightthreeprime)) {
			if(kmers.get(rightthreeprime).intValue() >= this.minKmerOccurrencesInProbes) return true;
		}
		String leftthreeprimeRC = Sequence.reverseSequence(leftthreeprime);
		if(kmers.containsKey(leftthreeprimeRC)) {
			if(kmers.get(leftthreeprimeRC).intValue() >= this.minKmerOccurrencesInProbes) return true;
		}
		String rightthreeprimeRC = Sequence.reverseSequence(rightthreeprime);
		if(kmers.containsKey(rightthreeprimeRC)) {
			if(kmers.get(rightthreeprimeRC).intValue() >= this.minKmerOccurrencesInProbes) return true;
		}
		return false;
	}

	/**
	 * Check if a left primer has double match with any probe when paired with any right primer in set
	 * @param a an ArrayDesigner object
	 * @param leftPrimer the left primers
	 * @param rightPrimers the set of right primers
	 * @param kmers the probe kmers
	 * @return whether there is a double match to a probe
	 */
	private static boolean leftPrimerHasDoubleMatchWithRightSet(ArrayDesigner a, String leftPrimer, Collection<String> rightPrimers, HashSet<StringPair> kmers) {
		for(String rightPrimer : rightPrimers) {
			if(a.hasDoubleMatch(a, leftPrimer, rightPrimer, kmers)) return true;
		}
		return false;
	}

	/**
	 * Check if a left primer stops amplification of a probe in a set
	 * @param primer the primer
	 * @param probes the probe set
	 * @param length the length of match
	 * @return
	 */
	private static boolean leftPrimerHasMatchInProbeSet(String primer, Collection<Probe> probes, int length) {
		String threeprime = primer.substring(primer.length() - length);
		for(Probe probe : probes) {
			if(probe.getSequenceBases().contains(threeprime)) return true;
		}
		return false;
	}

	/**
	 * Check if a primer cross primes a primer that has already been kept
	 * @param primer the primer to check against existing primers
	 * @return whether a match is found among existing primers
	 */
	private boolean primerCrossPrimesExistingPrimer(String primer) {
		if(this.keptPrimerKmers.isEmpty()) return false;
		String threeprime = primer.substring(primer.length() - this.betweenGenePrimerDimerMatchLength);
		String threeprimeRC = Sequence.reverseSequence(threeprime);
		if(this.keptPrimerKmers.contains(threeprime)) return true;
		if(this.keptPrimerKmers.contains(threeprimeRC)) return true;
		return false;
	}

	/**
	 * Check if a primer forms dimer with another set of primers
	 * @param primer the primer
	 * @param otherPrimers the other primer
	 * @param length minimum match length
	 * @return
	 */
	private static boolean primerCrossPrimesInSet(String primer, Collection<String> otherPrimers, int length) {
		String threeprime = primer.substring(primer.length() - length);
		String threeprimeRC = Sequence.reverseSequence(threeprime);
		for(String otherPrimer : otherPrimers) {
			if(otherPrimer.contains(threeprime) || otherPrimer.contains(threeprimeRC)) return true;
		}
		return false;
	}

	/**
	 * Check if a right primer has double match with any probe when paired with any left primer in set
	 * @param a an ArrayDesigner object
	 * @param leftPrimers the left primers
	 * @param rightPrimer the right primer
	 * @param kmers the probe kmers
	 * @return whether there is a double match to a probe
	 */
	private static boolean rightPrimerHasDoubleMatchWithLeftSet(ArrayDesigner a, Collection<String> leftPrimers, String rightPrimer, HashSet<StringPair> kmers) {
		for(String leftPrimer : leftPrimers) {
			if(a.hasDoubleMatch(a, leftPrimer, rightPrimer, kmers)) return true;
		}
		return false;
	}

	/**
	 * Check if a right primer stops amplification of a probe in a set
	 * @param primer the primer
	 * @param probes the probe set
	 * @param length the length of match
	 * @return
	 */
	private static boolean rightPrimerHasMatchInProbeSet(String primer, Collection<Probe> probes, int length) {
		String threeprime = primer.substring(primer.length() - length);
		String threeprimeRC = Sequence.reverseSequence(threeprime);
		for(Probe probe : probes) {
			if(probe.getSequenceBases().contains(threeprimeRC)) return true;
		}
		return false;
	}

	/**
	 * Main method designs array(s) and outputs file of oligos and primers.
	 * If command arguments are not correct, prints help menu with arguments.
	 * @param args The command line arguments which are printed as a help menu if entered incorrectly
	 * @throws IOException 
	 * @throws SearchException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, SearchException, InterruptedException {
		
		System.out.println("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		System.out.println("        ARRAY DESIGNER");
		System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		
		ArrayDesigner arrayDesigner = new ArrayDesigner();
	
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Check if primer files are already available
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		File pulldownPrimers = new File("ThirdPrimers_pulldown.out");
		File depletionPrimers = new File("SecondPrimers_depletion.out");
		if(pulldownPrimers.exists()) {
			System.out.println("Warning: reading pulldown primers from file ThirdPrimers_pulldown.out.");
		}
		
		if(arrayDesigner.mainArrayHasDepletionSequences && depletionPrimers.exists()) {
			System.out.println("Warning: reading depletion primers from file SecondPrimers_depletion.out.");
		}
		
		
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Parse command line and load input data
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		// Set argument values from command line
		arrayDesigner.loadCommandArgs(args);		
		
		System.out.println("\nPrimer filtering parameters:");
		System.out.println("Single primer match to any probe on array: " + arrayDesigner.singlePrimerMatchToProbe);
		System.out.println("Max allowable probes with " + arrayDesigner.singlePrimerMatchToProbe + "-mer: " + arrayDesigner.minKmerOccurrencesInProbes);
		System.out.println("Double primer match to any probe on array: " + arrayDesigner.doublePrimerMatchToProbe);
		System.out.println("Single primer match to probe within gene: " + arrayDesigner.singlePrimerMatchToProbeWithinGene);
		System.out.println("Match length for within gene primer dimer: " + arrayDesigner.withinGenePrimerDimerMatchLength);
		System.out.println("Match length for between gene primer dimer: " + arrayDesigner.betweenGenePrimerDimerMatchLength + "\n");
		
		// Load sequences
		arrayDesigner.loadSequences();		
		// Load genes
		arrayDesigner.loadGenes();
	
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Classify transcripts by gene and species
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
		arrayDesigner.loadSequencesBySpeciesGene();			
			
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Make qPCR primers for pulldown genes and add to qPCR primer input, then write to new file
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
		arrayDesigner.makeAndCombineQpcrPrimers();					
		
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Create initial tiling paths
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		arrayDesigner.makeIsoformProbes();
		
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Get final set of filtered, unique probes
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		// For each gene:
		// Write the initial isoform probes to a file
		// Filter repetetive probes
		// Collapse to genes and eliminate redundancy
		// Remove probes that cross-hybridize in transcriptome
		// Identify probes containing a qPCR primer
		// Assign probes to other paths (even/odd, domains)
		// Read data back in
		arrayDesigner.makeFinalGeneProbes();
		
		// Write filtered probes to file
		if(arrayDesigner.outprobesfile != null) arrayDesigner.writeAllProbes();
		
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Design PCR tails, assign to probes and remove cross-amplifying primers
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		arrayDesigner.designPrimers();
		if(arrayDesigner.outprimersfile != null) arrayDesigner.writeAllPrimers();
		
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Write output
		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		arrayDesigner.writeFullDesign();		
		System.out.println();
		System.out.println("All done.");
		
		
	}

	/**
	 * The RNA classes of interest
	 * NOTE: the number of primers designed depends on the number of RNA classes. Comment some classes out to save time if applicable.
	 * @author prussell
	 */
	private enum RnaClass {
		//MRNA,
		LINCRNA;
		//SNRNA,
		//SNORNA,
		//RRNA,
		//RRNA_PRECURSOR;
		//MICRORNA,
		//GUIDERNA,
		//ANTISENSE_RNA,
		//RNASE_P_RNA,
		//TELOMERASE_RNA,
		//SMALL_CYTOPLASMIC_RNA,
		//RNASE_MRP_RNA,
		//MISC_SEQUENCE,
		//NOT_SET;
		
		/**
		 * Get the number of RNA classes
		 * @return the number of classes
		 */
		public static int numberOfClasses() {
			return values().length;
		}
		
		/**
		 * Parse from a string
		 * @param classname
		 * @return the RnaClass corresponding to the string name
		 */
		public static RnaClass getRnaClass(String classname) {
			//if(classname.equals("mRNA")) return MRNA;
			if(classname.equals("lincRNA")) return LINCRNA;
			//if(classname.equals("snRNA")) return SNRNA;
			//if(classname.equals("snoRNA")) return SNORNA;
			//if(classname.equals("microRNA")) return MICRORNA;
			//if(classname.equals("guideRNA")) return GUIDERNA;
			//if(classname.equals("antisenseRNA")) return ANTISENSE_RNA;
			//if(classname.equals("RNasePRNA")) return RNASE_P_RNA;
			//if(classname.equals("telomeraseRNA")) return TELOMERASE_RNA;
			//if(classname.equals("smallCytoplasmicRNA")) return SMALL_CYTOPLASMIC_RNA;
			//if(classname.equals("RNaseMRPRNA")) return RNASE_MRP_RNA;
			//if(classname.equals("rRNA")) return RRNA;
			//if(classname.equals("rRNAprecursor")) return RRNA_PRECURSOR;
			throw new IllegalArgumentException("Rna class name " + classname + " not recognized.");
		}
		
		@Override
		public String toString() {
			switch(this) {
			//case MRNA: return "mRNA";
			case LINCRNA: return "lincRNA";
			//case SNRNA: return "snRNA";
			//case SNORNA: return "snoRNA";
			//case RRNA: return "rRNA";
			//case RRNA_PRECURSOR: return "rRNAprecursor";
			//case MICRORNA: return "microRNA";
			//case GUIDERNA: return "guideRNA";
			//case ANTISENSE_RNA: return "antisenseRNA";
			//case RNASE_P_RNA: return "RNasePRNA";
			//case TELOMERASE_RNA: return "telomeraseRNA";
			//case SMALL_CYTOPLASMIC_RNA: return "smallCytoplasmicRNA";
			//case RNASE_MRP_RNA: return "RNaseMRPRNA";
			//case MISC_SEQUENCE: return "misc_sequence";
			//case NOT_SET: return "NotSet";
			}
			throw new IllegalArgumentException("Unknown RNA class.");
		}
	}
	
	/**
	 * The species of interest
	 * NOTE: the number of primers designed depends on the number of species. Comment some species out to save time if applicable.
	 * @author prussell
	 */
	private enum Species {
		MOUSE;
		//HUMAN,
		//DOG,
		//COW,
		//RAT,
		//OPOSSUM,
		//PIG,
		//MACAQUE;
		
		/**
		 * Parse string to species object
		 * @param name the string name
		 * @return the corresponding Species
		 */
		public static Species getSpecies(String name) {
			if(name.equals("Mus musculus") || name.equals("Mouse") || name.equals("mouse")) return MOUSE;
			//if(name.equals("Bos taurus") || name.equals("Cow")) return COW;
			//if(name.equals("Canis lupus familiaris") || name.equals("Dog") || name.equals("Canis lupus") || name.equals("Canis familiaris")) return DOG;
			//if(name.equals("Macaca mulatta") || name.equals("Macaque")) return MACAQUE;
			//if(name.equals("Homo sapiens") || name.equals("Human") || name.equals("human")) return HUMAN;
			//if(name.equals("Monodelphis domestica") || name.equals("Opossum")) return OPOSSUM;
			//if(name.equals("Rattus norvegicus") || name.equals("Rattus rattus") || name.equals("Rat")) return RAT;
			//if(name.equals("Sus scrofa") || name.equals("Pig") || name.equals("pig")) return PIG;
			throw new IllegalArgumentException("Species name " + name + " not recognized.");
		}
		
		/**
		 * Get the number of species
		 * @return the number of species
		 */
		public static int numberOfSpecies() {
			return values().length;
		}

		
		
		/**
		 * Get string representation of species name
		 */
		@Override
		public String toString() {
			switch(this) {
			case MOUSE: return "Mouse";
			//case HUMAN: return "Human";
			//case DOG: return "Dog";
			//case COW: return "Cow";
			//case RAT: return "Rat";
			//case OPOSSUM: return "Opossum";
			//case MACAQUE: return "Macaque";
			//case PIG: return "Pig";
			}
			throw new IllegalArgumentException("Unknown species.");
		}
				
	}
	
	/**
	 * The experimental goal for a probe
	 * @author prussell
	 */
	protected enum ProbePurpose {
		PULLDOWN,
		DEPLETE_MAIN_ARRAY,
		DEPLETE_ALL_MRNA,
		NOT_SET;
		
		@Override
		public String toString() {
			switch(this) {
			case PULLDOWN: return "Pulldown";
			case DEPLETE_MAIN_ARRAY: return "Depletion(MainArray)";
			case DEPLETE_ALL_MRNA: return "Depletion(mRNA)";
			case NOT_SET: return "NotSet";
			}
			throw new IllegalArgumentException("Unknown purpose.");
		}
		
		/**
		 * Parse from a string
		 * @param purposename Probe purpose
		 * @return the RnaClass corresponding to the string name
		 */
		public static ProbePurpose getProbePurpose(String purposename) {
			if(purposename.equals("Pulldown")) return PULLDOWN;
			if(purposename.equals("Depletion(MainArray)")) return DEPLETE_MAIN_ARRAY;
			if(purposename.equals("Depletion(mRNA)")) return DEPLETE_ALL_MRNA;
			if(purposename.equals("NotSet")) return NOT_SET;
			throw new IllegalArgumentException("Probe purpose name " + purposename + " not recognized.");
		}
	
	}
	
	/**
	 * A Probe object consists of a Sequence object (the probe sequence) plus coordinates on the parent transcript and other path info.
	 * Probe is naturally tied to a parent Sequence object.
	 * Note: this class has a natural ordering that is inconsistent with equals.
	 * The ordering is by start coordinate along the implied parent sequence.
	 * @author prussell
	 */
	protected final class Probe implements Comparable<Probe> {
		
		public Probe() {
			this.sequence = new Sequence("");
			this.startPos = -1;
			this.endPos = -1;
			this.tilingPath = -1;
			this.evenOdd = -1;
			this.domain = -1;
			this.qpcr = 0;
		}
		
		
		/**
		 * Ordering by start position then end position then alphabetical
		 */
		@Override
		public int compareTo(Probe p) {
			// if this starts earlier return less than
			if(this.startPos < p.startPos) {
				return -1;
			}
			// if this starts later return greater than
			if(this.startPos > p.startPos) {
				return 1;
			}
			// if this ends earlier return less than
			if(this.endPos < p.endPos) {
				return -1;
			}
			// if this ends later return greater than
			if(this.endPos > p.endPos) {
				return 1;
			}
			// alphabetical
				return this.getSequenceBases().compareTo(p.getSequenceBases());
			
			
		}
		
		/**
		 * Get sequence of probe
		 * @return the Sequence object
		 */
		public Sequence getSequence() {
			return this.sequence;
		}
		
		/**
		 * Get sequence bases
		 * @return The sequence as String
		 */
		public String getSequenceBases() {
			return this.sequence.getSequenceBases();
		}
		
		/**
		 * Parse the string representation of a probe to a new Probe object
		 * @param data a whitespace-separated line of data of the form start position, end position, tiling path, even/odd, domain, qpcr, sequence bases
		 */
		public void readFromString(String data) {
			StringParser p = new StringParser();
			p.parse(data);
			if(p.getFieldCount() != 7) {
				throw new IllegalArgumentException("String must be a whitespace-separated line of data of the form start position, end position, tiling path, even/odd, domain, qpcr, sequence bases");
			}
			Sequence seq = new Sequence("");
			seq.setSequenceBases(p.asString(6));
			this.setSequence(seq);
			this.setStartPos(Integer.parseInt(p.asString(0)));
			this.setEndPos(Integer.parseInt(p.asString(1)));
			this.setTilingPath(Integer.parseInt(p.asString(2)));
			this.setEvenOdd(Integer.parseInt(p.asString(3)));
			this.setDomain(Integer.parseInt(p.asString(4)));
			this.setQpcrPrimer(Integer.parseInt(p.asString(5)));
		}
		
		/**
		 * String representation
		 * @return a tab separated line of data: start position, end position, tiling path, even/odd, domain, qpcr, sequence bases
		 */
		@Override
		public String toString() {
			return Integer.valueOf(this.startPos).toString() + "\t" + Integer.valueOf(this.endPos).toString()+ "\t" + Integer.valueOf(this.tilingPath).toString()+ "\t" + Integer.valueOf(this.evenOdd).toString()+ "\t" + Integer.valueOf(this.domain).toString()+ "\t" + Integer.valueOf(this.qpcr).toString() + "\t" + this.getSequenceBases();
		}
		
		/**
		 * Get 0-based start position of probe along parent sequence
		 * @return the start position
		 */
		public int getStartPos() {
			return this.startPos;
		}
		
		/**
		 * Get 0-based end position of probe along parent sequence
		 * @return the end position
		 */
		public int getEndPos() {
			return this.endPos;
		}
		
		/**
		 * Get tiling path number for probe
		 * @return the tiling path number
		 */
		public int getTilingPath() {
			return this.tilingPath;
		}
		
		/**
		 * Get even or odd path with respect to its tiling path
		 * @return 0 or 1 or -1 if has not been set
		 */
		public int getEvenOdd() {
			return this.evenOdd;
		}
		
		/**
		 * Get the domain with respect to the parent sequence
		 * @return the domain number or -1 if domain has not been set
		 */
		public int getDomain() {
			return this.domain;
		}
		
		/**
		 * the value of the qpcr primer flag
		 * @return 1 if yes, 0 if no, -1 if not set
		 */
		public int getContainsQpcrPrimer() {
			return this.qpcr;
		}
		
		/**
		 * Whether the probe contains a qpcr primer
		 * @return Whether the probe contains a qPCR primer
		 */
		public boolean containsQpcrPrimer() {
			if(this.qpcr == 1) return true;
			if(this.qpcr == 0 || this.qpcr == -1) return false;
			throw new IllegalStateException("qPCR state " + this.qpcr + " is illegal.");
		}
		
		/**
		 * Set probe sequence
		 * @param seq the Sequence object
		 */
		public void setSequence(Sequence seq) {
			this.sequence = seq;
		}
		
		/**
		 * Set 0-based start position of probe along parent sequence
		 * @param pos the start position
		 */
		public void setStartPos(int pos) {
			this.startPos = pos;
		}
		
		/**
		 * Set 0-based end position of probe along parent sequence
		 * @param pos the end position
		 */
		public void setEndPos(int pos) {
			this.endPos = pos;
		}
		
		/**
		 * Set tiling path number
		 * @param t the tiling path number
		 */
		public void setTilingPath(int t) {
			this.tilingPath = t;
		}
		
		/**
		 * Set the even/odd path number
		 * @param i must be 0, 1 or -1
		 */
		public void setEvenOdd(int i) {
			if(i!=1 && i!=0 && i!=-1) {
				throw new IllegalArgumentException("Not a valid even/odd path number");
			}
			this.evenOdd=i;
		}
		
		/**
		 * Set the domain number
		 * @param domain The domain number
		 */
		public void setDomain(int d) {
			this.domain = d;
		}
		
		/**
		 * Set whether probe contains a qpcr primer
		 * 1 for yes, 0 for no, -1 for not set
		 * @param p
		 */
		public void setQpcrPrimer(int p) {
			if(p != 1 && p != 0 && p != -1) {
				throw new IllegalArgumentException("qPCR primer must be 1, 0 or -1");
			}
			this.qpcr = p;
		}
		
		// Probe attributes
		// the probe sequence, not the parent sequence
		private Sequence sequence; 
		// the start position along the parent sequence
		private int startPos;
		// the end position along the parent sequence
		private int endPos;
		// the tiling path with respect to the parent sequence
		private int tilingPath;
		// even or odd position within the tiling path
		private int evenOdd;
		// domain with respect to the parent sequence
		private int domain;
		// whether the probe contains a qPCR primer
		private int qpcr;
		
	}
	
	
	protected final class StringPair {
		
		public StringPair(String a, String b) {
			this.first = a;
			this.second = b;
		}
		
		@Override
		public boolean equals(Object o) {
			if(o == null || o.getClass() != getClass()) return false;
			StringPair stringpair = (StringPair)o;
			if(this.first().equals(stringpair.first()) && this.second().equals(stringpair.second())) return true;
			return false;
		}
		
		@Override
		public int hashCode() {
			return((this.first + this.second).hashCode());
		}
		
		public String first() {return this.first;}
		public String second() {return this.second;}
		
		private String first;
		private String second;
	}

	
	/**
	 * A pair consisting of species name and gene name
	 * @author prussell
	 *
	 */
	protected final class SpeciesGene {
		
		public SpeciesGene() {
			this.species = null;
			this.gene = null;
			this.purpose = ProbePurpose.NOT_SET;
			//this.rnaClass = RnaClass.NOT_SET;
		}
		
		/**
		 * Parse the string representation that is generated by toString and get attributes
		 * @param stringRepresentation
		 */
		public void getFromStringRepresentation(String stringRepresentation) {
			this.setSpecies(getSpeciesFromStringRepresentation(stringRepresentation));
			this.setGene(getGeneNameFromStringRepresentation(stringRepresentation));
			this.setRnaClass(getRnaClassFromStringRepresentation(stringRepresentation));
			this.setPurpose(getPurposeFromStringRepresentation(stringRepresentation));
		}
		
		public SpeciesGene(String speciesName, String geneName) {
			this.species = Species.getSpecies(speciesName);
			this.gene = geneName;
			this.purpose = ProbePurpose.NOT_SET;
			//this.rnaClass = RnaClass.NOT_SET;
		}
		
		public SpeciesGene(Species s, String geneName) {
			this.species = s;
			this.gene = geneName;
			this.purpose = ProbePurpose.NOT_SET;
			//this.rnaClass = RnaClass.NOT_SET;
		}
		
		/**
		 * String representation consists of the fields separated by exclamation points
		 */
		@Override
		public String toString() {
		    return this.species.toString() + "!" + this.gene + "!" + this.purpose.toString() + "!" + this.rnaClass.toString();
		}
		
		/**
		 * Parse the string representation of a SpeciesGene and get the species
		 * @param representation
		 * @return the species
		 */
		public Species getSpeciesFromStringRepresentation(String representation) {
			StringParser parser = new StringParser();
			parser.parse(representation,"!");
			return Species.getSpecies(parser.asString(0));
		}
		
		/**
		 * Parse the string representation of a SpeciesGene and get the gene name
		 * @param representation
		 * @return the gene name
		 */
		public String getGeneNameFromStringRepresentation(String representation) {
			StringParser parser = new StringParser();
			parser.parse(representation,"!");
			return parser.asString(1);
		}
		
		/**
		 * Parse the string representation of a SpeciesGene and get the probe purpose
		 * @param representation
		 * @return the probe purpose
		 */
		public ProbePurpose getPurposeFromStringRepresentation(String representation) {
			StringParser parser = new StringParser();
			parser.parse(representation,"!");
			return ProbePurpose.getProbePurpose(parser.asString(2));
		}
		
		/**
		 * Parse the string representation of a SpeciesGene and get the RNA class
		 * @param representation
		 * @return the RNA class
		 */
		public RnaClass getRnaClassFromStringRepresentation(String representation) {
			StringParser parser = new StringParser();
			parser.parse(representation,"!");
			return RnaClass.getRnaClass(parser.asString(3));
		}
		
		
		/**
		 * Only check if species and gene are the same; class and purpose can be different
		 * @param speciesgene
		 * @return Whether species and gene are the same
		 */
		public boolean sameGeneAndSpecies(SpeciesGene speciesgene) {
			if(this.gene.length() != speciesgene.gene.length()) return false;
			if(this.species != speciesgene.species) return false;
			for(int i=0; i<this.gene.length(); i++) {
				if(this.gene.charAt(i) != speciesgene.gene.charAt(i)) return false;
			}
			return true;
		}
		
		/**
		 * Returns true if the two objects' species, gene and purpose fields represent the same values
		 * Ignores RNA class because this should be redundant, but purpose is important because the same gene could appear twice with different purposes
		 */
		@Override public boolean equals(Object o) {
			if(o == null || o.getClass() != getClass()) return false;
			
			SpeciesGene speciesgene = (SpeciesGene)o;
			
			if(this.gene.length() != speciesgene.gene.length()) return false;
			if(this.species != speciesgene.species) return false;
			if(this.purpose != speciesgene.purpose) return false;
			if(this.rnaClass != speciesgene.rnaClass) return false; 
			for(int i=0; i<this.gene.length(); i++) {
				if(this.gene.charAt(i) != speciesgene.gene.charAt(i)) return false;
			}
			return true;
		}
		
		/**
		 * Returns a hashcode of the concatenated string combining species and gene names and purpose
		 */
		@Override public int hashCode() {
			return this.species.toString().concat(this.gene).concat(this.purpose.toString()).concat(this.rnaClass.toString()).hashCode();
		}
		
		/**
		 * Set the species from a String
		 * @param speciesName
		 */
		public void setSpecies(String speciesName) {
			this.species = Species.getSpecies(speciesName);
		}
		
		/**
		 * Set the species from a Species
		 * @param s The species
		 */
		public void setSpecies(Species s) {
			this.species = s;
		}
		
		/**
		 * Set the gene name
		 * @param geneName
		 */
		public void setGene(String geneName) {
			this.gene = geneName;
		}
		
		/**
		 * Set the purpose for the gene
		 * @param p The probe purpose
		 */
		public void setPurpose(ProbePurpose p) {
			this.purpose = p;
		}
		
		/**
		 * Set the RNA class for the gene
		 * @param r The RNA class
		 */
		public void setRnaClass(RnaClass r) {
			this.rnaClass = r;
		}
		
		/**
		 * Get the species name as a String
		 * @return the species name
		 */
		public String getSpeciesName() {
			return this.species.toString();
		}
		
		/**
		 * Get the species as a Species
		 * @return the species
		 */
		public Species getSpecies() {
			return this.species;
		}
		
		/**
		 * Get the gene name
		 * @return the gene name
		 */
		public String getGeneName() {
			return this.gene;
		}
		
		/**
		 * Get the purpose for gene on array
		 * @return the purpose
		 */
		public ProbePurpose getPurpose() {
			return this.purpose;
		}
		
		/**
		 * Get the RNA class
		 * @return the RNA class
		 */
		public RnaClass getRnaClass() {
			return this.rnaClass;
		}
		
		// The species
		private Species species;
		// The gene name
		private String gene;
		// The purpose of the gene on array
		private ProbePurpose purpose;
		// The RNA class of the gene
		private RnaClass rnaClass;
	}
	
}


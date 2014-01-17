/**
 * 
 */
package nextgen.synbio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import nextgen.editing.RestrictionEnzymeFactory;
import nextgen.editing.TypeIISRestrictionEnzyme;
import nextgen.synbio.GibsonAssemblyOligoSet.FullOligo;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 * Oligo pool for Gibson assembly
 * Sets of sequences get their own primer and restriction enzyme
 */
public class GibsonAssemblyOligoPool {
	
	static Logger logger = Logger.getLogger(GibsonAssemblyOligoPool.class.getName());

	/**
	 * Read sequences from a list of fasta files and a fasta file of individual sequences
	 * @param fastaList List containing sequence set IDs and fasta files of the sequence sets, or null of not using
	 * @param fastaIndividualSeqs Fasta file of sequences that are their own sets, or null if not using
	 * @return Map of set ID to set of sequences under the ID
	 * @throws IOException
	 */
	private static Map<String, Collection<Sequence>> readSequencesAndAssignSetIDs(String fastaList, String fastaIndividualSeqs) throws IOException {
		Map<String, Collection<Sequence>> rtrn = new TreeMap<String, Collection<Sequence>>();
		if(fastaList != null) {
			logger.info("Reading list of fasta files from " + fastaList + " and loading sets of sequences...");
			FileReader r = new FileReader(fastaList);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			while(b.ready()) {
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() == 0) continue;
				if(s.getFieldCount() != 2) {
					b.close();
					throw new IllegalArgumentException("Line format: set_identifier   fasta_file");
				}
				String id = s.asString(0);
				String fasta = s.asString(1);
				FastaSequenceIO fsio = new FastaSequenceIO(fasta);
				rtrn.put(id, fsio.loadAll());
				logger.info(id + "\t" + fasta + "\t" + rtrn.get(id).size() + " sequences loaded.");
			}
			r.close();
			b.close();
			logger.info("Done loading fasta files from list.");
		}
		if(fastaIndividualSeqs != null) {
			logger.info("Reading individual sequences from " + fastaIndividualSeqs + "...");
			FastaSequenceIO fsio2 = new FastaSequenceIO(fastaIndividualSeqs);
			Collection<Sequence> indSeqs = fsio2.loadAll();
			int numLoaded = 0;
			for(Sequence seq : indSeqs) {
				String seqId = seq.getId();
				Collection<Sequence> thisSeq = new ArrayList<Sequence>();
				thisSeq.add(seq);
				rtrn.put(seqId, thisSeq);
				numLoaded++;
			}
			logger.info("Loaded " + numLoaded + " individual sequences.");
		}
		// Convert to upper case
		for(String id : rtrn.keySet()) {
			for(Sequence seq : rtrn.get(id)) {
				seq.setSequenceBases(seq.getSequenceBases().toUpperCase());
			}
		}
		return rtrn;
	}
	
	/**
	 * Design oligo sets for sequence sets and write output to table and fasta file
	 * @param outPrefix Output file prefix
	 * @param sequenceSets Sequence sets by set ID
	 * @param enzymes Restriction enzymes to try
	 * @param oligoSize Full oligo size
	 * @param overlapSize Gibson assembly overlap size
	 * @param primerLength Primer length
	 * @param primer3coreExecutable primer3_core executable
	 * @param divideSetsByCompatibleOligos If sequence sets do not have compatible enzymes, divide sets into subsets sharing a common compatible enzyme
	 * @throws IOException
	 */
	private static void designOligosAndWriteOutput(String outPrefix, Map<String, Collection<Sequence>> sequenceSets, Collection<TypeIISRestrictionEnzyme> enzymes, int oligoSize, int overlapSize, int primerLength, String primer3coreExecutable, boolean divideSetsByCompatibleOligos, double optimalTm) throws IOException {
		logger.info("Designing oligo sets and writing output...");
		boolean writeHeader = true;
		String errorFile = outPrefix + "_ERROR";
		logger.warn("Writing important error messages to " + errorFile + ".");
		FileWriter errorWriter = new FileWriter(errorFile);
		for(String setId : sequenceSets.keySet()) {
			logger.info("");
			logger.info("***** " + setId + " *****");
			if(divideSetsByCompatibleOligos) {
				Map<TypeIISRestrictionEnzyme, GibsonAssemblyOligoSet> subsetsByEnzyme = GibsonAssemblyOligoSet.divideByCompatibleEnzymes(sequenceSets.get(setId), enzymes, oligoSize, overlapSize, primerLength, primer3coreExecutable, errorWriter, optimalTm);
				for(TypeIISRestrictionEnzyme enzyme : subsetsByEnzyme.keySet()) {
					String prefix = setId + "_" + enzyme.getName();
					GibsonAssemblyOligoSet oligoSet = subsetsByEnzyme.get(enzyme);
					Collection<FullOligo> oligos = oligoSet.designOligoSet(errorWriter);
					GibsonAssemblyOligoSet.writeOutput(oligos, prefix, outPrefix, writeHeader, !writeHeader);
					writeHeader = false;
				}
			} else {
				GibsonAssemblyOligoSet oligoSet = new GibsonAssemblyOligoSet(sequenceSets.get(setId), enzymes, oligoSize, overlapSize, primerLength, primer3coreExecutable, optimalTm);
				Collection<FullOligo> oligos = oligoSet.designOligoSet(errorWriter);
				GibsonAssemblyOligoSet.writeOutput(oligos, setId, outPrefix, writeHeader, !writeHeader);
				writeHeader = false;
			}
		}
		errorWriter.close();
		logger.info("Done writing oligo pool.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-e", "File containing list of possible restriction enzymes, one per line (options: " + RestrictionEnzymeFactory.RestrictionEnzymeName.commaSeparatedList() + ")", true);
		p.addIntArg("-s", "Oligo size including primers, etc.", false, GibsonAssemblyOligoSet.DEFAULT_OLIGO_SIZE);
		p.addIntArg("-v", "Overlap size for Gibson assembly", false, GibsonAssemblyOligoSet.DEFAULT_OVERLAP_SIZE);
		p.addIntArg("-p", "Primer length", false, GibsonAssemblyOligoSet.DEFAULT_PRIMER_SIZE);
		p.addStringArg("-o", "Output file prefix", true);
		p.addStringArg("-p3", "Primer3core executable", true);
		p.addStringArg("-f", "Fasta file of sequences that get individual primers", false, null);
		p.addStringArg("-fl", "File containing list of fasta files. Each fasta file is a set of sequences that get one primer. Each line of list file: set_identifier   fasta_file", false, null);
		p.addBooleanArg("-d", "If sequence sets do not have compatible enzymes, divide sets into subsets sharing a common compatible enzyme", false, false);
		p.addBooleanArg("-debug", "Debug logging", false, false);
		p.addDoubleArg("-tm", "Optimal TM for primers", true);
		p.parse(args);
		if(p.getBooleanArg("-debug")) {
			logger.setLevel(Level.DEBUG);
			GibsonAssemblyOligoSet.logger.setLevel(Level.DEBUG);
		}
		Collection<TypeIISRestrictionEnzyme> enzymes = RestrictionEnzymeFactory.readFromFileAsTypeIIS(p.getStringArg("-e"));
		int oligoSize = p.getIntArg("-s");
		int overlapSize = p.getIntArg("-v");
		int primerLength = p.getIntArg("-p");
		String outPrefix = p.getStringArg("-o");
		String primer3core = p.getStringArg("-p3");
		String fastaFile = p.getStringArg("-f");
		String fastaList = p.getStringArg("-fl");
		boolean divide = p.getBooleanArg("-d");
		double optimalTm = p.getDoubleArg("-tm");
		
		Map<String, Collection<Sequence>> sequenceSets = readSequencesAndAssignSetIDs(fastaList, fastaFile);
		designOligosAndWriteOutput(outPrefix, sequenceSets, enzymes, oligoSize, overlapSize, primerLength, primer3core, divide, optimalTm);
		
		logger.info("");
		logger.info("All done.");


	}

}

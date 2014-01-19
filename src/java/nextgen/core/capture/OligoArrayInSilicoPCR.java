package nextgen.core.capture;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;


import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.datastructures.Pair;
import broad.core.motif.SearchException;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.isPCRLike;
import broad.core.sequence.Sequence;

/**
 * Perform in silico PCR on primers and corresponding oligo sets
 * Verify that primers amplify all the oligos they are associated with
 * Verify that primers do not amplify any oligos in other sets
 * Write all output to a file 
 * @author prussell
 *
 */
public class OligoArrayInSilicoPCR {
	
	private Map<Pair<String>, List<Sequence>> oligoSetsByPrimer;
	private FileWriter resultsWriter;
	private static Logger logger = Logger.getLogger(OligoArrayInSilicoPCR.class.getName());
	
	/**
	 * @param input Input file: either table of primers and oligos (format: oligo_ID  left_primer  right_primer  oligo_sequence) or file containing list of such table files
	 * @param output Output file for record of problems
	 * @throws IOException
	 */
	public OligoArrayInSilicoPCR(String input, String output) throws IOException {
		// Check if input file is a table of oligos or a list of table files
		FileReader r = new FileReader(input);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		s.parse(b.readLine());
		r.close();
		b.close();
		if(s.getFieldCount() == 1) {
			logger.info("Interpreting " + input + " as list of table files");
			oligoSetsByPrimer = readFromTables(input);
		} else if(s.getFieldCount() == 4) {
			logger.info("Interpreting " + input + " as table of oligos");
			oligoSetsByPrimer = readFromTable(input);
		} else {
			throw new IllegalArgumentException("Input file must be either table of oligos (format: oligo_ID  left_primer  right_primer  oligo_sequence) or list of such files");
		}
		logger.info("Writing output to " + output);
		resultsWriter = new FileWriter(output);
		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				try {
					resultsWriter.close();
				} catch (IOException e) {
					e.printStackTrace();
					throw new IllegalStateException("Couldn't close file writer");
				}
			}
		});
	}
	
	/**
	 * Check that the primer pair amplifies all the oligos
	 * Write any problems to output stream and keep going
	 * @param primerPair Primer pair
	 * @param oligos Oligos
	 * @throws SearchException
	 * @throws IOException
	 */
	private void checkSelf(Pair<String> primerPair, List<Sequence> oligos) throws SearchException, IOException {
		isPCRLike ispcr = new isPCRLike(primerPair, oligos);
		Collection<Sequence> amplicons = ispcr.getAllPossibleAmplicons();
		Collection<String> ampliconNames = new TreeSet<String>();
		for(Sequence a : amplicons) {
			ampliconNames.add(a.getId());
		}
		for(Sequence oligo : oligos) {
			if(!ampliconNames.contains(oligo.getId())) {
				resultsWriter.write("Oligo " + oligo.getId() + " not amplified by self primer pair " + primerPair.getValue1() + " " + primerPair.getValue2() + "\n");
			} else {
				logger.debug("Oligo " + oligo.getId() + " correctly amplified by self primer pair " + primerPair.getValue1() + " " + primerPair.getValue2());
			}
		}
	}
	
	/**
	 * Check that the primer pair does not amplify any part of any of the oligos
	 * Write any problems to output stream and keep going
	 * @param primerPair Primer pair
	 * @param oligos Oligos
	 * @throws SearchException
	 * @throws IOException
	 */
	private void checkOther(Pair<String> primerPair, List<Sequence> oligos) throws SearchException, IOException {
		isPCRLike ispcr = new isPCRLike(primerPair, oligos);
		Collection<Sequence> amplicons = ispcr.getAllPossibleAmplicons();
		if(amplicons.isEmpty()) {
			logger.debug("Primers " + primerPair.getValue1() + " " + primerPair.getValue2() + " correctly fail to amplify " + oligos.size() + " oligos");
		}
		for(Sequence amplicon : amplicons) {
			resultsWriter.write("Primers " + primerPair.getValue1() + " " + primerPair.getValue2() + " amplify sequence " + amplicon.getId() + " " + amplicon.getSequenceBases() + "\n");
		}
	}
	
	/**
	 * Check that all primer pairs amplify their associated oligos
	 * Check that primer pairs do not amplify any part of other oligos
	 * Write any problems to output stream and keep going
	 * @throws SearchException
	 * @throws IOException
	 */
	public void checkAll() throws SearchException, IOException {
		for(Pair<String> primer : oligoSetsByPrimer.keySet()) {
			logger.info("Checking primers " + primer.getValue1() + " " + primer.getValue2());
			checkSelf(primer, oligoSetsByPrimer.get(primer));
			for(Pair<String> otherPrimer : oligoSetsByPrimer.keySet()) {
				if(primer.equals(otherPrimer)) continue;
				checkOther(primer, oligoSetsByPrimer.get(otherPrimer));
			}
		}
	}
	
	/**
	 * 
	 * @param tableFileList
	 * @return
	 * @throws IOException
	 */
	private static Map<Pair<String>, List<Sequence>> readFromTables(String tableFileList) throws IOException {
		
		// Read the table file names from the list
		Collection<String> tableFiles = new ArrayList<String>();
		FileReader r = new FileReader(tableFileList);
		BufferedReader b = new BufferedReader(r);
		while(b.ready()) {
			tableFiles.add(b.readLine());
		}
		r.close();
		b.close();
		
		Map<Pair<String>, List<Sequence>> rtrn = new HashMap<Pair<String>, List<Sequence>>();
		
		for(String table : tableFiles) {
			Map<Pair<String>, List<Sequence>> array = readFromTable(table);
			for(Pair<String> primers : array.keySet()) {
				if(rtrn.containsKey(primers)) {
					throw new IllegalArgumentException("Primer pair " + primers.getValue1() + " " + primers.getValue2() + " is used twice");
				}
			}
			rtrn.putAll(array);
		}
		
		return rtrn;
		
	}
	
	/**
	 * Read oligos and primers from a table
	 * Line format: oligo_ID  left_primer  right_primer  oligo_sequence
	 * @param file Table file
	 * @return Map of pair of primers to collection of oligos they are supposed to amplify
	 * @throws IOException
	 */
	private static Map<Pair<String>, List<Sequence>> readFromTable(String file) throws IOException {
		
		Map<Pair<String>, List<Sequence>> rtrn = new HashMap<Pair<String>, List<Sequence>>();
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			if(s.getFieldCount() == 0) {
				continue;
			}
			if(s.getFieldCount() != 4) {
				r.close();
				b.close();
				throw new IllegalArgumentException("Line format: oligo_ID  left_primer  right_primer  oligo_sequence");
			}
			Sequence oligoSeq = new Sequence(s.asString(0));
			oligoSeq.setSequenceBases(s.asString(3));
			Pair<String> primers = new Pair<String>(s.asString(1), s.asString(2));
			if(!rtrn.containsKey(primers)) {
				rtrn.put(primers, new ArrayList<Sequence>());
			}
			rtrn.get(primers).add(oligoSeq);
		}
		r.close();
		b.close();
		
		return rtrn;
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws SearchException 
	 */
	public static void main(String[] args) throws IOException, SearchException {
		
		//logger.setLevel(Level.DEBUG);
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input file: either table of primers and oligos (format: oligo_ID  left_primer  right_primer  oligo_sequence) or file containing list of such table files", true);
		p.addStringArg("-o", "Output file", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		
		OligoArrayInSilicoPCR a = new OligoArrayInSilicoPCR(input, output);
		a.checkAll();
		
		logger.info("");
		logger.info("All done.");

	}

}

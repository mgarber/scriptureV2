package nextgen.core.capture;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.PSLRecord;
import nextgen.core.readers.PSLReader;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

/**
 * Filter probes produced by oligo array designer
 * Blat the probe sequences against the genome, then provide the blat alignments to this class
 * Filters probes that align with a certain length and % identity more than a specified number of times
 * @author prussell
 *
 */
public class OligoArrayBlatFilter {
	
	private Map<String, Map<Filter, Integer>> probeAlignments;
	private Collection<Filter> filters;
	private Collection<String> probesToRemove;
	private static Logger logger = Logger.getLogger(OligoArrayBlatFilter.class.getName());
	
	/**
	 * @param pslFile PSL file of probe sequences blatted to genome
	 * @param alignmentFilters Filters to apply
	 * @throws FileNotFoundException
	 */
	private OligoArrayBlatFilter(String pslFile, Collection<Filter> alignmentFilters) throws FileNotFoundException {
		filters = alignmentFilters;
		readFile(pslFile);
	}
	
	/**
	 * Dummy constructor do not use
	 */
	private OligoArrayBlatFilter() {}

	/**
	 * Read the PSL file and store counts of probe alignments with respect to filter criteria
	 * @param pslFile PSL file of probes blatted to genome
	 * @throws FileNotFoundException
	 */
	private void readFile(String pslFile) throws FileNotFoundException {
		logger.info("");
		logger.info("Reading from PSL file " + pslFile);
		probeAlignments = new TreeMap<String, Map<Filter, Integer>>();
		PSLReader pslReader = new PSLReader(pslFile);
		while(pslReader.hasNext()) {
			PSLRecord record = pslReader.next();
			String name = record.getName();
			for(Filter filter : filters) {
				if(filter.hasMinLengthAndPctId(record)) {
					if(!probeAlignments.containsKey(name)) {
						probeAlignments.put(name, new HashMap<Filter, Integer>());
					}
					if(!probeAlignments.get(name).containsKey(filter)) {
						probeAlignments.get(name).put(filter, Integer.valueOf(1));
					} else {
						probeAlignments.get(name).put(filter, Integer.valueOf(probeAlignments.get(name).get(filter).intValue() + 1));
					}
				}
			}
		}
		logger.info("Determining which probes to remove");
		probesToRemove = new TreeSet<String>();
		for(String probe : probeAlignments.keySet()) {
			for(Filter filter : probeAlignments.get(probe).keySet()) {
				int numAlignments = probeAlignments.get(probe).get(filter).intValue();
				if(numAlignments > filter.getMaxNumAlignments()) {
					probesToRemove.add(probe);
				}
			}
		}
	}
	
	/**
	 * Write new filtered oligo table using original oligo table created by oligo pool designer
	 * @param inTable Table of oligos with ID in first column
	 * @param outTable Output table
	 * @throws IOException
	 */
	private void writeFilteredFile(String inTable, String outTable) throws IOException {
		logger.info("");
		logger.info("Filtering " + inTable + " and writing to " + outTable);
		String outRemoved = outTable + ".removed_probes.out";
		FileWriter tableWriter = new FileWriter(outTable);
		FileWriter removedWriter = new FileWriter(outRemoved);
		FileReader r = new FileReader(inTable);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			String name = s.asString(0);
			if(probesToRemove.contains(name)) {
				removedWriter.write(name + "\n");
			} else {
				tableWriter.write(line + "\n");
			}
		}
		r.close();
		b.close();
		tableWriter.close();
		removedWriter.close();
		logger.info("Done writing filtered table.");
	}

	private static Filter fromString(String description) {
		StringParser p = new StringParser();
		p.parse(description,",");
		OligoArrayBlatFilter o = new OligoArrayBlatFilter();
		return o.new Filter(p.asInt(0), p.asFloat(1), p.asInt(2));
	}
	
	private static Collection<Filter> fromStrings(Collection<String> descriptions) {
		Collection<Filter> rtrn = new ArrayList<Filter>();
		for(String d : descriptions) {
			rtrn.add(fromString(d));
		}
		return rtrn;
	}

	/**
	 * A filter for probes based on blat alignments
	 * @author prussell
	 *
	 */
	private class Filter {
		
		private int minLen;
		private float minPctId;
		private int maxNumAlignments;
		
		/**
		 * @param len Minimum length of aligment
		 * @param pctId Minimum percent identity of aligment
		 * @param numAlignments Max number of times a probe can align to the genome with these criteria
		 */
		private Filter(int len, float pctId, int numAlignments) {
			if(minLen < 0) {
				throw new IllegalArgumentException("Min length must be >= 0");
			}
			if(pctId < 0 || pctId > 1) {
				throw new IllegalArgumentException("Pct identity must be between 0 and 1");
			}
			if(numAlignments < 0) {
				throw new IllegalArgumentException("Num alignments must be >= 0");
			}
			minLen = len;
			minPctId = pctId;
			maxNumAlignments = numAlignments;
			logger.info("Created filter with min length " + minLen + ", min percent identity " + minPctId + ", and max number of alignments " + maxNumAlignments + ".");
		}
		
		public int getMaxNumAlignments() {
			return maxNumAlignments;
		}
		
		/**
		 * @param alignment A PSL record of a probe aligned to genome with blat
		 * @return Whether the probe meets the length and percent identity criteria for the filter, meaning it may be filtered if there are enough alignment like this one
		 */
		public boolean hasMinLengthAndPctId(PSLRecord alignment) {
			return alignment.getSize() >= minLen && alignment.getPercentIdentity() >= minPctId;
		}
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("Filter probes produced by oligo array designer. First Blat the probe sequences against the genome, then provide the blat alignments to this program. The program filters probes that align with a certain length and % identity more than a specified number of times");
		p.addStringListArg("-f", "Filter to use (repeatable). Format: <min_alignment_length(int)>,<min_pct_identity(float)>,<max_num_alignments> e.g. 60,0.9,1", true);
		p.addStringArg("-i", "Oligo designer full output table with oligo ID in first column", true);
		p.addStringArg("-o", "Output filtered oligo design table", true);
		p.addStringArg("-p", "PSL file of probes aligned to genome", true);
		p.parse(args, true);
		String inputTable = p.getStringArg("-i");
		String outputTable = p.getStringArg("-o");
		String pslFile = p.getStringArg("-p");
		ArrayList<String> filters = p.getStringListArg("-f");
		
		OligoArrayBlatFilter oabf = new OligoArrayBlatFilter(pslFile, fromStrings(filters));
		oabf.writeFilteredFile(inputTable, outputTable);
		
		logger.info("");
		logger.info("All done.");
		
	}

}

/**
 *  Jesse Engreitz
 *  March 12, 2012
 *  
 *  Design probes for hybrid capture of methylation targets.  
 *  Designs three sets of probes, against: uncoverted sequence, bisulfite-converted fully unmethylated, and bisulfite-converted fully methylated DNA
 */

package broad.pda.methylation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;

import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.datastructures.Pair;
import broad.core.error.ParseException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.primer3.PrimerPair;

import broad.pda.methylation.MethylationUtils;
import broad.pda.methylation.OutputPadlock;

/**
 * @author engreitz
 *
 */
public class HybridCaptureDesigner {

	static public int REPEAT_CUTOFF = 40;
	static public int OVERLAP = 90;
	static public int PROBE_LENGTH = 120;
	static public int OVERHANG = (int)(Math.ceil(PROBE_LENGTH / 2));
	
	//private static final Pattern SOFT_MASKED_PAT = Pattern.compile("[Nacgt]+");
	public HybridCaptureDesigner() {}
	
	public HybridCaptureDesigner(String regionFile, String outFile, String reference, int arraySize) throws IOException, ParseException {
		Writer outputOrder = new BufferedWriter(new FileWriter(outFile + ".order.txt"));
		Writer outputTable = new BufferedWriter(new FileWriter(outFile + ".txt"));
		Writer outputFasta = new BufferedWriter(new FileWriter(outFile + ".fa"));
		Writer outputPrimer = new BufferedWriter(new FileWriter(outFile + ".primers.txt"));
		FastaSequenceIO fsio = new FastaSequenceIO(reference);
		BEDReader bed = new BEDReader(regionFile);
		List<BED> regions = bed.getAnnotationList();

		List<Probe> allProbes = new ArrayList<Probe>();
		
		getProbes(allProbes, regions, fsio, PROBE_LENGTH, OVERLAP, OVERHANG, REPEAT_CUTOFF);
		allProbes = filterProbes(allProbes, REPEAT_CUTOFF);
		assignPrimers(allProbes, outputPrimer);
		
		// Write output
		List<String> senseProbes = new ArrayList<String>();
		List<String> antisenseProbes = new ArrayList<String>();
		
		// Write the master probe list
		for (Probe probe : allProbes) {
			outputTable.write(probe.toString() + "\n");
			outputFasta.write(probe.toFasta() + "\n");
			senseProbes.add(probe.toString(false));
			antisenseProbes.add(probe.toString(true));
		}
		
		// Write the table of probes to order (fill array with sense + antisense)
		List<String> finalList = MethylationUtils.fillArray(senseProbes, antisenseProbes, arraySize);
		for (String line : finalList) 
			outputOrder.write(line + "\n");
		outputOrder.close();
		outputTable.close();
		outputPrimer.close();
	}
		
	private void getProbes(List<Probe> allProbes, List<BED> regions, FastaSequenceIO fsio, int probeLength, int overlap, int overhang, int repeatCutoff) throws IOException {
		// Get the genomic sequence for each region
		List<SequenceRegion> sequences = new ArrayList<SequenceRegion>();
		for (BED region : regions) {
			region.setStart(region.getStart() - overhang);
			region.setEnd(region.getEnd() + overhang);
			SequenceRegion newRegion = new SequenceRegion(region.getChromosome(), region);
			newRegion.setName(region.getName());
			sequences.add(newRegion);
		}
		System.out.printf("Extracting sequences ...");
		fsio.extractRegions(sequences);
		System.out.println(" done.");
		
		addTiledProbes(allProbes, sequences, false, false, probeLength, overlap);
		addTiledProbes(allProbes, sequences, true, false, probeLength, overlap);
		addTiledProbes(allProbes, sequences, true, true, probeLength, overlap);
	}
	
	
	private void addTiledProbes(List<Probe> probeList, List<SequenceRegion> maskedSequences, Boolean converted, Boolean methylated, int probeLength, int overlap) {
		for (SequenceRegion sr : maskedSequences) {
			addTiledProbesForSequence(probeList, sr, converted, methylated, probeLength, overlap);
		}
	}
	
	private void addTiledProbesForSequence(List<Probe> probeList, SequenceRegion sr, Boolean converted, Boolean methylated, int probeLength, int overlap) {
		String seq = new String(sr.getSequenceBases());
		String revComp = Sequence.reverseSequence(seq);
		
		if (converted) {
			//System.out.println(methylated + " " + seq);
			seq = MethylationUtils.bisulfiteConvert(seq, methylated);
			revComp = MethylationUtils.bisulfiteConvert(revComp, methylated);
			//System.out.println(methylated + " " + seq);
		}
		//System.out.println(sr.getName() + " " + seq.length());
		for (int i = 0; i <= sr.getLength() - probeLength; i += probeLength - overlap) {
			String newId = sr.getName() + "|" + i;
			
			// Create probe as the reverse complement of the target sequence
			Probe newProbe = new Probe(Sequence.reverseSequence(seq.substring(i, i+probeLength+1)), newId, sr.getChromosome(), sr.getStart() + i, sr.getStart() + i + probeLength, "+", converted, methylated);
			probeList.add(newProbe);
			
			// For bisulfite converted stuff, we need to design separately against + and - strands
			if (converted) {
				String segment = new StringBuffer(revComp).reverse().substring(i, i+probeLength+1);
				segment = new StringBuffer(segment).reverse().toString();
				//System.out.println(segment);
				Probe revProbe = new Probe(Sequence.reverseSequence(segment), newId, sr.getChromosome(), sr.getStart() + i, sr.getStart() + i + probeLength, "-", converted, methylated);
				probeList.add(revProbe);
				//System.out.println(Sequence.reverseSequence(segment));
			}
		}
		


	}
	
	
	private List<Probe> filterProbes(List<Probe> probeList, int repeatCutoff) {
		List<Probe> filtered = new ArrayList<Probe>();
		for (Probe probe : probeList) {
			if (probe.sequence.countSoftMaskedBases() > repeatCutoff) {
				//System.err.println("Removing " + probe.sequence.getSequenceBases());
				continue;
			}
			filtered.add(probe);
		}
		return filtered;
	}
	
	
	
	private void assignPrimers(List<Probe> probeList, Writer out) throws FileNotFoundException, IOException {
		// Get the list of primers from Pam
		List<PrimerPair> primers = loadPrimers();
		System.out.println(primers.size() + " primers loaded.");
		// Choose three primers: one for unconverted, one for converted methylated, and one for converted unmethylated
		
		// Check that primer pairs cannot form primer dimers
		int maxPrimerOverlap = 2;
		Iterator<PrimerPair> itr = primers.iterator();
		while (itr.hasNext()) {
			PrimerPair pp = itr.next();
			String leftOverlap = getKmerTail(pp.getLeftPrimer(), maxPrimerOverlap);
			String rightOverlap = getKmerTail(Sequence.reverseSequence(pp.getRightPrimer()), maxPrimerOverlap);
			if (leftOverlap.equals(rightOverlap)) {
				itr.remove();
			}
		}
		System.out.println(primers.size() + " primers remaining.");
				
		
		// Check that primer pairs cannot prime off of any non-specific place 
		// Generate a HashSet of kmers
		int k = 9;  // Chosen as the highest value where we still retain some of our list of ~9000 primer pairs
		HashSet<String> kmers = new HashSet<String>();
		for (Probe probe : probeList)
			addAllSubstringsHashSet(kmers, probe.sequence.getSequenceBases(), k);
		
		Iterator<PrimerPair> jtr = primers.iterator();
		while (jtr.hasNext()) {
			PrimerPair pp = jtr.next();
			String leftOverlap = getKmerTail(pp.getLeftPrimer(), k);
			String rightOverlap = getKmerTail(pp.getRightPrimer(), k);
			if (kmers.contains(leftOverlap) ||
				kmers.contains(rightOverlap) ||
				kmers.contains(Sequence.reverseSequence(leftOverlap)) ||
				kmers.contains(Sequence.reverseSequence(rightOverlap)))
				jtr.remove();
		}
		System.out.println(primers.size() + " primers remaining.");
		if (primers.size() < 3) {
			System.out.println("No primers assigned.");
		} else {
			PrimerPair ppUnconverted = primers.get(0);
			PrimerPair ppConvertedUnmethylated = primers.get(1);
			PrimerPair ppConvertedMethylated = primers.get(2);
			out.write("Unconverted\t" + ppUnconverted.toString() + "\n");
			out.write("ConvertedUnmethylated\t" + ppConvertedUnmethylated.toString() + "\n");
			out.write("ConvertedMethylated\t" + ppConvertedMethylated.toString() + "\n");
			
			for (Probe probe : probeList) {
				if (!probe.converted) {
					probe.setPrimers(ppUnconverted);
				} else if (probe.methylated) {
					probe.setPrimers(ppConvertedMethylated);
				} else {
					probe.setPrimers(ppConvertedUnmethylated);
				}
			}
		}
	}
	
	
	private String getKmerTail(String s, int k) {
		int l = s.length();
		return s.substring(l - k);
	}
	
	
	private String primerFile = "primers.txt";
	private List<PrimerPair> loadPrimers() throws FileNotFoundException, IOException {
		List<PrimerPair> primers = new ArrayList<PrimerPair>();
		File inFile = new File(primerFile);
		BufferedReader b = new BufferedReader(new FileReader(inFile));
		
		while (b.ready()) {
			PrimerPair pp = new PrimerPair(b.readLine().split("\t"));
			primers.add(pp);
		}
		return primers;
	}
	
	
	/**
	 * Add all substrings of a string for a specified length to a hashset
	 * @param s the string
	 * @param length the length of substrings
	 * @return set of all substrings
	 */
	private static void addAllSubstringsHashSet(HashSet<String> set, String s, int length) {
		for(int i=0; i <= s.length()-length; i++) {
			set.add(s.substring(i,i + length));
		}
	}
		
		/*
		List<Pair<Integer> > chunks = new ArrayList<Pair<Integer> >();
		// Split sequence into chunks if it is interrupted by long stretches of repetitive sequence
		int lastEnd = 0;
		Matcher m = SOFT_MASKED_PAT.matcher(seq);
		while (m.find()) {
			if ((m.start() > 0) & (m.end() - m.start() > REPEAT_CUTOFF)) {
				System.out.println("Repeat found at coordinates " + m.start() + "-" + m.end() + "\t" + (m.end() - m.start()));

				int start = lastEnd;
				int end = m.start() + REPEAT_CUTOFF;
				int hardStart = (start < PROBE_LENGTH) ? PROBE_LENGTH - start : 0;
				//int hardEnd = (end > )
					
				//System.out.println("Segment is: " + start + " " + end + " " + hardStart + " " + hardEnd);

				
				chunks.add(new Pair<Integer>(lastEnd, m.start())); 
				lastEnd = m.end();
			}
		}
		
		if (chunks.size() == 0)
			chunks.add(new Pair<Integer>(0, seq.length()));
	
		for (Pair<Integer> chunk : chunks) {
			String bases = seq.substring(chunk.getValue1(), chunk.getValue2());
			Integer hardStart = 
			//addTiledProbesForChunk(probeList, bases, chunk.getValue1() < PROBE_LENGTH, chunk.getValue2() > seq.length() - PROBE_LENGTH);

	}
	
	
	private static void addTiledProbesForChunk(List<Sequence> probeList) {
	
	}
	*/
	
	
	public class Probe {
		public Boolean converted, methylated;
		public String chr, strand;
		public int start, end;
		public Sequence sequence;
		public PrimerPair pp;
		
		public Probe(String seq, String name, String chr, int start, int end, String strand, Boolean converted, Boolean methylated) {
			this.sequence = new Sequence(name);
			this.sequence.setSequenceBases(seq);
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.converted = converted;
			this.methylated = methylated;
			this.strand = strand;
		}
		
		public void setPrimers(PrimerPair pp) {
			this.pp = pp;
		}
		
		public String getFullName() {
			String name = "";
			
			if (!converted)
				name = name + "Unconverted";
			else if (converted & !methylated)
				name = name + "ConvertedUnmethylated";
			else if (converted & methylated)
				name = name + "ConvertedMethylated";
			else
				name = name + "Unknown";
			
			name = name + "|" + sequence.getId() + "|" + strand;
			
			return name;
		}
		
		public String getFullSequence() {
			return getFullSequence(false);
		}
		
		public String getFullSequence(boolean antisense) {
			String rtrn = pp.getLeftPrimer() + sequence.getSequenceBases() + Sequence.reverseSequence(pp.getRightPrimer());
			if (antisense) rtrn = Sequence.reverseSequence(rtrn);
			return rtrn;
		}
		
		public String toString() {
			return toString(false);
		}
	
		public String toString(boolean antisense) {
			return String.format("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chr, start, end, getFullName(), OVERHANG, strand, converted.toString(), methylated.toString(), pp.getLeftPrimer(), pp.getRightPrimer(), sequence.getSequenceBases(), getFullSequence(antisense));
		}
		
		public String toFasta() {
			StringBuilder s = new StringBuilder(">");
			s.append(getFullName());
			s.append("\n");
			s.append(getFullSequence());
			return s.toString();
		}
	}
	
	
	
	public boolean systemTest(String reference) throws IOException {
		System.out.println("Conducting test ...");
		List<BED> regions = new ArrayList<BED>();
		BED testbed = new BED("testRegion", "chr10", 47487574, 47487605);
		regions.add(testbed);
		
		FastaSequenceIO fsio = new FastaSequenceIO(reference);
		List<Probe> allProbes = new ArrayList<Probe>();
		
		getProbes(allProbes, regions, fsio, testbed.length(), 1, 0, 100);
		
		//allProbes = filterProbes(allProbes);
		String original = "CAGGCCGGAAATGCCGGCTTGAAGCTGGCTTC";
		String unconverted = "GAAGCCAGCTTCAAGCCGGCATTTCCGGCCTG";
		String convertedMethylatedPlus = "AAAACCAACTTCAAACCGACATTTCCGACCTA";
		String convertedUnmethylatedPlus = "AAAACCAACTTCAAACCAACATTTCCAACCTA";
		String convertedMethylatedMinus = "CAAACCGAAAATACCGACTTAAAACTAACTTC";
		String convertedUnmethylatedMinus = "CAAACCAAAAATACCAACTTAAAACTAACTTC";
		
		System.out.println("The correct probes should be:");
		System.out.println("Unconverted:\t\t" + unconverted);
		System.out.println("ConMeth+:\t\t" + convertedMethylatedPlus);
		System.out.println("ConMeth-:\t\t" + convertedMethylatedMinus);
		System.out.println("ConUnmeth+:\t\t" + convertedUnmethylatedPlus);
		System.out.println("ConUnmeth-:\t\t" + convertedUnmethylatedMinus);
		
		for (Probe probe : allProbes) 
			System.out.println(probe);
		
		
		return true;
	}
	
	private static String USAGE = "HybridCaptureDesigner parameters:\n\n\t-in\tRegion file in BED format\n\t-ref\tBasename for reference FASTA file and index\n";
	
	/**
	 * @param args
	 * 
	 * Reference FASTA file should be repeat soft-masked.
	 */
	public static void main(String[] args) throws IOException, ParseException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "write");
		
		final String in = argmap.getInput();
		final String out = argmap.getOutput();
		final String ref = argmap.getMandatory("ref");
		final int arraySize = argmap.getInteger("arraySize");
		
		final boolean test = argmap.isPresent("test");
		if (test) {
			HybridCaptureDesigner hb = new HybridCaptureDesigner();
			hb.systemTest(ref);
		} else {
			new HybridCaptureDesigner(in, out, ref, arraySize);
		}
	}

}

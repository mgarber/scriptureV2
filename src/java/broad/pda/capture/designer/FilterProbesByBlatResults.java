package broad.pda.capture.designer;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class FilterProbesByBlatResults {

	public static boolean hasHyb(String matchLengths, int minLength) {
		StringParser p = new StringParser();
		p.parse(matchLengths,",");
		for(int i=0; i<p.getFieldCount(); i++) {
			if(p.asInt(i) >= minLength) return true;
		}
		return false;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-p","fasta file of all probes from main array",true);
		p.addStringArg("-rp", "psl file of blat results for ribosomal RNA against transcriptome, perfect matches", true);
		p.addStringArg("-ri", "psl file of blat results for ribosomal RNA against transcriptome, imperfect matches", true);
		p.addStringArg("-gp", "psl file of blat results for non-ribosomal pulldown probes to the genome, perfect matches",true);
		p.addStringArg("-gi", "psl file of blat results for non-ribosomal pulldown probes to the genome, imperfect matches", true);
		//p.addStringArg("-o", "output fasta file of probes passing filter", true);
		//p.addStringArg("-or", "output fasta file of removed probes", true);
		p.addStringArg("-idt", "full design table", true);
		p.addStringArg("-idf", "full design fasta", true);
		p.addStringArg("-odt", "output filtered design table", true);
		p.addStringArg("-odf", "output filtered design fasta", true);
		p.parse(args);
		String rPerfect = p.getStringArg("-rp");
		String rImperfect = p.getStringArg("-ri");
		String gPerfect = p.getStringArg("-gp");
		String gImperfect = p.getStringArg("-gi");
		String probeFasta = p.getStringArg("-p");
		//String outFasta = p.getStringArg("-o");
		//String outRemoved = p.getStringArg("-or");
		String designTable = p.getStringArg("-idt");
		String designFasta = p.getStringArg("-idf");
		String outDesignTable = p.getStringArg("-odt");
		String outDesignFasta = p.getStringArg("-odf");
		
		HashMap<String,Boolean> keepProbes = new HashMap<String,Boolean>();
		List<String> probesToFilter = new ArrayList<String>();
		
		
		// Check ribosomal RNA perfect matches to transcriptome
		FileReader readRp = new FileReader(rPerfect);
		System.out.println("Reading perfect rRNA matches to transcriptome from file " + rPerfect);
		StringParser stringparse = new StringParser();
		StringParser nameparse = new StringParser();
		BufferedReader brp = new BufferedReader(readRp);
		while(brp.ready()) {
			String line = brp.readLine();
			stringparse.parse(line);
			if(stringparse.getFieldCount() != 21) continue;
			String probeName = stringparse.asString(9);
			if(probesToFilter.contains(probeName)) continue;
			if(hasHyb(stringparse.asString(18),30)) {
				if(keepProbes.containsKey(probeName)) {
					keepProbes.put(probeName, Boolean.valueOf(false));
					probesToFilter.add(probeName);
					System.out.println("Filtering " + probeName);
				}
				else keepProbes.put(probeName, Boolean.valueOf(true));
			}
		}
		keepProbes.clear();
		
		// Check ribosomal RNA imperfect matches to transcriptome
		FileReader readRi = new FileReader(rImperfect);
		System.out.println("Reading imperfect rRNA matches to transcriptome from file " + rImperfect);
		BufferedReader bri = new BufferedReader(readRi);
		while(bri.ready()) {
			String line = bri.readLine();
			stringparse.parse(line);
			if(stringparse.getFieldCount() != 21) continue;
			String probeName = stringparse.asString(9);
			if(probesToFilter.contains(probeName)) continue;
			if(hasHyb(stringparse.asString(18),60)) {
				if(keepProbes.containsKey(probeName)) {
					keepProbes.put(probeName, Boolean.valueOf(false));
					probesToFilter.add(probeName);
					System.out.println("Filtering " + probeName);
				}
				else keepProbes.put(probeName, Boolean.valueOf(true));
			}
		}
		keepProbes.clear();
		
		// Check nonribosomal RNA perfect matches to genome
		FileReader readNp = new FileReader(gPerfect);
		System.out.println("Reading perfect non-rRNA matches to genome from file " + gPerfect);
		BufferedReader bnp = new BufferedReader(readNp);
		while(bnp.ready()) {
			String line = bnp.readLine();
			stringparse.parse(line);
			if(stringparse.getFieldCount() != 21) continue;
			String probeName = stringparse.asString(9);
			if(probesToFilter.contains(probeName)) continue;
			if(hasHyb(stringparse.asString(18),30)) {
				if(keepProbes.containsKey(probeName)) {
					keepProbes.put(probeName, Boolean.valueOf(false));
					probesToFilter.add(probeName);
					System.out.println("Filtering " + probeName);

				}
				else keepProbes.put(probeName, Boolean.valueOf(true));
			}
		}
		keepProbes.clear();
		
		// Check ribosomal RNA imperfect matches to transcriptome
		FileReader readNi = new FileReader(gImperfect);
		System.out.println("Reading imperfect non-rRNA matches to genome from file " + gImperfect);
		BufferedReader bni = new BufferedReader(readNi);
		while(bni.ready()) {
			String line = bni.readLine();
			stringparse.parse(line);
			if(stringparse.getFieldCount() != 21) continue;
			String probeName = stringparse.asString(9);
			if(probesToFilter.contains(probeName)) continue;
			if(hasHyb(stringparse.asString(18),60)) {
				if(keepProbes.containsKey(probeName)) {
					keepProbes.put(probeName, Boolean.valueOf(false));
					probesToFilter.add(probeName);
					System.out.println("Filtering " + probeName);

				}
				else keepProbes.put(probeName, Boolean.valueOf(true));
			}
		}
		keepProbes.clear();
		
		
		List<Sequence> probesToKeep = new ArrayList<Sequence>();
		List<String> probesToRemove = new ArrayList<String>();
		
		//FileWriter or = new FileWriter(outRemoved);
				
		// Write the kept probes to fasta file
		// Write the removed probes to fasta file
		//System.out.println("Writing unique probes to file " + outFasta);
		//System.out.println("Writing removed probes to file " + outRemoved);
		// First read in all the probes
		FastaSequenceIO fsio = new FastaSequenceIO(probeFasta);
		List<Sequence> allprobes = fsio.loadAll();
		// now check if they are marked for filtering
		for(Sequence probe : allprobes) {
			if(!probesToFilter.contains(probe.getId())) {
				probesToKeep.add(probe);
			} else {
				probesToRemove.add(probe.getSequenceBases());
				//or.write(probe.getSequenceBases() + "\n");
			}
		}
		//fsio.write(probesToKeep,outFasta);
	
		FileWriter odf = new FileWriter(outDesignFasta);
		FileWriter odt = new FileWriter(outDesignTable);
		
		// Filter the full design fasta and table
		System.out.println("Writing filtered full probe sequences to file " + outDesignFasta);
		FastaSequenceIO fsiodesign = new FastaSequenceIO(designFasta);
		List<Sequence> fullProbeSeqs = fsiodesign.loadAll();
		for(Sequence s : fullProbeSeqs) {
			boolean matches = false;
			for(String ss : probesToRemove) {
				if(s.contains(ss) && s.getId().contains("Pulldown")) matches = true;
			}
			if(!matches) {
				odf.write(">" + s.getId() + "\n" + s.getSequenceBases() + "\n");
			}
		}
		odf.close();
		
		StringParser lineparse = new StringParser();
		System.out.println("Writing filtered design to file " + outDesignTable);
		FileReader dr = new FileReader(designTable);
		BufferedReader br = new BufferedReader(dr);
		while(br.ready()) {
			boolean matches = false;
			String line = br.readLine();
			lineparse.parse(line);
			String fullprobe = lineparse.asString(19);
			for(String ss : probesToRemove) {
				if(fullprobe.contains(ss) && lineparse.asString(3).equals("Pulldown")) matches = true;
			}
			if(!matches) {
				odt.write(line + "\n");
			}
		}
		odt.close();
		System.out.println("All done.");
	
	}

}

/* 
 * Jesse Engreitz
 * March 7, 2012
 * Given an input BED file containing genomic regions, output a BED file containing the locations of 
 * all CpG nucleotides contained in these regions.
 */

package broad.pda.methylation;


import java.io.IOException;
import java.io.File;
import java.io.FileWriter;
import java.util.Collection;
import java.util.List;
import java.util.ArrayList;

import net.sf.picard.reference.*;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.error.ParseException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.SequenceRegion;
import broad.pda.annotation.BEDFileParser;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

import broad.pda.methylation.MethylationUtils;

public class WriteCpGBED {
	
	public WriteCpGBED(String regionFile, String reference) throws IOException, ParseException {
		FastaSequenceIO fsio = new FastaSequenceIO(reference);
		
		//FastaSequenceIndex index = new FastaSequenceIndex(new File(referenceIndex));
		//IndexedFastaSequenceFile ref = new IndexedFastaSequenceFile(new File(reference), index);

		// Load the region file
		BEDReader bed = new BEDReader(regionFile);
		List<BED> regions = bed.getAnnotationList();
		
		// Get the genomic sequence for each region
		List<SequenceRegion> sequences = new ArrayList<SequenceRegion>();
		for (BED region : regions) {
			SequenceRegion newRegion = new SequenceRegion(region.getChromosome(), region);
			newRegion.setName(region.getName());
			sequences.add(newRegion);
		}
		
		fsio.extractRegions(sequences);
		
		// Find CpG sites in each region
		List<BED> cpgBed = new ArrayList<BED>();
		for (SequenceRegion seq : sequences) {
			//System.out.println(seq.getSequenceBases());
			List<Integer> cpgIndices = MethylationUtils.findCpGs(seq.getSequenceBases());
			for (Integer startIndex : cpgIndices) {
				String cpgName = seq.getName() + "." + startIndex;
				cpgBed.add(new BED(cpgName, seq.getChromosome(), seq.getStart() + startIndex - 1, seq.getStart() + startIndex + 1));
			}
		}
		
		// Write the output BED file
		for (BED cpg : cpgBed) 
			System.out.println(cpg);
	}

	
	private static String USAGE = "WriteCpGBED parameters:\n\n\t-in\tRegion file in BED format\n\t-ref\tBasename for reference FASTA file and index\n";
	
	public static void main (String[] args) throws IOException, ParseException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "write");
		
		final String in = argmap.getInput();
		final String ref = argmap.get("ref");
		
		new WriteCpGBED(in, ref);
	}
	
	
}

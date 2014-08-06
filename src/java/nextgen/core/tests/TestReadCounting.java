package nextgen.core.tests;

import broad.core.parser.CommandLineParser;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import broad.pda.annotation.BEDFileParser;

public class TestReadCounting {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-a", "Bam file of alignments", true);
		p.addStringArg("-b", "Bed gene annotation", true);
		p.addStringArg("-g", "Gene to test", true);
		p.parse(args);
		String bamFile = p.getStringArg("-a");
		String bedFile = p.getStringArg("-b");
		String geneName = p.getStringArg("-g");
		
		File file = new File(bedFile);
		Map<String, Gene> genesByName = BEDFileParser.loadDataByName(file);
		TranscriptomeSpace t = new TranscriptomeSpace(BEDFileParser.loadDataByChr(file));
		
		AlignmentModel transcriptomeAlignments = new AlignmentModel(bamFile, t);
		transcriptomeAlignments.addFilter(new FragmentLengthFilter(t,2000));
		transcriptomeAlignments.addFilter(new GenomicSpanFilter(300000));
		
		System.err.println("Computing count for gene " + geneName);
		Gene gene = genesByName.get(geneName);
		long start = System.currentTimeMillis();
		double numReads = transcriptomeAlignments.getCount(gene);
		long stop = System.currentTimeMillis();
		System.err.println("Read count: " + numReads);
		long time = stop - start;
		System.err.println("Time: " + time);
		
	}

}

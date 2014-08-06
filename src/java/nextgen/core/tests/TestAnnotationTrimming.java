/**
 * 
 */
package nextgen.core.tests;

import broad.core.parser.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.protection.TwoSampleScanSkellamPeakCaller;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;

/**
 * @author prussell
 *
 */
public final class TestAnnotationTrimming {
	
	private TranscriptomeSpace transcriptomeSpace;
	private Map<String, Collection<Gene>> genesByChr;
	private static Logger logger = Logger.getLogger(TestAnnotationTrimming.class.getName());
	
	private TestAnnotationTrimming(String bedFile) throws IOException {
		genesByChr = BEDFileParser.loadDataByChr(new File(bedFile));
		transcriptomeSpace = new TranscriptomeSpace(genesByChr);
	}
	
	private void testTrimming(String outBedFile) throws IOException {
		FileWriter w = new FileWriter(outBedFile);
		for(String chr : genesByChr.keySet()) {
			for(Gene gene : genesByChr.get(chr)) {
				if(gene.getSize() < 70) continue;
				logger.info("Trimming gene " + gene.getName());
				Gene trimStart = new Gene(gene);
				Gene trimEnd = new Gene(gene);
				Gene trimBoth = new Gene(gene);
				trimStart.setName(gene.getName() + "_trim_start_50");
				trimEnd.setName(gene.getName() + "_trim_end_30");
				trimBoth.setName(gene.getName() + "_trim_start_40_end_20");
				trimStart = new Gene(trimStart.trim(50, 0));
				trimEnd.trim(0, 30);
				trimBoth.trim(40, 20);
				w.write(trimStart.toBED() + "\n");
				w.write(trimEnd.toBED() + "\n");
				w.write(trimBoth.toBED() + "\n");
			}
		}
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file", true);
		p.addStringArg("-o", "Output bed file", true);
		p.parse(args);
		String bedFile =  p.getStringArg("-b");
		String outBedFile = p.getStringArg("-o");
		TestAnnotationTrimming tat = new TestAnnotationTrimming(bedFile);
		tat.testTrimming(outBedFile);
		
	}

}

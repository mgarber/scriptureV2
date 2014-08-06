/**
 * 
 */
package nextgen.core.utils;

import broad.core.parser.CommandLineParser;

import java.io.File;
import java.io.IOException;

import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class AlignmentModelGlobalStats {

	private static Logger logger = Logger.getLogger(AlignmentModelGlobalStats.class.getName());
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-t", "Bed file for transcriptome space stats", false, null);
		p.addStringArg("-g", "Chromosome size file for genomic space stats", false, null);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-t");
		String sizeFile = p.getStringArg("-g");
		
		if(bedFile != null) {
			logger.info("Creating transcriptome space from annotation in " + bedFile);
			TranscriptomeSpace t = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(bedFile)));
			logger.info("Calculating transcriptome space stats for bam file " + bamFile);
			AlignmentModel a = new AlignmentModel(bamFile, t);
			a.computeGlobalStats();
			logger.info("Done calculating transcriptome space stats.");
		}
		
		if(sizeFile != null) {
			logger.info("Creating genomic space from chromosome sizes in " + sizeFile);
			GenomicSpace g = new GenomicSpace(sizeFile);
			logger.info("Calculating genomic space stats for bam file " + bamFile);
			AlignmentModel a = new AlignmentModel(bamFile, g);
			a.computeGlobalStats();
			logger.info("Done calculating genomic space stats.");
		}
		
	}

}

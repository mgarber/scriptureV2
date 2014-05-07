package nextgen.core.pipeline.util;

import broad.core.parser.CommandLineParser;

public class BamFileSplitter {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-i", "Input bam file", true);
		p.addIntArg("-n", "Number of files to split into", true);
		p.parse(args);
		String bam = p.getStringArg("-i");
		int n = p.getIntArg("-n");
		
		BamUtils.splitBam(bam, n);

	}

}

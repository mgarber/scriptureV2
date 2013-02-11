/**
 * 
 */
package broad.pda.capture.designer;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.PipelineUtils;

/**
 * @author prussell
 *
 */
public final class ArrayDesignerInSilicoPcrTest {

	/**
	 * Constructor stores location of array design produced by ArrayDesigner
	 * @param designFile the output design table produced by ArrayDesigner
	 * @param probeFile the output fasta file of full probe sequences produced by ArrayDesigner
	 */
	public ArrayDesignerInSilicoPcrTest(String designfile, String probefile) {
		this.designFile = designFile;
		this.probeFile = probefile;
		this.tempIsPcrOutFile = "isPcrOut.fa";
		this.tempIsPcrPrimerInput = "isPcrPrimerInput";
		this.okPrimerPairs = new HashSet<String>();
		this.noMatchPrimerPairs = new HashSet<String>();
	}
	
	/**
	 * 
	 * @param designFileLine a line from the output design file produced by ArrayDesigner
	 * @param leftPrimer the left primer
	 * @param rightPrimer the right primer
	 * @param outfile output fasta file for isPCR
	 * @param checkSpecies check the species name for non specificity
	 * @param checkGene check the gene name for non specificity
	 * @param checkClass check the RNA class for non specificity
	 * @param checkPurpose check the probe purpose for non specificity
	 * @param checkTilingPath check the tiling path for non specificity
	 * @param checkEvenOdd check the even/odd status for non specificity
	 * @param checkDomain check the domain number for non specificity
	 * @param checkContainsQpcrPrimer check the qPCR primer status for non specificity
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private void runIsPcr(String designFileLine, String leftPrimer, String rightPrimer, String outfile, boolean checkSpecies,
			boolean checkGene, boolean checkClass, boolean checkPurpose, boolean checkTilingPath, boolean checkEvenOdd,
			boolean checkDomain, boolean checkContainsQpcrPrimer) throws IOException, InterruptedException {

		String bothPrimers = leftPrimer + rightPrimer;
		if(this.okPrimerPairs.contains(bothPrimers)) {
			System.out.println("PRIMER_PAIR_OK\t" + leftPrimer + "\t" + rightPrimer + "\t" + designFileLine);
			return;
		}
		if(this.noMatchPrimerPairs.contains(bothPrimers)) {
			System.out.println("NO_MATCHES\tLEFT_PRIMER\t" + leftPrimer);
			System.out.println("NO_MATCHES\tRIGHT_PRIMER\t" + rightPrimer);
			System.out.println("NO_MATCHES\tDESIGN\t" + designFileLine);
			return;
		}		
		
		// Write the primer pair file for input to isPCR
		FileWriter writer = new FileWriter(this.tempIsPcrPrimerInput);
		writer.write("name\t" + leftPrimer + "\t" + rightPrimer + "\n");
		writer.close();
		
		String jobID = "job__" + Long.valueOf(System.currentTimeMillis()).toString();
		String cmmd = "/seq/mguttman/scripts/isPCR/isPcr " + this.probeFile + " " + this.tempIsPcrPrimerInput + " " + outfile;
		
		// Submit the job
		int exitCode = PipelineUtils.bsubSmallProcess(Runtime.getRuntime(), jobID , cmmd , "bsub_output_ispcr_test");
		PipelineUtils.waitForJobs(jobID, Runtime.getRuntime(), 10, false);
		
		boolean foundNonSpecificMatch = false;
		boolean noMatches = false;
		
		// Parse the output
		FastaSequenceIO fsio = new FastaSequenceIO(outfile);
		List<Sequence> matches = fsio.loadAll();
		if(matches.size() == 0) {
			System.out.println("NO_MATCHES\tLEFT_PRIMER\t" + leftPrimer);
			System.out.println("NO_MATCHES\tRIGHT_PRIMER\t" + rightPrimer);
			System.out.println("NO_MATCHES\tDESIGN\t" + designFileLine);
			this.noMatchPrimerPairs.add(bothPrimers);
			noMatches = true;
			
		}
		for(Sequence match : matches) {
			String name = match.getId();
			if(matchIsNonSpecific(designFileLine, name, checkSpecies, checkGene, checkClass, checkPurpose, checkTilingPath, checkEvenOdd, checkDomain, checkContainsQpcrPrimer)) {
				System.out.println("NON_SPECIFIC_MATCH\tLEFT_PRIMER\t" + leftPrimer);
				System.out.println("NON_SPECIFIC_MATCH\tRIGHT_PRIMER\t" + rightPrimer);
				System.out.println("NON_SPECIFIC_MATCH\tFROM_PROBE\t" + designFileLine);
				System.out.println("NON_SPECIFIC_MATCH\tMATCHES_PROBE\t" + name);
				foundNonSpecificMatch = true;
			} 
		}
		
		if(!foundNonSpecificMatch && !noMatches) {
			this.okPrimerPairs.add(bothPrimers);
			System.out.println("PRIMER_PAIR_OK\t" + leftPrimer + "\t" + rightPrimer + "\t" + designFileLine);
		}

	}
	
	/**
	 * Check primer pairs for a probe for non specific priming in full probe set
	 * @param designFileLine a line from the output design file produced by ArrayDesigner
	 * @param outfile a file to write in silico PCR output to
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	private void checkAllPossiblePrimerPairs(String designFileLine, String outfile) throws IOException, InterruptedException {
		
		StringParser stringparse = new StringParser();
		stringparse.parse(designFileLine);
		
		String probePurpose = stringparse.asString(3);
		String speciesClassLeftPrimer = stringparse.asString(10);
		String speciesClassRightPrimer = stringparse.asString(11);
		String geneLeftPrimer = stringparse.asString(12);
		String geneRightPrimer = stringparse.asString(13);
		String tilingPathLeftPrimer = stringparse.asString(14);
		String evenOddLeftPrimer = stringparse.asString(15);
		String qpcrRightPrimer = stringparse.asString(16);
		String domainRightPrimer = stringparse.asString(17);
		
		if(probePurpose.equals("Pulldown")) {
			
			// Gene left, gene right
			this.runIsPcr(designFileLine, geneLeftPrimer, geneRightPrimer, this.tempIsPcrOutFile, false, true, false, false, false, false, false, false);
			// Gene left, qpcr right
			this.runIsPcr(designFileLine, geneLeftPrimer, qpcrRightPrimer, this.tempIsPcrOutFile, false, true, false, false, false, false, false, true);
			// Gene left, domain right
			this.runIsPcr(designFileLine, geneLeftPrimer, domainRightPrimer, this.tempIsPcrOutFile, false, true, false, false, false, false, true, false);
			// Tiling path left, gene right
			this.runIsPcr(designFileLine, tilingPathLeftPrimer, geneRightPrimer, this.tempIsPcrOutFile, false, true, false, false, true, false, false, false);
			// Even odd left, gene right
			this.runIsPcr(designFileLine, evenOddLeftPrimer, geneRightPrimer, this.tempIsPcrOutFile, false, true, false, false, false, true, false, false);
			// Tiling path left, qpcr right
			this.runIsPcr(designFileLine, tilingPathLeftPrimer, qpcrRightPrimer, this.tempIsPcrOutFile, false, true, false, false, true, false, false, true);
			// Tiling path left, domain right
			this.runIsPcr(designFileLine, tilingPathLeftPrimer, domainRightPrimer, this.tempIsPcrOutFile, false, true, false, false, true, false, true, false);
			// Even odd left, qpcr right
			this.runIsPcr(designFileLine, evenOddLeftPrimer, qpcrRightPrimer, this.tempIsPcrOutFile, false, true, false, false, false, true, false, true);
			// Even odd left, domain right
			this.runIsPcr(designFileLine, evenOddLeftPrimer, domainRightPrimer, this.tempIsPcrOutFile, false, true, false, false, false, true, true, false);

		}
		
		if(probePurpose.contains("Depletion")) {
			
			// Species class left, species class right
			this.runIsPcr(designFileLine, speciesClassLeftPrimer, speciesClassRightPrimer, this.tempIsPcrOutFile, true, false, true, false, false, false, false, false);
			// Gene left, gene right
			this.runIsPcr(designFileLine, geneLeftPrimer, geneRightPrimer, this.tempIsPcrOutFile, false, true, false, false, false, false, false, false);
			
		}
		
	}
	
	/**
	 * Check if an in silico PCR match is non specific
	 * @param designFileLine a line from the output design file produced by ArrayDesigner
	 * @param ispcrFastaName the fasta header of a match written by isPCR
	 * @param checkSpecies check the species name for non specificity
	 * @param checkGene check the gene name for non specificity
	 * @param checkClass check the RNA class for non specificity
	 * @param checkPurpose check the probe purpose for non specificity
	 * @param checkTilingPath check the tiling path for non specificity
	 * @param checkEvenOdd check the even/odd status for non specificity
	 * @param checkDomain check the domain number for non specificity
	 * @param checkContainsQpcrPrimer check the qPCR primer status for non specificity
	 * @return
	 */
	private static boolean matchIsNonSpecific(String designFileLine, String ispcrFastaName, boolean checkSpecies,
			boolean checkGene, boolean checkClass, boolean checkPurpose, boolean checkTilingPath, boolean checkEvenOdd,
			boolean checkDomain, boolean checkContainsQpcrPrimer) {
		
		StringParser stringparse = new StringParser();
		
		// Probe attributes from array design output
		stringparse.parse(designFileLine);
		String designSpecies = stringparse.asString(0);
		String designGene = stringparse.asString(1);
		String designClass = stringparse.asString(2);
		String designPurpose = stringparse.asString(3);
		String designTilingPath = stringparse.asString(6);
		String designEvenOdd = stringparse.asString(7);
		String designDomain = stringparse.asString(8);
		String designContainsQpcr = stringparse.asString(9);
		
		// Probe attributes from isPCR match
		stringparse.parse(ispcrFastaName,"_");
		String ispcrSpecies = stringparse.asString(0);
		
		// Gene names containing underscore
		if(stringparse.asString(1).equals("NM") || stringparse.asString(1).equals("Human")) {
			String ispcrGene = stringparse.asString(1) + "_" + stringparse.asString(2);
			String ispcrPurpose = stringparse.asString(3);
			String ispcrClass = stringparse.asString(4);
			String ispcrTilingPathLong = stringparse.asString(7);
			String ispcrEvenOddLong = stringparse.asString(8);
			String ispcrDomainLong = stringparse.asString(9);
			String ispcrContainsQpcrLong = stringparse.asString(10);
			// Further parse fasta name from isPCR match
			stringparse.parse(ispcrTilingPathLong,":");
			String ispcrTilingPath = stringparse.asString(1);
			stringparse.parse(ispcrEvenOddLong,":");
			String ispcrEvenOdd = stringparse.asString(1);
			stringparse.parse(ispcrDomainLong,":");
			String ispcrDomain = stringparse.asString(1);
			stringparse.parse(ispcrContainsQpcrLong,":");
			String ispcrContainsQpcr = stringparse.asString(1);
			// Check requested attributes for specificity
			if(checkSpecies && !designSpecies.equals(ispcrSpecies)) return true;
			if(checkGene && !designGene.equals(ispcrGene)) return true;
			if(checkClass && !designClass.equals(ispcrClass)) return true;
			if(checkPurpose && !designPurpose.equals(ispcrPurpose)) return true;
			if(checkTilingPath && !designTilingPath.equals(ispcrTilingPath)) return true;
			if(checkEvenOdd && !designEvenOdd.equals(ispcrEvenOdd)) return true;
			if(checkDomain && !designDomain.equals(ispcrDomain)) return true;
			if(checkContainsQpcrPrimer && !designContainsQpcr.equals(ispcrContainsQpcr)) return true;
			
		} else {
			String ispcrGene = stringparse.asString(1);
			String ispcrPurpose = stringparse.asString(2);
			String ispcrClass = stringparse.asString(3);
			String ispcrTilingPathLong = stringparse.asString(6);
			String ispcrEvenOddLong = stringparse.asString(7);
			String ispcrDomainLong = stringparse.asString(8);
			String ispcrContainsQpcrLong = stringparse.asString(9);
			// Further parse fasta name from isPCR match
			stringparse.parse(ispcrTilingPathLong,":");
			String ispcrTilingPath = stringparse.asString(1);
			stringparse.parse(ispcrEvenOddLong,":");
			String ispcrEvenOdd = stringparse.asString(1);
			stringparse.parse(ispcrDomainLong,":");
			String ispcrDomain = stringparse.asString(1);
			stringparse.parse(ispcrContainsQpcrLong,":");
			String ispcrContainsQpcr = stringparse.asString(1);
			// Check requested attributes for specificity
			if(checkSpecies && !designSpecies.equals(ispcrSpecies)) return true;
			if(checkGene && !designGene.equals(ispcrGene)) return true;
			if(checkClass && !designClass.equals(ispcrClass)) return true;
			if(checkPurpose && !designPurpose.equals(ispcrPurpose)) return true;
			if(checkTilingPath && !designTilingPath.equals(ispcrTilingPath)) return true;
			if(checkEvenOdd && !designEvenOdd.equals(ispcrEvenOdd)) return true;
			if(checkDomain && !designDomain.equals(ispcrDomain)) return true;
			if(checkContainsQpcrPrimer && !designContainsQpcr.equals(ispcrContainsQpcr)) return true;
		}

		return false;
		
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {

		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("Use in silico PCR to check for non-specific priming.");
		p.addStringArg("-p", "Fasta file of full probe sequences generated by ArrayDesigner", true);
		p.addStringArg("-d", "Full design table generated by ArrayDesigner", true);
		p.parse(args);
		String probeFile = p.getStringArg("-p");
		String designFile = p.getStringArg("-d");
		
		ArrayDesignerInSilicoPcrTest tester = new ArrayDesignerInSilicoPcrTest(designFile, probeFile);
		
		FileReader r = new FileReader(designFile);
		BufferedReader b = new BufferedReader(r);
		
		while(b.ready()) {
			String line = b.readLine();
			tester.checkAllPossiblePrimerPairs(line, tester.tempIsPcrOutFile);
		}
		
		
		
	}

	private String designFile;
	private String probeFile;
	private String tempIsPcrOutFile;
	private String tempIsPcrPrimerInput;
	private HashSet<String> okPrimerPairs;
	private HashSet<String> noMatchPrimerPairs;
	
}

package broad.pda.capture.designer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.TimeoutException;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.Primer3Configuration;
import broad.core.primer3.Primer3ConfigurationFactory;
import broad.core.primer3.Primer3IO;
import broad.core.primer3.Primer3SequenceInputTags;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.PrimerUtils;
import broad.core.primer3.qPCRPrimerDesigner;
import broad.core.sequence.Sequence;
import broad.core.util.PipelineUtils;

public final class PcrTailDesigner {
	
	/*
	 * Process:
	 * 
	 * Store the length of the primers
	 * Store the substructure of the two primers
	 * Store the marginal number of different primers in each category
	 * 
	 * Do this for each of 5x the desired number of instances of the outermost primer category
	 * 		Choose a random sequence
	 * 		Ask primer3 for the best primer pair in the sequence (full length left and right primers)
	 * 		Mutate parts of the left and right primers that will encode the next level of information
	 * 		Get the primer pair penalty for each mutated primer pair
	 * 		Keep the best 5x the desired number of instances of the inner category
	 * 
	 * If there is a third category:
	 * 		Repeat the above, mutating the part of the sequences generated in previous step that is unique to the third category
	 * 
	 * Filter primers that cross-hybridize to each other
	 * 			
	 * If any outer category did not keep enough inner categories, make new outer primers and repeat
	 * 
	 * Return the desired number of outer+inner categories
	 * 
	 */
		
	/**
	 *  Constructor with no parameters is private
	 */
	private PcrTailDesigner() {
		this.primerSize = Integer.valueOf(-1);
		
		this.leftOuterLength = Integer.valueOf(-1);
		this.leftOuterStart = Integer.valueOf(0);
		this.leftOuterEnd = Integer.valueOf(-1);
		this.leftSecondLength = Integer.valueOf(-1);
		this.leftSecondStart = Integer.valueOf(-1);
		this.leftSecondEnd = Integer.valueOf(-1);
		this.leftThirdLength = Integer.valueOf(-1);
		this.leftThirdStart = Integer.valueOf(-1);
		this.leftThirdEnd = Integer.valueOf(-1);

		
		this.rightOuterLength = Integer.valueOf(-1);
		this.rightOuterStart = Integer.valueOf(0);
		this.rightOuterEnd = Integer.valueOf(-1);
		this.rightSecondLength = Integer.valueOf(-1);
		this.rightSecondStart = Integer.valueOf(-1);
		this.rightSecondEnd = Integer.valueOf(-1);
		this.rightThirdLength = Integer.valueOf(-1);
		this.rightThirdStart = Integer.valueOf(-1);
		this.rightThirdEnd = Integer.valueOf(-1);
		
		this.outerCount = Integer.valueOf(-1);
		this.secondCount = Integer.valueOf(-1);
		this.thirdCount = Integer.valueOf(-1);
		
		this.threeSegments = false;
		this.verbose = false;
		this.veryVerbose = false;
		
		initialOuterPrimerPairs = new ArrayList<PrimerPair>();
		this.secondPrimerPairs = new HashMap<PrimerPair, ArrayList<PrimerPair>>();
		this.thirdPrimerPairs = new HashMap< PrimerPair, Map< PrimerPair, ArrayList<PrimerPair> > >();
	
		this.projectName = "";
		
	}

	/**
	 *  Constructor when there are two overlapping primers
	 * @param primerSize
	 * @param leftOuterLength
	 * @param leftInnerLength
	 * @param rightOuterLength
	 * @param rightInnerLength
	 * @param outerCount
	 * @param innerCount
	 */
	public PcrTailDesigner(int primerSize, int leftOuterLength, int leftInnerLength, int rightOuterLength, int rightInnerLength, int outerCount, int innerCount, String projectName) {
		this.primerSize = Integer.valueOf(primerSize);
		this.leftOuterLength = Integer.valueOf(leftOuterLength);
		this.leftSecondLength = Integer.valueOf(leftInnerLength);
		this.rightOuterLength = Integer.valueOf(rightOuterLength);
		this.rightSecondLength = Integer.valueOf(rightInnerLength);
		this.outerCount = Integer.valueOf(outerCount);
		this.secondCount = Integer.valueOf(innerCount);
		this.leftOuterStart = Integer.valueOf(0);
		this.leftOuterEnd = Integer.valueOf(leftOuterLength - 1);
		this.leftSecondStart = Integer.valueOf(primerSize - leftInnerLength);
		this.leftSecondEnd = Integer.valueOf(primerSize - 1);
		this.rightOuterStart = Integer.valueOf(0);
		this.rightOuterEnd = Integer.valueOf(rightOuterLength - 1);
		this.rightSecondStart = Integer.valueOf(primerSize - rightInnerLength);
		this.rightSecondEnd = Integer.valueOf(primerSize - 1);
		this.threeSegments = Boolean.valueOf(false);
		
		// So methods will throw an exception if a method tries to use the nonexistent third nested primer
		this.leftThirdStart = Integer.valueOf(-1);
		this.leftThirdEnd = Integer.valueOf(-1);
		this.leftThirdLength = Integer.valueOf(-1);
		this.rightThirdStart = Integer.valueOf(-1);
		this.rightThirdEnd = Integer.valueOf(-1);
		this.rightThirdLength = Integer.valueOf(-1);
		
		this.threeSegments = false;
		this.verbose = false;
		this.veryVerbose = false;
		
		this.initialOuterPrimerPairs = new ArrayList<PrimerPair>();
		this.secondPrimerPairs = new HashMap<PrimerPair, ArrayList<PrimerPair>>();
		this.thirdPrimerPairs = new HashMap< PrimerPair, Map< PrimerPair, ArrayList<PrimerPair> > >();
		
		this.projectName = projectName;
		
		this.printParameters();
	}
	
	/**
	 * Constructor when there are three overlapping primers
	 * @param primerSize
	 * @param leftOuterLength
	 * @param leftMiddleLength
	 * @param leftMiddleStart
	 * @param leftInnerLength
	 * @param rightOuterLength
	 * @param rightMiddleLength
	 * @param rightMiddleStart
	 * @param rightInnerLength
	 * @param outerCount
	 * @param middleCount
	 * @param innerCount
	 */
	public PcrTailDesigner(int primerSize, int leftOuterLength, int leftMiddleLength, int leftMiddleStart, int leftInnerLength, int rightOuterLength, int rightMiddleLength, int rightMiddleStart, int rightInnerLength, int outerCount, int middleCount, int innerCount, String projectName) {
		this.primerSize = Integer.valueOf(primerSize);
		this.leftOuterLength = Integer.valueOf(leftOuterLength);
		this.leftOuterStart = Integer.valueOf(0);
		this.leftOuterEnd = Integer.valueOf(leftOuterLength - 1);
		this.leftSecondLength = Integer.valueOf(leftMiddleLength);
		this.leftSecondStart = Integer.valueOf(leftMiddleStart);
		this.leftSecondEnd = Integer.valueOf(leftMiddleStart + leftMiddleLength - 1);
		this.leftThirdLength = Integer.valueOf(leftInnerLength);
		this.leftThirdStart = Integer.valueOf(primerSize - leftInnerLength);
		this.leftThirdEnd = Integer.valueOf(primerSize - 1);
		this.rightOuterLength = Integer.valueOf(rightOuterLength);
		this.rightOuterStart = Integer.valueOf(0);
		this.rightOuterEnd = Integer.valueOf(rightOuterLength - 1);
		this.rightSecondLength = Integer.valueOf(rightMiddleLength);
		this.rightSecondStart = Integer.valueOf(rightMiddleStart);
		this.rightSecondEnd = Integer.valueOf(rightMiddleStart + rightMiddleLength - 1);
		this.rightThirdLength = Integer.valueOf(rightInnerLength);
		this.rightThirdStart = Integer.valueOf(primerSize - rightInnerLength);
		this.rightThirdEnd = Integer.valueOf(primerSize - 1);
		this.outerCount = Integer.valueOf(outerCount);
		this.secondCount = Integer.valueOf(middleCount);
		this.thirdCount = Integer.valueOf(innerCount);
		this.threeSegments = Boolean.valueOf(true);
		
		this.verbose = Boolean.valueOf(false);
		this.veryVerbose = Boolean.valueOf(false);
		
		this.initialOuterPrimerPairs = new ArrayList<PrimerPair>();
		this.secondPrimerPairs = new HashMap<PrimerPair, ArrayList<PrimerPair>>();
		this.thirdPrimerPairs = new HashMap< PrimerPair, Map< PrimerPair, ArrayList<PrimerPair> > >();
		
		this.printParameters();
		
		this.projectName = projectName;
		
	}
	
	/**
	 * Get command line arguments and check for validity, set attributes
	 * @param args the command line
	 */
	private void getCommandArgs(String[] args) {
		// Add command arguments
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("Generate primers and write to a file.");
		p.addIntegerArg("-size", "Primer size", true);
		p.addIntegerArg("-l1length", "Length of outermost segment of left primer", true);
		p.addIntegerArg("-l2length", "Length of second segment of left primer", true);
		p.addIntegerArg("-l3length", "Length of third segment of left primer", false, Integer.valueOf(-1));
		p.addIntegerArg("-l2start", "Zero-based start position of left second segment if there are three segments", false, Integer.valueOf(-1));
		p.addIntegerArg("-r1length", "Length of outermost segment of right primer", true);
		p.addIntegerArg("-r2length", "Length of second segment of right primer", true);
		p.addIntegerArg("-r3length", "Length of third segment of right primer", false, Integer.valueOf(-1));
		p.addIntegerArg("-r2start", "Zero-based start position of right second segment if there are three segments", false, Integer.valueOf(-1));
		p.addIntegerArg("-count1", "Desired number of outermost primer", true);
		p.addIntegerArg("-count2", "Desired number of second primer", true);
		p.addIntegerArg("-count3", "Desired number of third primer", false, Integer.valueOf(-1));
		p.addBooleanArg("-v", "Verbose", false, false);
		p.addBooleanArg("-vv", "Very verbose", false, false);
		
		p.parse(args);
		
		this.primerSize = p.getIntegerArg("-size");
		this.leftOuterLength = p.getIntegerArg("-l1length");
		this.outerCount = p.getIntegerArg("-count1");
		this.leftSecondLength = p.getIntegerArg("-l2length");
		this.secondCount = p.getIntegerArg("-count2");
		this.leftThirdLength = p.getIntegerArg("-l3length");
		this.thirdCount = p.getIntegerArg("-count3");
		this.leftSecondStart = p.getIntegerArg("-l2start");
		this.rightOuterLength = p.getIntegerArg("-r1length");
		this.rightSecondLength = p.getIntegerArg("-r2length");
		this.rightThirdLength = p.getIntegerArg("-r3length");
		this.rightSecondStart = p.getIntegerArg("-r2start");
		this.verbose = p.getBooleanArg("-v");
		this.veryVerbose = p.getBooleanArg("-vv");
		
		// If there are two segments on each side
		if(this.leftThirdLength.intValue() == -1 && this.rightThirdLength.intValue() == -1) {
			this.leftSecondStart = Integer.valueOf(this.primerSize.intValue() - this.leftSecondLength.intValue());
			this.leftSecondEnd = Integer.valueOf(this.primerSize.intValue() - 1);
			this.leftOuterEnd = Integer.valueOf(this.leftOuterLength.intValue() - 1);
			this.rightSecondStart = Integer.valueOf(this.primerSize.intValue() - this.rightSecondLength.intValue());
			this.rightSecondEnd = Integer.valueOf(this.primerSize.intValue() - 1);
			this.rightOuterEnd = Integer.valueOf(this.rightOuterLength.intValue() - 1);
		} else {
			// Check that all necessary information is present to make three segments
			// If not, throw exception
			// Set threeSegments to true
			
			
			
			
		}
		
		this.printParameters();

	}
	
	/**
	 * Print parameters to stdout
	 */
	private void printParameters() {
		if(!this.threeSegments) {
			System.out.println("There are two segments per primer.");
			System.out.println("Primer size: " + this.primerSize + "\n");
			System.out.println("Left outer segment length/start/end: " + this.leftOuterLength + "/" + this.leftOuterStart + "/" + this.leftOuterEnd);
			System.out.println("Left inner segment length/start/end: " + this.leftSecondLength + "/" + this.leftSecondStart + "/" + this.leftSecondEnd);
			System.out.println("\nRight outer segment length/start/end: " + this.rightOuterLength + "/" + this.rightOuterStart + "/" + this.rightOuterEnd);
			System.out.println("Right inner segment length/start/end: " + this.rightSecondLength + "/" + this.rightSecondStart + "/" + this.rightSecondEnd);
			System.out.println("\nNumber of primers (outer/inner): " + this.outerCount + "/" + this.secondCount);
		} else {
			System.out.println("There are three segments per primer.");
			System.out.println("Primer size: " + this.primerSize + "\n");
			System.out.println("Left outer segment length/start/end: " + this.leftOuterLength + "/" + this.leftOuterStart + "/" + this.leftOuterEnd);
			System.out.println("Left middle segment length/start/end: " + this.leftSecondLength + "/" + this.leftSecondStart + "/" + this.leftSecondEnd);
			System.out.println("Left inner segment length/start/end: " + this.leftThirdLength + "/" + this.leftThirdStart + "/" + this.leftThirdEnd);
			System.out.println("Right outer segment length/start/end: " + this.rightOuterLength + "/" + this.rightOuterStart + "/" + this.rightOuterEnd);
			System.out.println("Right middle segment length/start/end: " + this.rightSecondLength + "/" + this.rightSecondStart + "/" + this.rightSecondEnd);
			System.out.println("Right inner segment length/start/end: " + this.rightThirdLength + "/" + this.rightThirdStart + "/" + this.rightThirdEnd);
			System.out.println("\nNumber of primers (outer/middle/inner): " + this.outerCount + "/" + this.secondCount + "/" + this.thirdCount);
		}
		System.out.println();
	}
	
	/**
	 * Pick random sequences until a good primer pair is found then write the primer to the file
	 * This method is intended to be used in parallel to write many primers at once
	 * @param primerLength
	 * @param outFile
	 * @throws IOException
	 */
	private static void writeOneInitialPrimerPair(int primerLength, String outFile) throws IOException {
		
		// Get the name of the output file to also use as the primer ID
		StringParser stringparse = new StringParser();
		stringparse.parse(outFile,"/");
		String filename = stringparse.asString(stringparse.getFieldCount() - 1);
		
		// Repeat until a suitable primer pair is found
		boolean found = false;
		while(found == false) {
			String seq = Sequence.generateRandomSequence(5000);
			// Only ask for one primer pair
			Collection<PrimerPair> primers = qPCRPrimerDesigner.designSyntheticPrimers(seq, 1, primerLength, 5000);
			// A primer pair was found
			if(primers != null && !primers.isEmpty()) {
				PrimerPair primer = primers.iterator().next();
				if(primer.getPrimerPairPenalty() <= MAX_PRIMER_PENALTY) {
					found = true;
					// Write primer data to file to be read in by another method
					FileWriter writer = new FileWriter(outFile);
					// The primer pair ID is the name of the outfile
					String data = filename;
					data += " ";
					data += Integer.valueOf(primer.getLeftPrimerPosition()).toString();
					data += " ";
					data += primer.getLeftPrimer();
					data += " ";
					data += Integer.valueOf(primer.getRightPrimerPosition()).toString();
					data += " ";
					data += primer.getRightPrimer();
					data += " ";
					data += PrimerUtils.computeTM(primer.getLeftPrimer()); // recompute TM with our tool
					data += " ";
					data += PrimerUtils.computeTM(primer.getRightPrimer());;
					data += " ";
					data += Integer.valueOf(primer.getProductSize()).toString();
					data += " ";
					data += "no_comment";
					// Write the data to read into PrimerPair constructor later
					writer.write(data);
					writer.close();
				} 
			}
		}
	}
	
	
	/**
	 * Reads one primer pair from the first line of file
	 * Fields in the first line should be in the format of the String[] passed to the PrimerPair constructor
	 * @param fileName
	 * @return the PrimerPair object whose data is listed on the first line of file
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	PrimerPair readOnePrimerPairFromFile(String fileName) throws IOException, FileNotFoundException {
		FileReader reader = new FileReader(fileName);
		BufferedReader b = new BufferedReader(reader);
		String line = b.readLine();
		StringParser stringparse = new StringParser();
		stringparse.parse(line);
		// Check if line has the right number of fields to be a primer data line
		if(stringparse.getFieldCount() != 9) {
			throw new IllegalArgumentException("First line of file does not contain primer pair data.");
		}
		String[] data = stringparse.getStringArray();
		// Parse the array of data and create the corresponding primerPair object
		PrimerPair primerPair = new PrimerPair(data);
		return primerPair;
	}
	
	/**
	 * Get a set of primers from a data file
	 * Fields in each line should be in the format of the String[] passed to the PrimerPair constructor
	 * @param fileName
	 * @return the set of primer pairs
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	ArrayList<PrimerPair> readPrimerPairsFromFile(String fileName) throws FileNotFoundException, IOException {
		ArrayList<PrimerPair> rtrn = new ArrayList<PrimerPair>();
		FileReader reader = new FileReader(fileName);
		BufferedReader b = new BufferedReader(reader);
		// Read one line of data at a time and parse into primerPair
		while(b.ready()) {
			String line = b.readLine();
			StringParser stringparse = new StringParser();
			stringparse.parse(line);
			if(stringparse.getFieldCount() != 9) {
				throw new IllegalArgumentException("Line of file does not contain primer pair data.");
			}
			String[] data = stringparse.getStringArray();
			PrimerPair primerPair = new PrimerPair(data);
			rtrn.add(primerPair);
		}
		return rtrn;		
	}
	
	/**
	 * Get a new primer with extra tail at the end within suitable TM and GC of the original
	 * @param initPair initial primer pair
	 * @param extraLength number of bases to add
	 * @param outerSubprimerUniqueBases the starting position of the subprimer to extend
	 * @return the new longer primer
	 * @throws IOException
	 */
	private static PrimerPair getOneMutatedSubPrimer(PrimerPair initPair, int extraLength, int outerSubprimerUniqueBases) throws IOException {
		
		long numtries = 0L;
		DesignPCRTails dpt = new DesignPCRTails();
		
		// The initial outer primers
		String initLeftPrimer = initPair.getLeftPrimer();
		String initRightPrimer = initPair.getRightPrimer();
		
		// Recalculate TMs so the new primers will be calculated the same way
		double initLeftTM = PrimerUtils.computeTM(initLeftPrimer);
		double initRightTM = PrimerUtils.computeTM(initRightPrimer);
		double initLeftGC = PrimerUtils.percentGC(initLeftPrimer);
		double initRightGC = PrimerUtils.percentGC(initRightPrimer);
		
		// Get the part of initial primers that overlaps the new primer
		String initLeftOverlap = initLeftPrimer.substring(outerSubprimerUniqueBases);
		String initRightOverlap = initRightPrimer.substring(outerSubprimerUniqueBases);
		
		//String finalLeftPrimer = null;
		//String finalRightPrimer = null;
		//double finalLeftTM = 0;
		//double finalRightTM = 0;
		
		// Repeat until a suitable primer pair is found
		while(numtries < 1000000000000L) {
			
			// Generate a random extra tail
			String randLeft = Sequence.generateRandomSequence(extraLength);
			String randRight = Sequence.generateRandomSequence(extraLength);
			
			// The new primers to try
			String newLeftPrimer = initLeftOverlap + randLeft;
			String newRightPrimer = initRightOverlap + randRight;
			
			// Every 10,000 tries, print the reason a primer pair is being rejected
			//if(numtries % 10000 == 0) {
				//System.out.println("randleft=" + randLeft + " randright=" + randRight + " newleft=" + newLeftPrimer + " newright="+ newRightPrimer);
			//}
			
			// Get the TMs of the new inner primers
			double newLeftInnerTM = PrimerUtils.computeTM(newLeftPrimer);
			double newRightInnerTM = PrimerUtils.computeTM(newRightPrimer);
			
			// Compare TMs to original outer primers
			if(Math.abs(initLeftTM - newLeftInnerTM) > 1) {
				//if(numtries % 10000 == 0) System.out.println("Rejecting left primer because TM is too different " + newLeftPrimer + " oldTM=" + initLeftTM + " newTM=" + newLeftInnerTM);
				numtries++;
				randLeft = null;
				randRight = null;
				newLeftPrimer = null;
				newRightPrimer = null;
				continue;
			}
			
			if(Math.abs(initRightTM - newRightInnerTM) > 1) {
				//if(numtries % 10000 == 0) System.out.println("Rejecting right primer because TM is too different " + newRightPrimer + " oldTM=" + initRightTM + " newTM=" + newRightInnerTM);
				numtries++;
				randLeft = null;
				randRight = null;
				newLeftPrimer = null;
				newRightPrimer = null;
				continue;
			}
			
			// Get the GC of the new inner primers
			double newLeftInnerGC = PrimerUtils.percentGC(newLeftPrimer);
			double newRightInnerGC = PrimerUtils.percentGC(newRightPrimer);
			
			// Compare GC to original outer primers
			if(Math.abs(initLeftGC - newLeftInnerGC) > 5) {
				//if(numtries % 10000 == 0) System.out.println("Rejecting left primer because GC is too different " + newLeftPrimer + " oldGC=" + initLeftGC + " newGC=" + newLeftInnerGC);
				numtries++;
				randLeft = null;
				randRight = null;
				newLeftPrimer = null;
				newRightPrimer = null;
				continue;
			}
			
			if(Math.abs(initRightGC - newRightInnerGC) > 5) {
				//if(numtries % 10000 == 0) System.out.println("Rejecting right primer because GC is too different " + newRightPrimer + " oldGC=" + initRightGC + " newGC=" + newRightInnerGC);
				numtries++;
				randLeft = null;
				randRight = null;
				newLeftPrimer = null;
				newRightPrimer = null;
				continue;
			}
			
			// Check if the full length primers match over greater than dpt.minMatch bases
			if(dpt.complementary(newLeftPrimer, newRightPrimer)) {
				//if(numtries % 10000 == 0) System.out.println("Rejecting complementary primers " + newLeftPrimer + " " + newRightPrimer);
				numtries++;
				randLeft = null;
				randRight = null;
				newLeftPrimer = null;
				newRightPrimer = null;
				continue;
			}
			
			// If code gets to this point a primer has been found
			//System.out.println("Found new extended primer " + newLeftPrimer + " " + newRightPrimer + " " + newLeftInnerTM + " " + newRightInnerTM);
			
			String[] primerPairData = new String[9];
			
			primerPairData[0] = initPair.getPrimerPairId() + "_Subprimer"; // primer pair ID
			primerPairData[1] = Integer.toString(initPair.getLeftPrimerPosition()); // left primer position
			primerPairData[2] = newLeftPrimer; // left primer
			primerPairData[3] = Integer.toString(initPair.getRightPrimerPosition()); // right primer position
			primerPairData[4] = newRightPrimer; // right primer
			primerPairData[5] = Double.toString(newLeftInnerTM); // left primer TM
			primerPairData[6] = Double.toString(newRightInnerTM); // right primer TM
			primerPairData[7] = Integer.toString(initPair.getProductSize()); // product size
			primerPairData[8] = initPair.getComment(); // comment

			PrimerPair rtrn = new PrimerPair(primerPairData);
			return rtrn;
			
			
		}
	
		throw new IllegalArgumentException("Suitable extended primer could not be found.");
		
	}
	
	/**
	 * Get extended primers from a single primer pair and write to a file
	 * @param initPair
	 * @param outerSubprimerUniqueBases the starting position of the subprimer to extend
	 * @param extraLength
	 * @param numSubPrimers
	 * @param outFile
	 * @throws TimeoutException
	 * @throws IOException
	 */
	private static void writeGoodSubPrimers(PrimerPair initPair, int outerSubprimerUniqueBases, int extraLength, int numSubPrimers, String outFile) throws TimeoutException,IOException {
		
		//System.out.println("Trying to write " + numSubPrimers + " primer pairs to file " + outFile);
		
		ArrayList<PrimerPair> newMutations = new ArrayList<PrimerPair>();
		
		int numtries = 0;
		// First fill the set with mutated primers
		while(newMutations.size() < numSubPrimers) {
			// Throw exception if can't find subprimers after this many tries
			// Doesn't need to be very big because getOneMutatedSubPrimer tries a lot of random sequences
			if(numtries > 100 * numSubPrimers) {
				System.out.println("Too many tries: " + numtries);
				throw new TimeoutException("The primer pair does not yield enough suitable subprimers.");
			}
			try {
				PrimerPair newPrimer = getOneMutatedSubPrimer(initPair, extraLength, outerSubprimerUniqueBases);
				// Set the ID of the new primer to be the ID of the parent primer plus a number
				newPrimer.setPrimerPairId(initPair.getPrimerPairId() + "_" + numtries);
				newMutations.add(newPrimer);
				numtries++;
			} catch (IllegalArgumentException e) {
				// getOneMutatedSubPrimer throws this exception if no suitable subprimer is found
				numtries++;
			}
		}
		
		
		// Write the new primers to file
		FileWriter writer = new FileWriter(outFile);
		//System.out.println("Writing to file " + outFile);
		for(PrimerPair primer : newMutations) {
			
			String data = primer.getPrimerPairId();
			data += " ";
			data += Integer.valueOf(primer.getLeftPrimerPosition()).toString();
			data += " ";
			data +=  primer.getLeftPrimer();
			data += " ";
			data +=  Integer.valueOf(primer.getRightPrimerPosition()).toString();
			data += " ";
			data +=  primer.getRightPrimer();
			data += " ";
			data +=  Double.valueOf(primer.getLeftPrimerTM()).toString();
			data += " ";
			data += Double.valueOf(primer.getRightPrimerTM()).toString();
			data += " ";
			data +=  Integer.valueOf(primer.getProductSize()).toString();
			data += " ";
			data += "no_comment\n";
			// Write the data to read into PrimerPair constructor later
			writer.write(data);
		}
		writer.close();
		
	}
	
	/**
	 * Create initial full length primers to be mutated later to add more overlapping inner primers
	 * @throws IOException
	 */
	private void createInitialOuterPrimerPairs() throws IOException, InterruptedException {
		
		System.out.println("Creating initial primer pairs...");
		
		if(this.leftOuterLength != this.rightOuterLength) {
			throw new IllegalStateException("Outer primer lengths must be equal.");
		}
		
		File dir = new File("TmpInitialPrimers_" + this.projectName);
		System.out.println("Creating directory " + dir);
		boolean madedir = dir.mkdirs();
		if(!dir.exists()) {
			System.err.println("Could not create directory " + dir);
			//System.exit(-1);
		}
		System.out.println("Creating primer and writing a temp file for each initial primer to directory " + dir + "...");
		
		// Keep track of progress
		boolean tenPercent = false;
		boolean twentyPercent = false;
		boolean thirtyPercent = false;
		boolean fortyPercent = false;
		boolean fiftyPercent = false;
		boolean sixtyPercent = false;
		boolean seventyPercent = false;
		boolean eightyPercent = false;
		boolean ninetyPercent = false;
		
		this.initialOuterPrimerPairs.clear();
		int numtries = 0;
		ArrayList<String> jobIDs = new ArrayList<String>();
		
		int numToMake = this.outerCount.intValue();
		
		// Pick outer primers with acceptable scores
		while(numtries < numToMake) {
	
			// Submit 100 jobs at a time
			if(numtries > 0 && numtries % 100 == 0) {
				
				
				// Every 100 jobs, provide indication of progress
				if(!tenPercent) {
					if((float)numtries/((float)numToMake) > 0.1) {
						System.out.println("Finished 10% of primers");
						tenPercent = true;
					}
				}
				if(!twentyPercent) {
					if((float)numtries/((float)numToMake) > 0.2) {
						System.out.println("Finished 20% of primers");
						twentyPercent = true;
					}
				}
				if(!thirtyPercent) {
					if((float)numtries/((float)numToMake) > 0.3) {
						System.out.println("Finished 30% of primers");
						thirtyPercent = true;
					}
				}
				if(!fortyPercent) {
					if((float)numtries/((float)numToMake) > 0.4) {
						System.out.println("Finished 40% of primers");
						fortyPercent = true;
					}
				}
				if(!fiftyPercent) {
					if((float)numtries/((float)numToMake) > 0.5) {
						System.out.println("Finished 50% of primers");
						fiftyPercent = true;
					}
				}
				if(!sixtyPercent) {
					if((float)numtries/((float)numToMake) > 0.6) {
						System.out.println("Finished 60% of primers");
						sixtyPercent = true;
					}
				}
				if(!seventyPercent) {
					if((float)numtries/((float)numToMake) > 0.7) {
						System.out.println("Finished 70% of primers");
						seventyPercent = true;
					}
				}
				if(!eightyPercent) {
					if((float)numtries/((float)numToMake) > 0.8) {
						System.out.println("Finished 80% of primers");
						eightyPercent = true;
					}
				}
				if(!ninetyPercent) {
					if((float)numtries/((float)numToMake) > 0.9) {
						System.out.println("Finished 90% of primers");
						ninetyPercent = true;
					}
				}
				
				
				// Wait to submit more jobs
				Thread.sleep(1000);
			}
			
			// The file to write the primer to
			// Don't bother if the primer file already exists
			String outFile = "TmpInitialPrimers_" + this.projectName + "/Primer_" + Integer.valueOf(numtries).toString();
			File file = new File(outFile);
			if(file.exists()) {
				numtries++;
				continue;
			}
			
			// Parameters to submit job to LSF
			String argstring = this.leftOuterLength.toString() + " " + outFile;
			String jobID = "job_Primer_" + numtries + "_" + Long.valueOf(System.currentTimeMillis()).toString();
			jobIDs.add(jobID);
			// Use the main method which only does this task
			// For main array
			String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools/build:/seq/lincRNA/Pam/Software/gbtools/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
			// For mrna depletion array
			//String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
			// For short probes
			// cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
			
			// Submit the job
			try {
				int exitCode = PipelineUtils.bsubSmallProcess(Runtime.getRuntime(), jobID , cmmd , "TmpInitialPrimers_" + this.projectName + "/bsub_output_" + jobID);
			} catch (InterruptedException e) {
				System.err.println("Caught InterruptedException when trying to submit job " + jobID);
				e.printStackTrace();
			} catch (IOException e) {
				System.err.println("Caught IOException when trying to submit job " + jobID);
				e.printStackTrace();					
			}
			
			numtries++;
		}
		
		// Wait for the last job to finish
		// Throw exception if more than 10% of jobs failed
		System.out.println("Waiting for jobs to finish.");
		try {
			PipelineUtils.waitForAllJobs(jobIDs, Runtime.getRuntime());
			System.out.println("Done creating initial primer temp files.\n");
			
			// Read all the primers back in from the files
			System.out.println("Reading initial primers back in from temp files...");
			for(int p = 0; p < numToMake; p++) {
				String primerFile = dir + "/Primer_" + Integer.valueOf(p).toString();
				PrimerPair primerPair = this.readOnePrimerPairFromFile(primerFile);
				this.initialOuterPrimerPairs.add(primerPair);
			}
			System.out.println("Read and stored " + this.initialOuterPrimerPairs.size() + " initial primers.");
		} catch (IllegalArgumentException e) {
			// Catch the exception thrown when at least one job fails
			// Now check how bad the situation is
			System.out.println("Done creating initial primer temp files.\n");
			System.out.println("At least one LSF job failed. Checking if results are salvageable.");
			System.out.println("Reading initial primers back in from temp files...");
			int failedPrimers = 0;
			for(int p = 0; p < this.outerCount.intValue() * 2; p++) {
				String primerFile = dir + "/Primer_" + Integer.valueOf(p).toString();
				File tempfile = new File(primerFile);
				if(!tempfile.exists()) {
					failedPrimers++;
					continue;
				}
				PrimerPair primerPair = this.readOnePrimerPairFromFile(primerFile);
				this.initialOuterPrimerPairs.add(primerPair);
			}
			if(failedPrimers / (this.outerCount.intValue() * 2) > 0.1) {
				System.out.println(failedPrimers + " jobs failed. Exiting.");
				throw new IllegalArgumentException("Too many LSF jobs failed.");
			}
			System.out.println("Read and stored " + this.initialOuterPrimerPairs.size() + " initial primers.");
			
		}

		String outfile = "OuterPrimers_" + this.projectName + ".out";
		System.out.println("Writing initial outer primers to file " + outfile);
		
		FileWriter writer = new FileWriter(outfile);
		for(PrimerPair primer : this.initialOuterPrimerPairs) {
			writer.write(primer.getPrimerFieldsAsStringForConstructor() + "\n");
		}
		writer.close();
				
		
		
	}

	/**
	 * Mutate the outer primers to add more overlapping second primers
	 * Adds an initial primer and its subprimers to secondPrimerPairs if enough good subprimers can be found in time
	 */
	private void createSecondOverlappingPrimers() throws IOException, InterruptedException {
		
		if(this.leftOuterEnd != this.rightOuterEnd) {
			throw new IllegalStateException("Outer primer lengths must be equal.");
		}
		
		if(this.leftSecondEnd != this.rightSecondEnd) {
			throw new IllegalStateException("Second primers must end at same position.");
		}

		// How many extra bases to add
		int extraLength = this.leftSecondEnd.intValue() - this.leftOuterEnd.intValue();
		
		System.out.println("\nMutating initial primers to add second inner overlapping primer.\n");
		System.out.println("Writing a temp file of subprimers for each outer primer to directory TmpSecondPrimers_" + this.projectName + "/");
		
		File dir = new File("TmpSecondPrimers_" + this.projectName);
		System.out.println("Creating directory " + dir);
		boolean madedir = dir.mkdirs();
		if(!dir.exists()) {
			System.err.println("Could not create directory " + dir);
			//System.exit(-1);
		}
		System.out.println("Creating primers and writing a temp file for each second primer to directory " + dir + "...");

		
		this.secondPrimerPairs.clear();
		
		if(this.initialOuterPrimerPairs.isEmpty()) {
			throw new IllegalStateException("Can't try to make second level primers because outer primer set is empty. Try calling createInitialOuterPrimerPairs() first.");
		}
		
		// Keep track of jobs
		int numtries = 0;
		ArrayList<String> jobIDs = new ArrayList<String>();
		
		// Keep track of progress
		boolean tenPercent = false;
		boolean twentyPercent = false;
		boolean thirtyPercent = false;
		boolean fortyPercent = false;
		boolean fiftyPercent = false;
		boolean sixtyPercent = false;
		boolean seventyPercent = false;
		boolean eightyPercent = false;
		boolean ninetyPercent = false;
		int totalprimers = this.initialOuterPrimerPairs.size();

		
		// Do this for every initial outer primer
		for(PrimerPair initialpair : this.initialOuterPrimerPairs) {
			
			// Submit 100 jobs at a time
			if(numtries > 0 && numtries % 100 == 0) {
				
				// Every 100 jobs, provide indication of progress
				if(!tenPercent) {
					if((float)numtries/((float)totalprimers) > 0.1) {
						System.out.println("Finished 10% of primers");
						tenPercent = true;
					}
				}
				if(!twentyPercent) {
					if((float)numtries/((float)totalprimers) > 0.2) {
						System.out.println("Finished 20% of primers");
						twentyPercent = true;
					}
				}
				if(!thirtyPercent) {
					if((float)numtries/((float)totalprimers) > 0.3) {
						System.out.println("Finished 30% of primers");
						thirtyPercent = true;
					}
				}
				if(!fortyPercent) {
					if((float)numtries/((float)totalprimers) > 0.4) {
						System.out.println("Finished 40% of primers");
						fortyPercent = true;
					}
				}
				if(!fiftyPercent) {
					if((float)numtries/((float)totalprimers) > 0.5) {
						System.out.println("Finished 50% of primers");
						fiftyPercent = true;
					}
				}
				if(!sixtyPercent) {
					if((float)numtries/((float)totalprimers) > 0.6) {
						System.out.println("Finished 60% of primers");
						sixtyPercent = true;
					}
				}
				if(!seventyPercent) {
					if((float)numtries/((float)totalprimers) > 0.7) {
						System.out.println("Finished 70% of primers");
						seventyPercent = true;
					}
				}
				if(!eightyPercent) {
					if((float)numtries/((float)totalprimers) > 0.8) {
						System.out.println("Finished 80% of primers");
						eightyPercent = true;
					}
				}
				if(!ninetyPercent) {
					if((float)numtries/((float)totalprimers) > 0.9) {
						System.out.println("Finished 90% of primers");
						ninetyPercent = true;
					}
				}

				// Wait to submit more jobs
				Thread.sleep(1000);
			}
			
			// Don't bother if the primer file already exists
			String outFile = "TmpSecondPrimers_" + this.projectName + "/" + initialpair.getPrimerPairId();
			File file = new File(outFile);
			if(file.exists()) {
				numtries++;
				continue;
			}
			
			// Parameters to submit job to LSF
			String argstring = initialpair.getPrimerPairId();
			argstring += " ";
			argstring += initialpair.getLeftPrimerPosition();
			argstring += " ";
			argstring += initialpair.getLeftPrimer();
			argstring += " ";
			argstring += initialpair.getRightPrimerPosition();
			argstring += " ";
			argstring += initialpair.getRightPrimer();
			argstring += " ";
			argstring += initialpair.getLeftPrimerTM();
			argstring += " ";
			argstring += initialpair.getRightPrimerTM();
			argstring += " ";
			argstring += initialpair.getProductSize();
			argstring += " ";
			argstring += initialpair.getComment();
			argstring += " ";
			argstring += Integer.valueOf(extraLength).toString();
			argstring += " ";
			argstring += Integer.valueOf(this.secondCount.intValue()).toString();
			argstring += " ";
			argstring += outFile;
			argstring += " ";
			argstring += this.leftSecondStart.toString();
		
			
			String jobID = "job_" + initialpair.getPrimerPairId() + "_" + Long.valueOf(System.currentTimeMillis()).toString();
			jobIDs.add(jobID);
			// Use the main method which only does this task
			// For main array
			String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools/build:/seq/lincRNA/Pam/Software/gbtools/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
			// For mrna depletion array
			//String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
			// For short probes
			//String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
				
			if(this.verbose) System.out.println("Submitting command to LSF: " + cmmd);
				
			// Submit the job
			try {
				int exitCode = PipelineUtils.bsubSmallProcess(Runtime.getRuntime(), jobID , cmmd , "TmpSecondPrimers_" + this.projectName + "/bsub_output_" + jobID);
			} catch (InterruptedException e) {
				System.err.println("Caught InterruptedException when trying to submit job " + jobID);
				e.printStackTrace();
			} catch (IOException e) {
				System.err.println("Caught IOException when trying to submit job " + jobID);
				e.printStackTrace();					
			}
			
			numtries++;

		}
		
		// Wait for half of jobs to finish
		// This is enough because we originally made 10x the required number of outer primers
		System.out.println("Waiting for jobs to finish.");
		try {
			ArrayList<String> stillRunning = PipelineUtils.waitForEnoughJobs(jobIDs, jobIDs.size() / 2, Runtime.getRuntime());
			System.out.println("Done creating second subprimer temp files.\n");
			// Read all the primers back in from the files and add to secondPrimerPairs
			System.out.println("Reading second primers back in from temp files.");
			this.secondPrimerPairs.clear();
			int newPrimers = 0;
			for(PrimerPair p : this.initialOuterPrimerPairs) {
				String tempfile = "TmpSecondPrimers_" + this.projectName + "/" + p.getPrimerPairId();
				File f = new File(tempfile);
				if(f.exists()) {
					ArrayList<PrimerPair> newMutations = new ArrayList<PrimerPair>();
					newMutations = this.readPrimerPairsFromFile(tempfile);
					this.secondPrimerPairs.put(p, newMutations);
					newPrimers += this.secondPrimerPairs.get(p).size();
				}
			}
			PipelineUtils.bkillAll(stillRunning, Runtime.getRuntime());
			System.out.println("Done reading second primers. Added " + newPrimers + " new extended primers.");
		} catch (IllegalArgumentException e) {
			System.out.println("Not enough jobs were successful.");
			e.printStackTrace();
		}
				
		// Write the second level primers to a file
		String outfile = "SecondPrimers_" + this.projectName + ".out";
		System.out.println("Writing second primers to file " + outfile);
		
		FileWriter writer = new FileWriter(outfile);
		for(PrimerPair primer : this.initialOuterPrimerPairs) {
			if(!this.secondPrimerPairs.containsKey(primer)) continue;
			writer.write("OUTER_PRIMER\t" + primer.getPrimerFieldsAsStringForConstructor() + "\n");
			for(PrimerPair second : this.secondPrimerPairs.get(primer)) {
				writer.write("INNER_PRIMER\t" + second.getPrimerFieldsAsStringForConstructor() + "\n");
			}
		}
		writer.close();
		
	}
	
	/**
	 * Mutate the second primers to add more overlapping third primers
	 */
	private void createThirdOverlappingPrimers() throws IOException, InterruptedException {
		
		if(this.leftThirdLength != this.rightThirdLength) {
			throw new IllegalStateException("Third primer lengths must be equal.");
		}
		
		if(this.leftSecondEnd != this.rightSecondEnd) {
			throw new IllegalStateException("Second primers must end at same position.");
		}

		// How many extra bases to add
		int extraLength = this.primerSize.intValue() - 1 - this.leftSecondEnd.intValue();
		
		System.out.println("\nMutating second level primers to add third inner overlapping primer.\n");
		System.out.println("Writing a temp file of subprimers for each second primer to directory TmpThirdPrimers_" + this.projectName + "/");
		
		File dir = new File("TmpThirdPrimers_" + this.projectName);
		System.out.println("Creating directory " + dir);
		boolean madedir = dir.mkdirs();
		if(!dir.exists()) {
			System.err.println("Could not create directory " + dir);
			//System.exit(-1);
		}
		System.out.println("Creating primers and writing a temp file for each third primer to directory " + dir + "...");

		
		this.thirdPrimerPairs.clear();
		
		if(this.secondPrimerPairs.isEmpty()) {
			throw new IllegalStateException("Can't try to make third level primers because second primer set is empty. Try calling createSecondOverlappingPrimers() first.");
		}
		
		// Keep track of jobs
		int numtries = 0;
		ArrayList<String> jobIDs = new ArrayList<String>();
		
		// Keep track of progress
		boolean tenPercent = false;
		boolean twentyPercent = false;
		boolean thirtyPercent = false;
		boolean fortyPercent = false;
		boolean fiftyPercent = false;
		boolean sixtyPercent = false;
		boolean seventyPercent = false;
		boolean eightyPercent = false;
		boolean ninetyPercent = false;
		// Count second level primers
		int totalprimers = 0;
		
		for(PrimerPair p : this.initialOuterPrimerPairs) {
			if(this.secondPrimerPairs.containsKey(p)) totalprimers += this.secondPrimerPairs.get(p).size();
		}

		
		// Do this for every second level primer
		for(PrimerPair initialpair : this.initialOuterPrimerPairs) {
			
			if(this.secondPrimerPairs.containsKey(initialpair)) {
			
				for(PrimerPair secondpair : this.secondPrimerPairs.get(initialpair)) {
			
					// Submit 100 jobs at a time
					if(numtries > 0 && numtries % 100 == 0) {
						
						// Every 100 jobs, provide indication of progress
						if(!tenPercent) {
							if((float)numtries/((float)totalprimers) > 0.1) {
								System.out.println("Finished 10% of primers");
								tenPercent = true;
							}
						}
						if(!twentyPercent) {
							if((float)numtries/((float)totalprimers) > 0.2) {
								System.out.println("Finished 20% of primers");
								twentyPercent = true;
							}
						}
						if(!thirtyPercent) {
							if((float)numtries/((float)totalprimers) > 0.3) {
								System.out.println("Finished 30% of primers");
								thirtyPercent = true;
							}
						}
						if(!fortyPercent) {
							if((float)numtries/((float)totalprimers) > 0.4) {
								System.out.println("Finished 40% of primers");
								fortyPercent = true;
							}	
						}
						if(!fiftyPercent) {
							if((float)numtries/((float)totalprimers) > 0.5) {
								System.out.println("Finished 50% of primers");
								fiftyPercent = true;
							}
						}
						if(!sixtyPercent) {
							if((float)numtries/((float)totalprimers) > 0.6) {
								System.out.println("Finished 60% of primers");
								sixtyPercent = true;
							}
						}
						if(!seventyPercent) {
							if((float)numtries/((float)totalprimers) > 0.7) {
								System.out.println("Finished 70% of primers");
								seventyPercent = true;
							}
						}
						if(!eightyPercent) {
							if((float)numtries/((float)totalprimers) > 0.8) {
								System.out.println("Finished 80% of primers");
								eightyPercent = true;
							}
						}
						if(!ninetyPercent) {
							if((float)numtries/((float)totalprimers) > 0.9) {
								System.out.println("Finished 90% of primers");
								ninetyPercent = true;
							}
						}

						// Wait to submit more jobs
						Thread.sleep(20);
					}
			
					// Don't bother if the file exists
					String outFile = "TmpThirdPrimers_" + this.projectName + "/" + secondpair.getPrimerPairId();
					File file = new File(outFile);
					if(file.exists()) {
						numtries++;
						continue;
					}
					
					// Parameters to submit job to LSF
					String argstring = secondpair.getPrimerPairId();
					argstring += " ";
					argstring += secondpair.getLeftPrimerPosition();
					argstring += " ";
					argstring += secondpair.getLeftPrimer();
					argstring += " ";
					argstring += secondpair.getRightPrimerPosition();
					argstring += " ";
					argstring += secondpair.getRightPrimer();
					argstring += " ";
					argstring += secondpair.getLeftPrimerTM();
					argstring += " ";
					argstring += secondpair.getRightPrimerTM();
					argstring += " ";
					argstring += secondpair.getProductSize();
					argstring += " ";
					argstring += secondpair.getComment();
					argstring += " ";
					argstring += Integer.valueOf(extraLength).toString(); // length to add
					argstring += " ";
					argstring += Integer.valueOf(this.thirdCount.intValue()).toString();
					argstring += " ";
					argstring += outFile;
					argstring += " ";
					argstring += Integer.valueOf(this.leftThirdStart.intValue() - this.leftSecondStart.intValue()).toString();
		
			
					String jobID = "job_" + secondpair.getPrimerPairId() + "_" + Long.valueOf(System.currentTimeMillis()).toString();
					jobIDs.add(jobID);
					// Use the main method which only does this task
					// For main array
					String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools/build:/seq/lincRNA/Pam/Software/gbtools/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
					// For mrna depletion array
					//String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_mrna_depletion/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;
					// For short probes
					//String cmmd = "java -classpath /seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/build:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/igv.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/log4j-1.2.14.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/ATV-3.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/JConnector-5.0.6.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/a_sam-1.48.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/collections-generic-4.01.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/colt-1.2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/commons-math-1.2-SNAPSHOT.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/epsgraphics-1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jaligner-1.0-MG.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jgraphx.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-algorithms-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-api-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-graph-impl-2.0.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/jung-visualization-2.0.1.jar:/seq/lincRNA/Pam/Software/gbtools_copy_for_short_probes/lib/junit-4.4.jar broad.pda.capture.designer.PcrTailDesigner " + argstring;

					if(this.verbose) System.out.println("Submitting command to LSF: " + cmmd);
				
					// Submit the job
					try {
						int exitCode = PipelineUtils.bsubSmallProcess(Runtime.getRuntime(), jobID , cmmd , "TmpThirdPrimers_" + this.projectName + "/bsub_output_" + jobID);
					} catch (InterruptedException e) {
						System.err.println("Caught InterruptedException when trying to submit job " + jobID);
						e.printStackTrace();
					} catch (IOException e) {
						System.err.println("Caught IOException when trying to submit job " + jobID);
						e.printStackTrace();					
					}
			
					numtries++;

				}
			
			}
		}
		
		// Wait for half of jobs to finish
		System.out.println("Waiting for jobs to finish.");
		try {
			ArrayList<String> stillRunning = PipelineUtils.waitForEnoughJobs(jobIDs, jobIDs.size() / 2, Runtime.getRuntime());
			System.out.println("Done creating third subprimer temp files.\n");
			// Read all the primers back in from the files and add to thirdPrimerPairs
			System.out.println("Reading third primers back in from temp files.");
			this.thirdPrimerPairs.clear();
			int newPrimers = 0;
			for(PrimerPair p : this.initialOuterPrimerPairs) {
				if(!this.secondPrimerPairs.containsKey(p)) continue;
				// For each outer primer, create a map from its second primers to their third primers
				Map< PrimerPair, ArrayList<PrimerPair >> allForOuter = new HashMap< PrimerPair, ArrayList<PrimerPair > >();				
				for(PrimerPair s : this.secondPrimerPairs.get(p)) {
					String tempfile = "TmpThirdPrimers_" + this.projectName + "/" + s.getPrimerPairId();
					File f = new File(tempfile);
					if(f.exists()) {
						ArrayList<PrimerPair> newMutations = new ArrayList<PrimerPair>();
						newMutations = this.readPrimerPairsFromFile(tempfile);
						allForOuter.put(s, newMutations);
					}
				}
				// If some third primers were created for this outer primer, add to thirdPrimerPairs
				if(!allForOuter.isEmpty()) {
					this.thirdPrimerPairs.put(p, allForOuter);
					// Count the number of new primers added
					for(PrimerPair q : this.thirdPrimerPairs.get(p).keySet()) {
						newPrimers += this.thirdPrimerPairs.get(p).get(q).size();
					}
				}
			}
			PipelineUtils.bkillAll(stillRunning, Runtime.getRuntime());
			System.out.println("Done reading third primers. Added " + newPrimers + " new extended primers.");
		} catch (IllegalArgumentException e) {
			System.out.println("Not enough jobs were successful.");
			e.printStackTrace();
		}
				
		// Write the third level primers to a file
		String outfile = "ThirdPrimers_" + this.projectName + ".out";
		System.out.println("Writing third primers to file " + outfile);
		
		FileWriter writer = new FileWriter(outfile);
		for(PrimerPair primer : this.initialOuterPrimerPairs) {
			if(!this.thirdPrimerPairs.containsKey(primer)) continue;
			writer.write("OUTER_PRIMER\t" + primer.getPrimerFieldsAsStringForConstructor() + "\n");
			for(PrimerPair second : this.thirdPrimerPairs.get(primer).keySet()) {
				writer.write("SECOND_PRIMER\t" + second.getPrimerFieldsAsStringForConstructor() + "\n");
				for(PrimerPair third : this.thirdPrimerPairs.get(primer).get(second)) {
					writer.write("THIRD_PRIMER\t" + third.getPrimerFieldsAsStringForConstructor() + "\n");
				}
			}
		}
		writer.close();
		
		
	}
	
	
	/**
	 * Create the final full primer set
	 * @throws IOException
	 */
	public void createPrimers() throws IOException, InterruptedException {
		this.initialOuterPrimerPairs.clear();
		this.secondPrimerPairs.clear();
		this.thirdPrimerPairs.clear();
		this.createInitialOuterPrimerPairs();
		this.createSecondOverlappingPrimers();
		if(this.threeSegments) this.createThirdOverlappingPrimers();
	}
	
	/**
	 * Get the full set of outer primer pairs
	 * @return the primer set
	 */
	public Collection<PrimerPair> getOuterPrimers() {
		if(this.initialOuterPrimerPairs.isEmpty()) {
			throw new IllegalStateException("Outer primer set is empty. Try calling createPrimers() first.");
		}
		return this.initialOuterPrimerPairs;
	}
	
	/**
	 * Get the full set of second level primer pairs
	 * @return the primer set
	 */
	public Map<PrimerPair, ArrayList<PrimerPair>> getSecondPrimers() {
		if(this.secondPrimerPairs.isEmpty()) {
			throw new IllegalStateException("Second primer set is empty. Try calling createPrimers() first.");
		}
		return this.secondPrimerPairs;
	}
	
	/**
	 * Get the full set of third level primer pairs
	 * @return the primer set
	 */
	public Map<PrimerPair, Map<PrimerPair, ArrayList<PrimerPair>>> getThirdPrimers() {
		if(this.thirdPrimerPairs.isEmpty()) {
			throw new IllegalStateException("Third primer set is empty. Try calling createPrimers() first.");
		}
		return this.thirdPrimerPairs;
	}
	
	
	/**
	 * Just write subprimers to a file
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		// The situation where we are writing a single initial primer pair
		if(args.length == 2) {
			
			int primerSize = Integer.parseInt(args[0]);
			String outFile = args[1];
			writeOneInitialPrimerPair(primerSize,outFile);
			
		}
		
		// The situation where we are writing the set of suitable subprimers for an existing primer pair
		if(args.length == 13) {
			
			String[] primerPairData = new String[9];
		
			primerPairData[0] = args[0]; // primer pair ID
			primerPairData[1] = args[1]; // left primer position
			primerPairData[2] = args[2]; // left primer
			primerPairData[3] = args[3]; // right primer position
			primerPairData[4] = args[4]; // right primer
			primerPairData[5] = args[5]; // left primer TM
			primerPairData[6] = args[6]; // right primer TM
			primerPairData[7] = args[7]; // product size
			primerPairData[8] = args[8]; // comment
		
			PrimerPair primer = new PrimerPair(primerPairData);
		
			int extraLength = Integer.parseInt(args[9]);
			int numSubPrimers = Integer.parseInt(args[10]);
			String outFile = args[11];
			int outerSubprimerUniqueBases = Integer.parseInt(args[12]);
			
			try {
				writeGoodSubPrimers(primer, outerSubprimerUniqueBases, extraLength, numSubPrimers, outFile);
			} catch (TimeoutException e) {
				System.out.println("Could not find subprimers for primer " + args[2] + " " + args[4]);
			}
		}
	
	}

	

	private Integer primerSize;
	
	private Integer leftOuterLength;
	private Integer leftOuterStart;
	private Integer leftOuterEnd;
	private Integer leftSecondLength;
	private Integer leftSecondStart;
	private Integer leftSecondEnd;
	private Integer leftThirdLength;
	private Integer leftThirdStart;
	private Integer leftThirdEnd;
	
	private Integer rightOuterLength;
	private Integer rightOuterStart;
	private Integer rightOuterEnd;
	private Integer rightSecondLength;
	private Integer rightSecondStart;
	private Integer rightSecondEnd;
	private Integer rightThirdLength;
	private Integer rightThirdStart;
	private Integer rightThirdEnd;
	
	private Integer outerCount;
	private Integer secondCount;
	private Integer thirdCount;
	
	private boolean threeSegments;
	private boolean verbose;
	private boolean veryVerbose;
	
	private Collection<PrimerPair> initialOuterPrimerPairs;
	private Map< PrimerPair, ArrayList<PrimerPair> > secondPrimerPairs;
	private Map< PrimerPair, Map< PrimerPair, ArrayList<PrimerPair> > > thirdPrimerPairs;
	
	private String projectName;
	
	private static final float MAX_PRIMER_PENALTY = (float)0.05;
	
}



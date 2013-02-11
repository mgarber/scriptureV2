/**
 * 
 */
package broad.pda.methylation;

import java.io.*;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

//import org.apache.commons.lang3.StringUtils;
import net.sf.samtools.util.StringUtil;

import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

import net.sf.picard.util.*;

import broad.core.sequence.Sequence;
import broad.pda.methylation.MethylationUtils;

/**
 * @author engreitz
 *
 */
public class OutputPadlock {

	public static String COMMON_INSERT = "CTTCAGCTTCCCGATATCCGACGGTAGTGT";
	public static String LEFT_AMP = "AGGACCGGATCAACT";
	public static String RIGHT_AMP = "CATTGCGTGAACCGA";
	
	public static String AGILENT_COMMON_INSERT = "GTTGGAGGCTCATCGTTCCTATTCAGGCAGATGTTATCGAGGTCCGAC";
	public static String AGILENT_LEFT_AMP = "GTCATATCGGTCACTGTT";
	public static String AGILENT_RIGHT_AMP = "GATCAGGATACACACTACCC";
	
	public OutputPadlock(String in, int arraySize, boolean agilentFlag) throws FileNotFoundException {
		InputStream is = new FileInputStream(in);
		TabbedInputParser ip = new TabbedInputParser(false, is);
		Iterator<String []> itr = ip.iterator();
		
		List<String> senseOligos = new ArrayList<String>();
		List<String> antisenseOligos = new ArrayList<String>();
		
		while (itr.hasNext()) {
			List<String> oligos = new ArrayList<String>();
			String[] fields = itr.next();
			
			createPadlock(oligos, fields[5].toUpperCase(), fields[7].toUpperCase(), agilentFlag);
			
			if (oligos.size() > 4) {
				System.err.println("Found uber CpG " + fields[0] + "\t" + oligos.size() + "\t" + fields[5] + "\t" + fields[7]);
			}	
			
			// Add the oligo and its antisense to the list
			for (String oligo : oligos) {
				senseOligos.add(StringUtil.join("\t", fields) + "\t" + oligo);
				antisenseOligos.add(StringUtil.join("\t", fields) + "\t" + Sequence.reverseSequence(oligo));
				//System.out.println(StringUtils.join(fields, "\t") + "\t" + oligo);
			}	
		}
		
		// Print oligos
		List<String> toPrint = MethylationUtils.fillArray(senseOligos, antisenseOligos, arraySize);
		for (String line : toPrint)
			System.out.println(line);
	}
	
	
	/*
	 *  Creates primers for fully methylated and fully unmethylated sequence.
	 */
	private void createPadlock(List<String> oligos, String primer1, String primer2, boolean agilentFlag) {
		int i = primer1.indexOf("CG");
		int j = primer2.indexOf("CG");
		String oligo = makeOligo(primer1, primer2, false, agilentFlag);
		oligos.add(oligo);
		if ((i != -1) || (j != -1)) {
			// If a CG is present in either primer, create another methylated version
			String oligo2 = makeOligo(primer1, primer2, true, agilentFlag);
			oligos.add(oligo2);
		}
	}
	
	
	private String makeOligo(String primer1, String primer2, boolean methylated, boolean agilentFlag) {
		String bs1 = MethylationUtils.bisulfiteConvert(primer1, methylated);
		String bs2 = MethylationUtils.bisulfiteConvert(primer2, methylated);
		String leftAmp = agilentFlag ? AGILENT_LEFT_AMP : LEFT_AMP;
		String rightAmp = agilentFlag ? AGILENT_RIGHT_AMP : RIGHT_AMP;
		String commonInsert = agilentFlag ? AGILENT_COMMON_INSERT : COMMON_INSERT;
		
		return leftAmp + Sequence.reverseSequence(bs1) + commonInsert + Sequence.reverseSequence(bs2) + rightAmp;
	}
	
	
	
	/*
	 *  Creates primers with all combinations of CG conversion
	 
	private void createPadlockRecursive(List<String> oligos, String primer1, String primer2) {
		int i = primer1.indexOf("CG");
		if (i != -1) {
			createPadlockRecursive(oligos, primer1.subSequence(0,i) + "X" + primer1.subSequence(i+2, primer1.length()), primer2);
			createPadlockRecursive(oligos, primer1.subSequence(0,i) + "Y" + primer1.subSequence(i+2, primer1.length()), primer2);
		} else {
			int j = primer2.indexOf("CG");
			if (j != -1) {
				createPadlockRecursive(oligos, primer1, primer2.subSequence(0,j) + "X" + primer2.subSequence(j+2, primer2.length()));
				createPadlockRecursive(oligos, primer1, primer2.subSequence(0,j) + "Y" + primer2.subSequence(j+2, primer2.length()));
			}
		}
		
		primer1 = primer1.replaceAll("X", "CG");
		primer1 = primer1.replaceAll("Y", "TG");
		primer2 = primer2.replaceAll("X", "CG");
		primer2 = primer2.replaceAll("Y", "TG");
		String oligo = Sequence.reverseSequence(primer1) + COMMON_INSERT + Sequence.reverseSequence(primer2);
		oligos.add(oligo);
	}
	*/
	
	
	private static String USAGE = "OutputPadlock parameters:\n\n\t-in\tppDesigner output\n";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, ParseException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "write");
		
		final String in = argmap.getInput();
		final Integer arraySize = argmap.getInteger("arraySize");
		final boolean agilentFlag = argmap.containsKey("agilent");
		new OutputPadlock(in, arraySize, agilentFlag);
	}

}

package broad.core.primer3;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;

public class Primer3IO {
	//private static final String PRIMER3_CMD = "/seq/mgarber/tools/primer3-1.1.1/src/primer3_core";
	private static final String PRIMER3_CMD = "/seq/mguttman/scripts/Primer3/primer3-2.2.1-alpha/src/primer3_core";
	
	
	Process primer3Proc;
	BufferedReader primer3StdOut;
	BufferedReader primer3StdErr;
	BufferedWriter primer3StdIn;

	public Primer3IO() {
		super();
	}

	public void startPrimer3Communications() throws IOException {
		primer3Proc = Runtime.getRuntime().exec(PRIMER3_CMD);
		primer3StdIn = new BufferedWriter(new OutputStreamWriter(primer3Proc.getOutputStream()));
		primer3StdOut = new BufferedReader(new InputStreamReader(primer3Proc.getInputStream()));
		primer3StdErr = new BufferedReader(new InputStreamReader(primer3Proc.getErrorStream()));
	}
	
	public void endPrimer3Communications() throws IOException {
		primer3StdIn.close();
		primer3StdOut.close();
		primer3StdErr.close();
		primer3Proc.destroy();
	}
	
	public void write(Collection<PrimerPair> inputTags, Primer3Configuration configuration, String file) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		Iterator<PrimerPair> inputTagIt = inputTags.iterator();
		Primer3SequenceInputTags input = null;
		while(inputTagIt.hasNext()) {
			input = inputTagIt.next().getUsedInput();
			//System.out.println("Creating Record for "+input.getPrimerSequenceId());
			writeRecord(bw,input, configuration);
		}
		bw.close();
	}
	/*
	public PrimerPair parseReply(BufferedReader primer3StdOut) throws IOException {
		PrimerPair pp = new PrimerPair();
		String line = null;
		String [] lineInfo = null;
		while((line = primer3StdOut.readLine()) != null && !(line.equals("="))) {
			lineInfo = line.split("=");
			if(lineInfo.length != 2) {
				System.out.println("Primer 3 output was not a boulder record:<"+line+">");
				return null;
			} 
			if("PRIMER_ERROR".equals(lineInfo[0])) {
				System.out.println("Primer 3 Error<"+lineInfo[1]+">");
				return null;
			}
			//System.out.println(line);
			if("PRIMER_LEFT_SEQUENCE".equals(lineInfo[0])) {
				pp.setLeftPrimer(lineInfo[1]);
			} else if ("PRIMER_RIGHT_SEQUENCE".equals(lineInfo[0])) {
				pp.setRightPrimer(lineInfo[1]);
			} else if ("PRIMER_LEFT".equals(lineInfo[0])) {
				pp.setLeftPrimerPosition(Integer.parseInt(lineInfo[1].split(",")[0]));
			} else if ("PRIMER_RIGHT".equals(lineInfo[0])) {
				pp.setRightPrimerPosition(Integer.parseInt(lineInfo[1].split(",")[0]));
			} else if ("PRIMER_SEQUENCE_ID".equals(lineInfo[0])) {
				pp.setPrimerPairId(lineInfo[1]);
			} else if ("PRIMER_LEFT_TM".equals(lineInfo[0])) {
				pp.setLeftPrimerTM(Float.parseFloat(lineInfo[1]));
			} else if ("PRIMER_RIGHT_TM".equals(lineInfo[0])) {
				pp.setRightPrimerTM(Float.parseFloat(lineInfo[1]));
			} else if ("PRIMER_PRODUCT_SIZE".equals(lineInfo[0])) {
				pp.setProductSize(Integer.parseInt(lineInfo[1]));
			} else if ("PRIMER_PAIR_PENALTY".equals(lineInfo[0])) {
				pp.setPrimerPairPenalty(Integer.parseInt(lineInfo[0]));
			}
		}
		return pp;
	}*/
	
	/*public Collection<PrimerPair> parseReply(BufferedReader primer3StdOut) throws IOException {
		Stack<PrimerPair> primers=new Stack();
		String line = null;
		PrimerPair pp=new PrimerPair();
		primers.push(pp);
		String [] lineInfo = null;
		while((line = primer3StdOut.readLine()) != null && !(line.equals("="))) {
			lineInfo = line.split("=");
			if(lineInfo.length != 2) {
				System.out.println("Primer 3 output was not a boulder record:<"+line+">");
				return null;
			} 
			if("PRIMER_ERROR".equals(lineInfo[0])) {
				System.out.println("Primer 3 Error<"+lineInfo[1]+">");
				return null;
			}
			//System.out.println(line);
			
			if(lineInfo[0].startsWith("PRIMER_PAIR_PENALTY")){
				pp=new PrimerPair();
				if(!primers.isEmpty()){PrimerPair test=primers.pop(); if(test!=null && test.getLeftPrimer()!=null && test.getRightPrimer()!=null){primers.push(test);}}
				primers.push(pp);
			}
			
			pp=primers.pop();
			
			if(lineInfo[0].startsWith("PRIMER_LEFT") && lineInfo[0].endsWith("_SEQUENCE")){
				pp.setLeftPrimer(lineInfo[1]);
			}
			else if (lineInfo[0].startsWith("PRIMER_RIGHT") && lineInfo[0].endsWith("_SEQUENCE")) {
				pp.setRightPrimer(lineInfo[1]);
			} else if ("PRIMER_SEQUENCE_ID".equals(lineInfo[0])) {
				pp.setPrimerPairId(lineInfo[1]);
			} else if (lineInfo[0].startsWith("PRIMER_LEFT") && lineInfo[0].endsWith("_TM")) {
				pp.setLeftPrimerTM(Float.parseFloat(lineInfo[1]));
			} else if (lineInfo[0].startsWith("PRIMER_RIGHT") && lineInfo[0].endsWith("_TM")) {
				pp.setRightPrimerTM(Float.parseFloat(lineInfo[1]));
			} 
			
			else if (lineInfo[0].startsWith("PRIMER_PRODUCT_SIZE") && !(lineInfo[0].equalsIgnoreCase("PRIMER_PRODUCT_SIZE_RANGE"))) {
				pp.setProductSize(Integer.parseInt(lineInfo[1]));
			} else if (lineInfo[0].startsWith("PRIMER_PAIR_PENALTY")) {
				pp.setPrimerPairPenalty(Float.parseFloat(lineInfo[1]));
			}
			primers.push(pp);
		}
		return primers;
	}*/
	
	/**Updated for Primer3 upgrade version 2.0**/
	//TODO Update to new version of parsing
	public Collection<PrimerPair> parseReply(BufferedReader primer3StdOut) throws IOException {
		Map<Integer, PrimerPair> primerMap=new TreeMap();
		Stack<PrimerPair> primers=new Stack();
		String line = null;
		PrimerPair pp=new PrimerPair();
		primers.push(pp);
		String [] lineInfo = null;
		while((line = primer3StdOut.readLine()) != null && !(line.equals("="))) {
			lineInfo = line.split("=");
			if(lineInfo.length != 2) {
				System.out.println("Primer 3 output was not a boulder record:<"+line+">");
				return null;
			} 
			if("PRIMER_ERROR".equals(lineInfo[0])) {
				System.out.println("Primer 3 Error<"+lineInfo[1]+">");
				return null;
			}
			
			if(lineInfo[0].startsWith("PRIMER_")){
				String[] tokens=lineInfo[0].split("_");
				try{
					if(tokens.length>2){
					int num=Integer.parseInt(tokens[2]);
					PrimerPair pair=new PrimerPair();
					if(primerMap.containsKey(num)){pair=primerMap.get(num);}
					String primer=tokens[1];
					String tag="";
					for(int i=3; i<tokens.length; i++){
						tag+=tokens[i]+"_";
					}
					pair.setFields(primer, tag, lineInfo[1]);
					primerMap.put(num, pair);
					}else{System.err.println(lineInfo[0]+" "+lineInfo[1]);}
				}catch(NumberFormatException ex){}
			}
			
			
		}

		return primerMap.values();
	}
	
	public Collection<PrimerPair> findPrimerPair(Primer3SequenceInputTags input, Primer3Configuration configuration) throws IOException {
		//writeRecord(new BufferedWriter(new OutputStreamWriter(System.out)), input, configuration);
		writeRecord(primer3StdIn, input, configuration);
		
		try {
			Collection<PrimerPair> result = parseReply(primer3StdOut);		
			for(PrimerPair pp: result){
				pp.setUsedInput(input);
				pp.setConfigurationUsed(configuration);
			}
			return result;		
		} catch (NullPointerException e) {
			return null;
		}
		
	}
	
	public void writeRecord(BufferedWriter bw, Primer3SequenceInputTags input, Primer3Configuration conf) throws IOException {
		//writeRecordItem(bw, "PRIMER_SEQUENCE_ID",input.getPrimerSequenceId());
		writeRecordItem(bw, "SEQUENCE_ID",input.getPrimerSequenceId()); //For upgraded version of Primer3
		//writeRecordItem(bw,"SEQUENCE",input.getSequence().getSequenceBases());
		writeRecordItem(bw,"SEQUENCE_TEMPLATE",input.getSequence().getSequenceBases());//For upgraded version of Primer3
		//writeRecordItem(bw,"EXCLUDED_REGION",input.regionListToString(input.getExcludedRegions()));
		writeRecordItem(bw,"SEQUENCE_EXCLUDED_REGION",input.regionListToString(input.getExcludedRegions()));//For upgraded version of Primer3
		//writeRecordItem(bw,"INCLUDED_REGION",input.regionListToString(input.getIncludedRegions()));
		writeRecordItem(bw,"SEQUENCE_INCLUDED_REGION ",input.regionListToString(input.getIncludedRegions()));//For upgraded version of Primer3
		//writeRecordItem(bw,"TARGET",input.regionListToString(input.getTargets()));
		writeRecordItem(bw,"SEQUENCE_TARGET",input.regionListToString(input.getTargets(), 200));//For upgraded version of Primer3
		writeRecordItem(bw, "SEQUENCE_OVERLAP_JUNCTION_LIST", input.getJunctions());
		
		if(input.getPrimerStartCodonPosition() > -1000000) {
			//writeRecordItem(bw,"PRIMER_START_CODON_POSITION",input.getPrimerStartCodonPosition());
			writeRecordItem(bw,"SEQUENCE_START_CODON_POSITION",input.getPrimerStartCodonPosition());//For upgraded version of Primer3
		}
		
		//writeRecordItem(bw, "PRIMER_LEFT_INPUT",input.getPrimerLeftInput());
		//writeRecordItem(bw, "PRIMER_RIGHT_INPUT",input.getPrimerRightInput());
		
		writeRecordItem(bw, "SEQUENCE_PRIMER",input.getPrimerLeftInput());//For upgraded version of Primer3
		writeRecordItem(bw, "SEQUENCE_PRIMER_REVCOMP",input.getPrimerRightInput());//For upgraded version of Primer3

		writeRecordItem(bw, "PRIMER_OPT_TM",conf.optimalMeltingTemp);
		writeRecordItem(bw, "PRIMER_MIN_TM",conf.minMeltingTemp);
		
		writeRecordItem(bw, "PRIMER_MAX_TM",conf.maxMeltingTemp);
		
		writeRecordItem(bw, "PRIMER_LIBERAL_BASE" ,conf.interpretBasesLiberally ? 1 : 0);
		
		writeRecordItem(bw, "PRIMER_MIN_QUALITY",conf.minQualityScore);
		
		writeRecordItem(bw, "PRIMER_MIN_GC",conf.minGCContent);
		
		writeRecordItem(bw, "PRIMER_MAX_GC",conf.maxGCContent);
		
		writeRecordItem(bw, "PRIMER_GC_CLAMP", conf.useGCClamp ? String.valueOf(1) : String.valueOf(0));
		
		//writeRecordItem(bw, "PRIMER_SALT_CONC", conf.saltConcentration);
		writeRecordItem(bw, "PRIMER_SALT_MONOVALENT", conf.saltConcentration);
		
		writeRecordItem(bw, "PRIMER_DNA_CONC",conf.DNAConcentration);
		
		//writeRecordItem(bw, "PRIMER_NUM_NS_ACCEPTED", conf.numNBasesAccepted);
		writeRecordItem(bw, "PRIMER_MAX_NS_ACCEPTED", conf.numNBasesAccepted);
		
		//writeRecordItem(bw, "PRIMER_SELF_ANY",conf.selfAnyAlignScore);
		writeRecordItem(bw, "PRIMER_MAX_SELF_ANY",conf.selfAnyAlignScore);
		
		//writeRecordItem(bw, "PRIMER_SELF_END", conf.selfEndAlignScore);
		writeRecordItem(bw, "PRIMER_MAX_SELF_END", conf.selfEndAlignScore);
		
		writeRecordItem(bw, "PRIMER_OPT_SIZE",conf.optimalPrimerSize);
		
		writeRecordItem(bw, "PRIMER_MIN_SIZE", conf.minPrimerSize);
		
		writeRecordItem(bw, "PRIMER_MAX_SIZE", conf.maxPrimerSize);
		
		writeRecordItem(bw, "PRIMER_NUM_RETURN", conf.maxNumPrimersToReturn);
		
		writeRecordItem(bw, "PRIMER_PRODUCT_SIZE_RANGE",conf.minProductSize + "-" + conf.maxProductSize);
		
		writeRecordItem(bw, "PRIMER_PRODUCT_OPT_SIZE", conf.optimalProductSize);
		
		//writeRecordItem(bw, "PRIMER_MAX_DIFF_TM", conf.maxMeltingTempDiff);
		writeRecordItem(bw, "PRIMER_PAIR_MAX_DIFF_TM", conf.maxMeltingTempDiff);
		
		writeRecordItem(bw, "PRIMER_WT_TM_GT", conf.overMeltingTempPenaltyWeight);
		
		writeRecordItem(bw, "PRIMER_WT_TM_LT", conf.underMeltingTempPenaltyWeight);
		
		writeRecordItem(bw, "PRIMER_WT_GC_PERCENT_GT", conf.overGCContentPenaltyWeight);
		
		writeRecordItem(bw, "PRIMER_WT_GC_PERCENT_LT", conf.underGCContentPenaltyWeight);
		
	    writeRecordItem(bw, "PRIMER_PAIR_WT_PRODUCT_SIZE_GT", conf.underProductSizePenaltyWeight);
	    
	    writeRecordItem(bw, "PRIMER_PAIR_WT_PRODUCT_SIZE_LT", conf.overProductSizePenaltyWeight);
		
		writeRecordItem(bw, "PRIMER_WT_SIZE_GT", conf.overLengthPenaltyWeight);
		
		writeRecordItem(bw, "PRIMER_WT_SIZE_lT",conf.underLengthPenaltyWeight);
		
		writeRecordItem(bw, "PRIMER_PICK_ANYWAY", conf.canViolateConstraints ? 1 : 0);
		
		writeRecordItem(bw, "PRIMER_MISPRIMING_LIBRARY", conf.missprimingLibraryFile == null ? "" : conf.missprimingLibraryFile);
		
		writeRecordItem(bw, "PRIMER_MAX_POLY_X", conf.maxPolyX);
		
		if(input.getJunctions()!=null && !input.getJunctions().isEmpty()){
			writeRecordItem(bw, "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION", conf.primerMin3PrimeOverlapOfJunction);
			writeRecordItem(bw, "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION", conf.primerMin5PrimerOverlapOfJunction);
		}
		
		if(conf.primerNumReturn >0) {
			writeRecordItem(bw, "PRIMER_NUM_RETURN", conf.primerNumReturn);
		}
		
		
		if(conf.maxMisspriming > 0) {
			//writeRecordItem(bw, "PRIMER_MAX_MISPRIMING",String.valueOf(conf.maxMisspriming));
			writeRecordItem(bw, "PRIMER_MAX_LIBRARY_MISPRIMING",String.valueOf(conf.maxMisspriming));
		}
		
		writeRecordItem(bw, "PRIMER_EXPLAIN_FLAG", String.valueOf(1));
		
		//writeRecordItem(bw, conf.primerWindowSize;
		
		//writeRecordItem(bw, conf.maxPrimersToRetain);
		
		bw.write("=\n");
		bw.flush();
	}
	
	private void writeRecordItem(BufferedWriter bw, String key, String value) throws IOException {
		//System.out.print(" Printing key<"+key+"> value<"+value+"> ...");
		if((value != null) && (value.length() > 0)) {
			bw.write(key);
			bw.write("=");
			bw.write(value);
			bw.newLine();
		}
		//System.out.print("  Done with " + key);
	}
	
	private void writeRecordItem(BufferedWriter bw, String key, double value) throws IOException {
		//System.out.print(" Printing key<"+key+"> value<"+value+"> ...");
		if(value > 0) {
			bw.write(key);
			bw.write("=");
			bw.write(String.valueOf(value));
			bw.newLine();
		}
		//System.out.print("  Done with " + key);
	}
	private void writeRecordItem(BufferedWriter bw, String key, int value) throws IOException {
		//System.out.print(" Printing key<"+key+"> value<"+value+"> ...");
		if(value > 0) {
			bw.write(key);
			bw.write("=");
			bw.write(String.valueOf(value));
			bw.newLine();
		}
		//System.out.print("  Done with " + key);
	}

}

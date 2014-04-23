package nextgen.core.programs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import org.apache.log4j.Logger;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.fastq.FastqSequence;

public class FastqConsolidateBarcode {
	public static String USAGE = "Takes a molecular barcode and uses it to rename the read on the read fastq file " +
			"-task twoFile <Requires user to supply two fastq files - one with barcode and other without." +
			"\n\t\t-bcfile <Barcode fastq file> " +
			"\n\t\t-in <Target fastq file> " +
			"\n\t\t-out <Output file name or stdout if not specified.>"+
			"-task oneFile <Assumes barcode information is part of the same fastq file" +
			"\n\t\t-in <Target fastq file> " +
			"\n\t\t-length <Length of barcode.>"+
			"\n\t\t-offset <Distance from 5' end of the read where the barcode starts. Default = 0>"+
			"\n\t\t-out <Output file name or stdout if not specified.>";
	
	
	static Logger logger = Logger.getLogger(FastqConsolidateBarcode.class.getName());
	
	public static void main (String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "adjustbc");
		
		if(argMap.getTask().equalsIgnoreCase("twoFile")){
			runTwoFile(argMap);
		}
		else if(argMap.getTask().equalsIgnoreCase("oneFile")){
			runOneFile(argMap);
		}
		else{
			logger.error("Enter correct task"+"\n"+USAGE);
		}
		
	}
	
	private static void runTwoFile(ArgumentMap argMap) throws IOException{
		String readFqFile = argMap.getInput();
		String bcFqFile   = argMap.getMandatory("bcfile");
		
		FastqParser readFQParser = new FastqParser();
		FastqParser bcFQParser   = new FastqParser();
		
		readFQParser.start(new File(readFqFile));
		bcFQParser.start(new File(bcFqFile));
		BufferedWriter bw = argMap.getOutputWriter();
		try {
			while (readFQParser.hasNext()) {
				FastqSequence read = readFQParser.next();
				FastqSequence bc   = bcFQParser.next();
				if(read != null){
					String name = read.getName().split("\\s")[0];
					String bcName = bc.getName().split("\\s")[0];
				
				
					if(!name.equals(bcName)){
						throw new IllegalStateException("Read name and barcode name are not the same for " + name + " and " + bcName);
					}
				
					String [] sampleBCInfo = read.getName().split(":");
					String sampleBC = sampleBCInfo[sampleBCInfo.length-1];
					String barcode = bc.getSequence();
					name = name + "_#_" + sampleBC + barcode;
					//logger.debug("new readname: " + name);
					read.setName(name);
					read.write(bw);
					bw.newLine();
				}	
			}
		} finally {		
			readFQParser.close();
			bcFQParser.close();
			bw.close();
		}
	}
	
	private static void runOneFile(ArgumentMap argMap) throws IOException{
		String readFqFile = argMap.getInput();
		int length = new Integer(argMap.getMandatory("length"));
		int offset = argMap.getInteger("offset", 0);
		FastqParser readFQParser = new FastqParser();
		
		readFQParser.start(new File(readFqFile));
		BufferedWriter bw = argMap.getOutputWriter();

		try {
			while (readFQParser.hasNext()) {
				FastqSequence read = readFQParser.next();
				
				if(read != null){
					String name = read.getName().split("\\s")[0];
					String sequence = read.getSequence();
								
					String [] sampleBCInfo = read.getName().split(":");
					String sampleBC = sampleBCInfo[sampleBCInfo.length-1];
					String barcode = sequence.substring(offset,length+offset);
					name = name + "_#_" + sampleBC + barcode;
					//logger.debug("new readname: " + name);
					read.setName(name);
					read.trimStartBases(length+offset);
					read.write(bw);
					bw.newLine();
				}	
			}
		} finally {		
			readFQParser.close();
			bw.close();
		}
	}
}

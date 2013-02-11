package broad.pda.seq.fastq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

/**
 * Reads through a fastq file and writes records to different files depending on 
 * the barcode found within the start or end of the read. Reads are trimmed to remove
 * the barcode.
 * @author mgarber
 *
 */
public class FastqSplitterByBarcode {
	static Logger logger = Logger.getLogger(FastqSplitterByBarcode.class.getName());
	static private Pattern nameBarcodePattern = Pattern.compile("#[A,C,G,T,N]+");
	//static private Pattern endOfReadMamePairInfo = Pattern.compile("\\(/[0-9]+$\\)");
	static private Pattern endOfReadMamePairInfo = Pattern.compile("(/[0-9]+$)");
	
	static final String USAGE = "FastqSplitterByBarcode. Reads a fastq file and writes reads to \ndifferent files depending on their barcode: " +
			"\n\t-in <Path to the fastq file (NO standard input is supported must input a file paht)>" +
			"\n\t-outdir <Ouput directory for split files>" +
			"\n\t-bcsize <Size, in bases, of the expected barcode>" +
			"\n\t-pair2File <If paired end this specifies the path to the second read file. The barcode is assumed to be on the first file>" +
			"\n\t-atEnd <add this flag is the barcode is expected at the end, not the begginging of the read>" +
			"\n\t-toTrimSecondRead <number of bases to trim the second read>" +
			"\n\t-atEndOfSecondRead <Add this flag if the sequence to trim is at the end of the second read>"+
			"\n\t-addToReadName <Add this flag if the trimmed sequence from the second read should added to the read name's end after a _BC_ separator>" +
			"\n\t-outdir2 <Ouput directory for split files of paired read if not the same as outdir>" +
			"\n\t-minReadsToOutput <Minimum number of reads associated to a barcode in order to output its reads> " +
			"\n";
	
	public static void main (String [] args) throws Exception {
		Globals.setHeadless(true);
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "splitfile");
		
		boolean atEnd = argMap.containsKey("atEnd");
		int toTrim    = argMap.getInteger("bcsize");
		boolean isPaired = argMap.containsKey("pair2File");
		File pair1File = new File(argMap.getInput());
		int toTrim2  = argMap.containsKey("toTrimSecondRead") ?argMap.getInteger("toTrimSecondRead") :0;
		boolean atEnd2 = argMap.containsKey("atEndOfSecondRead");
		int minReads = argMap.getInteger("minReadsToOutput");
		boolean addTrimm2ToReadName =  argMap.containsKey("addToReadName");
		String outdir = argMap.getOutputDir();
		String outdir2 = argMap.containsKey("outdir2") ? argMap.get("outdir2") : outdir;
		FastqParser fqparser = new FastqParser();
		FastqParser pairFQParser = new FastqParser();
		
		BarcodeInformation bcInfo = generateBCReport(pair1File, toTrim, atEnd, minReads);
		bcInfo.write(pair1File.getAbsolutePath() + ".bc.report");
		
		String baseOutName = pair1File.getName().replace(".fq", "");
		
		if(isPaired) {
			pairFQParser.start(new File(argMap.getMandatory("pair2File")));
		}
		
		fqparser.start(pair1File);
		Map<String, BufferedWriter> bcOutMap = new HashMap<String, BufferedWriter>(); 
		Map<String, BufferedWriter> bcOutMapP2 = new HashMap<String, BufferedWriter>(); 
		int numReads = 0;
		while(fqparser.hasNext()) {
			FastqSequence fs = fqparser.next();
			FastqSequence pairedFs = null;
			if(isPaired) {
				pairedFs = pairFQParser.next();
			}
			if(fs == null) {
				logger.info("Null fastqsequence at record " + numReads);
				continue;
			}
			
			String nameBC = getNameBarcode(fs);
			numReads++;
			FastqSequence trimmed = atEnd ? fs.trimEndBases(toTrim) : fs.trimStartBases(toTrim);
			String bc = trimmed.getSequence();
			
			if(bcInfo.getInnerBarcodeTotalReads(nameBC, bc) < minReads) {
				//logger.debug("read with barcodes " + nameBC + " --- " + bc + " had only "+ bcInfo.getInnerBarcodeTotalReads(nameBC, bc) + " reads. IGNORING");
				continue;
			}
			String jointBC = toTrim > 0 ? nameBC+"_"+bc : nameBC;
			BufferedWriter bcBW = bcOutMap.get(jointBC);
			BufferedWriter bcBWP2 = bcOutMapP2.get(jointBC); 
			if(bcBW == null) {
				String fileName = "/"+ baseOutName + "."+jointBC + ".fq";
				String out1Name =  outdir + fileName;
				logger.info("new barcode found, opening " + out1Name + " for writting");
				bcBW = new BufferedWriter(new FileWriter(out1Name));
				bcOutMap.put(jointBC, bcBW);
				if(isPaired) {
					String out2Name = outdir2 + fileName;
					logger.info("opening " + out2Name + " for writting paired reads");
					bcBWP2 = new BufferedWriter(new FileWriter(out2Name));
					bcOutMapP2.put(jointBC, bcBWP2);					
				}
			}
			
			if(isPaired) {
				FastqSequence trimmed2 = atEnd2 ? pairedFs.trimEndBases(toTrim2) : pairedFs.trimStartBases(toTrim2);
				if(addTrimm2ToReadName) {
					String prefix = "_BC_".intern();
					Matcher m = endOfReadMamePairInfo.matcher(fs.getName());
					//logger.debug("Could find pattern? " + m.find());
					String newP1Name = m.replaceAll(prefix+trimmed2.getSequence()+"$1");
					fs.setName(newP1Name);

					m = endOfReadMamePairInfo.matcher(fs.getDescription());
					String newP1Descr = m.replaceAll(prefix+trimmed2.getSequence()+"$1");
					fs.setDescription(newP1Descr);
					
					m = endOfReadMamePairInfo.matcher(pairedFs.getName());
					String newP2Name = m.replaceAll(prefix+trimmed2.getSequence()+"$1");
					pairedFs.setName(newP2Name);
					
					m = endOfReadMamePairInfo.matcher(pairedFs.getDescription());
					String newP2Descr = m.replaceAll(prefix+trimmed2.getSequence()+"$1");
					pairedFs.setDescription(newP2Descr);
				}
				pairedFs.write(bcBWP2);
				bcBWP2.newLine();
			}
			fs.write(bcBW);
			bcBW.newLine();
			
			if(numReads % 500000 == 0) {
				logger.info(numReads + " processed");
			}
		}
		
		try {
			logger.debug("closing source fastq file");
			fqparser.close();
		} catch (IOException e){
			logger.error("Could not close source fastq file ",e);
		}
		
		if(isPaired) {
			try {
				logger.debug("closing source paired fastq file");
				pairFQParser.close();
			} catch (IOException e){
				logger.error("Could not close source paired source file ",e);
			}
		}
		
		for(String bc : bcOutMap.keySet()) {
			try {
				logger.debug("closing "+bc+" destination  fastq file");
				bcOutMap.get(bc).close();
				if(isPaired) {
					logger.debug("closing "+bc+" second pair destination  fastq file");
					bcOutMapP2.get(bc).close();
				}
			} catch (IOException e) {
				logger.error("Could not close output file for bc: " + bc ,e);
			}
		}
	}

	private static BarcodeInformation generateBCReport(File fqf, int innerBarcodeSize, boolean innerBCAtEnd, int minReadsToReport) throws IOException {
		FastqParser fqparser = new FastqParser();
		BarcodeInformation bi = new BarcodeInformation(minReadsToReport);
		fqparser.start(fqf);
		int numReads = 0;
		while(fqparser.hasNext()) {
			FastqSequence fs = fqparser.next();
			if(fs == null) {
				logger.info("Null fastqsequence at record " + bi.totalReads());
				continue;
			}
			numReads++;
			String nameBC = getNameBarcode(fs);
			String bc = "NA";
			if(innerBarcodeSize > 0 ) {
				FastqSequence bcSeq = innerBCAtEnd ? fs.trimEndBases(innerBarcodeSize) : fs.trimStartBases(innerBarcodeSize);
				bc = bcSeq.getSequence();
			}
			bi.add(nameBC, bc);
			if(numReads % 5000000 == 0) {
				logger.debug(numReads + " processed for report");
			}
			
		}		
		fqparser.close();
		return bi;
	}

	private static String getNameBarcode(FastqSequence fs) {
		String name = fs.getName();
		String bc = "NA";
		Matcher m = nameBarcodePattern.matcher(name);
		//logger.debug("read name: " + name);
		if(m.find()) {
			//logger.debug("read matched, extracting (" +m.start()+1 +","+m.end()+")");
			bc = name.substring(m.start()+1, m.end());
			//logger.debug("bc: " + bc);
		} else {
			//logger.debug("read did not match");
		}
		return bc;
	}
	
	public static class BarcodeInformation {
		private int minReadsToReport = 0;
		private Map<String, Integer> nameBCNumbers;
		private Map<String, Map<String, Integer>> nameInnerBCNumbers;
		private int totalReads = 0;
		
		public BarcodeInformation(int minReadsToReport) {
			nameBCNumbers = new HashMap<String, Integer>();
			nameInnerBCNumbers = new HashMap<String, Map<String, Integer>>();
			this.minReadsToReport = minReadsToReport;
		}
		
		public int totalReads() {
			return totalReads;
		}

		public int getOuterBarcodeTotalReads(String nameBC) {
			return 0;
		}

		public void add(String nameBC, String bc) {
			totalReads++;
			//logger.debug("adding: " + nameBC + ", " + bc);
			if (nameBCNumbers.containsKey(nameBC)) {
				nameBCNumbers.put(nameBC, nameBCNumbers.get(nameBC)+1);
				Map<String, Integer> innerBCMap = nameInnerBCNumbers.get(nameBC);
				if(innerBCMap == null) {
					innerBCMap = new HashMap<String, Integer>();
					innerBCMap.put(bc, 1);
					nameInnerBCNumbers.put(nameBC, innerBCMap);
				} else {
					if(innerBCMap.containsKey(bc)) {
						innerBCMap.put(bc, innerBCMap.get(bc)+1);
					} else {
						innerBCMap.put(bc, 1);
					}
				}				
				
			} else {
				nameBCNumbers.put(nameBC, 1);
				Map<String, Integer> innerMap = new HashMap<String, Integer>();
				innerMap.put(bc, 1);
				nameInnerBCNumbers.put(nameBC, innerMap);
			}
			
			//logger.debug("added counts to " + nameBC + " it has now " + nameBCNumbers.get(nameBC) + " also added inner bc " + bc + " map for " + nameBC + ": " + nameInnerBCNumbers.get(nameBC).get(bc));
			
		}

		public void write(String file) throws IOException{
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			
			bw.write("Total reads: " + totalReads);
			bw.newLine();
			
			bw.write("High count barcodes (count > " + minReadsToReport + ")");
			bw.newLine();
			List<String> sortedNameBCs = new ArrayList<String>();
			for(String nameBC : nameBCNumbers.keySet()) {
				if(nameBCNumbers.get(nameBC) > minReadsToReport) {
					sortedNameBCs.add(nameBC);
				}
			}
			
			Collections.sort(sortedNameBCs, new Comparator<String>() {
				public int compare(String o1, String o2) {
					return nameBCNumbers.get(o2) - nameBCNumbers.get(o1);
				}
			});
			
			int totalReadsReported = 0;
			for( String nBC : sortedNameBCs) {
				bw.write(nBC + "\t" + nameBCNumbers.get(nBC));
				bw.newLine();
				totalReadsReported += nameBCNumbers.get(nBC);
				List<String> sortedInnerBCs = new ArrayList<String>();
				final Map<String, Integer> innerBCMap = this.nameInnerBCNumbers.get(nBC);
				
				//logger.debug("inner map for " + nBC + "with counts " + nameBCNumbers.get(nBC)+ ": "+innerBCMap);
				for(String innerBC : innerBCMap.keySet()) {
					if(innerBCMap.get(innerBC) > minReadsToReport) {
						sortedInnerBCs.add(innerBC);
					}
				}
				
				Collections.sort(sortedInnerBCs, new Comparator<String>() {
					public int compare(String o1, String o2) {
						return innerBCMap.get(o2) - innerBCMap.get(o1);
					}
				});
				
				int totalInnerBarcodeReads = 0;
				for( String iBC : sortedInnerBCs) {
					bw.write("\t"+iBC + "\t" + innerBCMap.get(iBC));
					bw.newLine();
					totalInnerBarcodeReads += innerBCMap.get(iBC);
				}
				
				bw.write("\t#scattered reads: " + (nameBCNumbers.get(nBC) - totalInnerBarcodeReads));
				bw.newLine();
				
			}
			
			bw.write("#reported reads: "+ totalReadsReported + " #unreported: " + (totalReads - totalReadsReported));
			bw.newLine();
			
			
			bw.close();
			
		}

		public int getInnerBarcodeTotalReads(String nameBC, String bc) {
			int total = 0;
			if(this.nameBCNumbers.containsKey(nameBC) ) {
				if(nameInnerBCNumbers.get(nameBC).containsKey(bc)) {
					total = nameInnerBCNumbers.get(nameBC).get(bc);
				}
			}
			return total;
		}
		
	}

}

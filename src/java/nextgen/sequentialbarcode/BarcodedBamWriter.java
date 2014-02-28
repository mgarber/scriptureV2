package nextgen.sequentialbarcode;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.utils.CountLogger;
import nextgen.sequentialbarcode.fragmentgroup.FragmentGroup;
import nextgen.sequentialbarcode.fragmentgroup.NamedBarcodedFragmentGroup;

/**
 * Tools to add barcode tag to a bam file and write new bam file
 * @author prussell
 *
 */
public class BarcodedBamWriter {
	
	/**
	 * Sam tag for barcode attribute
	 */
	public static String BARCODES_SAM_TAG = "XB";
	
	/**
	 * Sam tag for fragment group attribute
	 */
	public static String FRAGMENT_GROUP_SAM_TAG = "XF";
	
	private static Logger logger = Logger.getLogger(BarcodedBamWriter.class.getName());
	
	/**
	 * Set the barcode attribute of a SAM record
	 * @param record SAM record
	 * @param barcodes Barcodes to set
	 */
	private static void setBarcodes(SAMRecord record, BarcodeSequence barcodes) {
		record.setAttribute(BARCODES_SAM_TAG, barcodes.toSamAttributeString());
	}
	
	/**
	 * Set the fragment group attribute of a SAM record
	 * @param record SAM record
	 * @param fragments Fragment group to set
	 */
	private static void setFragmentGroup(SAMRecord record, FragmentGroup fragments) {
		record.setAttribute(FRAGMENT_GROUP_SAM_TAG, fragments.toSamAttributeString());
	}
	
	/**
	 * Read 2 column table into a map
	 * @param tableFile Table file
	 * @return Map of first column to second column
	 * @throws IOException
	 */
	private static Map<String, String> fromTable(String tableFile) throws IOException{
		logger.info("");
		logger.info("Reading information from table " + tableFile);
		FileReader r = new FileReader(tableFile);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		Map<String, String> rtrn = new HashMap<String, String>();
		while(b.ready()) {
			s.parse(b.readLine());
			rtrn.put(s.asString(0), s.asString(1));
		}
		r.close();
		b.close();
		logger.info("Done reading table.");
		return rtrn;
	}
	
	/**
	 * Using barcode mapping, add barcode and fragment group attributes to a bam file
	 * @param inputBam Input bam file
	 * @param barcodesByReadId Map of read ID to barcode sequence
	 * @throws IOException
	 */
	private static void writeBarcodedGroupedBam(String inputBam, Map<String, BarcodeSequence> barcodesByReadId) throws IOException {
		
		String barcodedBam = getBarcodedBamFileName(inputBam);
		File barcodedBamFile = new File(barcodedBam);
		if(!barcodedBamFile.exists()) {
			writeBarcodedBam(inputBam, barcodesByReadId);
		}
		
		String outputBam = getBarcodedGroupedBamFileName(inputBam);
		
		Map<BarcodeSequence, FragmentGroup> fragmentGroupsByBarcode = getFragmentGroupsByBarcode(barcodedBam, getNameMappingFileName(barcodedBam));
		
		logger.info("");
		logger.info("Writing grouped version of " + barcodedBam + " to " + outputBam + "...");
		
		BAMFileWriter w = new BAMFileWriter(new File(outputBam));
		
		// Get header from original bam file
		SAMFileReader r = new SAMFileReader(new File(barcodedBam));
		w.setHeader(r.getFileHeader());
		
		int numDone = 0;
		// Go through the bam file in order
		SAMRecordIterator iter = r.iterator();
		while(iter.hasNext()) {
			if(numDone != 0 && numDone % 10000 == 0) {
				logger.info("Finished " + numDone + " reads.");
			}
			SAMRecord record = iter.next();
			BarcodeSequence barcodes = BarcodeSequence.fromSamRecord(record);
			if(fragmentGroupsByBarcode.containsKey(barcodes)) {
				setFragmentGroup(record, fragmentGroupsByBarcode.get(barcodes));
			}
			w.addAlignment(record);
			numDone++;
		}
		
		r.close();		
		w.close();
		
		logger.info("Done writing file.");
		

	}
	
	private static String BARCODED_BAM_SUFFIX = ".barcode.bam";
	private static String GROUPED_BAM_SUFFIX = ".grouped.bam";
	private static String NAME_MAPPING_SUFFIX = ".shortNameMapping";
	
	/**
	 * Get name of barcoded grouped bam file
	 * @param inputBam Regular or barcoded bam file
	 * @return Name of barcoded grouped bam file
	 */
	private static String getBarcodedGroupedBamFileName(String inputBam) {
		return getBarcodedBamFileName(inputBam).replaceAll(".bam", "") + GROUPED_BAM_SUFFIX;
	}
	
	/**
	 * Get name of barcoded bam file
	 * @param inputBam Regular bam file
	 * @return Name of barcoded bam file
	 */
	public static String getBarcodedBamFileName(String inputBam) {
		return inputBam.replaceAll(".bam", "") + BARCODED_BAM_SUFFIX;
	}
	
	/**
	 * Find out if a barcoded bam file exists for the bam file
	 * @param regularBam Regular bam file
	 * @return True iff corresponding barcoded bam file exists in same directory
	 */
	public static boolean barcodedBamExists(String regularBam) {
		return new File(getBarcodedBamFileName(regularBam)).exists();
	}
	
	/**
	 * Get name of names mapping file, a table with mapping from full read name to shortened read name used in the grouped bam file
	 * @param inputBam Regular or barcoded bam file
	 * @return Names mapping file name
	 */
	private static String getNameMappingFileName(String inputBam) {
		String s1 = inputBam.replaceAll(BARCODED_BAM_SUFFIX, "");
		return s1.replaceAll(".bam", "") + NAME_MAPPING_SUFFIX;
	}
	
	
	/**
	 * Add barcode attribute to bam file entries using a barcode mapping
	 * @param inputBam Regular bam file
	 * @param barcodesByReadId Table of read name and barcode sequence
	 * @throws IOException
	 */
	public static void writeBarcodedBam(String inputBam, String barcodeTable) throws IOException {
		writeBarcodedBam(inputBam, readBarcodesFromTable(barcodeTable));
	}
	
	/**
	 * Add barcode attribute to bam file entries using a barcode mapping
	 * @param inputBam Regular bam file
	 * @param barcodesByReadId Map of read name to barcode sequence
	 * @throws IOException
	 */
	private static void writeBarcodedBam(String inputBam, Map<String, BarcodeSequence> barcodesByReadId) throws IOException {
		
		String outputBam = getBarcodedBamFileName(inputBam);
		String outputNames = getNameMappingFileName(outputBam);
		
		logger.info("");
		logger.info("Writing barcoded version of " + inputBam + " to " + outputBam + "...");
		logger.info("Writing mapping of read names to shortened IDs to " + outputNames + "...");
		
		BAMFileWriter w = new BAMFileWriter(new File(outputBam));
		FileWriter nw = new FileWriter(outputNames);
		
		// Get header from original bam file
		SAMFileReader r = new SAMFileReader(new File(inputBam));
		w.setHeader(r.getFileHeader());
		
		int numDone = 0;
		int skipped = 0;
		int unmapped = 0;
		// Go through the bam file in order
		SAMRecordIterator iter = r.iterator();
		int recordNum = 0;
		while(iter.hasNext()) {
			numDone++;
			if(numDone != 0 && numDone % 100000 == 0) {
				logger.info("Finished " + numDone + " reads. Skipped " + unmapped + " unmapped reads and " + skipped + " reads with barcodes not in map.");
			}
			SAMRecord record = iter.next();
			String oldName = record.getReadName();
			String newName = "r" + recordNum;
			nw.write(oldName + "\t" + newName + "\n");
			recordNum++;
			
			if(record.getReadUnmappedFlag()) {
				unmapped++;
				continue;
			}
			if(!barcodesByReadId.containsKey(oldName)) {
				logger.debug("READ_NOT_FOUND\t" + oldName);
				skipped++;
				continue;
			}
			BarcodeSequence barcodes = barcodesByReadId.get(oldName);
			setBarcodes(record, barcodes);
			w.addAlignment(record);
		}
		
		r.close();		
		w.close();
		nw.close();
		
		logger.info("Done writing file.");
				
	}
	
	/**
	 * Get mapping of read ID to barcode sequence from a table
	 * @param tableFile Table file with line format: read_ID   barcode_sequence_as_SAM_attribute
	 * @return Mapping of read ID to barcode in the read
	 * @throws IOException
	 */
	private static Map<String, BarcodeSequence> readBarcodesFromTable(String tableFile) throws IOException {
		
		logger.info("");
		logger.info("Reading fragment barcodes from file " + tableFile + "...");
		
		int numDone = 0;
		Map<String, BarcodeSequence> rtrn = new HashMap<String, BarcodeSequence>();
		FileReader r = new FileReader(tableFile);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			numDone++;
			if(numDone % 100000 == 0) {
				logger.info("Finished " + numDone);
			}
			s.parse(b.readLine());
			String name = s.asString(0);
			rtrn.put(name, BarcodeSequence.fromSamAttributeString(s.asString(1)));
			logger.debug("GOT_READ\t" + name);
		}
		r.close();
		b.close();
		logger.info("Done reading from table.");
		return rtrn;
	}
	
	/**
	 * Get name of fragment group table based on barcode table file name
	 * @param barcodeTable Barcode table
	 * @return Name of fragment group table
	 */
	private static String getFragmentGroupTableName(String barcodeTable) {
		return barcodeTable.replaceAll(".bam", "") + ".fragment_groups";
	}
	
	/**
	 * Read barcoded bam file and determine fragment groups
	 * @param barcodedBam Barcoded bam file
	 * @param shortenedReadNameTable Table with shortened read names. Line format: read_ID   shortened_read_name
	 * @return Map of barcode sequence to the group of fragments sharing those barcodes
	 * @throws IOException
	 */
	private static Map<BarcodeSequence, FragmentGroup> getFragmentGroupsByBarcode(String barcodedBam, String shortenedReadNameTable) throws IOException {
		logger.info("");
		logger.info("Reading fragments by barcode from " + barcodedBam);
		
		String tableFile = getFragmentGroupTableName(barcodedBam);
		if(new File(tableFile).exists()) {
			return readFragmentGroupsFromTable(tableFile);
		}
		
		Map<String,String> shortNames = fromTable(shortenedReadNameTable);
		
		Map<BarcodeSequence, FragmentGroup> rtrn = new HashMap<BarcodeSequence, FragmentGroup>();
		SAMFileReader r = new SAMFileReader(new File(barcodedBam));
		SAMRecordIterator iter = r.iterator();
		int numDone = 0;
		while(iter.hasNext()) {
			numDone++;
			if(numDone % 100000 == 0) {
				logger.info("Read " + numDone + " reads. There are " + rtrn.size() + " different barcode signatures.");
			}
			SAMRecord record = iter.next();
			@SuppressWarnings("unused")
			String dummy = record.getSAMString(); // This somehow prevents a picard bug from surfacing. Do not remove.
			record.setReadName(shortNames.get(record.getReadName()));
			BarcodedFragment fragment = new BarcodedFragmentImpl(record);
			BarcodeSequence b = fragment.getBarcodes();
			if(!rtrn.containsKey(b)) {
				rtrn.put(b, new NamedBarcodedFragmentGroup(b));
			}
			rtrn.get(b).addFragment(fragment);
			record = null;
		}
		logger.info("Got fragments for " + rtrn.size() + " different barcode signatures.");
		r.close();
		
		logger.info("");
		logger.info("Writing fragments by barcode in " + barcodedBam + " to file " + tableFile);
		FileWriter w = new FileWriter(tableFile);
		CountLogger cl = new CountLogger(rtrn.size(), 100);
		for(BarcodeSequence barcodes : rtrn.keySet()) {
			cl.advance();
			w.write(barcodes.toSamAttributeString() + "\t" + rtrn.get(barcodes).toSamAttributeString() + "\n");
		}
		w.close();
		logger.info("Done writing file.");
		
		return rtrn;
	}
	
	/**
	 * Read fragment groups from table
	 * @param fragmentGroupsTable Table with line format: barcode_sequence_sam_attribute_string  fragment_group_sam_attribute_string
	 * @return Map of barcode sequence to the group of fragments with the barcodes
	 * @throws IOException
	 */
	private static Map<BarcodeSequence, FragmentGroup> readFragmentGroupsFromTable(String fragmentGroupsTable) throws IOException {
		logger.info("");
		logger.info("Reading fragment groups from table " + fragmentGroupsTable);
		Map<BarcodeSequence, FragmentGroup> rtrn = new HashMap<BarcodeSequence, FragmentGroup>();
		FileReader r = new FileReader(fragmentGroupsTable);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		int numDone = 0;
		while(b.ready()) {
			numDone++;
			if(numDone % 100000 == 0) {
				logger.info("Finished " + numDone + " barcodes.");
			}
			s.parse(b.readLine());
			if(s.getFieldCount() != 2) {
				b.close();
				throw new IllegalArgumentException("Line format: barcode_signature   fragment_group");
			}
			BarcodeSequence barcodes = BarcodeSequence.fromSamAttributeString(s.asString(0));
			FragmentGroup fragments = NamedBarcodedFragmentGroup.fromSamAttributeStrings(s.asString(0), s.asString(1));
			rtrn.put(barcodes, fragments);
		}
		r.close();
		b.close();
		logger.info("Got fragment groups for " + rtrn.size() + " barcodes.");
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-ib", "Input bam file", true);
		p.addStringArg("-bt", "Table of barcodes by read name (required if barcoded bam file does not exist)", false, null);
		p.addBooleanArg("-d", "Debug logging on", false, false);
		p.parse(args);
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
		}
		
		String inputBam = p.getStringArg("-ib");
		String barcodedBam = getBarcodedBamFileName(inputBam);
		String barcodeTable = p.getStringArg("-bt");
		File tableFile = new File(barcodeTable);
		File barcodeBamFile = new File(barcodedBam);
		if(!barcodeBamFile.exists() && barcodeTable == null) {
			throw new IllegalArgumentException("Must provide barcode table");
		}
		
		writeBarcodedGroupedBam(inputBam, tableFile.exists() ? readBarcodesFromTable(barcodeTable) : null);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
	
	
}

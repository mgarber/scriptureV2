package nextgen.sequentialbarcode.database;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.sequentialbarcode.BarcodedBamWriter;
import nextgen.sequentialbarcode.BarcodedFragmentImpl;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;

public class DatabaseWriter {
	
	private static Logger logger = Logger.getLogger(DatabaseWriter.class.getName());
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-ib", "Input bam file", true);
		p.addStringArg("-bt", "Table of barcodes by read name (required if barcoded bam file does not exist)", false, null);
		p.addStringArg("-dbh", "Berkeley DB environment home directory", true);
		p.addStringArg("-dbs", "Database entity store name", true);
		p.addBooleanArg("-d", "Debug logging on", false, false);
		p.parse(args);
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
		}
		
		// Write barcoded bam file if necessary
		String inputBam = p.getStringArg("-ib");
		String barcodeTable = p.getStringArg("-bt");
		if(!BarcodedBamWriter.barcodedBamExists(inputBam)) {
			if(barcodeTable == null) {
				throw new IllegalArgumentException("Must provide table of barcodes by read name with -bt option");
			}
			BarcodedBamWriter.writeBarcodedBam(inputBam, barcodeTable);
		}
		
		// Write to database
		String envHome = p.getStringArg("-dbh");
		String storeName = p.getStringArg("-dbs");
		BarcodedFragmentImpl.DataAccessor dataAccessor = BarcodedFragmentImpl.getDataAccessor(envHome, storeName, false);
		
		// Iterate through bam file and enter into database
		logger.info("");
		logger.info("Writing barcoded SAM records to database.");
		String barcodedBam = BarcodedBamWriter.getBarcodedBamFileName(inputBam);
		SAMFileReader samReader = new SAMFileReader(new File(barcodedBam));
		SAMRecordIterator iter = samReader.iterator();
		int numDone = 0;
		
		while(iter.hasNext()) {
			numDone++;
			if(numDone % 100000 == 0) {
				logger.info("Finished entering " + numDone + " records.");
			}
			SAMRecord record = iter.next();
			BarcodedFragmentImpl fragment = new BarcodedFragmentImpl(record);
			dataAccessor.put(fragment);
		}
		
		// Close data accessor and sam reader
		dataAccessor.close();
		samReader.close();
		
		logger.info("");
		logger.info("All done.");


	}

}

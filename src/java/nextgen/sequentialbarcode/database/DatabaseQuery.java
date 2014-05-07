package nextgen.sequentialbarcode.database;


import org.apache.log4j.Logger;

import com.sleepycat.persist.EntityCursor;

import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.berkeleydb.JoinedEntityCursor;
import nextgen.sequentialbarcode.BarcodedFragmentImpl;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

public class DatabaseQuery {
	
	private static Logger logger = Logger.getLogger(DatabaseQuery.class.getName());
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-e", "Database environment home directory", true);
		p.addStringArg("-s", "Database store name", true);
		p.addStringArg("-b", "Barcode sequence", false, null);
		p.addStringArg("-in", "Interval (chr:start-end)", false, null);
		p.addStringArg("-c", "Entire chromosome", false, null);
		p.addStringArg("-bb", "Barcoded bam file", true);
		p.addIntArg("-mb", "Min number of barcodes to print fragment", false, 0);
		p.addBooleanArg("-es", "For interval or chromosome, exclude reads in the region queried", false, false);
		p.parse(args);
		String envHome = p.getStringArg("-e");
		int minBarcodes = p.getIntArg("-mb");
		String storeName = p.getStringArg("-s");
		String barcodes = p.getStringArg("-b");
		String interval = p.getStringArg("-in");
		String bam = p.getStringArg("-bb");
		String wholeChr = p.getStringArg("-c");
		boolean excludeSelf = p.getBooleanArg("-es");
		
		if(interval != null && barcodes != null) {
			throw new IllegalArgumentException("Choose one: -b or -in");
		}
		
		BarcodedFragmentImpl.DataAccessor dataAccessor = BarcodedFragmentImpl.getDataAccessor(envHome, storeName, true, false);
		
		// Single barcode sequence
		if(barcodes != null) {
			EntityCursor<BarcodedFragmentImpl> fragments = dataAccessor.getAllFragmentsWithBarcodes(barcodes);
			logger.info("");
			logger.info("All fragments with barcodes " + barcodes + ":");
			System.out.println("");
			System.out.println("******************************************************");
			System.out.println("");
			for(BarcodedFragmentImpl fragment : fragments) {
				System.out.println(fragment.getId() + "\t" + fragment.getMappedLocation().toUCSC() + "\t" + fragment.getBarcodes().toSamAttributeString() + "\t" + fragment.getUnpairedSequence());			
			}
			fragments.close();
			System.out.println("");
			System.out.println("******************************************************");
			System.out.println("");
		}
		
		// Single window
		if(interval != null) {
			StringParser s = new StringParser();
			s.parse(interval,":");
			if(s.getFieldCount() != 2) {
				throw new IllegalArgumentException("Interval must be in format chr:start-end");
			}
			String chr = s.asString(0);
			String coords = s.asString(1);
			s.parse(coords,"-");
			if(s.getFieldCount() != 2) {
				throw new IllegalArgumentException("Interval must be in format chr:start-end");
			}
			int start = s.asInt(0);
			int end = s.asInt(1);
			JoinedEntityCursor<BarcodedFragmentImpl> iter = dataAccessor.getAllFragmentsWithBarcodesMatchingFragmentInWindow(bam, chr, start, end, false);
			logger.info("");
			logger.info("All fragments with barcodes matching a read mapped to " + interval + " (at least " + minBarcodes + " barcodes):");
			System.out.println("");
			System.out.println("******************************************************");
			System.out.println("");
			while(iter.hasNext()) {
				BarcodedFragmentImpl fragment = iter.next();
				if(excludeSelf) {
					if(fragment.getMappedLocation().overlaps(new BasicAnnotation(chr, start, end))) {
						continue;
					}
				}
				if(fragment.getNumBarcodes() >= minBarcodes) {
					System.out.println(fragment.getId() + "\t" + fragment.getMappedLocation().toUCSC() + "\t" + fragment.getBarcodes().toSamAttributeString() + "\t" + fragment.getUnpairedSequence());
				}
			}
			System.out.println("");
			System.out.println("******************************************************");
			System.out.println("");
			iter.close();
		}
		
		// Entire chromsome
		if(wholeChr != null) {
			JoinedEntityCursor<BarcodedFragmentImpl> iter = dataAccessor.getAllFragmentsWithBarcodesMatchingFragmentInChr(bam, wholeChr);
			logger.info("");
			logger.info("All fragments with barcodes matching a read mapped to " + wholeChr + " (at least " + minBarcodes + " barcodes):");
			System.out.println("");
			System.out.println("******************************************************");
			System.out.println("");
			while(iter.hasNext()) {
				BarcodedFragmentImpl fragment = iter.next();
				if(excludeSelf) {
					if(fragment.getMappedLocation().getChr().equals(wholeChr)) {
						continue;
					}
				}
				if(fragment.getNumBarcodes() >= minBarcodes) {
					System.out.println(fragment.getId() + "\t" + fragment.getMappedLocation().toUCSC() + "\t" + fragment.getBarcodes().toSamAttributeString() + "\t" + fragment.getUnpairedSequence());
				}
			}
			System.out.println("");
			System.out.println("******************************************************");
			System.out.println("");
			iter.close();
		}

		
		dataAccessor.close();
		
		logger.info("");
		logger.info("All done.");
		
	}

}

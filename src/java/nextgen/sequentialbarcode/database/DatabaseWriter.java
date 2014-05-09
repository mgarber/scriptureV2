package nextgen.sequentialbarcode.database;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.job.OGSJob;
import nextgen.core.pipeline.util.BamUtils;
import nextgen.core.pipeline.util.OGSUtils;
import nextgen.sequentialbarcode.BarcodedBamWriter;
import nextgen.sequentialbarcode.BarcodedFragmentImpl;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import broad.core.parser.CommandLineParser;

public class DatabaseWriter {
	
	private class WriterThread implements Runnable {
		
		private BarcodedFragmentImpl.DataAccessor dataAccessor;
		private String barcodedBam;
		
		protected WriterThread(BarcodedFragmentImpl.DataAccessor accessor, String barcodeBam) {
			dataAccessor = accessor;
			barcodedBam = barcodeBam;
		}
		
		@Override
		public void run() {
			SAMFileReader samReader = new SAMFileReader(new File(barcodedBam));
			SAMRecordIterator iter = samReader.iterator();
			int numDone = 0;
			
			logger.info("");
			logger.info("Writing barcoded SAM records to database.");
			
			long time = System.currentTimeMillis();
			
			while(iter.hasNext()) {
				numDone++;
				if(numDone % 100000 == 0) {
					long newTime = System.currentTimeMillis();
					long diff = (newTime - time) / 1000;
					time = newTime;
					logger.info("Finished entering " + numDone + " records (" + diff + " seconds).");
				}
				SAMRecord record = iter.next();
				BarcodedFragmentImpl fragment = new BarcodedFragmentImpl(record);
				dataAccessor.put(fragment);
			}
			
			// Close data accessor and sam reader
			dataAccessor.close();
			samReader.close();			
		}
		
	}
	
	private static Logger logger = Logger.getLogger(DatabaseWriter.class.getName());
	private static Session drmaaSession;

	/**
	 * @param args
	 * @throws IOException 
	 * @throws DrmaaException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, DrmaaException, InterruptedException {
		
		drmaaSession = OGSUtils.getDrmaaSession();
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-ib", "Input bam file", true);
		p.addStringArg("-bt", "Table of barcodes by read name (required if barcoded bam file does not exist)", false, null);
		p.addStringArg("-dbh", "Berkeley DB environment home directory", true);
		p.addStringArg("-dbs", "Database entity store name", true);
		p.addBooleanArg("-d", "Debug logging on", false, false);
		p.addIntArg("-rj", "When batching out writing of barcoded bam file, number of reads per job", false, 10000000);
		p.addStringArg("-bbj", "Barcoded bam writer jar file (needed to write barcoded bam file)", false, null);
		p.addStringArg("-pj", "Picard jar directory (needed to write barcoded bam file)", false, null);
		p.addBooleanArg("-tw", "Multithread writing of database", false, false);
		p.addIntArg("-twn", "When multithreading writing of database, number of bam files to split into", false, 10);
		p.addStringArg("-bsj", "Bam splitter jar file for multithreading on smaller bam files, required if smaller files do not already exist", false, null);
		p.parse(args);
		
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
		}
		
		String inputBam = p.getStringArg("-ib");

		logger.info("");
		logger.info("Writing berkeleydb database for " + inputBam + ".");
		
		// Write barcoded bam file if necessary
		String barcodeTable = p.getStringArg("-bt");
		int readsPerJob = p.getIntArg("-rj");
		String barcodedBamWriterJar = p.getStringArg("-bbj");
		String picardJarDir = p.getStringArg("-pj");
		if(!BarcodedBamWriter.barcodedBamExists(inputBam)) {
			if(barcodeTable == null) {
				throw new IllegalArgumentException("Must provide table of barcodes by read name with -bt option");
			}
			if(barcodedBamWriterJar == null) {
				throw new IllegalArgumentException("Must provide BarcodedBamWriter.jar file with -bbj option");
			}
			if(picardJarDir == null) {
				throw new IllegalArgumentException("Must provide Picard jar directory with -pj option");
			}
			drmaaSession = OGSUtils.getDrmaaSession();
			BarcodedBamWriter.batchWriteBarcodedBam(inputBam, barcodeTable, drmaaSession, readsPerJob, barcodedBamWriterJar, picardJarDir);
		}
		
		// Write to database
		String envHome = p.getStringArg("-dbh");
		String storeName = p.getStringArg("-dbs");
		File e = new File(envHome);
		e.mkdir();

		// Make smaller bam files if requested
		if(p.getBooleanArg("-tw")) {
			int numBams = p.getIntArg("-twn");
			String barcodedBam = BarcodedBamWriter.getBarcodedBamFileName(inputBam);
			String splitterJar = p.getStringArg("-bsj");
			// Split bam file
			// Get names of split files
			Collection<String> smallBams = BamUtils.splitBam(BarcodedBamWriter.getBarcodedBamFileName(inputBam), numBams, true);
			// Check if all small bam files already exist
			boolean allExist = true;
			for(String smallBam : smallBams) {
				String barcodedSmallBam = BarcodedBamWriter.getBarcodedBamFileName(smallBam);
				File f = new File(barcodedSmallBam);
				if(!f.exists()) {
					allExist = false;
					break;
				}
			}
			if(allExist) {
				logger.warn("");
				logger.warn("All small bam files already exist. Not rewriting.");
			} else {
				// Submit bam splitter to cluster
				if(splitterJar == null) {
					throw new IllegalArgumentException("Must provide bam splitter jar file to batch out database writing");
				}
				String splitterCmmd = "java -jar -Xmx29g -Xms15g -Xmn10g " + splitterJar + " -i " + barcodedBam + " -n " + numBams;
				OGSJob splitterJob = new OGSJob(drmaaSession, splitterCmmd);
				splitterJob.submit();
				splitterJob.waitFor();
			}
			// Create threads
			Collection<Thread> threads = new ArrayList<Thread>();
			for(String smallBam : smallBams) {
				File oldFile = new File(smallBam);
				String barcodedSmallBam = BarcodedBamWriter.getBarcodedBamFileName(smallBam);
				File newFile = new File(barcodedSmallBam);
				if(!newFile.exists()) {
					boolean success = oldFile.renameTo(newFile);
					if(!success) {
						throw new IllegalStateException("Could not rename " + smallBam + " to " + barcodedSmallBam);
					}
				}
				BarcodedFragmentImpl.DataAccessor dataAccessor = BarcodedFragmentImpl.getDataAccessor(envHome, storeName, false, true);
				dataAccessor.setCachePercent(57 / smallBams.size());
				WriterThread wt = new DatabaseWriter().new WriterThread(dataAccessor, barcodedSmallBam);
				Thread thread = new Thread(wt);
				thread.start();
				threads.add(thread);
			}
			for(Thread t : threads) {
				t.join();
			}
		} else {
			BarcodedFragmentImpl.DataAccessor dataAccessor = BarcodedFragmentImpl.getDataAccessor(envHome, storeName, false, false);
			String barcodedBam = BarcodedBamWriter.getBarcodedBamFileName(inputBam);
			WriterThread wt = new DatabaseWriter().new WriterThread(dataAccessor, barcodedBam);
			Thread thread = new Thread(wt);
			thread.start();
			thread.join();
		}
		
		logger.info("");
		logger.info("All done.");


	}

}

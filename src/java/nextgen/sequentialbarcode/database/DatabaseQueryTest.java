package nextgen.sequentialbarcode.database;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.apache.log4j.Logger;

import com.sleepycat.persist.EntityCursor;

import nextgen.core.annotation.Gene;
import nextgen.core.berkeleydb.JoinedEntityCursor;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.feature.Window;
import nextgen.core.model.AlignmentModel;
import nextgen.sequentialbarcode.BarcodedFragmentImpl;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

public class DatabaseQueryTest {
	
	private static Logger logger = Logger.getLogger(DatabaseQueryTest.class.getName());
	
	/**
	 * Concurrent query to individual genes
	 * @author prussell
	 *
	 */
	private class GeneQuery implements Runnable {
		
		private ConcurrentLinkedQueue<Gene> geneQueue;
		private BarcodedFragmentImpl.DataAccessor dataAccessor;
		private AlignmentModel alignmentModel;
		
		/**
		 * @param name Identifier for this object
		 * @param data Alignment data
		 * @param genes Concurrent linked queue populated with the genes
		 * @param environmentHome Database environment home
		 * @param storeName Entity store name
		 * @param readOnly Database is read only
		 */
		public GeneQuery(AlignmentModel data, ConcurrentLinkedQueue<Gene> genes, String environmentHome, String storeName, boolean readOnly) {
			geneQueue = genes;
			alignmentModel = data;
			dataAccessor = BarcodedFragmentImpl.getDataAccessor(environmentHome, storeName, true);
		}
		
		@Override
		public void run() {
			while(!geneQueue.isEmpty()) {
				Gene gene = geneQueue.poll();
				JoinedEntityCursor<BarcodedFragmentImpl> iter = dataAccessor.getAllFragmentsWithBarcodesMatchingFragmentInWindow(alignmentModel, gene, false);
				if(iter == null) {
					continue;
				}
				int count = 0;
				while(iter.hasNext()) {
					@SuppressWarnings("unused")
					BarcodedFragmentImpl fragment = iter.next();
					count++;
				}
				logger.info(gene.getName() + "\t" + gene.toUCSC() + "\t" + count);
			}
			dataAccessor.close();
		}
		
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-e", "Database environment home", true);
		p.addStringArg("-s", "Database store name", true);
		p.addStringArg("-id", "Fragment ID to get", false, null);
		p.addStringArg("-b", "Barcodes to get (SAM attribute format)", false, null);
		p.addStringArg("-chr", "Interval chromosome", false, null);
		p.addIntArg("-start", "Interval start", false, -1);
		p.addIntArg("-end", "Interval end", false, -1);
		p.addStringArg("-bb", "Barcoded bam file to query interval or genome wide windows", false, null);
		p.addIntArg("-w", "Window size. Scan non-overlapping windows over the whole genome, get all reads with same barcodes as reads in window, print summary. Leave out to query whole genes.", false, -1);
		p.addStringArg("-cs", "Chromosome size file for scanning windows in genomic space", false, null);
		p.addStringArg("-bd", "Bed file for transcriptome space", false, null);
		p.addIntArg("-nt", "Number of threads for whole genes", false, 1);
		p.parse(args);
		String envHome = p.getStringArg("-e");
		String storeName = p.getStringArg("-s");
		String id = p.getStringArg("-id");
		String barcodes = p.getStringArg("-b");
		String chr = p.getStringArg("-chr");
		int start = p.getIntArg("-start");
		int end = p.getIntArg("-end");
		String bam = p.getStringArg("-bb");
		int window = p.getIntArg("-w");
		String chrSizeFile = p.getStringArg("-cs");
		String bedFile = p.getStringArg("-bd");
		int numThreads = p.getIntArg("-nt");
		
		BarcodedFragmentImpl.DataAccessor dataAccessor = BarcodedFragmentImpl.getDataAccessor(envHome, storeName, true);
		
		CoordinateSpace coordSpace = null;
		if(chrSizeFile != null && bedFile != null) {
			throw new IllegalArgumentException("Must choose genomic space or transcriptome space.");
		}
		if(chrSizeFile != null) {
			coordSpace = new GenomicSpace(chrSizeFile);
		}
		if(bedFile != null) {
			coordSpace = new TranscriptomeSpace(BEDFileParser.loadDataByChr(bedFile));
		}
		
		AlignmentModel alignmentModel = null;
		if(bam != null) {
			alignmentModel = new AlignmentModel(bam, coordSpace, false);
		}
		
		// Whole genes with multithreading
		if(bedFile != null && window <= 0 && bam != null) {
			
			logger.info("");
			logger.info("Scanning all individual genes and counting number of genome-wide fragments sharing barcodes with fragments in each gene...");
			
			Collection<Gene> genes = BEDFileParser.loadData(new File(bedFile));
			ConcurrentLinkedQueue<Gene> geneQueue = new ConcurrentLinkedQueue<Gene>(genes);
			
			Collection<Thread> threads = new ArrayList<Thread>();
			
			for(int i = 0; i < numThreads; i++) {
				AlignmentModel model = new AlignmentModel(bam, coordSpace, false);
				GeneQuery q = (new DatabaseQueryTest()).new GeneQuery(model, geneQueue, envHome, storeName, true);
				Thread t = new Thread(q);
				threads.add(t);
				t.start();
			}
			
			for(Thread t : threads) {
				t.join();
			}
			
		}
		
		// Windows
		if((chrSizeFile != null || bedFile != null) && window > 0 && bam != null) {
			
			logger.info("");
			logger.info("Scanning " + window + "bp non-overlapping windows and counting number of genome-wide fragments sharing barcodes with fragments in each window...");
			
			Iterator<? extends Window> windowIterator = coordSpace.getWindowIterator(window, 0);
			while(windowIterator.hasNext()) {
				Window w = windowIterator.next();
				JoinedEntityCursor<BarcodedFragmentImpl> iter = dataAccessor.getAllFragmentsWithBarcodesMatchingFragmentInWindow(alignmentModel, w, false);
				if(iter == null) {
					continue;
				}
				int count = 0;
				while(iter.hasNext()) {
					@SuppressWarnings("unused")
					BarcodedFragmentImpl fragment = iter.next();
					count++;
				}
				try {
					GeneWindow gw = (GeneWindow)w;
					logger.info(gw.getSourceAnnotations().iterator().next().getName() + "\t" + w.toUCSC() + "\t" + count);
				} catch(Exception e) {
					logger.info(w.toUCSC() + "\t" + count);
				}
			}
			
		}
		
		// Single gene
		if(id != null) {
			BarcodedFragmentImpl fragment = dataAccessor.getFragmentByID(id);
			logger.info("");
			logger.info("Fragment with ID " + id + ":");
			logger.info(fragment.getId() + "\t" + fragment.getBarcodes().toSamAttributeString());
		}
		
		// Single barcode sequence
		if(barcodes != null) {
			EntityCursor<BarcodedFragmentImpl> fragments = dataAccessor.getAllFragmentsWithBarcodes(barcodes);
			logger.info("");
			logger.info("All fragments with barcodes " + barcodes + ":");
			for(BarcodedFragmentImpl fragment : fragments) {
				logger.info(fragment.getId() + "\t" + fragment.getBarcodes().toSamAttributeString());
			}
			fragments.close();
		}
		
		// Single window
		if(chr != null && bam != null && start >= 0 && end > 0) {
			JoinedEntityCursor<BarcodedFragmentImpl> iter = dataAccessor.getAllFragmentsWithBarcodesMatchingFragmentInWindow(bam, chr, start, end, false);
			while(iter.hasNext()) {
				BarcodedFragmentImpl fragment = iter.next();
				logger.info(fragment.getId() + "\t" + fragment.getMappedLocation().toUCSC() + "\t" + fragment.getBarcodes().toSamAttributeString());
			}
			iter.close();
		}
		
		dataAccessor.close();
		
		logger.info("");
		logger.info("All done.");
		
	}

}

package nextgen.sequentialbarcode;


import static com.sleepycat.persist.model.Relationship.MANY_TO_ONE;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.log4j.Logger;

import com.sleepycat.je.DatabaseException;
import com.sleepycat.persist.EntityCursor;
import com.sleepycat.persist.PrimaryIndex;
import com.sleepycat.persist.SecondaryIndex;
import com.sleepycat.persist.model.Entity;
import com.sleepycat.persist.model.PrimaryKey;
import com.sleepycat.persist.model.SecondaryKey;

import broad.core.parser.StringParser;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.berkeleydb.DatabaseEnvironment;
import nextgen.core.berkeleydb.DatabaseStore;
import nextgen.core.berkeleydb.JoinedEntityCursor;
import nextgen.core.model.AlignmentModel;
import nextgen.sequentialbarcode.fragmentgroup.FragmentGroup;
import nextgen.sequentialbarcode.fragmentgroup.NamedBarcodedFragmentGroup;
import nextgen.sequentialbarcode.readlayout.Barcode;
import nextgen.sequentialbarcode.readlayout.BarcodeSet;
import nextgen.sequentialbarcode.readlayout.ReadLayout;
import nextgen.sequentialbarcode.readlayout.ReadSequenceElement;

/**
 * A basic implementation of a barcoded fragment
 * @author prussell
 *
 */
@Entity(version=0)
public class BarcodedFragmentImpl implements BarcodedFragment {
	@PrimaryKey
	protected String id;
	protected String read1sequence;
	protected String read2sequence;
	protected String unpairedSequence;
	protected Annotation location;
	protected ReadLayout read1layout;
	protected ReadLayout read2layout;
	protected BarcodeSequence barcodes;
	@SecondaryKey(relate=MANY_TO_ONE)
	private String barcodeString;
	public static Logger logger = Logger.getLogger(BarcodedFragmentImpl.class.getName());
	protected int barcodeMaxMismatches;
	protected FragmentGroup fragmentGroup;
	
	/**
	 * @param fragmentId Fragment ID
	 * @param barcodeSignature Barcodes for the fragment
	 * @param mappedChr Mapped chromosome for the fragment
	 * @param mappedStart Mapped start
	 * @param mappedEnd Mapped end
	 */
	public BarcodedFragmentImpl(String fragmentId, BarcodeSequence barcodeSignature, String mappedChr, int mappedStart, int mappedEnd) {
		this(fragmentId, barcodeSignature, new BasicAnnotation(mappedChr, mappedStart, mappedEnd));
	}
	
	/**
	 * @param fragmentId Fragment ID
	 * @param barcodeSignature Barcodes for the fragment
	 * @param mappedLocation Mapped location of the fragment
	 */
	public BarcodedFragmentImpl(String fragmentId, BarcodeSequence barcodeSignature, Annotation mappedLocation) {
		id = StringParser.firstField(fragmentId);
		setBarcodes(barcodeSignature);
		location = mappedLocation;
		fragmentGroup = new NamedBarcodedFragmentGroup(barcodes);
	}
	
	/**
	 * Get number of barcodes
	 * @return Number of barcodes
	 */
	public int getNumBarcodes() {
		return barcodes.getNumBarcodes();
	}
	
	private void setBarcodes(BarcodeSequence bs) {
		barcodes = bs;
		barcodeString = barcodes.toString();
	}
	
	/**
	 * Instantiate from a SAM record by reading location and attributes
	 * @param samRecord SAM record
	 */
	public BarcodedFragmentImpl(SAMRecord samRecord) {
		
		String fragmentId = StringParser.firstField(samRecord.getReadName());
		read1sequence = null;
		read2sequence = null;
		String seq = samRecord.getReadString();
		if(samRecord.getReadPairedFlag()) {
			if(samRecord.getFirstOfPairFlag()) {
				read1sequence = seq;
			} else if(samRecord.getSecondOfPairFlag()) {
				read2sequence = seq;
			} 
		} else {
			unpairedSequence = seq;
		}
		String barcodeString = samRecord.getStringAttribute(BarcodedBamWriter.BARCODES_SAM_TAG);
		BarcodeSequence barcodeSignature = BarcodeSequence.fromSamAttributeString(barcodeString);
		
		id = fragmentId;
		setBarcodes(barcodeSignature);
		location = new BasicAnnotation(samRecord);
		fragmentGroup = NamedBarcodedFragmentGroup.fromSAMRecord(samRecord);
		
	}
	
	/**
	 * @param fragmentId Fragment ID
	 * @param read1seq Read1 sequence
	 * @param read2seq Read2 sequence
	 * @param barcodeSignature Barcodes for the fragment
	 */
	public BarcodedFragmentImpl(String fragmentId, String read1seq, String read2seq, BarcodeSequence barcodeSignature) {
		this(fragmentId, read1seq, read2seq, barcodeSignature, null);
	}

	/**
	 * @param fragmentId Fragment ID
	 * @param read1seq Read1 sequence
	 * @param read2seq Read2 sequence
	 * @param barcodeSignature Barcodes for the fragment
	 * @param mappedLocation Mapped location of the fragment
	 */
	public BarcodedFragmentImpl(String fragmentId, String read1seq, String read2seq, BarcodeSequence barcodeSignature, Annotation mappedLocation) {
		id = StringParser.firstField(fragmentId);
		read1sequence = read1seq;
		read2sequence = read2seq;
		setBarcodes(barcodeSignature);
		location = mappedLocation;
		fragmentGroup = new NamedBarcodedFragmentGroup(barcodes);
	}

	/**
	 * @param fragmentId Fragment ID
	 * @param read1seq Read1 sequence
	 * @param read2seq Read2 sequence
	 * @param layoutRead1 Read1 layout or null if not specified
	 * @param layoutRead2 Read2 layout or null if not specified
	 * @param maxMismatchesBarcode Max number of mismatches in each barcode when matching to reads
	 */
	public BarcodedFragmentImpl(String fragmentId, String read1seq, String read2seq, ReadLayout layoutRead1, ReadLayout layoutRead2, int maxMismatchesBarcode) {
		id = StringParser.firstField(fragmentId);
		read1sequence = read1seq;
		read2sequence = read2seq;
		read1layout = layoutRead1;
		read2layout = layoutRead2;
		barcodeMaxMismatches = maxMismatchesBarcode;
		fragmentGroup = new NamedBarcodedFragmentGroup(barcodes);
	}
	
	private BarcodedFragmentImpl() {}

	public BarcodeSequence getBarcodes() {
		if(barcodes == null) {
			findBarcodes();
		}
		return barcodes;
	}
	
	public void findBarcodes() {
		barcodes = new BarcodeSequence();
		if(read1layout != null && read1sequence != null) {
			List<List<ReadSequenceElement>> read1elements = read1layout.getMatchedElements(read1sequence);
			if(read1elements != null) {
				for(int i = 0; i < read1elements.size(); i++) {
					ReadSequenceElement parentElement = read1layout.getElements().get(i);
					if(parentElement.getClass().equals(Barcode.class) || parentElement.getClass().equals(BarcodeSet.class)) {
						for(ReadSequenceElement elt : read1elements.get(i)) {
							barcodes.appendBarcode((Barcode)elt);
						}
						continue;
					}
				}
			}
			read1elements = null;
		}
		if(read2layout != null && read2sequence != null) {
			List<List<ReadSequenceElement>> read2elements = read2layout.getMatchedElements(read2sequence);
			if(read2elements != null) {
				for(int i = 0; i < read2elements.size(); i++) {
					ReadSequenceElement parentElement = read2layout.getElements().get(i);
					if(parentElement.getClass().equals(Barcode.class) || parentElement.getClass().equals(BarcodeSet.class)) {
						for(ReadSequenceElement elt : read2elements.get(i)) {
							barcodes.appendBarcode((Barcode)elt);
						}
						continue;
					}
				}
			}
			read2elements = null;
		}
		setBarcodes(barcodes);
	}
	
	public String getId() {
		return id;
	}
	
	public String getUnpairedSequence() {
		return unpairedSequence;
	}
	
	public String getRead1Sequence() {
		return read1sequence;
	}
	
	public String getRead2Sequence() {
		return read2sequence;
	}
	
	public ReadLayout getRead1Layout() {
		return read1layout;
	}
	
	public ReadLayout getRead2Layout() {
		return read2layout;
	}
	
	public Annotation getMappedLocation() {
		return location;
	}
	
	/**
	 * Set mapped location of fragment
	 * @param mappedLocation Mapped location
	 */
	public void setMappedLocation(Annotation mappedLocation) {
		location = mappedLocation;
	}
	
	public int compareTo(BarcodedFragment other) {
		if(location != null && other.getMappedLocation() != null) {
			int l = location.compareTo(other.getMappedLocation());
			if(l != 0) return l;
		}
		return id.compareTo(other.getId());
	}

	@Override
	public FragmentGroup getFragmentGroup() {
		return fragmentGroup;
	}

	@Override
	public void addFragmentWithSameBarcodes(BarcodedFragment fragment) {
		fragmentGroup.addFragment(fragment);
	}

	@Override
	public String getFullInfoString() {
		return location.getFullInfoString();
	}
	
	/**
	 * Get data accessor for Berkeley DB for this entity type
	 * @param environmentHome Database environment home
	 * @param storeName Database entity store name
	 * @param readOnly Database is read only
	 * @return Data accessor
	 */
	public static DataAccessor getDataAccessor(String environmentHome, String storeName, boolean readOnly) {
		logger.info("");
		logger.info("Getting data accessor for database environment " + environmentHome + " and entity store " + storeName + ".");
		BarcodedFragmentImpl b = new BarcodedFragmentImpl();
		return b.new DataAccessor(environmentHome, storeName, readOnly);
	}
	
	/*
	 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/simpleda.html
	 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persist_index.html
	 */
	/**
	 * Data accessor which provides access to database
	 * Opens a database environment
	 * Must call close() when finished with this object
	 * @author prussell
	 *
	 */
	public class DataAccessor {
		
		PrimaryIndex<String, BarcodedFragmentImpl> primaryIndex;
		SecondaryIndex<String, String, BarcodedFragmentImpl> secondaryIndexBarcodeString;
		private DatabaseEnvironment environment;
		private DatabaseStore entityStore;
		
		/**
		 * Data accessor object. MUST BE CLOSED WHEN DONE.
		 * @param environmentHome Database environment home directory
		 * @param storeName Name of database store
		 * @param readOnly Whether environment should be read only
		 */
		public DataAccessor(String environmentHome, String storeName, boolean readOnly) {
		
			/*
			 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persist_first.html
			 */
			environment = new DatabaseEnvironment();
			try {
				environment.setup(new File(environmentHome), readOnly);
			} catch(DatabaseException e) {
				logger.error("Problem with environment setup: " + e.toString());
				System.exit(-1);
			}
			
			entityStore = new DatabaseStore();
			try {
				entityStore.setup(environment.getEnvironment(), storeName, readOnly);
			} catch(DatabaseException e) {
				logger.error("Problem with store setup: " + e.toString());
				System.exit(-1);
			}
			
			primaryIndex = entityStore.getStore().getPrimaryIndex(String.class, BarcodedFragmentImpl.class);
			secondaryIndexBarcodeString = entityStore.getStore().getSecondaryIndex(primaryIndex, String.class, "barcodeString");

			// Attach shutdown hook to close
			/*Runtime.getRuntime().addShutdownHook(new Thread() {
				@Override
				public void run() {
					close();
				}
			});*/
			
		}
		
		/*
		 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/simpleput.html
		 */
		public void put(Collection<BarcodedFragmentImpl> fragments) {
			for(BarcodedFragmentImpl fragment : fragments) {
				primaryIndex.put(fragment);
			}
		}
		
		public void put(BarcodedFragmentImpl fragment) {
			primaryIndex.put(fragment);
		}

		
		/**
		 * Close the store and the environment
		 */
		public void close() {
			logger.info("");
			logger.info("Closing database environment and entity store.");
			entityStore.close();
			environment.close();
		}
		
		/**
		 * Get a cursor over all fragments
		 * Must close when finished
		 * Can either use an iterator over the cursor or treat like a collection
		 * @return Cursor over all fragments
		 */
		public EntityCursor<BarcodedFragmentImpl> getAllFragments() {
			return primaryIndex.entities();
		}
		
		
		/**
		 * Get an iterator over all fragments sharing barcodes with fragments mapping to the window
		 * @param barcodedBam Barcoded bam file
		 * @param chr Window chromsome
		 * @param start Window start
		 * @param end Window end
		 * @param contained Fully contained reads only
		 * @return Iterator over all fragments sharing barcodes with fragments mapping to the window, or null if the window contains no mappings
		 * @throws Exception
		 */
		public JoinedEntityCursor<BarcodedFragmentImpl> getAllFragmentsWithBarcodesMatchingFragmentInChr(String barcodedBam, String chr) throws Exception {
			return getAllFragmentsWithBarcodesMatchingFragmentInWindow(barcodedBam, chr, 0, Integer.MAX_VALUE, true);
		}
		

		
		/**
		 * Get an iterator over all fragments sharing barcodes with fragments mapping to the window
		 * @param barcodedBam Barcoded bam file
		 * @param chr Window chromsome
		 * @param start Window start
		 * @param end Window end
		 * @param contained Fully contained reads only
		 * @return Iterator over all fragments sharing barcodes with fragments mapping to the window, or null if the window contains no mappings
		 * @throws Exception
		 */
		public JoinedEntityCursor<BarcodedFragmentImpl> getAllFragmentsWithBarcodesMatchingFragmentInWindow(String barcodedBam, String chr, int start, int end, boolean contained) throws Exception {
			SAMFileReader reader = new SAMFileReader(new File(barcodedBam));
			SAMRecordIterator iter = reader.query(chr, start, end, contained);
			List<EntityCursor<BarcodedFragmentImpl>> cursors = new ArrayList<EntityCursor<BarcodedFragmentImpl>>();
			if(!iter.hasNext()) {
				reader.close();
				return null;
			}
			while(iter.hasNext()) {
				try {
					SAMRecord record = iter.next();
					String b = record.getStringAttribute(BarcodedBamWriter.BARCODES_SAM_TAG);
					cursors.add(getAllFragmentsWithBarcodes(b));
				} catch (Exception e) {
					for(EntityCursor<BarcodedFragmentImpl> c : cursors) {
						c.close();
					}
					reader.close();
					throw e;
				}
			}
			reader.close();
			return new JoinedEntityCursor<BarcodedFragmentImpl>(cursors);
		}
		
		/**
		 * Get an iterator over all fragments sharing barcodes with fragments mapping to the region
		 * @param alignmentModel Alignment model loaded with barcoded bam file
		 * @param region Region to query
		 * @param contained Fully contained reads only
		 * @return Iterator over all fragments sharing barcodes with fragments mapping to the region, or null if the region contains no mappings
		 * @throws Exception
		 */
		public JoinedEntityCursor<BarcodedFragmentImpl> getAllFragmentsWithBarcodesMatchingFragmentInWindow(AlignmentModel alignmentModel, Annotation region, boolean contained) {
			CloseableIterator<Alignment> iter = alignmentModel.getOverlappingReads(region, contained);
			List<EntityCursor<BarcodedFragmentImpl>> cursors = new ArrayList<EntityCursor<BarcodedFragmentImpl>>();
			if(!iter.hasNext()) {
				iter.close();
				return null;
			}
			while(iter.hasNext()) {
				try {
					Alignment align = iter.next();
					SAMRecord record = align.toSAMRecord();
					String b = record.getStringAttribute(BarcodedBamWriter.BARCODES_SAM_TAG);
					cursors.add(getAllFragmentsWithBarcodes(b));
				} catch (Exception e) {
					for(EntityCursor<BarcodedFragmentImpl> c : cursors) {
						c.close();
					}
					e.printStackTrace();
				}
			}
			iter.close();
			return new JoinedEntityCursor<BarcodedFragmentImpl>(cursors);
		}
		
		/**
		 * Get a cursor over all fragments with a particular barcode string
		 * Must close when finished
		 * Can either use an iterator over the cursor or treat like a collection
		 * @param barcodeString Barcode string
		 * @return Cursor over all fragments with this barcode string
		 */
		public EntityCursor<BarcodedFragmentImpl> getAllFragmentsWithBarcodes(String barcodeString) {
			return secondaryIndexBarcodeString.subIndex(barcodeString).entities();
		}
		
		/**
		 * Get a fragment by its ID
		 * @param id ID
		 * @return The fragment with this ID
		 */
		public BarcodedFragmentImpl getFragmentByID(String id) {
			return primaryIndex.get(id);
		}
		
	}
	
}

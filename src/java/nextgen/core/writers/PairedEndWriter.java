package nextgen.core.writers;

import java.io.File;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.BAMIndex;
import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMTag;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AlignmentPair;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;

import org.apache.commons.io.output.NullOutputStream;
import org.apache.log4j.Logger;

import broad.core.datastructures.Pair;

public class PairedEndWriter {


	public static String PAIRED_END_EXTENSION = ".PairedEnd.bam";
	public static String getDefaultFile(String file) {
		return file + PAIRED_END_EXTENSION;
	}
	
	/**
	 * Writes a file to represent a BAM file where the paired end reads are concatenated
	 * We will represent the basic alignment start and end for the whole fragment
	 * The attributes will have information for the component parts
	 */
	static public final String mateSequenceFlag="ms";
	static public final String mateCigarFlag="mc";
	//static final String mateEndFlag="me";
	static public final String readStartFlag="rs";
	//static final String readEndFlag="re";
	static public final String readCigarFlag="rc";
	static public final String mateLineFlag="mateLine";
	static Logger logger = Logger.getLogger(PairedEndWriter.class.getName());
		
	private final String output;
	private BAMFileWriter writer;
	private SAMFileReader reader;
	private SAMFileHeader header;
	private BAMRecordCodec testCodec;
	private int maxAllowableInsert=5000000;
	
	
	/**
	 * @param bamFile Input SAM or BAM file to extract header and/or reads from.
	 */
	public PairedEndWriter(File bamFile) {
		this(bamFile, getDefaultFile(bamFile.getAbsolutePath()));
	}
			
	/**
	 * @param bamFile Input SAM or BAM file to extract header and/or reads from.
	 * @param output Output path
	 */
	public PairedEndWriter(File bamFile, String output) {
		this.output=output;
		
		reader = new SAMFileReader(bamFile);
		header = reader.getFileHeader();
		
		//set header to pairedEnd
		header.setAttribute(mateLineFlag, "mergedPairedEndFormat");
		
		//We are going to write a BAM File directly
		File outFile = new File(this.output);
		//if (outFile.exists()) outFile.delete();
		writer=new BAMFileWriter(outFile);
		writer.setSortOrder(SortOrder.coordinate, false);
		writer.setHeader(header);
		
		testCodec = new BAMRecordCodec(header);
		testCodec.setOutputStream(new NullOutputStream());
	}
	
	public void setMaxAllowableInsert(int x) {
		maxAllowableInsert = x;
	}
	
	/**
	 * Convert the bamFile provided in the constructor to paired end format.
	 * If no transcription read is supplied, set to unstranded
	 */

	public void convertInputToPairedEnd(){
		this.convertInputToPairedEnd(TranscriptionRead.UNSTRANDED);
	}

	/**
	 * Convert the bamFile provided in the constructor to paired end format.
	 * FOR STRANDED DATA
	 */
	public void convertInputToPairedEnd(TranscriptionRead txnRead) {
		SAMRecordIterator iter = reader.iterator();		
		Map<String, AlignmentPair> tempCollection=new TreeMap<String, AlignmentPair>();
		int numRead = 0;
		while(iter.hasNext()) {
			SAMRecord record=iter.next();
			String name=record.getReadName();
			//If the read is unmapped, skip
			if(record.getReadUnmappedFlag()) continue;
			//If the read is not paired or the mate is unmapped, write it as it is
			if(!record.getReadPairedFlag() || record.getMateUnmappedFlag()){
				//mate unmapped so just write it
				if(record.getReadPairedFlag()) {record.setMateUnmappedFlag(true);} //revised for single end @zhuxp
				//If first read is in direction of trasncription
				if(txnRead.equals(TranscriptionRead.FIRST_OF_PAIR)){
					//This is the first read
					if(record.getFirstOfPairFlag()){
						//Orientation of fragment is same as that of read
					}
					//This is the other mate
					//Reverse its orientation
					else{
						record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
					}
				}
				//Second read is the transcription read
				else if(txnRead.equals(TranscriptionRead.SECOND_OF_PAIR)){
					//This is the first read
					//Reverse orientation
					if(record.getFirstOfPairFlag()){
						record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
					}
					//This is the other mate
					else{
						//NOTHING
					}
				}//UNSTRANDED
				else{
					//NOTHING
				}
				try {
					writer.addAlignment(record);
				} catch (Exception e) {
					logger.error(e.getMessage());
					logger.error("Skipping read " + name );
					continue;
				}
			}
			// read is paired && mate is mapped	
			else{
					
				//create or get the existing pair
				AlignmentPair pair = tempCollection.containsKey(name) ? pair=tempCollection.get(name) :  new AlignmentPair();
				
				//add to pair
				pair.add(record);
					
				//add to Collection
				tempCollection.put(name, pair);

				//If so
				if(pair.isComplete()){
					//Remove from collection
					tempCollection.remove(name);
					//Make paired line for each combo
					Collection<SAMRecord> fragmentRecords = pair.makePairs();
					//write to output
					writeAll(fragmentRecords);
				}
			}		
			numRead++;
			if(numRead % 1000000 == 0) {
				logger.info("Processed " + numRead + " reads, free mem: " + Runtime.getRuntime().freeMemory() + " tempCollection size : " + tempCollection.size() );
			}
		}
		
		//Write remainder
		writeRemainder(tempCollection);
		
		close();
	}
	
	
	private void writeRemainder(Map<String, AlignmentPair> tempCollection) {
		logger.warn("WARNING Remainder: "+tempCollection.size()+" writing as single end reads");
		for(String name: tempCollection.keySet()){
			Pair<Collection<SAMRecord>> pair=tempCollection.get(name);
			System.out.println(name);
			Collection<SAMRecord> records;
			
			if(pair.hasValue1() && pair.hasValue2()){
				//throw new IllegalStateException("There are samples in both pairs that are unaccounted for: "+name);
				logger.error("There are samples in both pairs that are unaccounted for: "+name);

			}
			else{
				if(pair.hasValue1()){records=pair.getValue1();}
				else{records=pair.getValue2();}
				
				for(SAMRecord record: records) {
					// If mate is unpaired, fix SAMRecord settings accordingly
		            record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
		            record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
		            record.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
		            record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
		            record.setMateNegativeStrandFlag(!record.getReadNegativeStrandFlag());
		            record.setMateUnmappedFlag(true);
		            record.setAttribute(SAMTag.MQ.name(), null);
		            record.setInferredInsertSize(0);
				}
			}
		}
	}
	
	


	private void writeAll(Collection<SAMRecord> fragmentRecords) {
		int numHits=fragmentRecords.size();
				
		for(SAMRecord fragment: fragmentRecords){
			fragment.setAttribute("NH", numHits);
			fragment.setMateUnmappedFlag(false);
			addRecord(fragment);
		} 
	}


	public String getOutputFile(){return this.output;}
	
	

	/**
	 * Write a single Alignment
	 * @param alignment
	 */
	public void addAlignment(final Alignment alignment) {
		addRecord(alignment.toSAMRecord());
	}
	
	
	public void addRecord(final SAMRecord record) {
		// Since the BAMWriter stores records in cache and does not encode until the
		// buffer is full, we have to use a second Codec to try to convert each read
		// in order to figure out which one is failing
		boolean encoded = true;
		
		//if distance between pairs is greater than the max allowable we will skip it
		if (record.getInferredInsertSize() > maxAllowableInsert) {
			logger.warn("Skipping read " + record.toString() + " because insert size is greater than " + maxAllowableInsert);
			encoded = false;
		} else {
			try {
				testCodec.encode(record);
			} catch (RuntimeException e) {
				encoded = false;
				logger.error(e.getMessage());
				if (e.getMessage().indexOf("operator maps off end") >= 0) {
					logger.error("Known issue: skipping read " + record.toString() );
					logger.error("(This can happen when reads map greater than ~20Mb away from each other)");
					// TODO I think this is a bug with samtools
				} 
				else {
					logger.error("Failing on record: " + record);
					throw e;
				}
			}
		}

		if (encoded) writer.addAlignment(record);
	}
	
	
	public void close() {
		writer.close();
		
		//Now build a BAM index
		File transcriptomeBamIdxFile = new File( this.output + BAMIndex.BAMIndexSuffix);
		if(transcriptomeBamIdxFile.exists()) { transcriptomeBamIdxFile.delete();}
		SAMFileReader reader2 = new SAMFileReader(new File(this.output));
		BuildBamIndex.createIndex(reader2,transcriptomeBamIdxFile);
		reader2.close();
	}
	
}

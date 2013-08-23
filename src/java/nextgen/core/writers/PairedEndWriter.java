package nextgen.core.writers;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
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
import broad.core.math.Statistics;

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
	private int maxAllowableInsert=500000;
	
	
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

	public void convertInputToPairedEnd(TranscriptionRead txnRead) {
		SAMRecordIterator iter = reader.iterator();		
		Map<String, AlignmentPair> tempCollection=new TreeMap<String, AlignmentPair>();
		int numRead = 0;
		int single = 0;
		int paired = 0;
		String prevChr=null;
		int lastDistanceChecked=0;
		int lastPurgeAt = 0;
		long lastTime = System.currentTimeMillis();
		long writeTime = 0;
		long writeAllTime = 0;
		long mapFinding = 0;
		long pairMaking = 0;
		long iterTime  = 0;
		long removeTime = 0;
		List<Double>  pairMakingTimes = new ArrayList<Double>(); 
		Collection<String> cc = new ArrayList<String>();
		while(iter.hasNext()) {
			long  start = System.nanoTime();
			SAMRecord record=iter.next();
			iterTime = System.nanoTime() - start;
			
			String name=record.getReadName();
			
			if(prevChr==null){
				prevChr=record.getReferenceName();
			}
			//If the read is unmapped, skip
			if(record.getReadUnmappedFlag()) continue;
			
			//If the read is not paired or the mate is unmapped, write it as it is
			if(!record.getReadPairedFlag() || record.getMateUnmappedFlag() || !passesDistanceChecks(record)){ //MG: if we check passes distance checks later, we never write them.
				//mate unmapped so just write it
				if(record.getReadPairedFlag()) {record.setMateUnmappedFlag(true);} //revised for single end @zhuxp
				
				//MG: The implementation of this was so verbose it was hard to understand what was being done. 8/6/13
				//If first read is in direction of transcription change orientation of second read
				if(txnRead.equals(TranscriptionRead.FIRST_OF_PAIR) && !record.getFirstOfPairFlag()){
					record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
				}
				//Second read is the transcription read change the orientation of the first read
				else if(txnRead.equals(TranscriptionRead.SECOND_OF_PAIR) && record.getFirstOfPairFlag()){
						record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
				}//UNSTRANDED: DO nothing
				
				//writer.addAlignment(record); MG: We need to be consistent how we write alignments
				start = System.nanoTime();
				addRecord(record);
				writeTime += System.nanoTime() - start;
				single++;
			}
			// read is paired && mate is mapped	
			else{
				//create or get the existing pair
				
				start = System.nanoTime();
				AlignmentPair pair = tempCollection.get(name);
				if(pair == null) {
					pair = new AlignmentPair();
					tempCollection.put(name, pair);	
				}
				
				pair.add(record);
				mapFinding += System.nanoTime() - start;
				//MG: moved passesDistanceCheck to previous clause so here we only consider hopeful pairs
								

				if(pair.isComplete()){
					paired++;
					//Remove from collection
					start = System.nanoTime();
					tempCollection.remove(name);
					removeTime += System.nanoTime() - start;
					//Make paired line for each combo
					
					start = System.nanoTime();
					Collection<SAMRecord> fragmentRecords = pair.makePairs();
					long pmTime = System.nanoTime() - start;
					pairMaking += pmTime;
					pairMakingTimes.add((double)pmTime);
					//write to output
					long t =  Math.round(pmTime/(double)1000000);
					if(t > 50) {
						logger.debug("Making pairs for " + name + " took " + t + " and has " + fragmentRecords.size() + " frags ");
					}
					
					start = System.nanoTime();
					writeAll(fragmentRecords); //Buffering now.
					writeAllTime += System.nanoTime() - start;
				}
			}	
			if(!record.getReferenceName().equals(prevChr)){
				
				//PURGE TEMP COLLECTION
				writeRemainder(tempCollection);
				tempCollection=new TreeMap<String, AlignmentPair>();
				prevChr = record.getReferenceName();
				lastDistanceChecked=0;
			}
			//If the size of the tempCollection is more than 100,000 and the last distance checked was more than 100kB ago
/*			if(tempCollection.size()>100000 && (record.getAlignmentStart()-lastDistanceChecked)>100000){
				lastDistanceChecked = record.getAlignmentStart();
				//Purge the single alignments that are too far away
				tempCollection = purgeByDistance(tempCollection,record.getAlignmentEnd());
				
			}*/
			numRead++;
			if(numRead % 1000000 == 0) {
				long t = System.currentTimeMillis();
				logger.info("Processed " + numRead + " reads, this batch took: " + (t - lastTime) + " free mem: " + Runtime.getRuntime().freeMemory() + " tempCollection size : " + tempCollection.size()+" on "+record.getReferenceName()+" Single alignments : "+single+ " Paired alignments "+paired);
				logger.info(String.format("Times: it: %d\tws: %d\twp: %d\tpm: %d\tmf: %d\tmr: %d", Math.round(iterTime/(double)1000000), Math.round(writeTime/(double)1000000),  Math.round(writeAllTime/(double)1000000),   Math.round(pairMaking/(double)1000000),   Math.round(mapFinding/(double)1000000),  Math.round(removeTime/(double)1000000)));
				logger.info(String.format("Pairmaling stats max:%f\t mean:%f\tmin: %f\tmedian: %f " , Statistics.max(pairMakingTimes) , Statistics.mean(pairMakingTimes), Statistics.min(pairMakingTimes), Statistics.median(pairMakingTimes)));
				iterTime = 0; writeTime = 0; writeAllTime = 0; pairMaking = 0; mapFinding = 0; removeTime = 0;
				lastTime = t;
			}
		}
		
		//Write remainder
		writeRemainder(tempCollection);
		
		close();
	}
	/**
	 * Convert the bamFile provided in the constructor to paired end format.
	 */
	public void convertInputToPairedEndOld(TranscriptionRead txnRead) {
		SAMRecordIterator iter = reader.iterator();		
		Map<String, AlignmentPair> tempCollection=new LinkedHashMap<String, AlignmentPair>();
		int numRead = 0;
		int single = 0;
		int paired = 0;
		int temp = 0;
		String prevChr=null;
		int lastDistanceChecked=0;
		int lastPurgeAt = 0;
		
		Collection<String> cc = new ArrayList<String>();
		while(iter.hasNext()) {
			SAMRecord record=iter.next();
			String name=record.getReadName();
			
			if(prevChr==null){
				prevChr=record.getReferenceName();
			}
			//If the read is unmapped, skip
			if(record.getReadUnmappedFlag()) continue;
			
			//If the read is not paired or the mate is unmapped, write it as it is
			if(!record.getReadPairedFlag() || record.getMateUnmappedFlag()){
				//mate unmapped so just write it
				if(record.getReadPairedFlag()) {record.setMateUnmappedFlag(true);} //revised for single end @zhuxp
				
				//MG: The implementation of this was so verbose it was hard to understand what was being done. 8/6/13
				//If first read is in direction of transcription change orientation of second read
				if(txnRead.equals(TranscriptionRead.FIRST_OF_PAIR) && !record.getFirstOfPairFlag()){
					record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
				}
				//Second read is the transcription read change the orientation of the first read
				else if(txnRead.equals(TranscriptionRead.SECOND_OF_PAIR) && record.getFirstOfPairFlag()){
						record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
				}//UNSTRANDED: DO nothing
				
				try {
					writer.addAlignment(record);
					single++;
				} finally {
					//MG: I removed the catch, I think it is the wrong thing to do to keep going after a read could not be written to file.
					//logger.error("Error writing read to file " + record.getSAMString() );
				}
			}
			// read is paired && mate is mapped	
			else{
				//create or get the existing pair
				AlignmentPair pair = tempCollection.containsKey(name) ? pair=tempCollection.get(name) :  new AlignmentPair();
				
				//add to pair
				pair.add(record);
					
				//add to Collection
				if(passesDistanceChecks(record)){
					
					if(!tempCollection.containsKey(name))
						temp++;
					tempCollection.put(name, pair);					
					//If so
					if(pair.isComplete()){
						paired++;
						//Remove from collection
						tempCollection.remove(name);
						temp--;
						//Make paired line for each combo
						Collection<SAMRecord> fragmentRecords = pair.makePairs();
						//write to output
						writeAll(fragmentRecords);
					}
				}
			}	
			if(!record.getReferenceName().equals(prevChr)){
				
				//PURGE TEMP COLLECTION
				writeRemainder(tempCollection);
				tempCollection=new TreeMap<String, AlignmentPair>();
				prevChr = record.getReferenceName();
				lastDistanceChecked=0;
			}
			//If the size of the tempCollection is more than 100,000 and the last distance checked was more than 100kB ago
			//TODO: Implement purging so that already available proper pairs are written out.
			if((record.getAlignmentStart()-lastDistanceChecked) > maxAllowableInsert){
				lastDistanceChecked = record.getAlignmentStart();
				//Purge the single alignments that are too far away
				tempCollection = purgeByDistance(tempCollection,record.getAlignmentEnd());
				
			}
			numRead++;
			if(numRead % 1000000 == 0) {
				logger.info("Processed " + numRead + " reads, free mem: " + Runtime.getRuntime().freeMemory() + " tempCollection size : " + tempCollection.size()+" on "+record.getReferenceName()+" Single alignments : "+single+ " Paired alignments "+paired+" In temp : "+temp);
			}
		}
		
		//Write remainder
		writeRemainder(tempCollection);
		
		close();
	}
	
	/**
	 * This function checks whether the given record and its mate:
	 * 			- are on the same chromosome
	 * 			- are less than the maximum allowable distance apart
	 * @param record
	 * @return
	 */
	private boolean passesDistanceChecks(SAMRecord record){
	
		if(!record.getReferenceName().equals(record.getMateReferenceName())){
			return false;
		}
		else{
			if (record.getInferredInsertSize() > maxAllowableInsert) {
				return false;
			}
			else{		
				return true;
			}
		}
	}
	
	
	private void writeRemainder(Map<String, AlignmentPair> tempCollection) {
		if(tempCollection.size() > 0) {
			logger.warn("WARNING Remainder: "+tempCollection.size()+" writing as single end reads");
		}
		for(String name: tempCollection.keySet()){
			Pair<Collection<SAMRecord>> pair=tempCollection.get(name);
			Collection<SAMRecord> records;
			
			if(pair.hasValue1() && pair.hasValue2()){
				//throw new IllegalStateException("There are samples in both pairs that are unaccounted for: "+name);
				logger.error("There are samples in both pairs that are unaccounted for: "+name);
				
				Collection<SAMRecord> fragmentRecords = tempCollection.get(name).makePairs();
				//write to output
				writeAll(fragmentRecords);
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
	
	
	/**
	 * Removes all pairs from the collection such that one of the mates is farther than the max insert size allowed
	 * @param tempCollection
	 */
	private Map<String, AlignmentPair> purgeByDistanceOld(Map<String, AlignmentPair> tempCollection,int currentPosition) {
		
		Map<String, AlignmentPair> newCollection = new TreeMap<String, AlignmentPair>();
		for(String key : tempCollection.keySet()){
			AlignmentPair pair=tempCollection.get(key);
			boolean toAdd=true;
			Iterator<SAMRecord> recordIter;
			//For alignments where we are still waiting for the second one
			if(!pair.hasValue1() || !pair.hasValue2()){
				if(pair.hasValue1()){
					recordIter=pair.getValue1().iterator();
				}
				else{recordIter=pair.getValue2().iterator();}
				
				//pair.removeRecordsFartherThan(currentPosition);
				
				while(recordIter.hasNext()) {
					SAMRecord record = recordIter.next();
					//if(we are farther than the current 
					if((currentPosition-record.getAlignmentEnd()) > maxAllowableInsert){
						recordIter.remove();
						if(pair.valuesAreEmpty())
							toAdd=false;
					}
				}
			}
			if(toAdd)
				newCollection.put(key, pair);
		}
		return newCollection;
	}

	
	/**
	 * Removes all pairs from the collection such that one of the mates is farther than the max insert size allowed
	 * TODO: Implement an inteligent way to deal with multiple hits
	 * @param tempCollection
	 */
	private Collection<AlignmentPair> purgeByDistance(Map<String, AlignmentPair> tempCollection,int currentPosition) {
		
		ArrayList <AlignmentPair> toRemove = new ArrayList< AlignmentPair>();
		Iterator<String> pairNameIt = tempCollection.keySet().iterator();
		boolean done = false;
		while(!done && pairNameIt.hasNext()){
			String key = pairNameIt.next();
			AlignmentPair pair=tempCollection.get(key);
			boolean toAdd=true;
			Iterator<SAMRecord> recordIter;
			//For alignments where we are still waiting for the second one
			if(!pair.hasValue1() || !pair.hasValue2()){
				if(pair.hasValue1()){
					recordIter=pair.getValue1().iterator();
				}
				else{recordIter=pair.getValue2().iterator();}
				
				//pair.removeRecordsFartherThan(currentPosition);
				
				while(recordIter.hasNext()) {
					SAMRecord record = recordIter.next();
					//if(we are farther than the current 
					if((currentPosition-record.getAlignmentEnd()) > maxAllowableInsert){
						recordIter.remove();
						if(pair.valuesAreEmpty())
							toAdd=false;
					}
				}
			}
			if(toAdd)
				newCollection.put(key, pair);
		}
		return newCollection;
	}

	private void writeAll(Collection<SAMRecord> fragmentRecords) {
		for(SAMRecord fragment: fragmentRecords){
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
/*		if (record.getInferredInsertSize() > maxAllowableInsert) {
			logger.warn("Skipping read " + record.toString() + " because insert size is greater than " + maxAllowableInsert);
			encoded = false;
		} else {
*/			try {
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
//		}

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
	
	public static void main (String args []) {
		
	}
	
}

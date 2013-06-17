package nextgen.core.writers;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;


import broad.core.datastructures.Pair;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.BAMIndex;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.FragmentAlignment;
import nextgen.core.alignment.SingleEndAlignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.readFilters.FirstNucleotideFilter;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.ProperPairFilter;

/**
 * @author prussell
 *
 */
public class FilteredBamWriter {

	static Logger logger = Logger.getLogger(FilteredBamWriter.class.getName());	
	private Collection<Predicate<Alignment>> readFilters;
	private File bamFile;
	private TranscriptionRead transcriptionRead;
		
	/**
	 * @param inputBam Input bam file to filter
	 * @param txnRead Read in direction of transcription
	 */
	public FilteredBamWriter(String inputBam, TranscriptionRead txnRead){
		readFilters = new ArrayList<Predicate<Alignment>>();
		bamFile = new File(inputBam);
		transcriptionRead = txnRead;
	}

	/**
	 * Add a new read filter
	 * @param filter Read filter
	 */
	public void addReadFilter(Predicate<Alignment> filter) {
		readFilters.add(filter);
	}
	
	/**
	 * Add a new fragment length filter with max length
	 * @param bedAnnotation Bed file of annotations to use when computing fragment sizes
	 * @param maxLen Max fragment length inclusive
	 * @throws IOException
	 */
	public void addFragmentLengthFilter(String bedAnnotation, int maxLen) throws IOException {
		addFragmentLengthFilter(bedAnnotation, 0, maxLen);
	}
	
	/**
	 * Add a new fragment length filter with min and max lengths
	 * @param bedAnnotation Bed file of annotations to use when computing fragment sizes
	 * @param minLen Min fragment length inclusive
	 * @param maxLen Max fragment length inclusive
	 * @throws IOException
	 */
	public void addFragmentLengthFilter(String bedAnnotation, int minLen, int maxLen) throws IOException {
		logger.info("Adding fragment length filter. Coordinate space from " + bedAnnotation + ". Min length = " + minLen + ". Max length = " + maxLen + ".");
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(bedAnnotation));
		TranscriptomeSpace coord = new TranscriptomeSpace(genes);
		addReadFilter(new FragmentLengthFilter(coord, minLen, maxLen));
	}
	
	
	/**
	 * Add a filter for first nucleotide of fragment
	 * @param geneBedFile Bed file of genes
	 * @param genomeFasta Genome fasta file
	 * @param txnRead Transcription read
	 * @param firstNucleotide Nucleotide to require at first fragment position
	 * @param offsetAlongParent Offset along the parent transcript from the beginning position
	 * @throws IOException
	 */
	public void addFirstNucleotideFilter(String geneBedFile, String genomeFasta, TranscriptionRead txnRead, char firstNucleotide, int offsetAlongParent) throws IOException {
		logger.info("Adding first nucleotide filter for nucleotide " + firstNucleotide + "...");
		addReadFilter(new FirstNucleotideFilter(bamFile.getName(), geneBedFile, genomeFasta, txnRead, firstNucleotide, offsetAlongParent));
	}
	
	/**
	 * Add a new genomic span filter
	 * @param maxSpan Max genomic span inclusive
	 */
	public void addGenomicSpanFilter(int maxSpan) {
		logger.info("Adding genomic span filter. Max span = " + maxSpan + ".");
		addReadFilter(new GenomicSpanFilter(maxSpan));
	}
	
	/**
	 * Add a new proper pair filter
	 */
	public void addProperPairFilter() {
		logger.info("Adding proper pair filter.");
		addReadFilter(new ProperPairFilter());
	}
	
	/**
	 * Filter the reads and write filtered file
	 * @param output Output filtered file
	 */
	public void writeFilteredFile(String output) {
		logger.info("Writing to file " + output + "...");
		SAMFileReader reader = new SAMFileReader(bamFile);
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iter = reader.iterator();
		
		BAMFileWriter writer=new BAMFileWriter(new File(output));
		writer.setSortOrder(SortOrder.coordinate, false);
		writer.setHeader(header);
		
		Map<String, Pair<Collection<SAMRecord>>> tempCollection=new TreeMap<String, Pair<Collection<SAMRecord>>>();
		Map<String, Pair<Integer>> numHitsMap=new TreeMap<String, Pair<Integer>>();
		int numRecords = 0;
		int numWritten = 0;
		
		while(iter.hasNext()) {
			SAMRecord record=iter.next();
			logger.debug("RECORD\t" + record.getReadName());
			if(record.getReadUnmappedFlag()) {
				logger.debug("RECORD_UNMAPPED\t" + record.getReadName());
				continue;
			}
			numRecords++;
			if(numRecords % 100000 == 0) {
				logger.info("Got " + numRecords + " records. Wrote " + numWritten + ".");
			}
			String name=record.getReadName();

			
			if(!record.getReadPairedFlag() || record.getMateUnmappedFlag()){
				Alignment align = new SingleEndAlignment(record);
				logger.debug("NOT_PAIRED_OR_MATE_UNMAPPED\t" + align.getReadName());
				// APPLY FILTERS
				if(!isValid(align)) continue;
				logger.debug("WRITING_ALIGNMENT\t" + record.getReadName());
				writer.addAlignment(record);
				numWritten++;
			} 	else {
				logger.debug("BOTH_MATES_MAPPED\t" + name);
				Pair<Integer> numHits=new Pair<Integer>(); //get from NH flag
				if(numHitsMap.containsKey(name)){
					numHits=numHitsMap.get(name);
				}
				numHits=getNumHits(record, numHits);
				numHitsMap.put(name, numHits);
					
				//create or get the existing pair
				Pair<Collection<SAMRecord>> pair=new Pair<Collection<SAMRecord>>();
				if(tempCollection.containsKey(name)){
					pair=tempCollection.get(name);
				}
					
				//add to pair
				add(record, pair);
				logger.debug("Added to pair\t" + record.getReadName());
					
				//add to Collection
				tempCollection.put(name, pair);
					
				//if pair has both first and second complete then make pair and write
				boolean isComplete=isCompletePair(pair, numHits);
					
				//If so
				if(isComplete){
					logger.debug("PAIR_IS_COMPLETE\t" + record.getReadName());
					//Remove from collection
					tempCollection.remove(name);
					//Make paired line for each combo
					// APPLIES FILTERS
					logger.debug("MAKING_PAIRS\t" + record.getReadName());
					Collection<Pair<SAMRecord>> fragmentRecords=makePairs(pair);
					//write to output
					if(!fragmentRecords.isEmpty()) {
						logger.debug("FRAGMENT_HAS_RECORDS\t" + fragmentRecords.size() + " records\t" + record.getReadName());
						numWritten += fragmentRecords.size();
						writeAll(fragmentRecords, writer);
					}
				} else {
					logger.debug("NOT_COMPLETE_PAIR\t" + record.getReadName());
				}
				numHitsMap.put(name, numHits);
			}		
		}
		
		//Write remainder
		writeRemainder(tempCollection, writer);
		
		//Close the writer
		writer.close();
		
		
		logger.info("Done writing filtered bam file.");
		logger.info("Building bam index.");
		
		//Now build a BAM index
		File transcriptomeBamIdxFile = new File( output + BAMIndex.BAMIndexSuffix);
		if(transcriptomeBamIdxFile.exists()) { transcriptomeBamIdxFile.delete();}
		SAMFileReader reader2 = new SAMFileReader(new File(output));
		BuildBamIndex.createIndex(reader2,transcriptomeBamIdxFile);
		reader2.close();
		
		logger.info("Done building bam index.");

	}
	
	private void writeRemainder(Map<String, Pair<Collection<SAMRecord>>> tempCollection, BAMFileWriter writer) {
		System.err.println("WARNING Remainder: "+tempCollection.size()+" writing as single end reads");
		for(String name: tempCollection.keySet()){
			Pair<Collection<SAMRecord>> pair=tempCollection.get(name);
			Collection<SAMRecord> records;
			
			if(pair.hasValue1() && pair.hasValue2()){
				throw new IllegalArgumentException("There are samples in both pairs that are unaccounted for");
			}
			
			if(pair.hasValue1()){
				records=pair.getValue1();}
			else{
				records=pair.getValue2();}
			
			for(SAMRecord record: records){
				Alignment align = new SingleEndAlignment(record);
				// APPLY FILTERS
				if(!isValid(align)) continue;
				record.setMateUnmappedFlag(true);
				writer.addAlignment(record);
			}
		}
	}


	private static void writeAll(Collection<Pair<SAMRecord>> fragmentRecords, BAMFileWriter writer) {
		
		for(Pair<SAMRecord> fragment: fragmentRecords){
			logger.debug("WRITING_ALIGNMENT\t" + fragment.getValue1().getReadName());
			writer.addAlignment(fragment.getValue1());
			logger.debug("WRITING_ALIGNMENT\t" + fragment.getValue2().getReadName());
			writer.addAlignment(fragment.getValue2());
		}
		
	}



	private static Pair<Integer> getNumHits(SAMRecord record, Pair<Integer> numHits) {
		Object nh=record.getAttribute("NH");
		
		if(nh!=null){
			int num=new Integer(nh.toString()).intValue();
			if(record.getFirstOfPairFlag()){
				numHits.setValue1(Integer.valueOf(num));
				logger.debug("FIRST_OF_PAIR\tNUM_HITS\t" + numHits.getValue1());
			}
			else{
				numHits.setValue2(Integer.valueOf(num));
				logger.debug("SECOND_OF_PAIR\tNUM_HITS\t" + numHits.getValue2());
			}
		}
		else{
			//The NH flag is not set so we will default to 1,1
			logger.debug("NH_FLAG_NOT_SET");
			return new Pair<Integer>(Integer.valueOf(1),Integer.valueOf(1));
		}
				
		return numHits;
	}


	private static boolean isCompletePair(Pair<Collection<SAMRecord>> pair, Pair<Integer> numHits) {
		//First check that we have both ends
		if(pair.hasValue1() && pair.hasValue2()){
			logger.debug("PAIR_HAS_BOTH_VALUES");
			//if so, check that the size matches the expected numHits
			int size1=pair.getValue1().size();
			int size2=pair.getValue2().size();
			
			int expectedSize1=numHits.getValue1().intValue();
			int expectedSize2=numHits.getValue2().intValue();
			
			if(size1==expectedSize1 && size2==expectedSize2){
				logger.debug("SIZES_OK_PAIR_IS_COMPLETE");
				return true;
			}
			logger.debug("SIZES_UNEXPECTED");
		} else {
			logger.debug("PAIR_MISSING_ONE_VALUE");
		}
	
		return false;
	}
	
	private boolean isValid(Alignment read){
		for(Predicate<Alignment> filter: readFilters){
			try {
				boolean passes=filter.evaluate(read);
				if(!passes){
					logger.debug("FAILS_FILTER\t" + read.getReadName() + "\t" + filter.getClass().getName());
					return false;
				}
			} catch (NullPointerException e) {
				logger.debug("CAUGHT_EXCEPTION\t" + read.getReadName() + "\t" + filter.getClass().getName());
				return false;
			}
		}
		logger.debug("IS_VALID\t" + read.getReadName());
		return true;
	}

	
	private Collection<Pair<SAMRecord>> makePairs(Pair<Collection<SAMRecord>> pair) {
		Collection<Pair<SAMRecord>> rtrn=new ArrayList<Pair<SAMRecord>>();
		
		Collection<SAMRecord> pair1=pair.getValue1();
		Collection<SAMRecord> pair2=pair.getValue2();
		
		//match up all pair1 and pair2
		for(SAMRecord r1: pair1){
			for(SAMRecord r2: pair2){
				if(isCompatiblePair(r1, r2)){
					logger.debug("PAIR_COMPATIBLE\t" + r1.getReadName() + "\t" + r2.getReadName());
					Pair<SAMRecord> p=new Pair<SAMRecord>(r1, r2);
					Alignment align = new FragmentAlignment(new SingleEndAlignment(r1), new SingleEndAlignment(r2), transcriptionRead);
					logger.debug("FRAGMENT_ALIGNMENT\t" + align.getReadName());
					if(isValid(align)) {
						logger.debug("ADDING\t" + align.getReadName());
						rtrn.add(p);
					} else {
						logger.debug("NOT_VALID\t" + align.getReadName());
					}
				}
				logger.debug("PAIR_NOT_COMPATIBLE\t" + r1.getReadName() + "\t" + r2.getReadName());
			}
		}
		
		return rtrn;
	}
	
	
	private static boolean isCompatiblePair(SAMRecord r1, SAMRecord r2) {
		if(r1.getReferenceName()==r2.getReferenceName()){
			if((r1.getAlignmentStart()==r2.getMateAlignmentStart())&& (r1.getMateAlignmentStart()==r2.getAlignmentStart())){
				if(r1.getMateNegativeStrandFlag()==r2.getReadNegativeStrandFlag() && r1.getReadNegativeStrandFlag()==r2.getMateNegativeStrandFlag()){
					if(r1.getMateReferenceName().equalsIgnoreCase(r2.getReferenceName()) && r2.getMateReferenceName().equalsIgnoreCase(r1.getReferenceName())){
						return true;
					}
				}
			}
		}
		return false;
	}



	private static void add(SAMRecord record, Pair<Collection<SAMRecord>> pair) {
		if(record.getFirstOfPairFlag()){
			Collection<SAMRecord> set=new ArrayList<SAMRecord>();
			if(pair.hasValue1()){set=pair.getValue1();}
			set.add(record);
			pair.setValue1(set);
		}
		else{
			Collection<SAMRecord> set=new ArrayList<SAMRecord>();
			if(pair.hasValue2()){set=pair.getValue2();}
			set.add(record);
			pair.setValue2(set);
		}
	}
	
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Input bam file", true);
		p.addIntArg("-maxg", "Max genomic span", false, -1);
		p.addIntArg("-minf", "Min fragment size", false, -1);
		p.addIntArg("-maxf", "Max fragment size", false, -1);
		p.addStringArg("-fn", "Nucleotide for first fragment position (A, C, G or T)", false, null);
		p.addStringArg("-a", "Bed gene annotation", false);
		p.addStringArg("-g", "Genome fasta file required for -fn", false, null);
		p.addIntArg("-of", "Offset for -fn", false, 0);
		p.addStringArg("-o", "Output bam file", true);
		p.addBooleanArg("-ft", "First read is transcription strand", false, false);
		p.addBooleanArg("-d","Debug logging", false, false);
		p.parse(args);
		String inputBam = p.getStringArg("-b");
		int maxGenomicSpan = p.getIntArg("-maxg");
		int minFragmentSize = p.getIntArg("-minf");
		int maxFragmentSize = p.getIntArg("-maxf");
		String bedFile = p.getStringArg("-a");
		String outFile = p.getStringArg("-o");
		String firstNuc = p.getStringArg("-fn");
		String genomeFasta = p.getStringArg("-g");
		boolean firstReadTranscriptionStrand = p.getBooleanArg("-ft");
		int offset = p.getIntArg("-of");
		boolean debug = p.getBooleanArg("-d");
		
		if(debug) {
			logger.setLevel(Level.DEBUG);
		}
		
		FilteredBamWriter fbw = new FilteredBamWriter(inputBam, firstReadTranscriptionStrand ? TranscriptionRead.FIRST_OF_PAIR : TranscriptionRead.SECOND_OF_PAIR);
		
		if(maxGenomicSpan >= 0) {
			fbw.addGenomicSpanFilter(maxGenomicSpan);
		}
		
		if(maxFragmentSize >= 0) {
			if(bedFile == null) {
				throw new IllegalArgumentException("For fragment size filter must provide bed annotation with option -a.");
			}
			int minSize = minFragmentSize < 0 ? 0 : minFragmentSize;
			fbw.addFragmentLengthFilter(bedFile, minSize, maxFragmentSize);
		}
		
		if(firstNuc != null) {
			if(!firstNuc.equals("A") && !firstNuc.equals("C") && !firstNuc.equals("G") && !firstNuc.equals("T")) {
				throw new IllegalArgumentException(firstNuc + " is not a valid first nucleotide. Must be A, C, G or T.");
			}
			if(bedFile == null) {
				throw new IllegalArgumentException("For first nucleotide filter must provide bed annotation with option -a.");
			}
			if(genomeFasta == null) {
				throw new IllegalArgumentException("For first nucleotide filter must provide genome fasta file with option -g.");
			}
			
			
			fbw.addFirstNucleotideFilter(bedFile, genomeFasta, fbw.transcriptionRead, firstNuc.charAt(0), offset); 
			
			
		}
		
		fbw.writeFilteredFile(outFile);
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}

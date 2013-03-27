package nextgen.core.writers;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.collections15.Predicate;
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
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
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
	private boolean skipFirstReads;
	private boolean skipSecondReads;
		
	/**
	 * @param inputBam Input bam file to filter
	 */
	public FilteredBamWriter(String inputBam){
		this.readFilters = new ArrayList<Predicate<Alignment>>();
		this.bamFile = new File(inputBam);
		this.skipFirstReads = false;
		this.skipSecondReads = false;
	}

	/**
	 * Add a new read filter
	 * @param filter Read filter
	 */
	public void addReadFilter(Predicate<Alignment> filter) {
		this.readFilters.add(filter);
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
	 * Don't write any second reads
	 */
	public void skipSecondReads() {
		this.skipSecondReads = true;
	}
	
	/**
	 * Don't write any first reads
	 */
	public void skipFirstReads() {
		this.skipFirstReads = true;
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
		SAMFileReader reader = new SAMFileReader(this.bamFile);
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iter = reader.iterator();
		
		BAMFileWriter writer=new BAMFileWriter(new File(output));
		writer.setSortOrder(SortOrder.coordinate, false);
		writer.setHeader(header);
		
		Map<String, Pair<Collection<SAMRecord>>> tempCollection=new TreeMap<String, Pair<Collection<SAMRecord>>>();
		Map<String, Pair<Integer>> numHitsMap=new TreeMap<String, Pair<Integer>>();
		int numRecords = 0;
		
		while(iter.hasNext()) {
			SAMRecord record=iter.next();
			if(record.getReadUnmappedFlag()) {
				continue;
			}
			numRecords++;
			if(numRecords % 100000 == 0) {
				logger.info("Got " + numRecords + " records.");
			}
			String name=record.getReadName();
			
			if(!record.getReadPairedFlag()){
				Alignment align = new SingleEndAlignment(record);
				// APPLY FILTERS
				if(!isValid(align)) continue;
				writer.addAlignment(record);
			} else if (record.getMateUnmappedFlag()) {
				if(this.skipFirstReads && record.getFirstOfPairFlag()) continue;
				if(this.skipSecondReads && record.getSecondOfPairFlag()) continue;
				Alignment align = new SingleEndAlignment(record);
				// APPLY FILTERS
				if(!isValid(align)) continue;
				writer.addAlignment(record);
			}
							
			else{
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
					
				//add to Collection
				tempCollection.put(name, pair);
					
				//if pair has both first and second complete then make pair and write
				boolean isComplete=isCompletePair(pair, numHits);
					
				//If so
				if(isComplete){
					//Remove from collection
					tempCollection.remove(name);
					//Make paired line for each combo
					// APPLIES FILTERS
					Collection<Pair<SAMRecord>> fragmentRecords=makePairs(pair);
					//write to output
					if(!fragmentRecords.isEmpty()) writeAll(fragmentRecords, writer, this.skipFirstReads, this.skipSecondReads);
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
				if(this.skipFirstReads) continue;
				records=pair.getValue1();}
			else{
				if(this.skipSecondReads) continue;
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
		writeAll(fragmentRecords, writer, false, false);
	}
	
	private static void writeAll(Collection<Pair<SAMRecord>> fragmentRecords, BAMFileWriter writer, boolean skipFirstReads, boolean skipSecondReads) {
		
		for(Pair<SAMRecord> fragment: fragmentRecords){
			if(!skipFirstReads) writer.addAlignment(fragment.getValue1());
			if(!skipSecondReads) writer.addAlignment(fragment.getValue2());
		}
		
	}



	private static Pair<Integer> getNumHits(SAMRecord record, Pair<Integer> numHits) {
		Object nh=record.getAttribute("NH");
		
		if(nh!=null){
			int num=new Integer(nh.toString()).intValue();
			if(record.getFirstOfPairFlag()){
				numHits.setValue1(Integer.valueOf(num));
			}
			else{
				numHits.setValue2(Integer.valueOf(num));
			}
		}
		else{
			//The NH flag is not set so we will default to 1,1
			return new Pair<Integer>(Integer.valueOf(1),Integer.valueOf(1));
		}
		
		return numHits;
	}


	private static boolean isCompletePair(Pair<Collection<SAMRecord>> pair, Pair<Integer> numHits) {
		//First check that we have both ends
		if(pair.hasValue1() && pair.hasValue2()){
			//if so, check that the size matches the expected numHits
			int size1=pair.getValue1().size();
			int size2=pair.getValue2().size();
			
			int expectedSize1=numHits.getValue1().intValue();
			int expectedSize2=numHits.getValue2().intValue();
			
			if(size1==expectedSize1 && size2==expectedSize2){
				return true;
			}
		}
		
		return false;
	}
	
	private boolean isValid(Alignment read){
		for(Predicate<Alignment> filter: this.readFilters){
			try {
				boolean passes=filter.evaluate(read);
				if(!passes){return false;}
			} catch (NullPointerException e) {
				return false;
			}
		}
		
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
					Pair<SAMRecord> p=new Pair<SAMRecord>(r1, r2);
					Alignment align = new FragmentAlignment(new SingleEndAlignment(r1), new SingleEndAlignment(r2));
					if(isValid(align)) {
						rtrn.add(p);
					}
				}
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
		p.addIntArg("-maxg", "Max genomic span", false, Integer.valueOf(-1));
		p.addIntArg("-minf", "Min fragment size", false, Integer.valueOf(-1));
		p.addIntArg("-maxf", "Max fragment size", false, Integer.valueOf(-1));
		p.addBooleanArg("-sf", "Skip first reads", false, Boolean.valueOf(false));
		p.addBooleanArg("-ss", "Skip second reads", false, Boolean.valueOf(false));
		p.addStringArg("-a", "Bed annotation for fragment lengths", false);
		p.addStringArg("-o", "Output bam file", true);
		p.parse(args);
		String inputBam = p.getStringArg("-b");
		Integer maxGenomicSpan = p.getIntArg("-maxg");
		Integer minFragmentSize = p.getIntArg("-minf");
		Integer maxFragmentSize = p.getIntArg("-maxf");
		String bedFile = p.getStringArg("-a");
		String outFile = p.getStringArg("-o");
		boolean skipFirstReads = p.getBooleanArg("-sf");
		boolean skipSecondReads = p.getBooleanArg("-ss");
		
		FilteredBamWriter fbw = new FilteredBamWriter(inputBam);
		if(skipFirstReads) fbw.skipFirstReads();
		if(skipSecondReads) fbw.skipSecondReads();
		
		if(maxGenomicSpan.intValue() >= 0) {
			fbw.addGenomicSpanFilter(maxGenomicSpan.intValue());
		}
		
		if(maxFragmentSize.intValue() >= 0) {
			if(bedFile == null) {
				throw new IllegalArgumentException("For fragment size filter must provide bed annotation with option -a.");
			}
			int minSize = minFragmentSize.intValue() < 0 ? 0 : minFragmentSize.intValue();
			fbw.addFragmentLengthFilter(bedFile, minSize, maxFragmentSize.intValue());
		}
		
		fbw.writeFilteredFile(outFile);
		
	}
	
}

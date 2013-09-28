/**
 * 
 */
package nextgen.editing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.utils.CountLogger;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.PcrPrimerDesigner;
import broad.core.primer3.Primer3Configuration;
import broad.core.primer3.PrimerPair;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class DeletionPlasmidUtils {
	
	private static Logger logger = Logger.getLogger(DeletionPlasmidUtils.class.getName());
	
	/**
	 * Get a plasmid created by deleting a region
	 * @param chromosome The full chromosome containing the gene of interest
	 * @param transcript The gene of interest as an annotation on the chromosome
	 * @param regionToDelete The deletion site as an annotation on the chromosome
	 * @return The plasmid beginning immediately 3' of the deletion site, circling around, and ending immediately 5' of the deletion site
	 */
	public static DeletionPlasmid getDeletionPlasmid(Sequence chromosome, Gene transcript, Gene regionToDelete) {
		
		if(!transcript.contains(regionToDelete)) {
			throw new IllegalArgumentException("Transcript must contain region to delete");
		}
		
		Sequence transcriptSequence = chromosome.getSubsequence(transcript);
		transcriptSequence.setId(transcript.getName());
		
		int firstDeletedPosition = transcript.genomicToTranscriptPosition(regionToDelete.transcriptToGenomicPosition(0));
		int lastDeletedPosition = transcript.genomicToTranscriptPosition(regionToDelete.transcriptToGenomicPosition(regionToDelete.size() - 1));
		
		return new DeletionPlasmid(transcriptSequence, firstDeletedPosition, lastDeletedPosition, transcript.getName() + "_deletion_" + regionToDelete.toUCSC());
		
	}
	
	/**
	 * Get genomic position of a position on the plasmid
	 * @param transcript Original transcript
	 * @param deletionRegion Deletion region
	 * @param positionOnPlasmidSequenceAfterDeletion Position within plasmid after deletion
	 * @return Corresponding genomic position
	 */
	public static int getGenomicPosition(Gene transcript, Gene deletionRegion, int positionOnPlasmidSequenceAfterDeletion) {
		
		return transcript.transcriptToGenomicPosition(getPositionOnOriginalTranscript(transcript, deletionRegion, positionOnPlasmidSequenceAfterDeletion));
		
	}
	
	/**
	 * Get corresponding position on original transcript of a position on the plasmid
	 * @param transcript Original transcript
	 * @param deletionRegion Deletion region
	 * @param positionOnPlasmidSequenceAfterDeletion Position within plasmid after deletion
	 * @return Corresponding transcript position
	 */
	public static int getPositionOnOriginalTranscript(Gene transcript, Gene deletionRegion, int positionOnPlasmidSequenceAfterDeletion) {

		int deletionSize = deletionRegion.size();
		int firstDeletedPosition = transcript.genomicToTranscriptPosition(deletionRegion.transcriptToGenomicPosition(0));
		int lastDeletedPosition = transcript.genomicToTranscriptPosition(deletionRegion.transcriptToGenomicPosition(deletionSize - 1));
		int numBasesBeforeDeletion = firstDeletedPosition;
		int numBasesAfterDeletion = transcript.size() - lastDeletedPosition - 1;
		
		if(positionOnPlasmidSequenceAfterDeletion < numBasesAfterDeletion) {
			return numBasesBeforeDeletion + deletionSize + positionOnPlasmidSequenceAfterDeletion;
		}
		
		return positionOnPlasmidSequenceAfterDeletion - numBasesAfterDeletion;
		
	}
	
	private static Gene getSingleGeneFromBedFile(String singleGeneBedFile) throws IOException {
		Map<String, Collection<Gene>> genesByChromosome = BEDFileParser.loadDataByChr(new File(singleGeneBedFile));
		int numGenes = 0;
		for(String chr : genesByChromosome.keySet()) {
			numGenes += genesByChromosome.get(chr).size();
		}
		if(numGenes != 1) {
			throw new IllegalArgumentException("Gene bed file must contain a single gene.");
		}
		Gene gene = genesByChromosome.values().iterator().next().iterator().next();
		return gene;
	}
	
	/**
	 * Get plasmids created by deleting a region
	 * Plasmid sequence begins immediately 3' of the deletion site, circles around, and ends immediately 5' of the deletion site
	 * @param chromosomeFastaFile Fasta file containing the chromosome the gene is on
	 * @param singleGeneBedFile Bed file containing the single gene of interest
	 * @param deletionsBedFile Bed file of deletion sites within the gene
	 * @return The plasmids
	 * @throws IOException
	 */
	public static Collection<DeletionPlasmid> getDeletionPlasmids(String chromosomeFastaFile, String singleGeneBedFile, String deletionsBedFile) throws IOException {
		logger.info("Constructing deletion plasmids from chromosome sequences...");
		logger.info("Loading chromosome sequences...");
		FastaSequenceIO fastaReader = new FastaSequenceIO(chromosomeFastaFile);
		List<Sequence> chromosomes = fastaReader.loadAll();
		Gene gene = getSingleGeneFromBedFile(singleGeneBedFile);
		String chr = gene.getChr();
		Sequence chrSequence = null;
		for(Sequence chromosome : chromosomes) {
			if(chromosome.getId().equals(chr)) {
				chrSequence = chromosome;
				break;
			}
		}
		if(chrSequence == null) {
			throw new IllegalArgumentException("Chromosome fasta file must include the chromosome the gene is on");
		}
		logger.info("Using chromosome sequence " + chrSequence.getId() + ".");
		
		logger.info("Loading deletions...");
		List<DeletionPlasmid> plasmids = new ArrayList<DeletionPlasmid>();
		Map<String, Collection<Gene>> deletions = BEDFileParser.loadDataByChr(new File(deletionsBedFile));
		for(String c : deletions.keySet()) {
			for(Gene deletion : deletions.get(c)) {
				logger.info("Getting plasmid sequence for deletion " + deletion.getName());
				plasmids.add(getDeletionPlasmid(chrSequence, gene, deletion));
			}
		}
		logger.info("Got all deletions.");
		return plasmids;
	}
	
	/**
	 * Get sequences of plasmids created by deleting a region and write to a fasta file.
	 * Output sequence begins immediately 3' of the deletion site, circles around, and ends immediately 5' of the deletion site
	 * @param chromosomeFastaFile Fasta file containing the chromosome the gene is on
	 * @param singleGeneBedFile Bed file containing the single gene of interest
	 * @param deletionsBedFile Bed file of deletion sites within the gene
	 * @param outFastaFile Output fasta file for plasmid sequences
	 * @throws IOException
	 */
	public static void writeDeletionPlasmidSequences(String chromosomeFastaFile, String singleGeneBedFile, String deletionsBedFile, String outFastaFile) throws IOException {
		writeDeletionPlasmidSequences(getDeletionPlasmids(chromosomeFastaFile, singleGeneBedFile, deletionsBedFile), outFastaFile);
	}
	
	/**
	 * Write sequences to a file
	 * Output sequence begins immediately 3' of the deletion site, circles around, and ends immediately 5' of the deletion site
	 * @param plasmids Plasmids to write
	 * @param outFastaFile Output fasta file for plasmid sequences
	 * @throws IOException
	 */
	public static void writeDeletionPlasmidSequences(Collection<DeletionPlasmid> plasmids, String outFastaFile) throws IOException	{
		logger.info("Writing deletions to file " + outFastaFile + "...");
		FastaSequenceIO fastaWriter = new FastaSequenceIO(outFastaFile);
		List<Sequence> plasmidSeqs = new ArrayList<Sequence>();
		for(DeletionPlasmid plasmid : plasmids) {
			plasmidSeqs.add(plasmid.getPlasmidSequence());
		}
		fastaWriter.write(plasmidSeqs);
		logger.info("Done writing fasta file.");

	}
	
	/**
	 * If there are regions we want to delete and we designed PCR primers facing out from the regions
	 * to copy the rest of the construct, more than just the desired region will be deleted. The actual
	 * plasmid generated depends on the locations of the primers. This method takes output from the
	 * primer creation step, as well as the full desired plasmids (with only the specific region deleted)
	 * and returns the actual plasmids that will have been constructed by these primers. The parent
	 * transcript field remains the full original transcript without any deletion.
	 * @param fullDesiredPlasmids The full desired plasmids (with only the specific region deleted)
	 * @param pcrPrimerDesignerOutput PcrPrimerDesigner output against the full desired plasmids
	 * Indicates the actual locations of the primers used to construct the plasmid with deletion
	 * @return The actual plasmids that will have been constructed using the given primers
	 * @throws IOException
	 */
	private static Collection<DeletionPlasmid> getActualPCRProducts(Collection<DeletionPlasmid> fullDesiredPlasmids, String pcrPrimerDesignerOutput) throws IOException {
		logger.info("Getting actual PCR products: plasmids constructed using existing primers...");
		Map<String, DeletionPlasmid> fullPlasmidsByName = new TreeMap<String, DeletionPlasmid>();
		for(DeletionPlasmid p : fullDesiredPlasmids) {
			fullPlasmidsByName.put(p.getId(), p);
		}
		
		Collection<DeletionPlasmid> rtrn = new ArrayList<DeletionPlasmid>();
				
		FileReader r = new FileReader(pcrPrimerDesignerOutput);
		BufferedReader b = new BufferedReader(r);
		String header = b.readLine();
		StringParser p = new StringParser();
		p.parse(header);
		int sequenceColumn = p.getIndexFor(PrimerPair.SEQUENCE_FIELD_NAME);
		int leftPrimerPosColumn = p.getIndexFor(PrimerPair.LEFT_PRIMER_POS_FIELD_NAME);
		int rightPrimerPosColumn = p.getIndexFor(PrimerPair.RIGHT_PRIMER_POS_FIELD_NAME);
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			String sequenceName = p.asString(sequenceColumn);
			if(!fullPlasmidsByName.containsKey(sequenceName)) {
				throw new IllegalArgumentException("Sequence " + sequenceName + " not in fasta file.");
			}
			if(line.contains(PrimerPair.NO_PRIMERS_NAME)) {
				logger.info("Skipping sequence " + sequenceName + " which has no primers.");
				continue;
			}
			DeletionPlasmid fullPlasmid = fullPlasmidsByName.get(sequenceName);
			int fullPlasmidSize = fullPlasmid.getPlasmidSize();
			int origSeqSize = fullPlasmid.getOriginalSequence().getLength();
			int leftPrimerPosOnFullPlasmid = p.asInt(leftPrimerPosColumn);
			int lastDeletedPosOnOriginalSequence = (fullPlasmid.getLastDeletedPosition() + leftPrimerPosOnFullPlasmid) % origSeqSize;
			int rightPrimerPosOnFullPlasmid = p.asInt(rightPrimerPosColumn);
			int firstDeletedPosOnOriginalSequence = (fullPlasmid.getFirstDeletedPosition() - (fullPlasmidSize - rightPrimerPosOnFullPlasmid - 1));
			if(firstDeletedPosOnOriginalSequence < 0) firstDeletedPosOnOriginalSequence += origSeqSize;
			String actualProductID = fullPlasmid.getId();
			DeletionPlasmid actualProduct = new DeletionPlasmid(fullPlasmid.getOriginalSequence(), firstDeletedPosOnOriginalSequence, lastDeletedPosOnOriginalSequence, actualProductID);
			
			logger.info("Finished plasmid " + actualProduct.toString());
			rtrn.add(actualProduct);
		}
		r.close();
		b.close();
		logger.info("Done getting actual products.");
		return rtrn;
	}
	
	private static String ORIG_SEQUENCE_ID_COLUMN = "OriginalTranscriptID";
	private static String PLASMID_ID_COLUMN = "PlasmidID";
	private static String FIRST_DELETED_POSITION_COLUMN = "FirstDeletedPosition";
	private static String LAST_DELETED_POSITION_COLUMN = "LastDeletedPosition";
	private static String ORIG_SEQUENCE_COLUMN = "OriginalTranscriptSequence";
	private static String PLASMID_SEQUENCE_COLUMN = "PlasmidSequence";
	
	/**
	 * @param plasmids Plasmids to write to table
	 * @param outFile Output file
	 * @throws IOException
	 */
	public static void writeToTable(Collection<DeletionPlasmid> plasmids, String outFile) throws IOException {
		logger.info("Writing " + plasmids.size() + " plasmids to table " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		String header = ORIG_SEQUENCE_ID_COLUMN + "\t";
		header += PLASMID_ID_COLUMN + "\t";
		header += FIRST_DELETED_POSITION_COLUMN + "\t";
		header += LAST_DELETED_POSITION_COLUMN + "\t";
		header += ORIG_SEQUENCE_COLUMN + "\t";
		header += PLASMID_SEQUENCE_COLUMN + "\t";
		w.write(header + "\n");
		for(DeletionPlasmid plasmid : plasmids) {
			String line = plasmid.getOriginalSequence().getId() + "\t";
			line += plasmid.getId() + "\t";
			line += plasmid.getFirstDeletedPosition() + "\t";
			line += plasmid.getLastDeletedPosition() + "\t";
			line += plasmid.getOriginalSequence().getSequenceBases() + "\t";
			line += plasmid.getPlasmidSequence().getSequenceBases() + "\t";
			w.write(line + "\n");
		}
		w.close();
		logger.info("Done writing table.");
	}
	
	/**
	 * @param tableFile Table written by this class
	 * @return Plasmids described in table
	 * @throws IOException
	 */
	public static Collection<DeletionPlasmid> readFromTable(String tableFile) throws IOException {
		logger.info("Reading plasmids from file " + tableFile + "...");
		FileReader r = new FileReader(tableFile);
		BufferedReader b = new BufferedReader(r);
		String header = b.readLine();
		StringParser p = new StringParser();
		p.parse(header);
		int origSeqIdCol = p.getIndexFor(ORIG_SEQUENCE_ID_COLUMN);
		int plasmidIdCol = p.getIndexFor(PLASMID_ID_COLUMN);
		int firstDelPosCol = p.getIndexFor(FIRST_DELETED_POSITION_COLUMN);
		int lastDelPosCol = p.getIndexFor(LAST_DELETED_POSITION_COLUMN);
		int origSeqCol = p.getIndexFor(ORIG_SEQUENCE_COLUMN);
		Collection<DeletionPlasmid> rtrn = new ArrayList<DeletionPlasmid>();
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			String origName = p.asString(origSeqIdCol);
			String origSeq = p.asString(origSeqCol);
			Sequence orig = new Sequence(origName);
			orig.setSequenceBases(origSeq);
			DeletionPlasmid plasmid = new DeletionPlasmid(orig, p.asInt(firstDelPosCol), p.asInt(lastDelPosCol), p.asString(plasmidIdCol));
			rtrn.add(plasmid);
		}
		logger.info("Read " + rtrn.size() + " plasmids.");
		return rtrn;
	}
	
	/**
	 * Write bed file of primer pairs produced by PcrPrimerDesigner on the plasmid sequences with deletions
	 * @param pcrPrimerDesignerOutputFile Output file from PcrPrimerDesigner
	 * @param singleGeneBedFile Bed file of single gene of interest
	 * @param deletionsBedFile Bed file of deletions within gene
	 * @param outBedFile Bed file of genomic positions of primer pairs
	 * @throws IOException
	 */
	private static void writePrimerPairBedFile(String pcrPrimerDesignerOutputFile, String singleGeneBedFile, String deletionsBedFile, String outBedFile) throws IOException {
		
		logger.info("Writing primer pair annotation to file " + outBedFile + "...");
		
		logger.info("Reading deletions from file " + deletionsBedFile + "...");
		Map<String, Collection<Gene>> deletionsByChr = BEDFileParser.loadDataByChr(new File(deletionsBedFile));
		Map<String, Gene> deletionsByCoord = new TreeMap<String, Gene>();
		for(String chr : deletionsByChr.keySet()) {
			for(Gene deletion : deletionsByChr.get(chr)) {
				deletionsByCoord.put(deletion.toUCSC(), deletion);
			}
		}
		
		Gene gene = getSingleGeneFromBedFile(singleGeneBedFile);
		
		FileWriter w = new FileWriter(outBedFile);
		
		logger.info("Reading primer pair information from file " + pcrPrimerDesignerOutputFile + "...");
		FileReader r = new FileReader(pcrPrimerDesignerOutputFile);
		BufferedReader b = new BufferedReader(r);
		String header = b.readLine();
		StringParser p = new StringParser();
		p.parse(header);
		int sequenceColumn = p.getIndexFor(PrimerPair.SEQUENCE_FIELD_NAME);
		int leftPrimerColumn = p.getIndexFor(PrimerPair.LEFT_PRIMER_FIELD_NAME);
		int rightPrimerColumn = p.getIndexFor(PrimerPair.RIGHT_PRIMER_FIELD_NAME);
		int leftPrimerPosColumn = p.getIndexFor(PrimerPair.LEFT_PRIMER_POS_FIELD_NAME);
		int rightPrimerPosColumn = p.getIndexFor(PrimerPair.RIGHT_PRIMER_POS_FIELD_NAME);
		while(b.ready()) {
			String line = b.readLine();
			p.parse(line);
			String sequenceName = p.asString(sequenceColumn);
			if(line.contains(PrimerPair.NO_PRIMERS_NAME)) {
				logger.info("Skipping sequence " + sequenceName + " which has no primers.");
				continue;
			}
			Gene deletion = null;
			for(String ucsc : deletionsByCoord.keySet()) {
				if(sequenceName.subSequence(sequenceName.length() - ucsc.length(), sequenceName.length()).equals(ucsc)) {
					deletion = deletionsByCoord.get(ucsc);
					break;
				}
			}
			if(deletion == null) {
				throw new IllegalArgumentException("Sequence referred to in " + sequenceName + " not found in bed file");
			}
			int leftPrimerStartPosOnPlasmid = p.asInt(leftPrimerPosColumn);
			int rightPrimerEndPosOnPlasmid = p.asInt(rightPrimerPosColumn);
			int leftPrimerSize = p.asString(leftPrimerColumn).length();
			int rightPrimerSize = p.asString(rightPrimerColumn).length();
			int leftPrimerEndPosOnPlasmid = leftPrimerStartPosOnPlasmid + leftPrimerSize - 1;
			int rightPrimerStartPosOnPlasmid = rightPrimerEndPosOnPlasmid - rightPrimerSize + 1;
			int leftPrimerStartGenomic = getGenomicPosition(gene, deletion, leftPrimerStartPosOnPlasmid);
			int leftPrimerEndGenomic = getGenomicPosition(gene, deletion, leftPrimerEndPosOnPlasmid);
			int rightPrimerStartGenomic = getGenomicPosition(gene, deletion, rightPrimerStartPosOnPlasmid);
			int rightPrimerEndGenomic = getGenomicPosition(gene, deletion, rightPrimerEndPosOnPlasmid);
			int leftPrimerGenomicStart = Math.min(leftPrimerStartGenomic, leftPrimerEndGenomic);
			int leftPrimerGenomicEnd = gene.transcriptToGenomicPosition(gene.genomicToTranscriptPosition(Math.max(leftPrimerStartGenomic, leftPrimerEndGenomic))) + 1;
			int rightPrimerGenomicStart = Math.min(rightPrimerStartGenomic, rightPrimerEndGenomic);
			int rightPrimerGenomicEnd = gene.transcriptToGenomicPosition(gene.genomicToTranscriptPosition(Math.max(rightPrimerStartGenomic, rightPrimerEndGenomic))) + 1;
			Gene leftPrimer = gene.trimAbsolute(leftPrimerGenomicStart, leftPrimerGenomicEnd);
			Gene rightPrimer = gene.trimAbsolute(rightPrimerGenomicStart, rightPrimerGenomicEnd);
			Collection<Annotation> primerPairExons = new TreeSet<Annotation>();
			primerPairExons.addAll(leftPrimer.getBlocks());
			primerPairExons.addAll(rightPrimer.getBlocks());
			Gene primerPair = new Gene(primerPairExons);
			primerPair.setName(sequenceName);
			w.write(primerPair.toBED() + "\n");
		}
		r.close();
		b.close();
		w.close();
		
		logger.info("Done writing bed file.");
		
	}
		
	private static int DEFAULT_MIN_PCR_PRODUCT_SIZE = 60;
	private static int DEFAULT_MAX_PCR_PRODUCT_SIZE = 1000;
	
	/**
	 * Get a PCR primer flanking the deleted region
	 * @param plasmid The plasmid
	 * @param maxRatioWtProductToDeletion Maximum ratio of total wild type product size (includes the deleted region and flanking regions) to the size of the deleted region itself
	 * @param pathPrimer3core primer3core executable
	 * @return The primer pair flanking the deleted region with the lowest penalty
	 * @throws IllegalArgumentException
	 * @throws IllegalAccessException
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private static PrimerPair getPrimerToCheckDeletion(DeletionPlasmid plasmid, double maxRatioWtProductToDeletion, String pathPrimer3core) throws IllegalArgumentException, IllegalAccessException, IOException {
		return getPrimerToCheckDeletion(plasmid, DEFAULT_MIN_PCR_PRODUCT_SIZE, DEFAULT_MAX_PCR_PRODUCT_SIZE, maxRatioWtProductToDeletion, pathPrimer3core);
	}
	
	/**
	 * Get PCR primers flanking the deleted regions and write to a file
	 * @param plasmids The plasmids
	 * @param minProductSize Minimum desired PCR product size
	 * @param maxProductSize Maximum desired PCR product size
	 * @param maxRatioWtProductToDeletion Maximum ratio of total wild type product size (includes the deleted region and flanking regions) to the size of the deleted region itself
	 * @param pathPrimer3core primer3core executable
	 * @param outPrimerFile Output table of primers
	 * @throws IOException
	 * @throws IllegalArgumentException
	 * @throws IllegalAccessException
	 */
	private static void writePrimersToCheckDeletion(Collection<DeletionPlasmid> plasmids, int minProductSize, int maxProductSize, double maxRatioWtProductToDeletion, String pathPrimer3core, String outPrimerFile) throws IOException, IllegalArgumentException, IllegalAccessException {
		logger.info("Writing PCR primers to check for deletions to file " + outPrimerFile + "...");
		FileWriter w = new FileWriter(outPrimerFile);
		String header = PrimerPair.SEQUENCE_FIELD_NAME + "\t";
		header += "DeletionSize\t";
		header += "DeletionProductSize\t";
		header += "WildTypeProductSize\t";
		header += PrimerPair.getPrimerPairInformationFieldNames() + "\t";
		w.write(header + "\n");
		CountLogger countLogger = new CountLogger(plasmids.size());
		for(DeletionPlasmid plasmid : plasmids) {
			PrimerPair primer = getPrimerToCheckDeletion(plasmid, minProductSize, maxProductSize, maxRatioWtProductToDeletion, pathPrimer3core);
			if(primer != null) {
				String line = primer.getPrimerPairId() + "\t";
				line += plasmid.getDeletionSize() + "\t";
				line += Integer.valueOf(primer.getProductSize() - plasmid.getDeletionSize()).toString() + "\t";
				line += primer.getProductSize() + "\t";
				line += primer.getPrimerPairInformation() + "\t";
				w.write(line + "\n");
			} else {
				w.write(plasmid.getId() + "\tNO_PRIMERS\n");
			}
			countLogger.advance();
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * Get a PCR primer flanking the deleted region
	 * @param plasmid The plasmid
	 * @param minProductSize Minimum desired PCR product size
	 * @param maxProductSize Maximum desired PCR product size
	 * @param maxRatioWtProductToDeletion Maximum ratio of total wild type product size (includes the deleted region and flanking regions) to the size of the deleted region itself
	 * @param pathPrimer3core primer3core executable
	 * @return The primer pair flanking the deleted region with the lowest penalty
	 * @throws IllegalArgumentException
	 * @throws IllegalAccessException
	 * @throws IOException
	 */
	private static PrimerPair getPrimerToCheckDeletion(DeletionPlasmid plasmid, int minProductSize, int maxProductSize, double maxRatioWtProductToDeletion, String pathPrimer3core) throws IllegalArgumentException, IllegalAccessException, IOException {
		
		logger.info("Designing PCR primers to check for deletion in plasmid " + plasmid.getId());
		
		// Create primer3 configuration
        int deletionSize = plasmid.getDeletionSize();
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 62.0;
        config.minMeltingTemp = 58.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 40.0;
        config.maxGCContent = 60.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 5.0;
        config.selfEndAlignScore = 2.0;
        config.optimalPrimerSize = 23;
        config.minPrimerSize = 20;
        config.maxNumPrimersToReturn = 10000;
        config.canViolateConstraints = false;
        config.primerWindowSize = 0;
        config.maxPolyX = 3;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.minProductSize = minProductSize;
        config.maxProductSize = Math.min(maxProductSize, (int)(deletionSize * maxRatioWtProductToDeletion));
        
        logger.debug("Got primer3 configuration. Min product size = " + config.minProductSize + ". Max product size = " + config.maxProductSize + ".");
        
		// Check that bounds are consistent
		if(minProductSize > maxProductSize) {
			throw new IllegalArgumentException("Min product size must be <= max product size");
		}
		if(minProductSize <= 0 || maxProductSize <= 0 || config.minPrimerSize <= 0 || maxRatioWtProductToDeletion <= 0) {
			throw new IllegalArgumentException("All parameters must be >0");
		}
		if(2*config.minPrimerSize > maxProductSize) {
			throw new IllegalArgumentException("Max product size must be at least 2 * min primer size");
		}
        if(deletionSize + 2*config.minPrimerSize > maxProductSize) {
        	logger.warn("Can't create primer to test for deletion because desired max product size (" + maxProductSize + ") is smaller than the smallest possible product.");
        	return null;
        }
        int maxProduct = (int)(deletionSize * maxRatioWtProductToDeletion);
        if(maxProduct < minProductSize) {
        	logger.warn("Can't create primer to test for deletion because desired min product size (" + minProductSize + ") is larger than the largest desired product size based on ratio to deletion size (" + maxProduct + ").");
        	return null;
         }
                
        // make all possible primers
        Sequence flankingSequence = plasmid.getSequenceFlankingDeletion(maxProductSize, maxProductSize, true);
        logger.debug("Flanking sequence " + flankingSequence.getId());
        Collection<PrimerPair> allPrimers = PcrPrimerDesigner.designPrimers(config, flankingSequence, pathPrimer3core);
        
        if(allPrimers.isEmpty()) {
        	logger.warn("Couldn't find any acceptable primers flanking deletion for plasmid " + plasmid.getId());
        	return null;
        } 
        logger.debug("Got " + allPrimers.size() + " total primers.");
   
        
        // Get the primer with lowest penalty that straddles the deletion
        float minPenalty = Float.MAX_VALUE;
        PrimerPair rtrn = null;
        for(PrimerPair primer : allPrimers) {
        	if(primer.getPrimerPairPenalty() < minPenalty && primer.getLeftPrimerPosition() < maxProductSize && primer.getRightPrimerPosition() > maxProductSize + deletionSize) {
        		minPenalty = primer.getPrimerPairPenalty();
        		rtrn = primer;
        	}
        }
        
        if(rtrn == null) {
        	logger.warn("Couldn't find any acceptable primers flanking deletion for plasmid " + plasmid.getId() + ". Min and max product size: " + config.minProductSize + "," + config.maxProductSize);
        	return null;
        }
        
        logger.info("Primer with lowest penalty: " + rtrn.getPrimerPairId() + "\t" + rtrn.getPrimerPairInformation());
        
        return rtrn;
        
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws IllegalAccessException 
	 * @throws IllegalArgumentException 
	 */
	public static void main(String[] args) throws IOException, IllegalArgumentException, IllegalAccessException {
		
		CommandLineParser p = new CommandLineParser();
		String programDescription = "\n";
		programDescription += "Task 1: write sequences with deletions\n";
		programDescription += "Task 2: write bed file of primer pairs created for sequences from task 1\n";
		programDescription += "Task 3: write table of actual plasmids created from these primers\n";
		programDescription += "Task 4: create PCR primers to test for deletion\n";
		p.setProgramDescription(programDescription);
		p.addStringArg("-c", "For task 1 or 3: fasta file of chromosomes", false, null);
		p.addStringArg("-g", "For task 1, 2 or 3: single gene bed file", false, null);
		p.addStringArg("-d", "For task 1, 2 or 3: existing bed file of deletions", false, null);
		p.addStringArg("-of", "For task 1: output fasta file of plasmid sequences (requires -c, -d and -g)", false, null);
		p.addStringArg("-p", "For task 2 or 3: existing output file from PcrPrimerDesigner", false, null);
		p.addStringArg("-op", "For task 2: output bed file of primer pairs (requires -d, -g and -p)", false, null);
		p.addStringArg("-ot", "For task 3: output table of actual PCR products (requires -c, -d, -p and -g)", false, null);
		p.addStringArg("-odp", "For task 4: output table of PCR primers to check for deletion (requires -pt and -p3)", false, null);
		p.addStringArg("-pt", "For task 4: table of actual PCR products written by task 3", false, null);
		p.addIntArg("-minp", "For task 4: minimum product size", false, DEFAULT_MIN_PCR_PRODUCT_SIZE);
		p.addIntArg("-maxp", "For task 4: maximum product size", false, DEFAULT_MAX_PCR_PRODUCT_SIZE);
		p.addDoubleArg("-maxr", "For task 4: maximum ratio of total wild type product size (includes the deleted region and flanking regions) to the size of the deleted region itself", false, 2);
		p.addStringArg("-p3", "For task 4: primer3_core executable", false, null);
		p.addBooleanArg("-debug", "Debug logging", false, false);
		p.parse(args);
		String chrFasta = p.getStringArg("-c");
		String geneBed = p.getStringArg("-g");
		String deletionBed = p.getStringArg("-d");
		String outFasta = p.getStringArg("-of");
		String primerFile = p.getStringArg("-p");
		String outPrimerBed = p.getStringArg("-op");
		String outTableActualProducts = p.getStringArg("-ot");
		String outPrimersToCheckDeletion = p.getStringArg("-odp");
		String tableActualProducts = p.getStringArg("-pt");
		int minProdSize = p.getIntArg("-minp");
		int maxProdSize = p.getIntArg("-maxp");
		double maxRatio = p.getDoubleArg("-maxr");
		String primer3core = p.getStringArg("-p3");
		boolean debug = p.getBooleanArg("-debug");
		
		if(debug) {
			logger.setLevel(Level.DEBUG);
		}
		
		if(args.length == 0) {
			p.printHelpMessage();
		}
		
		if(outFasta != null) {
			if(chrFasta == null) {
				throw new IllegalArgumentException("To write fasta file must provide -c option");
			}
			if(geneBed == null) {
				throw new IllegalArgumentException("To write fasta file must provide -g option");
			}
			if(deletionBed == null) {
				throw new IllegalArgumentException("To write fasta file must provide -d option");
			}
			writeDeletionPlasmidSequences(chrFasta, geneBed, deletionBed, outFasta);
		}
		
		if(outPrimerBed != null) {
			if(primerFile == null) {
				throw new IllegalArgumentException("To write bed file of primers must provide -p option");
			}
			if(geneBed == null) {
				throw new IllegalArgumentException("To write bed file of primers must provide -g option");
			}
			if(deletionBed == null) {
				throw new IllegalArgumentException("To write bed file of primers must provide -d option");
			}
			writePrimerPairBedFile(primerFile, geneBed, deletionBed, outPrimerBed);
		}
		
		if(outTableActualProducts != null) {
			if(chrFasta == null) {
				throw new IllegalArgumentException("To write table of actual plasmids must provide -c option");
			}
			if(primerFile == null) {
				throw new IllegalArgumentException("To write table of actual plasmids must provide -p option");
			}
			if(geneBed == null) {
				throw new IllegalArgumentException("To write table of actual plasmids must provide -g option");
			}
			if(deletionBed == null) {
				throw new IllegalArgumentException("To write table of actual plasmids must provide -d option");
			}

			Collection<DeletionPlasmid> fullDesiredPlasmids = getDeletionPlasmids(chrFasta, geneBed, deletionBed);
			Collection<DeletionPlasmid> actualPlasmids = getActualPCRProducts(fullDesiredPlasmids, primerFile);
			writeToTable(actualPlasmids, outTableActualProducts);
		}
		
		if(outPrimersToCheckDeletion != null) {
			if(tableActualProducts == null) {
				throw new IllegalArgumentException("To write PCR primers to check for deletion, must provide -pt option");
			}
			if(primer3core == null) {
				throw new IllegalArgumentException("To write PCR primers to check for deletion, must provide -p3 option");
			}
			Collection<DeletionPlasmid> plasmids = readFromTable(tableActualProducts);
			writePrimersToCheckDeletion(plasmids, minProdSize, maxProdSize, maxRatio, primer3core, outPrimersToCheckDeletion);
		}
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}

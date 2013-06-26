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

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.PrimerPair;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class PlasmidUtils {
	
	private static Logger logger = Logger.getLogger(PlasmidUtils.class.getName());
	
	/**
	 * Get the sequence of a plasmid created by deleting a region
	 * @param chromosome The full chromosome containing the gene of interest
	 * @param transcript The gene of interest as an annotation on the chromosome
	 * @param regionToDelete The deletion site as an annotation on the chromosome
	 * @return The sequence beginning immediately 3' of the deletion site, circling around, and ending immediately 5' of the deletion site
	 */
	public static Sequence getPlasmidSequenceAfterDeletion(Sequence chromosome, Gene transcript, Gene regionToDelete) {
		
		if(!transcript.contains(regionToDelete)) {
			throw new IllegalArgumentException("Transcript must contain region to delete");
		}
		
		Sequence transcriptSequence = chromosome.getSubsequence(transcript);
		
		int firstDeletedPosition = transcript.genomicToTranscriptPosition(regionToDelete.transcriptToGenomicPosition(0));
		int lastDeletedPosition = transcript.genomicToTranscriptPosition(regionToDelete.transcriptToGenomicPosition(regionToDelete.size() - 1));
		
		return getPlasmidSequenceAfterDeletion(transcriptSequence, firstDeletedPosition, lastDeletedPosition, transcript.getName() + "_deletion_" + regionToDelete.toUCSC());
		
	}
	
	/**
	 * Get the sequence of a plasmid created by deleting a region
	 * @param fullTranscriptSequence The original transcript sequence
	 * @param firstDeletedPosition The first position of deletion site in transcript coordinates
	 * @param lastDeletedPosition The last position of deletion site (inclusive) in transcript coordinates
	 * @param returnSequenceName Name of returned sequence or null for default
	 * @return The sequence beginning immediately 3' of the deletion site, circling around, and ending immediately 5' of the deletion site
	 */
	public static Sequence getPlasmidSequenceAfterDeletion(Sequence fullTranscriptSequence, int firstDeletedPosition, int lastDeletedPosition, String returnSequenceName) {
		
		String name = returnSequenceName == null ? fullTranscriptSequence.getId() + "_deletion_" + firstDeletedPosition + "_" + lastDeletedPosition : returnSequenceName;
		Sequence beforeDeletion = fullTranscriptSequence.getSubSequence(null, 0, firstDeletedPosition);
		Sequence afterDeletion = fullTranscriptSequence.getSubSequence(name, lastDeletedPosition + 1, fullTranscriptSequence.getLength());
		afterDeletion.append(beforeDeletion.getSequenceBases());
		afterDeletion.uppercase();
		return afterDeletion;
		
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
	 * Get sequences of plasmids created by deleting a region and write to a fasta file.
	 * Output sequence begins immediately 3' of the deletion site, circles around, and ends immediately 5' of the deletion site
	 * @param chromosomeFastaFile Fasta file containing the chromosome the gene is on
	 * @param singleGeneBedFile Bed file containing the single gene of interest
	 * @param deletionsBedFile Bed file of deletion sites within the gene
	 * @param outFastaFile Output fasta file for plasmid sequences
	 * @throws IOException
	 */
	public static void writePlasmidSequencesAfterDeletion(String chromosomeFastaFile, String singleGeneBedFile, String deletionsBedFile, String outFastaFile) throws IOException {
		
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
		List<Sequence> plasmidSeqs = new ArrayList<Sequence>();
		Map<String, Collection<Gene>> deletions = BEDFileParser.loadDataByChr(new File(deletionsBedFile));
		for(String c : deletions.keySet()) {
			for(Gene deletion : deletions.get(c)) {
				logger.info("Getting plasmid sequence for deletion " + deletion.getName());
				plasmidSeqs.add(getPlasmidSequenceAfterDeletion(chrSequence, gene, deletion));
			}
		}
		logger.info("Got all deletions.");
		
		logger.info("Writing deletions to file " + outFastaFile + "...");
		FastaSequenceIO fastaWriter = new FastaSequenceIO(outFastaFile);
		fastaWriter.write(plasmidSeqs);
		logger.info("Done writing fasta file.");
		
	}
	
	/**
	 * Write bed file of primer pairs produced by PcrPrimerDesigner on the plasmid sequences with deletions
	 * @param pcrPrimerDesignerOutputFile Output file from PcrPrimerDesigner
	 * @param singleGeneBedFile Bed file of single gene of interest
	 * @param deletionsBedFile Bed file of deletions within gene
	 * @param outBedFile Bed file of genomic positions of primer pairs
	 * @throws IOException
	 */
	public static void writePrimerPairBedFile(String pcrPrimerDesignerOutputFile, String singleGeneBedFile, String deletionsBedFile, String outBedFile) throws IOException {
		
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
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-c", "Fasta file of chromosomes (required for -of)", false, null);
		p.addStringArg("-g", "Single gene bed file", true);
		p.addStringArg("-d", "Bed file of deletions within gene", true);
		p.addStringArg("-of", "Output fasta file of plasmid sequences (requires -c)", false, null);
		p.addStringArg("-p", "Output file from PcrPrimerDesigner (required for -op)", false, null);
		p.addStringArg("-op", "Output bed file of primer pairs (requires -p)", false, null);
		p.parse(args);
		String chrFasta = p.getStringArg("-c");
		String geneBed = p.getStringArg("-g");
		String deletionBed = p.getStringArg("-d");
		String outFasta = p.getStringArg("-of");
		String primerFile = p.getStringArg("-p");
		String outPrimerBed = p.getStringArg("-op");
		
		if(outFasta != null) {
			if(chrFasta == null) {
				throw new IllegalArgumentException("To write fasta file must provide -c option");
			}
			writePlasmidSequencesAfterDeletion(chrFasta, geneBed, deletionBed, outFasta);
		}
		
		if(outPrimerBed != null) {
			if(primerFile == null) {
				throw new IllegalArgumentException("To write bed file of primers must provide -p option");
			}
			writePrimerPairBedFile(primerFile, geneBed, deletionBed, outPrimerBed);
		}
		
		logger.info("");
		logger.info("All done.");
		
	}
	
}

/**
 * 
 */
package broad.core.overlaputils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;



/**
 * @author prussell
 *
 */
public class GeneSetIntersect {

	private static Logger logger = Logger.getLogger(GeneSetIntersect.class.getName());
	
	/**
	 * Get overlap between two gene sets
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes) {
		return getOverlap(genes, intersectGenes, false);
	}
	
	/**
	 * Get overlap between two gene sets
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @param verbose Print messages
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, boolean verbose) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		return getOverlap(genes, intersectGenes, excludeGenes, verbose);
	}
	
	/**
	 * Get overlap between two gene sets and automatically exclude all genes in a set
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @param excludeGenes Exclude genes from first gene set if they overlap a gene from this set
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, Map<String, Collection<Gene>> excludeGenes) {
		return getOverlap(genes, intersectGenes, excludeGenes, false);
	}
	
	/**
	 * Get overlap between two gene sets and automatically exclude all genes in a set
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @param excludeGenes Exclude genes from first gene set if they overlap a gene from this set
	 * @param verbose Print messages
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, Map<String, Collection<Gene>> excludeGenes, boolean verbose) {
		if(excludeGenes.isEmpty()) {
			for(String chr : genes.keySet()) {
				Collection<Gene> empty = new TreeSet<Gene>();
				excludeGenes.put(chr, empty);
			}
		}
		Map<String, Collection<Gene>> rtrn = new TreeMap<String, Collection<Gene>>();
		for(String chr : genes.keySet()) {
			if(verbose) {
				logger.info(chr);
			}
			Collection<Gene> thisChrOverlaps = new TreeSet<Gene>();
			for(Gene gene : genes.get(chr)) {
				for(Gene other : excludeGenes.get(chr)) {
					if(gene.overlaps(other)) continue;
				}
				if(!intersectGenes.containsKey(chr)) continue;
				Gene overlap = gene.getOverlap(intersectGenes.get(chr));
				if(overlap == null) continue;
				overlap.setName(gene.getName());
				thisChrOverlaps.add(overlap);
			}
			rtrn.put(chr, thisChrOverlaps);
		}
		return rtrn;
	}
	
	/**
	 * Intersect a set of genes with another set and randomize positions within each gene in the second set
	 * @param genes The genes of interest
	 * @param intersectGenes The genes to intersect with and randomize within
	 * @return Genes that have been intersected with the intersect set and then randomized within each intersect gene
	 */
	public static Map<String, Collection<Gene>> getOverlapWithRandomizedPositions(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		for(String chr : genes.keySet()) {
			Collection<Gene> empty = new TreeSet<Gene>();
			excludeGenes.put(chr, empty);
		}
		return getOverlapWithRandomizedPositions(genes, intersectGenes, excludeGenes);
	}
	
	/**
	 * Intersect a set of genes with another set and randomize positions within each gene in the second set
	 * @param genes The genes of interest
	 * @param intersectGenes The genes to intersect with and randomize within
	 * @param excludeFeatures Do not consider genes that overlap one of these
	 * @return Genes that have been intersected with the intersect set and then randomized within each intersect gene
	 */
	public static Map<String, Collection<Gene>> getOverlapWithRandomizedPositions(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, Map<String, Collection<Gene>> excludeFeatures) {
		Map<String, Collection<Gene>> rtrn = new TreeMap<String, Collection<Gene>>();
		logger.info("Getting overlaps with randomized positions...");
		Random rand = new Random();
		for(String chr : intersectGenes.keySet()) {
			logger.info(chr);
			Collection<Gene> randomizedOverlapsThisChr = new TreeSet<Gene>();
			if(!genes.containsKey(chr)) continue;
			if(genes.get(chr).isEmpty()) continue;
			Map<String, Collection<Gene>> genesThisChr = new HashMap<String, Collection<Gene>>();
			genesThisChr.put(chr, genes.get(chr));
			for(Gene feature : intersectGenes.get(chr)) {
				for(Gene gene : genesThisChr.get(chr)) {
					if(feature.getStart() > gene.getEnd()) continue;
					if(feature.getEnd() < gene.getStart()) break;
					if(feature.overlaps(gene) && feature.getStart() >= gene.getStart() && feature.getEnd() <= gene.getEnd()) {
						int overlapperStart = gene.genomicToTranscriptPosition(feature.getStart());
						int overlapperEndInclusive = gene.genomicToTranscriptPosition(feature.getEnd() - 1);
						int overlapperStartOnTranscript = Math.min(overlapperStart,overlapperEndInclusive);
						int overlapperEndOnTranscript = Math.max(overlapperStart, overlapperEndInclusive) + 1;
						int overlapperLengthOnTranscript = overlapperEndOnTranscript - overlapperStartOnTranscript;
						int lastPossibleRandomStartPos = gene.getSize() - overlapperLengthOnTranscript;
						int randomStartOnTranscript = rand.nextInt(lastPossibleRandomStartPos + 1);
						int randomEndOnTranscript = randomStartOnTranscript + overlapperLengthOnTranscript + 1;
						Gene randomOverlap = gene.copy();
						try {
							randomOverlap.trim(randomStartOnTranscript, gene.getSize() - randomEndOnTranscript);
						} catch (IllegalArgumentException e) {
							logger.warn("Caught exception when randomizing position of feature " + feature.getName() + " within gene " + gene.getName() + ". Skipping.");
							continue;
						}
						randomOverlap.setName(gene.getName() + "_" + feature.getName() + "_random");
						randomOverlap.setOrientation(feature.getOrientation());
						randomizedOverlapsThisChr.add(randomOverlap);
					}
				}
			}
			rtrn.put(chr, randomizedOverlapsThisChr);
		}
		return rtrn;
	}
	
	/**
	 * Get sum of gene sizes
	 * @param genes The gene set by chromosome
	 * @return Total transcribed bases including double counts if genes overlap
	 */
	private static int totalBases(Map<String, Collection<Gene>> genes) {
		int rtrn = 0;
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				rtrn += gene.getSize();
			}
		}
		return rtrn;
	}
	
	/**
	 * Get total size of gene set
	 * @param genes The gene set by chromosome
	 * @return Total number of genes
	 */
	private static int totalGenes(Map<String, Collection<Gene>> genes) {
		int rtrn = 0;
		for(String chr : genes.keySet()) rtrn += genes.get(chr).size();
		return rtrn;
	}


	/**
	 * Get number of overlapping bases between two gene sets
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @return The total number of overlapping bases including double counts if multiple genes overlap the same gene
	 */
	public static int numOverlappingBases(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		return numOverlappingBases(genes, intersectGenes, excludeGenes);
	}
	
	/**
	 * Get number of overlapping bases between two gene sets and automatically exclude all genes in a set
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @param excludeGenes Exclude genes from first gene set if they overlap a gene from this set
	 * @return The total number of overlapping bases including double counts if multiple genes overlap the same gene
	 */
	public static int numOverlappingBases(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, Map<String, Collection<Gene>> excludeGenes) {
		Map<String, Collection<Gene>> overlap = getOverlap(genes, intersectGenes, excludeGenes);
		return totalBases(overlap);
	}

	/**
	 * Get number of overlapping genes between two gene sets
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @return The total number of overlapping bases including double counts if multiple genes overlap the same gene
	 */
	public static int numOverlappingGenes(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		return numOverlappingGenes(genes, intersectGenes, excludeGenes);
	}
	
	/**
	 * Get number of overlapping genes between two gene sets and automatically exclude all genes in a set
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @param excludeGenes Exclude genes from first gene set if they overlap a gene from this set
	 * @return The total number of overlapping bases including double counts if multiple genes overlap the same gene
	 */
	public static int numOverlappingGenes(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, Map<String, Collection<Gene>> excludeGenes) {
		Map<String, Collection<Gene>> overlap = getOverlap(genes, intersectGenes, excludeGenes);
		return totalGenes(overlap);
	}


	/**
	 * Get overlap between two gene sets specified in bed files and write to a file
	 * @param geneFile Bed file of first gene set
	 * @param intersectFile Bed file of genes to intersect with
	 * @param outFile Output bed file to write Genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(String geneFile, String intersectFile, String outFile, boolean randomizePositionsWithinEachGene) throws IOException {
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneFile));
		Map<String, Collection<Gene>> intersectGenes = BEDFileParser.loadDataByChr(new File(intersectFile));
		writeOverlap(genes, intersectGenes, outFile, randomizePositionsWithinEachGene);
	}

	/**
	 * Get overlap between two gene sets specified in bed files, automatically excluding all genes in a set, and write to a file
	 * @param geneFile Bed file of first gene set
	 * @param intersectFile Bed file of genes to intersect with
	 * @param excludeFile Exclude genes from first gene set if they overlap a gene from this set
	 * @param outFile Output bed file to write Genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(String geneFile, String intersectFile, String excludeFile, String outFile, boolean randomizePositionsWithinEachGene) throws IOException {
		logger.info("Loading genes...");
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneFile));
		logger.info("Loading annotations to intersect with...");
		Map<String, Collection<Gene>> intersectGenes = BEDFileParser.loadDataByChr(new File(intersectFile));
		logger.info("Loading annotations to exclude...");
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		if(excludeFile != null) excludeGenes.putAll(BEDFileParser.loadDataByChr(new File(excludeFile)));
		logger.info("Done loading all annotations.");
		writeOverlap(genes, intersectGenes, excludeGenes, outFile, randomizePositionsWithinEachGene);
	}

	/**
	 * Get overlap between two gene sets and write to a file
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @param outFile Output bed file to write Genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, String outFile, boolean randomizePositionsWithinEachGene) throws IOException { 
		Map<String,Collection<Gene>> excludeGenes = new TreeMap<String,Collection<Gene>>();
		writeOverlap(genes, intersectGenes, excludeGenes, outFile, randomizePositionsWithinEachGene);
	}
	
	/**
	 * Get overlap between two gene sets, automatically excluding all genes in a set, and write to a file
	 * @param genes First gene set by chromosome
	 * @param intersectGenes Genes to intersect with by chromosome
	 * @param excludeGenes Exclude genes from first gene set if they overlap a gene from this set
	 * @param outFile Output bed file to write Genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectGenes, Map<String, Collection<Gene>> excludeGenes, String outFile, boolean randomizePositionsWithinEachGene) throws IOException {
		FileWriter w = new FileWriter(outFile);
		Map<String, Collection<Gene>> overlaps = new TreeMap<String, Collection<Gene>>();
		if(!randomizePositionsWithinEachGene) {
			logger.info("Getting overlaps...");
			overlaps.putAll(getOverlap(genes, intersectGenes, excludeGenes));
		}
		else overlaps.putAll(getOverlapWithRandomizedPositions(genes, intersectGenes, excludeGenes));
		logger.info("Writing overlaps to file " + outFile + "...");
		for(String chr : overlaps.keySet()) {
			for(Gene gene : overlaps.get(chr)) {
				w.write(gene.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file containing genes of interest", true);
		p.addStringArg("-i1", "Bed file containing genes to intersect with first (set1)", true);
		p.addStringArg("-i2", "Bed file containing genes to intersect with second (set2)", false, null);
		p.addBooleanArg("-r", "Only has effect if -ob is provided. Intersect genes with set1 and write bed file. Ignore set2. If true, randomize position of overlap within set1 gene. ", false, false);
		p.addIntArg("-nr", "Randomize overlap positions within set1 this many times, intersect with set2, and print counts only.", false, 0);
		p.addStringArg("-e", "Bed file containing genes to exclude - exclude any gene in input file that overlaps one of these", false, null);
		p.addStringArg("-ob", "Output bed file for intersection with set1 (randomized or not)", false, null);
		p.addStringArg("-oc", "Output counts file for intersections with set1, randomized and intersected with set2. Requires -i2 and -nr.", false, null);
		p.parse(args);
		String geneFile = p.getStringArg("-b");
		String intersectFile1 = p.getStringArg("-i1");
		String intersectFile2 = p.getStringArg("-i2");
		String excludeFile = p.getStringArg("-e");
		int numRandom = p.getIntArg("-nr");
		boolean randomize = p.getBooleanArg("-r");
		String outBed = p.getStringArg("-ob");
		String outRandCounts = p.getStringArg("-oc");
		
		if(outBed != null) {
			logger.info("Writing overlap between " + geneFile + " and " + intersectFile1 + " to file " + outBed + "...");
			writeOverlap(geneFile, intersectFile1, excludeFile, outBed, randomize);
		}
		
		if(outRandCounts != null) {
			
			if(intersectFile2 == null) {
				throw new IllegalArgumentException("-i2 is required for randomized overlap counts.");
			}
			
			FileWriter w = new FileWriter(outRandCounts);
			w.write("test\ttotal_bases_overlap\n");
			Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneFile));
			Map<String, Collection<Gene>> intersectFeatures1 = BEDFileParser.loadDataByChr(new File(intersectFile1));
			Map<String, Collection<Gene>> intersectFeatures2 = BEDFileParser.loadDataByChr(new File(intersectFile2));
			Map<String, Collection<Gene>> excludeFeatures = new TreeMap<String, Collection<Gene>>();
			if(excludeFile != null) excludeFeatures.putAll(BEDFileParser.loadDataByChr(new File(excludeFile)));
			
			logger.info("");
			logger.info("Getting overlap between genes and set1...");
			Map<String, Collection<Gene>> overlapSet1 = getOverlap(genes, intersectFeatures1, excludeFeatures, true);
			logger.info("Overlap is " + totalGenes(overlapSet1) + " genes and " + totalBases(overlapSet1) + " bases.");
			
			logger.info("");
			logger.info("Getting overlap between genes and set1+set2...");
			Map<String, Collection<Gene>> overlapSet2 = getOverlap(overlapSet1, intersectFeatures2, true);
			logger.info("Overlap is " + totalGenes(overlapSet2) + " genes and " + totalBases(overlapSet2) + " bases.");
			logger.info("COUNT\treal_overlap\t" + totalBases(overlapSet2));
			w.write("real_overlap\t" + totalBases(overlapSet2) + "\n");
			
			logger.info("");
			logger.info("Randomizing overlap positions of set1 within genes...");
			for(int i=0; i < numRandom; i++) {
				logger.info("Random iteration " + Integer.valueOf(i+1).toString());
				Map<String, Collection<Gene>> randomizedOverlap = getOverlapWithRandomizedPositions(genes, intersectFeatures1);
				Map<String, Collection<Gene>> randomizedWithSet2 = getOverlap(randomizedOverlap, intersectFeatures2, true);
				//writeOverlap(randomizedOverlap, intersectFeatures2, "test_" + i + "_" + outRandCounts + ".bed", false);
				logger.info("COUNT\toverlap_randomized_" + i + "\t" + totalBases(randomizedWithSet2));
				w.write("overlap_randomized_" + i + "\t" + totalBases(randomizedWithSet2) + "\n");
			}
			
			logger.info("");
			logger.info("All done. Wrote counts to file " + outRandCounts + ".");
			
			w.close();
		}

	}

}

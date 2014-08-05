/**
 * 
 */
package nextgen.core.programs;

import general.CommandLineParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class FastaAnnotationExtractor {

	private Map<String, Collection<Gene>> genes;
	private Map<String, Sequence> chromosomes;
	private static Logger logger = Logger.getLogger(FastaAnnotationExtractor.class.getName());
	
	/**
	 * @param genomeFasta Fasta file of sequences
	 * @param bedFile Bed file
	 * @throws IOException
	 */
	public FastaAnnotationExtractor(String genomeFasta, String bedFile) throws IOException {
		this(FastaSequenceIO.getChrSequencesFromFasta(genomeFasta), bedFile);
	}
	
	/**
	 * Load sequences and annotations
	 * @param chrsByName Chromosomes by name
	 * @param bedFile Genes in bed format
	 * @throws IOException
	 */
	public FastaAnnotationExtractor(Map<String, Sequence> chrsByName, String bedFile) throws IOException {
		
		logger.info("Loading genes from file " + bedFile + "...");
		
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		int numGenes = 0;
		for(String chr : genes.keySet()) {
			numGenes += genes.get(chr).size();
		}
		logger.info("Loaded " + numGenes + " genes.");
		
		chromosomes = chrsByName;
		
	}
	
	private void maskChromosomes(String bedRegions, boolean softmask) throws IOException {
		logger.info("Masking sequences in " + bedRegions + ". Softmask: " + softmask);
		Map<String, Collection<Gene>> regionsToMask = BEDFileParser.loadDataByChr(new File(bedRegions));
		for(String chr : regionsToMask.keySet()) {
			Sequence chrSeq = chromosomes.get(chr);
			for(Gene gene : regionsToMask.get(chr)) {
				chrSeq.mask(gene, softmask);
			}
		}
	}
	
	/**
	 * Get the spliced sequence of the gene
	 * @param gene The gene
	 * @return The sequence in 5' to 3' orientation
	 */
	private Sequence getGeneSequence(Gene gene) {
		String chr = gene.getReferenceName();
		return chromosomes.get(chr).getSubsequence(gene);
	}
	
	/**
	 * 
	 * @param outFile
	 * @throws IOException
	 */
	public void writeFasta(String outFile) throws IOException {
		logger.info("Writing gene sequences to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		BufferedWriter b = new BufferedWriter(w);
		FastaSequenceIO fsio = new FastaSequenceIO();
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				Sequence seq = getGeneSequence(gene);
				fsio.write(seq, b);
			}
		}
		b.close();
		logger.info("Done writing sequences.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of regions to extract", true);
		p.addStringArg("-f", "Fasta file of chromosomes", true);
		p.addStringArg("-m", "Bed file of regions to mask", false, null);
		p.addBooleanArg("-s", "If masking, mask to lower case (else Ns)", false, true);
		p.addStringArg("-o", "Output fasta file of sequences", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String fastaFile = p.getStringArg("-f");
		String outFile = p.getStringArg("-o");
		String maskBed = p.getStringArg("-m");
		boolean softmask = p.getBooleanArg("-s");
		
		FastaAnnotationExtractor f = new FastaAnnotationExtractor(fastaFile, bedFile);
		
		if(maskBed != null) {
			f.maskChromosomes(maskBed, softmask);
		}
		
		f.writeFasta(outFile);

	}

}

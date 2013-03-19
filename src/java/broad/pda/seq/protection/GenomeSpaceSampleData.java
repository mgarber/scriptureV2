/**
 * 
 */
package broad.pda.seq.protection;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;

/**
 * @author prussell
 *
 */
public class GenomeSpaceSampleData extends SampleData {

	private AlignmentModel genomeData;
	private Map<Gene, ScanStatisticScore> genomeScores;

	/**
	 * @param bamFile Bam alignment file
	 * @param chrSizeFile Chromosome size file
	 * @param genes Genes by chromosome
	 * @param window Window size
	 * @param step Step size
	 * @param expressionCutoff Genome wide scan P value cutoff for expression
	 * @throws IOException
	 */
	public GenomeSpaceSampleData(String bamFile, String chrSizeFile, Map<String, Collection<Gene>> genes, int window, int step, double expressionCutoff) throws IOException {
		super(bamFile, genes, window, step, expressionCutoff, true);
		genomeData = new AlignmentModel(bamFile, new GenomicSpace(chrSizeFile));
		genomeScores = new TreeMap<Gene, ScanStatisticScore>();
	}

	/**
	 * Whether the gene is significantly expressed genome wide
	 * @param gene The gene
	 * @return Whether the gene is expressed at the given significance level
	 */
	@Override
	public boolean isExpressed(Gene gene) {
		if(!genomeScores.containsKey(gene)) {
			ScanStatisticScore score = new ScanStatisticScore(genomeData, gene);
			logger.debug(gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + " count " + score.getCount() + " global lambda " + score.getGlobalLambda() + " size " + score.getCoordinateSpace().getSize(gene) + " global length " + score.getGlobalLength() + " pval " + score.getScanPvalue());
			genomeScores.put(gene, score);			
		}
		ScanStatisticScore score = genomeScores.get(gene);
		return score.getScanPvalue() <= expressionCutoffValue;
	}


	/**
	 * Get genome wide scan P value of number of fragments mapping to the gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The scan P value of the number of fragments mapping to the gene
	 */
	@Override
	public double getGeneScanPval(Gene gene) {
		if(!genomeScores.containsKey(gene)) {
			ScanStatisticScore score = new ScanStatisticScore(genomeData, gene);
			logger.debug(gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + " count " + score.getCount() + " global lambda " + score.getGlobalLambda() + " size " + score.getCoordinateSpace().getSize(gene) + " global length " + score.getGlobalLength() + " pval " + score.getScanPvalue());
			genomeScores.put(gene, score);			
		}
		ScanStatisticScore score = genomeScores.get(gene);
		return score.getScanPvalue();
	}

	
}

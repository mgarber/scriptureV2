/**
 * 
 */
package nextgen.core.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;

/**
 * @author prussell
 *
 */
public class WigWriter {
	
	private TranscriptomeSpace transcriptomeSpace;
	private GenomicSpace genomeSpace;
	private Map<String, Collection<Gene>> genesByChr;
	private AlignmentModel data;
	private Collection<String> chrNames;
	private boolean readBeginningPositionOnly;
	private boolean isTranscriptomeSpace;
	private static boolean DEFAULT_USE_FRAGMENTS = true;
	private static int DEFAULT_MAX_FRAGMENT_LENGTH = 2000;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	static Logger logger = Logger.getLogger(WigWriter.class.getName());

	/**
	 * Construct with a transcriptome space
	 * @param bamFile Bam alignments
	 * @param genesByChrName Genes by reference sequence name
	 * @param firstPositionOnly Only count the first position of each read
	 * @param useFragments Whether to convert reads to fragments
	 */
	public WigWriter(String bamFile, Map<String, Collection<Gene>> genesByChrName, boolean firstPositionOnly, boolean useFragments) {
		readBeginningPositionOnly = firstPositionOnly;
		genesByChr = genesByChrName;
		chrNames = new TreeSet<String>();
		chrNames.addAll(genesByChr.keySet());
		transcriptomeSpace = new TranscriptomeSpace(genesByChr);
		data = new AlignmentModel(bamFile, transcriptomeSpace, useFragments);
		isTranscriptomeSpace = true;
	}
	
	/**
	 * Construct with a genome space
	 * @param bamFile Bam alignments
	 * @param chrSizeFile Chromosome size file
	 * @param firstPositionOnly Only count the first position of each read
	 * @param useFragments Whether to convert reads to fragments
	 */
	public WigWriter(String bamFile, String chrSizeFile, boolean firstPositionOnly, boolean useFragments) {
		readBeginningPositionOnly = firstPositionOnly;
		genomeSpace = new GenomicSpace(chrSizeFile);
		chrNames = new TreeSet<String>();
		chrNames.addAll(genomeSpace.getReferenceNames());
		data = new AlignmentModel(bamFile, genomeSpace, useFragments);
		isTranscriptomeSpace = false;
	}
	
	
	
	/**
	 * Construct with either a genomic space or a transcriptome space
	 * @param bamFile Bam alignments
	 * @param geneBedFile Bed file of annotations for transcriptome space
	 * @param chrSizeFile File of chromosome sizes for genomic space
	 * @param firstPositionOnly Only count the first position of each read
	 * @param useFragments Whether to convert reads to fragments
	 * @throws IOException
	 */
	public WigWriter(String bamFile, String geneBedFile, String chrSizeFile, boolean firstPositionOnly, boolean useFragments) throws IOException {
		if((geneBedFile != null && chrSizeFile != null) || (geneBedFile == null && chrSizeFile == null)) {
			throw new IllegalArgumentException("Choose one: gene bed file (for transcriptome space) or chromosome size file (for genomic space)");
		}
		readBeginningPositionOnly = firstPositionOnly;
		if(geneBedFile != null) {
			genesByChr = BEDFileParser.loadDataByChr(new File(geneBedFile));
			chrNames = new TreeSet<String>();
			chrNames.addAll(genesByChr.keySet());
			transcriptomeSpace = new TranscriptomeSpace(genesByChr);
			data = new AlignmentModel(bamFile, transcriptomeSpace, useFragments);
			isTranscriptomeSpace = true;
			return;
		}
		if(chrSizeFile != null) {
			genesByChr = null;
			genomeSpace = new GenomicSpace(chrSizeFile);
			chrNames = new TreeSet<String>();
			chrNames.addAll(genomeSpace.getReferenceNames());
			data = new AlignmentModel(bamFile, genomeSpace, useFragments);
			isTranscriptomeSpace = false;
			return;
		}
	}
	
	/**
	 * Add a read filter to alignment model
	 * @param filter The filter
	 */
	public void addReadFilter(Predicate<Alignment> filter) {
		data.addFilter(filter);
	}
	
	/**
	 * Get counts by position
	 * @param region The region
	 * @return Map of position to count, includes only positions with at least one read
	 * @throws IOException
	 */
	private TreeMap<Integer, Double> getCounts(Annotation region) throws IOException {
		
		// Include entire read/fragment in transcriptome space
		if(!readBeginningPositionOnly && isTranscriptomeSpace) {
			/*
			 *  TODO:
			 *  Will be way too slow
			 *  Use AlignmentModel.scan() instead
			 */
			return transcriptomeSpace.getPositionCountMap(new Gene(region), data);
		}
		
		// Include entire read/fragment in genomic space
		if(!readBeginningPositionOnly && !isTranscriptomeSpace) {
			/*
			 * TODO:
			 * Implement for full read/fragment in genome space
			 */
			throw new IllegalArgumentException("Method not implemented for full reads/fragments in genomic space.");
		}
		
		// Count last read position only
		CloseableIterator<Alignment> iter = data.getOverlappingReads(region, false);
		TreeMap<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		while(iter.hasNext()) {
			Alignment read = iter.next();
			int lastPos = read.getFirstFragmentPositionStranded();
			read = null;
			Integer i = Integer.valueOf(lastPos);
			if(rtrn.containsKey(i)) {
				double oldVal = rtrn.get(i).doubleValue();
				double newVal = oldVal + 1;
				rtrn.put(i, Double.valueOf(newVal));
			} else {
				rtrn.put(i, Double.valueOf(1));
			}
			i = null;
		}
		iter.close();
		return rtrn;
		
	}
	
	/**
	 * Get counts across a whole chromosome for positions in transcriptome space
	 * @param chr Chromosome name
	 * @return Map of position to count, includes only positions with at least one read
	 * @throws IOException
	 */
	private TreeMap<Integer, Double> getCountsInTranscriptomeSpace(String chr) throws IOException {
		if(!isTranscriptomeSpace) {
			throw new IllegalStateException("Must instantiate alignment model with a transcriptome space");
		}
		TreeMap<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		for(Gene gene : genesByChr.get(chr)) {
			rtrn.putAll(getCounts(gene));
		}
		return rtrn;
	}
	
	/**
	 * Get counts across a whole chromosome
	 * @param chr Chromosome name
	 * @return Map of position to count, includes only positions with at least one read
	 * @throws IOException
	 */
	private TreeMap<Integer, Double> getCountsInGenomeSpace(String chr) throws IOException {
		if(isTranscriptomeSpace) {
			throw new IllegalStateException("Must instantiate alignment model with a genome space");
		}
		return getCounts(genomeSpace.getEntireChromosome(chr));
	}
	
	
	
	/**
	 * Write all counts to a wig file
	 * @param outFile Output file
	 * @throws IOException
	 */
	public void writeWig(String outFile) throws IOException {
		FileWriter w = new FileWriter(outFile);
		for(String chrName : chrNames) {
			logger.info("Writing counts for chromosome " + chrName + "...");
			TreeMap<Integer, Double> counts = new TreeMap<Integer, Double>();
			if(isTranscriptomeSpace) {
				if(genesByChr.get(chrName).isEmpty()) continue;
				counts = getCountsInTranscriptomeSpace(chrName);
			}
			else {
				counts = getCountsInGenomeSpace(chrName);
				// End position of chromosome is off the end - not valid position for wig format
				counts.remove(Integer.valueOf(Long.valueOf(genomeSpace.getLength(chrName)).intValue()));
			}
			w.write("variableStep chrom=" + chrName + "\n");
			for(Integer i : counts.keySet()) {
				int pos = i.intValue() + 1;
				w.write(pos + "\t" + counts.get(i).toString() + "\n");
			}
		}
		w.close();
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Gene bed file", false, null);
		p.addStringArg("-c", "Chromosome size file", false, null);
		p.addStringArg("-o", "Output file ending in .wig", true);
		p.addIntArg("-mf", "Max fragment length for paired reads", false, DEFAULT_MAX_FRAGMENT_LENGTH);
		p.addIntArg("-mg", "Max genomic span for paired reads", false, DEFAULT_MAX_GENOMIC_SPAN);
		p.addBooleanArg("-f", "Count beginning position of each read only", false, false);
		p.addBooleanArg("-pe", "Convert paired ends to fragments", false, DEFAULT_USE_FRAGMENTS);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-c");
		String outFile = p.getStringArg("-o");
		int maxFragmentLength = p.getIntArg("-mf");
		boolean fragments = p.getBooleanArg("-pe");
		int maxGenomicSpan = p.getIntArg("-mg");
		boolean firstPositionOnly = p.getBooleanArg("-f");
		
		WigWriter ww = new WigWriter(bamFile, bedFile, chrSizeFile, firstPositionOnly, fragments);
		ww.addReadFilter(new FragmentLengthFilter(ww.data.getCoordinateSpace(), maxFragmentLength));
		ww.addReadFilter(new GenomicSpanFilter(maxGenomicSpan));
		ww.writeWig(outFile);
		
	}

}

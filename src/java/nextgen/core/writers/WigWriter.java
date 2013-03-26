/**
 * 
 */
package nextgen.core.writers;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowScoreIterator;
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
	private AlignmentModel normData;
	private Collection<String> chrNames;
	private boolean readBeginningPositionOnly;
	private boolean isTranscriptomeSpace;
	private boolean normalize;
	private boolean bothFiles = false;
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
	public WigWriter(String bamFile, Map<String, Collection<Gene>> genesByChrName, boolean firstPositionOnly, boolean useFragments, boolean nor) {
		readBeginningPositionOnly = firstPositionOnly;
		genesByChr = genesByChrName;
		chrNames = new TreeSet<String>();
		chrNames.addAll(genesByChr.keySet());
		transcriptomeSpace = new TranscriptomeSpace(genesByChr);
		data = new AlignmentModel(bamFile, transcriptomeSpace, useFragments);
		isTranscriptomeSpace = true;
		normalize = nor;
	}
	
	/**
	 * Construct with a genome space
	 * @param bamFile Bam alignments
	 * @param chrSizeFile Chromosome size file
	 * @param firstPositionOnly Only count the first position of each read
	 * @param useFragments Whether to convert reads to fragments
	 */
	public WigWriter(String bamFile, String chrSizeFile, boolean firstPositionOnly, boolean useFragments, boolean nor) {
		readBeginningPositionOnly = firstPositionOnly;
		genomeSpace = new GenomicSpace(chrSizeFile);
		chrNames = new TreeSet<String>();
		chrNames.addAll(genomeSpace.getReferenceNames());
		data = new AlignmentModel(bamFile, genomeSpace, useFragments);
		isTranscriptomeSpace = false;
		normalize = nor;
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
	public WigWriter(String bamFile, String geneBedFile, String chrSizeFile, boolean firstPositionOnly, boolean useFragments, boolean nor) throws IOException {
		readBeginningPositionOnly = firstPositionOnly;
		normalize = nor;
		if (geneBedFile == null && chrSizeFile == null) {
			throw new IllegalArgumentException("Choose one or both: gene bed file or chromosome size file");
		}
		if((geneBedFile != null && chrSizeFile != null)) {
			//logger.info("constructing using both bed and size files.");
			genomeSpace = new GenomicSpace(chrSizeFile);
			genesByChr = BEDFileParser.loadDataByChr(new File(geneBedFile));
			chrNames = new TreeSet<String>();
			chrNames.addAll(genesByChr.keySet());
			transcriptomeSpace = new TranscriptomeSpace(genesByChr);
			normData = new AlignmentModel(bamFile,transcriptomeSpace, useFragments);
			isTranscriptomeSpace = false;
			bothFiles = true;
			data = new AlignmentModel(bamFile, genomeSpace, useFragments);
			return;
			
		}
		
		if(geneBedFile != null) {
			//logger.info("Entered proper constructor for use of transcriptome space only");
			genesByChr = BEDFileParser.loadDataByChr(new File(geneBedFile));
			chrNames = new TreeSet<String>();
			chrNames.addAll(genesByChr.keySet());
			if (!firstPositionOnly) {
				TreeMap<String,Collection<Gene>> collGenesByChr = new TreeMap<String,Collection<Gene>>();
				Collection<Gene> collGenes;
				for (String chrom : genesByChr.keySet()) {
					collGenes = collapseGenes(genesByChr.get(chrom));
					collGenesByChr.put(chrom, collGenes);
				}
				this.genesByChr = collGenesByChr;
			}
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
			//logger.info("Calculating counts using bed file only");
			double norm = data.getCount(region)/region.length(); 
			WindowScoreIterator<CountScore> wIter = data.scan(region,1,0);
			TreeMap<Integer,Double> rtrn = new TreeMap<Integer,Double>();
			//logger.info("Calculating scores for gene " + region + "...");
			while (wIter.hasNext()){
				CountScore score = wIter.next();
				Annotation a = score.getAnnotation();
				int pos = a.getStart();
				if (normalize) {
					if (score.getCount()>0) {
						rtrn.put(pos, score.getCount()/norm); 
					}
				}else {
					if (score.getCount()>0) {
						rtrn.put(pos, score.getCount());
					}
				}
			}
			return rtrn;
		}
		
		// Include entire read/fragment in genomic space
		if(!readBeginningPositionOnly && !isTranscriptomeSpace) {
			
			if (bothFiles) {
				//logger.info("Using both files to write counts...");
				Collection<Gene> collapsedGenes = collapseGenes(genesByChr.get(region.getChr()));
				Iterator<Gene> winItr = collapsedGenes.iterator();
				Gene currGene;

				double geneCount;
				CountScore score;

				//Map<String,Double> normValues = getNormalizationValues(collapsedGenes);
				TreeMap<Integer,Double> rtrn = new TreeMap<Integer,Double>();
				while (winItr.hasNext()) {
					currGene = winItr.next();
					//logger.info("Current Gene:\t" + currGene);
					if (normalize) {
						//logger.info("getting count for normalization");
						geneCount = normData.getCount(currGene,true);
						//logger.info("got count");
					} else {
						geneCount = currGene.getSize();
					}
					
					Iterator<? extends Annotation> exonItr = currGene.getExonSet().iterator();
					while (exonItr.hasNext()) {
					
						WindowScoreIterator<CountScore> wIter = data.scan(exonItr.next(), 1,0);
						Annotation a;
						int pos;
						double count;
						
						while (wIter.hasNext()){
							score = wIter.next();
							//logger.info("Checking score: " + score);
							a = score.getAnnotation();
							pos = a.getStart();
							count = score.getCount();
							if (normalize) {
								if (count>0) {
									if (geneCount>0 && currGene.getSize()>0)
										rtrn.put(pos, count/(geneCount/currGene.getSize())); 
								}
							} else {
								if (count>0) {
									rtrn.put(pos, count);
								}
							}
							score = null;
						}
					}
					
				}
				return rtrn;
			} else {
				//logger.info("Using chromosome size file to find counts...");
				double norm = 1;
				if (normalize) {
					norm = (double) data.getCount(region)/(genomeSpace.getLength(region.getChr()));
				}
				WindowScoreIterator<CountScore> wIter = data.scan(region, 1,0);
				TreeMap<Integer,Double> rtrn = new TreeMap<Integer,Double>();
				while (wIter.hasNext()){
					CountScore score = wIter.next();
					//logger.info("Checking score: " + score);
					Annotation a = score.getAnnotation();
					int pos = a.getStart();
					double count = score.getCount();
					if (normalize) {
						if (count>0) {
							rtrn.put(pos, (count/norm));
						}
					} else {
						if (count>0) {
							rtrn.put(pos, count);
						}
					}
				}
				return rtrn;
			}
		}
		
		// Count last read position only
		double norm = 1;
		if (normalize) {
			norm = data.getCount(region)/region.length();
		}
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
		if (normalize) {
			for (int key : rtrn.keySet()) {
				double val = rtrn.get(key);
				rtrn.put(key, val/norm);
			}
		}
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
	 * Get the one based wig format position for a zero based coordinate
	 * @param zeroBasedCoordinate The zero based coordinate
	 * @return The wig position
	 */
	public static int coordinateToWigPosition(int zeroBasedCoordinate) {
		return zeroBasedCoordinate + 1;
	}
	 	
	/**
	 * Get the zero based coordinate for a wig position
	 * @param wigPosition The wig position
	 * @return The zero based coordinate
	 */
	public static int wigPositionToCoordinate(int wigPosition) {
		return wigPosition - 1;
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
				int pos = coordinateToWigPosition(i.intValue());
				w.write(pos + "\t" + counts.get(i).toString() + "\n");
			}
		}
		w.close();
	}
	
	/**
	 * Collapse all overlapping genes within a collection into a non-overlapping set (considering strand)
	 * @param genes Collection of genes you wish to collapse
	 * @return Collection of collapsed genes
	 */
	private Collection<Gene> collapseGenes(Collection<Gene> genes) {
		TreeSet<Gene> rtrn = new TreeSet<Gene>();
		TreeSet<Gene> tmpGenes = new TreeSet<Gene>();
		tmpGenes.addAll(genes);
		IntervalTree<Gene> tree = makeTree(genes);
		Collection<Gene> skip = new TreeSet<Gene>();
		String name;
		
		for (Gene rtrnGene : genes) {
			name = rtrnGene.getName();
			
			if (tmpGenes.contains(rtrnGene)) {
				skip.add(rtrnGene);
				Iterator<Node<Gene>> overItr = tree.overlappers(rtrnGene.getStart(),rtrnGene.getEnd());
				
				while (overItr.hasNext()) {
					Gene g2 = overItr.next().getValue();
					if (rtrnGene.overlapsStranded(g2)) {
						rtrnGene = rtrnGene.takeUnion(g2);
						//logger.info("taking union...");
						skip.add(g2);
					}
				}
				

				rtrnGene.setName(name);
				//logger.info("finished collapsing gene into: " + rtrnGene);
				rtrn.add(rtrnGene);
				tmpGenes.removeAll(skip);
				Iterator<Gene> skipItr = skip.iterator();
				while (skipItr.hasNext()) {
					Gene rmv = skipItr.next();
					tree.remove(rmv.getStart(),rmv.getEnd(),rmv);
				}
				skip.clear();
			}
		}
		return rtrn;
	}
	
	/**
	 * Creates an interval tree using a collection of genes; speeds up the collapsing process
	 * @param genes Collection of genes used to make the tree
	 * @return Finished interval tree
	 */
	private IntervalTree<Gene> makeTree(Collection<Gene> genes) {
		IntervalTree<Gene> rtrn=new IntervalTree<Gene>();
		
		for(Gene gene: genes){
			int start=gene.getStart();
			int end=gene.getEnd();
			
			//if not already in then add
			Node<Gene> found=rtrn.find(start, end);
			if(found==null){
				rtrn.put(start, end, gene);
			}
			//else add isoform
			else{
				Gene ref=found.getValue();
				ref.addIsoform(gene);
				rtrn.put(start, end, ref);
			}
			
		}
		return rtrn;
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
		p.addBooleanArg("-n",  "Normalize position counts by average counts over region", false, false);

		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-c");
		String outFile = p.getStringArg("-o");

		int maxFragmentLength = p.getIntArg("-mf");
		boolean fragments = p.getBooleanArg("-pe");
		int maxGenomicSpan = p.getIntArg("-mg");
		boolean firstPositionOnly = p.getBooleanArg("-f");
		boolean normalize = p.getBooleanArg("-n");

		
		WigWriter ww = new WigWriter(bamFile, bedFile, chrSizeFile, firstPositionOnly, fragments, normalize);
		//ww.addReadFilter(new FragmentLengthFilter(ww.data.getCoordinateSpace(), maxFragmentLength));
		ww.addReadFilter(new GenomicSpanFilter(maxGenomicSpan));
		ww.writeWig(outFile);
		
	}

}

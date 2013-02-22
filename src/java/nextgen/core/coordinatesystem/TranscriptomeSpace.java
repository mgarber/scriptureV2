package nextgen.core.coordinatesystem;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import org.apache.log4j.Logger;
import nextgen.core.annotation.*;
import nextgen.core.feature.GeneWindow;
import nextgen.core.feature.Window;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScore;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.scripture.BuildScriptureCoordinateSpace;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.datastructures.Alignments;


public class TranscriptomeSpace implements CoordinateSpace{


	GeneTree geneTree;
	static Logger logger = Logger.getLogger(TranscriptomeSpace.class.getName());

	public TranscriptomeSpace(Map<String,Collection<Gene>> chrToGenesMap){
		this(null, chrToGenesMap,null);
	}
	
	/**
	 * Builds a transcriptome space from a bam file directly
	 * This will collapse all reads into an assembly of transcripts
	 * @param bamFile Alignment file
	 * @param genomeSeq A fasta file with the reference sequence 
	 */
	public TranscriptomeSpace(File bamFile,String genomeSeq){
		this(bamFile, new TreeMap<String, Collection<Gene>>(),genomeSeq);
	}
	
	/**
	 * Builds a transcriptome space from a bam file plus a known annotation file
	 * @param bamFile The reads to build the assembly from
	 * @param chrToGenesMap The known transcriptome to use
	 */
	public TranscriptomeSpace(File bamFile, Map<String, Collection<Gene>> chrToGenesMap,String genomeSeq){
		if(bamFile!=null){
			Map<String, Collection<Gene>> genes=assemble(bamFile,genomeSeq);
			merge(chrToGenesMap, genes);
		}
		
		this.geneTree=new GeneTree(chrToGenesMap);
	}
	
	/**
	 * This method will assemble reads into genes
	 * @param bamFile The reads to assemble
	 * @return
	 */
	private Map<String, Collection<Gene>> assemble(File bamFile,String genomeSeq) {
		BuildScriptureCoordinateSpace scriptureSpace = null;
		scriptureSpace = new BuildScriptureCoordinateSpace(bamFile,genomeSeq);
		return scriptureSpace.getPaths();
	}

	private void merge(Map<String, Collection<Gene>> chrToGenesMap,	Map<String, Collection<Gene>> genes) {
		for(String chr: genes.keySet()){
			Collection<Gene> mergedSet=genes.get(chr);
			if(chrToGenesMap.containsKey(chr)){
				Collection<Gene> set=chrToGenesMap.get(chr);
				mergedSet.addAll(set);
			}
			chrToGenesMap.put(chr, mergedSet);
		}
	}

	/**
	 * Return a fragment between start and end on chr 
	 */
	@Override
	public Collection<? extends Window> getFragment(String chr, int start, int end) {
		Collection<? extends Window> regions=this.geneTree.getRegion(new Alignments(chr, start,end));
		return regions;
	}

	public Collection<? extends Window> getFragment(Annotation annotation) {
		Collection<? extends Window> regions=this.geneTree.getRegion(annotation);
		return regions;
	}
	
	/**
	 * Returns an interval tree over the specified transcriptome region
	 */
	@Override
	public Collection<? extends Window> getOverlappingRegion(String chr,int start,int end) {
		Collection<? extends Window> fragments=this.getFragment(chr, start, end);
		return fragments;
	}

	/*
	 * Get a window iterator over a collection of GeneWindows
	 */
	@Override
	public Iterator<Window> getWindowIterator(Collection<Gene> baseGenes, int windowSize, int overlap) {
		return new WindowIterator(baseGenes, windowSize, overlap);
	}
	
	@Override
	public Iterator<Window> getWindowIterator(int windowSize, String chr, int start, int end, int overlap) {
		//Get regions
		Collection<GeneWindow> baseGenes=this.geneTree.getRegion(new Alignments(chr, start, end));
		
		//Pass to function to get iterator
		return new WindowIterator(baseGenes, windowSize, overlap);
		
		//Pass to function to get iterator
		//return getWindows(chr, start,end, windowSize, overlap).iterator();
	}
	
	@Override
	public Iterator<? extends Window> getWindowIterator(Annotation window, int windowSize, int overlap) {
		return getWindowIterator(windowSize, window.getChr(), window.getStart(), window.getEnd(), overlap);
	}

	/**
	 * Get a window iterator from start to finish on chromosome 
	 */
	@Override
	public Iterator<Window> getWindowIterator(String chr, int windowSize, int overlap) {
		//Get regions
		Collection<GeneWindow> baseGenes=this.geneTree.getRegion(chr);
		
		logger.info("here");
		
		//Pass to function to get iterator
		return new WindowIterator(baseGenes, windowSize, overlap);
	}

	@Override
	public Iterator<Window> getWindowIterator(int windowSize, int overlap) {
		//Get all genes
		Collection<GeneWindow> baseGenes=this.geneTree.getRegion();
		
		//Pass to function to get iterator
		return new WindowIterator(baseGenes, windowSize, overlap);
	}
	
	/**
	 * Scan windows over a Collection<GeneWindow>
	 * @param baseGenes The collection of regions to scan
	 * @param windowSize Window size
	 * @param overlap Overlap size
	 * @param processor The window processor
	 * @return A score iterator over windows
	 */
	public <T extends WindowScore> WindowScoreIterator<T> scan(Collection<Gene> baseGenes, int windowSize, int overlap, WindowProcessor<T> processor) {
		Gene gene = new Gene(baseGenes);
		Iterator<? extends Window> windowIterator = getWindowIterator(baseGenes, windowSize, overlap);
		return new WindowScoreIterator<T>(windowIterator, processor, gene);
	}
	
	/**
	 * Get list of position level counts in mature transcript
	 * @param region The region
	 * @param data Alignment data
	 * @return List of position counts
	 * @throws IOException 
	 */
	public List<Double> getPositionCountList(Gene region, AlignmentModel data) throws IOException {
		List<Double> rtrn = new ArrayList<Double>();
		Collection<Gene> baseGenes = new ArrayList<Gene>();
		baseGenes.add(region);
		WindowProcessor<CountScore> processor = new CountScore.Processor(data);
		WindowScoreIterator<CountScore> scoreIter = scan(baseGenes, 1, 0, processor);
		while(scoreIter.hasNext()) {
			rtrn.add(Double.valueOf(scoreIter.next().getCount()));
		}
		return rtrn;
	}
	
	/**
	 * Get map of position to of position level counts in mature transcript
	 * @param region The region
	 * @param data Alignment data
	 * @return Map of position to count
	 * @throws IOException 
	 */
	public TreeMap<Integer,Double> getPositionCountMap(Gene region, AlignmentModel data) throws IOException {
		TreeMap<Integer,Double> rtrn = new TreeMap<Integer,Double>();
		Collection<Gene> baseGenes = new ArrayList<Gene>();
		baseGenes.add(region);
		WindowProcessor<CountScore> processor = new CountScore.Processor(data);
		WindowScoreIterator<CountScore> scoreIter = scan(baseGenes, 1, 0, processor);
		int processed = 0;
		while(scoreIter.hasNext()) {
			rtrn.put(Integer.valueOf(processed), Double.valueOf(scoreIter.next().getCount()));
			processed++;
		}
		return rtrn;
	}



	private class WindowIterator implements Iterator<Window> {

		int windowSize;
		int step;
		Iterator<? extends Gene> genes;
		Iterator<? extends Window> currentGeneWindows;
		
		WindowIterator(Collection<? extends Gene> genes, int windowSize, int overlap){
			this.windowSize=windowSize;
			this.step=windowSize-overlap;
			this.genes=genes.iterator();
		}
		
		@Override
		public boolean hasNext() {
			//first, see if there are windows in the current gene
			if(currentGeneWindows!=null && currentGeneWindows.hasNext()){
				return true;
			}
			//else, see if there are still genes
			else if(genes.hasNext()){
				return true;
			}
			return false;
		}

		@Override
		public Window next() {
			//if there are windows still in currentGeneWindows send these
			if(currentGeneWindows!=null && currentGeneWindows.hasNext()){
				return currentGeneWindows.next();
			}
			else{
				//else make a new currentGeneWindows
				//make new currentGeneWindows
				Gene nextGene=genes.next();
				//logger.info(nextGene);
				this.currentGeneWindows=makeWindows(nextGene, this.windowSize, this.step);
				return next();
			}
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("TODO");
			
		}
		
		private Iterator<? extends Window> makeWindows(Gene gene, int windowSize, int step){
			Collection<? extends Window> rtrn=gene.getWindows(windowSize, step, 0);
			return rtrn.iterator();
			
			/*Collection<Window> temp=new TreeSet<Window>();
			//for each gene trim from relative 0 to relative end in increments of windowSize
			for(int windowStart=0; windowStart<(gene.getSize()+1)-(windowSize); windowStart+=step){
				int windowEnd=windowStart+windowSize;
				GeneWindow trimmed=gene.trimGene(windowStart, windowEnd);
				logger.info("iterating "+trimmed.getChr()+":"+trimmed.getStart()+"-"+trimmed.getEnd());
				//if(trimmed!=null && trimmed.getSize()==windowSize){
					temp.add(trimmed);
				//}
			}
			return temp.iterator();*/
		}
		
	}
	
	
	/**
	 * Returns an iterator over the subset of genes in the coordinate space
	 * @param genes
	 * @return
	 */
	Iterator<Collection<Window>> getWindowIterator(Collection<Gene> genes) {
		return null;
	}

	
protected class GeneTree {
		
		Map<String, IntervalTree<Gene>> tree;
		Map<String, Annotation> metaAnnotation; //This is the collapsed exons of the each chromosome	
		Map<String, Gene> genesByName;
		
		GeneTree(Map<String, Collection<Gene>> geneMap){
			tree=new TreeMap<String, IntervalTree<Gene>>();
			metaAnnotation=new TreeMap<String, Annotation>();
			genesByName = new TreeMap<String, Gene>();
			
			for(String chr: geneMap.keySet()){
				Collection<Gene> genes=geneMap.get(chr);
				for(Gene gene : genes) {
					genesByName.put(gene.getName(),gene);
				}
				IntervalTree<Gene> t=makeTree(genes);
				tree.put(chr, t);
				Annotation meta=makeMetaAnnotation(genes);
				metaAnnotation.put(chr, meta);
			}
			
		}
		
		Collection<String> getGeneNames() {
			return genesByName.keySet();
		}
		
		Map<String, Gene> getGenesByName() {
			return genesByName;
		}
		
		

		private Annotation makeMetaAnnotation(Collection<Gene> genes) {
			BasicAnnotation rtrn=null;
			
			for(Gene gene: genes){
				if(rtrn==null){
					rtrn=new BasicAnnotation(gene.getBlocks());
				}
				else{
					rtrn.addBlocks(gene.getBlocks());}
			}
			
			return rtrn;
		}



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
		 * Returns the trimmed region contained in this interval
		 * @param chr
		 * @param start
		 * @param end
		 * @return
		 */
		public Collection<GeneWindow> getRegion(Annotation region) {
			if(!tree.containsKey(region.getChr())){
				logger.error(region.getChr()+" is not in the genomic space");
				return null;
			}
			
			Collection<GeneWindow> rtrn=new TreeSet<GeneWindow>();
						
			//Get overlapping regions.
			// Note: the blocks in the overlappers do not necessarily intersect with the region
			Iterator<Node<Gene>> iter=tree.get(region.getChr()).overlappers(region.getStart(), region.getEnd());
						
			//iterate through and trim ends
			while(iter.hasNext()){
				Gene g=iter.next().getValue();
				Collection<Gene> iso=g.getIsoforms();
				for(Gene gene: iso){
					//if start, end exceeds the gene length then skip
					//if(start>=gene.getStart() && end<= gene.getEnd()){
						GeneWindow trimmed=gene.trimAbsolute(region.getStart(), region.getEnd());
						if(trimmed!=null && trimmed.size() > 0){
							if(trimmed.overlaps(region)){
								rtrn.add(trimmed);
							}
						}
					//}
					 //else{logger.error(g.getName()+" was null after trimming isoform " + gene.getName() + ":" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + " to absolute coordinates " + g.getChr() + ":" + region.getStart() + "-" + region.getEnd());}
				}
			}

			return rtrn;
		}
		
		/**
		 * Returns all genes that have chr
		 * @param chr
		 * @return
		 */
		public Collection<GeneWindow> getRegion(String chr) {
			if(!tree.containsKey(chr)){
				logger.error(chr+" is not in the genomic space");
				return null;
			}
			
			Collection<GeneWindow> rtrn=new TreeSet<GeneWindow>();
			//Get overlapping regions
			Iterator<Node<Gene>> iter=tree.get(chr).iterator();
			//iterate through and trim ends
			while(iter.hasNext()){
				Gene g=iter.next().getValue();
				Collection<Gene> iso=g.getIsoforms();
				for(Gene gene: iso){
					rtrn.add(new GeneWindow(gene));
				}
			}
			
			//return collection
			return rtrn;
		}
		
		/**
		 * Returns all genes
		 * @return
		 */
		public Collection<GeneWindow> getRegion() {
			Collection<GeneWindow> rtrn=new TreeSet<GeneWindow>();
			
			for(String chr: this.tree.keySet()){
				rtrn.addAll(getRegion(chr));
			}
			
			return rtrn;
		}



		public long getLength() {
			long rtrn=0;
			for(String chr: this.metaAnnotation.keySet()){
				rtrn+=metaAnnotation.get(chr).getSize();
			}
			return rtrn;
		}
		
		public long getLength(String chr) {
			return this.metaAnnotation.get(chr).getSize();
		}
	}


	@Override
	/**
	 * This will sum up the blocks ignoring the overlaps
	 */
	public long getLength() {
		return geneTree.getLength();
	}

	@Override
	public long getLength(String chr) {
		return geneTree.getLength(chr);
	}

	@Override
	public void permuteAnnotation(Annotation a) {
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public void permuteAnnotation(Annotation a, Annotation bounds) {
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public Collection<String> getReferenceNames() {
		return geneTree.getGeneNames();
	}
	
	@Override
	public Collection<? extends Annotation> getReferenceAnnotations() {
		return geneTree.getGenesByName().values();
	}
	
	@Override
	public Annotation getReferenceAnnotation(String name) {
		//System.out.println("Gene tree contains "+name+" : "+geneTree.getGenesByName().containsKey(name));

		if(geneTree.getGenesByName().containsKey(name)){
			return geneTree.getGenesByName().get(name);
		}
		else{
			return null;
		}
	}
	
}

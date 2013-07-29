package nextgen.core.scripture;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.jgrapht.GraphPath;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.Pair;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CollapseByIntersection;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.Path;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.JCSAlignmentModel;
import nextgen.core.readFilters.CanonicalSpliceFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.IndelFilter;
import nextgen.core.readFilters.PairedAndProperFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.readFilters.SameOrientationFilter;
import nextgen.core.readFilters.SplicedReadFilter;
import nextgen.core.scripture.OrientedChromosomeTranscriptGraph.TranscriptGraphEdge;
import nextgen.core.scripture.statistics.ConnectDisconnectedTranscripts;

public class BuildScriptureCoordinateSpace {

	private static final TranscriptionRead DEFAULT_TXN_READ =  TranscriptionRead.UNSTRANDED;
	static Logger logger = Logger.getLogger(BuildScriptureCoordinateSpace.class.getName());
	private JCSAlignmentModel model;
	String genomeSeq = null;
	int windowSize=20000000;
	private CoordinateSpace space;
	private Map<String, ChromosomeTranscriptGraph> graphs;
	private static boolean forceStrandSpecificity=true; //TODO This should be passed or at least determined from data
	private static double DEFAULT_MIN_COV_THRESHOLD = 0.2;
	private static double MIN_SPLICE_PERCENT = 0.05;
	private double coveragePercentThreshold = DEFAULT_MIN_COV_THRESHOLD;
	//private static TranscriptionRead DEFAULT_TXN_READ = TranscriptionRead.UNSTRANDED;
	String outName = null;
	private double DEFAULT_ALPHA = 0.01;
	private double alpha=DEFAULT_ALPHA;
	private static double MIN_SPLICE_READS = 3.0;
	private double THRESHOLD_SPURIOUS = 0.95;
	private double minSpliceReads = MIN_SPLICE_READS;
	private double minSplicePercent = MIN_SPLICE_PERCENT;
	int counter = 1000;
	int globalCounter = 1000;
	File bamFileName;
	private int constant = 10000;
	private double globalPairedLambda=0.0;
	File bamfile;
	//double globalFragments;
	
/*	public BuildScriptureCoordinateSpace(File bamFile){
		this(bamFile,DEFAULT_MIN_COV_THRESHOLD,null,null,true,DEFAULT_TXN_READ);
		logger.info("Genome sequence has not been provided");
	}
	
	
	public BuildScriptureCoordinateSpace(File bamFile,double threshold,String genomeDir,String outputName){
		//By default first read is transcription read
		this(bamFile,threshold,genomeDir,outputName,true,DEFAULT_TXN_READ);
	}
	*/
	
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir){
		this(bamFile,genomeDir,bamFile.getName()+".reconstructions",true,DEFAULT_TXN_READ,null);
	}
	/**
	 * 
	 * @param bamFile
	 * @param threshold
	 * @param genomeDir
	 * @param outputName
	 * @param forceStrandedness
	 */
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead strand,ArgumentMap argMap){
			
		bamfile=bamFile;
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		bamFileName = bamFile;
		forceStrandSpecificity = forceStrandedness;
		outName = outputName;
		model=new JCSAlignmentModel(bamfile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		space=model.getCoordinateSpace();

		setThresholds(argMap);
//		globalFragments = calculateGlobalFragments();
		logger.info("Parameters used: " +
				"\nIntron intention Filter: "+THRESHOLD_SPURIOUS+
				"\nPremature Assembly Filter: "+coveragePercentThreshold+
				"\nSplice junction Filter : "+
				"\n\tNumber of spliced reads : "+minSpliceReads+
				"\n\tPercentage of total spliced reads: "+minSplicePercent+
				"\nAlpha for single exon assemblies : "+alpha
				);
		
		assemble(strand);
		
		Map<String, Collection<Gene>> rtrn=getPaths();
/*		try {
			FileWriter writer=new FileWriter(outName+".08graph.all.paths.bed");
			for(String chr: rtrn.keySet()){
				for(Gene g:rtrn.get(chr)){
					writer.write(g.toBED()+"\n");
				}
			}			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}*/
		
		try {
			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			Map<String, Collection<Gene>> genes = setFPKMScores(rtrn,writer,writer2);
			writer.close();
			writer2.close();

//			write(outName+".10final.paths.bed",genes);
			
			postProcess(genes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	} 
	
	/**
	 * 
	 * @param bamFile
	 * @param threshold
	 * @param genomeDir
	 * @param outputName
	 * @param forceStrandedness
	 * @param chr
	 */
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead strand,String chr,ArgumentMap argMap){
			
		bamfile=bamFile;
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		bamFileName = bamFile;
		forceStrandSpecificity = forceStrandedness;

		outName = outputName;
		
		setThresholds(argMap);
		//assemble(strand);
		logger.info("Parameters used: " +
				"\nIntron Retention Filter: "+THRESHOLD_SPURIOUS+
				"\nPremature Assembly Filter: "+coveragePercentThreshold+
				"\nSplice junction Filter : "+
				"\n\tNumber of spliced reads : "+minSpliceReads+
				"\n\tPercentage of total spliced reads: "+minSplicePercent+
				"\nAlpha for single exon assemblies : "+alpha
				);
		
		ChromosomeTranscriptGraph graph=assemble(chr,strand);
		graphs.put(chr, graph);
		
		Map<String, Collection<Gene>> rtrn=getPaths();
/*		try {
			FileWriter writer=new FileWriter(outName+".08graph.all.paths.bed");
			for(Gene g:rtrn.get(chr)){
				writer.write(g.toBED()+"\n");
			}			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
*/		
		try {
			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			Map<String, Collection<Gene>> genes = setFPKMScores(rtrn,writer,writer2);
			writer.close();
			writer2.close();

//			write(outName+".10final.paths.bed",genes);
			
			postProcess(genes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	} 
	
	private void setThresholds(ArgumentMap argMap){ 
		
		if(argMap!=null){
			coveragePercentThreshold = argMap.getDouble("coverage", DEFAULT_MIN_COV_THRESHOLD);
			alpha = argMap.getDouble("alpha", DEFAULT_ALPHA);
			minSpliceReads = argMap.getDouble("minSpliceReads", MIN_SPLICE_READS);
			minSplicePercent = argMap.getDouble("percentSpliceReads", MIN_SPLICE_PERCENT);
		}
	}

	private void assemble(TranscriptionRead strand) {
		//Iterate over all chromosomes
		for(String chr: space.getReferenceNames()){
			logger.info("Reference name: "+chr);
			if(model.getRefSequenceLambda(chr)==0.0){
				logger.info(chr+" is not expressed in the alignment file");
			}
			else{
				ChromosomeTranscriptGraph graph=assemble(chr,strand);
				this.graphs.put(chr, graph);
			}
		}
	}

	private void postProcess(Map<String, Collection<Gene>> oldGenes) throws IOException{
		
		Map<String, Collection<Gene>> newGenes = new HashMap<String, Collection<Gene>>();
		for(String chr:oldGenes.keySet()){
			//For all paths on this chromosome
			
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:oldGenes.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			/*
			 * MERGE THE GENES
			 */
			//Step 1: Iterate through the working assembly and find paths that overlap
			//if overlaps but is incompatible, see if we can branch
			//iterate
			Iterator<Gene> iter=tree.toCollection().iterator();			

			Collection<Gene> considered = new HashSet<Gene>();
	 		logger.info("Enter merging");
			while(iter.hasNext()){
				Gene branch1=iter.next();
				considered.add(branch1);
				//get overlapping branches
				Iterator<Gene> overlappers=tree.overlappingValueIterator(branch1.getStart(), branch1.getEnd());
				//Collection<Assembly> toRemove=new TreeSet<Assembly>();
				while(overlappers.hasNext()){
					Gene branch2=overlappers.next();
					if(!considered.contains(branch2)){
						if(!branch1.equals(branch2) && compatible(branch1, branch2)){
							logger.debug("Merging: "+ branch1.getName()+ " and "+branch2.getName());
							Collection<Annotation> rtrn=new TreeSet<Annotation>();
							rtrn.addAll(branch1.getBlocks());
							rtrn.addAll(branch2.getBlocks());
							Gene merged=new Gene(rtrn);
							merged.setName(branch1.getName());
							//SET THE SCORE TO THE MAX FPKM
							merged.setBedScore(Math.max(branch1.getBedScore(),branch2.getBedScore()));
							//remove annotation1 and annotation2
							tree.remove(branch1.getStart(), branch1.getEnd(), branch1);
							tree.remove(branch2.getStart(), branch2.getEnd(), branch2);
							//add merged
							tree.put(merged.getStart(), merged.getEnd(), merged);
						}
					}
				}
			}
			
			/*
			 * FPKM THRESHOLD
			 */
			//MAKE A MAP OF GENE TO ISOFORMS
			Map<Gene,Set<Gene>> isoformMap = getIsoformMap(tree.toCollection());		
			Set<Gene> finalSet = new HashSet<Gene>();
			//For each gene
			for(Gene gene:isoformMap.keySet()){
				//For each transcript
				double max = Double.MIN_VALUE;
				for(Gene isoform:isoformMap.get(gene)){					
					double score = isoform.getBedScore();
					if(score>max){
						max = score;
					}
				}
				//CHeck if the FPKM of each is at least 10% of the mac FPKM
				for(Gene isoform:isoformMap.get(gene)){		
					double score = isoform.getBedScore();
					if(score>0.1*max && score>0.0){
						//filtered.remove(isoform.getStart(), isoform.getEnd(), isoform);
						finalSet.add(isoform);
					}
				}
			}
			newGenes.put(chr, finalSet);
		}
	
		//write(outName+".11postprocessed.paths.bed",newGenes);
		write(outName+".scripture.paths.bed",newGenes);
		connectDisconnectedTranscripts(newGenes);
	}
	
	
	private ChromosomeTranscriptGraph assemble(String chr,TranscriptionRead strand) {
		
		//Option 1: Scan through the space and collapse into compatible edges and nodes
		ChromosomeTranscriptGraph graph=assembleDirectly(chr,strand);
		
		return graph;
	}

	/**
	 * This function connects disconnected reconstructions using paired end reads.
	 * @param annotations
	 * @throws IOException 
	 */
	private void connectDisconnectedTranscripts(Map<String, Collection<Gene>> annotations) throws IOException{
		
		model.addFilter(new PairedAndProperFilter());
		int loop=0;
		boolean somethingWasConnected = true;
		double medianInsertSize=0.0; 
		Map<String,Collection<Gene>> temp = new HashMap<String,Collection<Gene>>();
		for(String chr:annotations.keySet()){
			if(!annotations.get(chr).isEmpty()){
				temp.put(chr, annotations.get(chr));
			}
		}
		medianInsertSize += model.getReadSizeDistribution(new TranscriptomeSpace(temp), 800, 100).getMedianOfAllDataValues();
		//double medianInsertSize = 600;
		logger.info("Median size = "+medianInsertSize);

		Map<String,Collection<Gene>> conn = null;
		
		while(somethingWasConnected || loop<10){
			somethingWasConnected =false;
			loop++;
			
			logger.debug("Connected disconnected transcripts: Loop "+loop);
			conn = new HashMap<String,Collection<Gene>>();
			for(String chr:annotations.keySet()){
				conn.put(chr, new TreeSet<Gene>());			
			}
			
			for(String chr:annotations.keySet()){
				//For all genes on this chromosome
				logger.debug("Connecting, Processing "+chr);
				Collection<Gene> newGenes = new TreeSet<Gene>();
				//MAKE AN INTERVAL TREE OF THE GENES on this chr
				IntervalTree<Gene> tree = new IntervalTree<Gene>();
				for(Gene g:annotations.get(chr)){
					conn.get(chr).add(g);
					tree.put(g.getStart(), g.getEnd(), g);
				}
				//For each transcript
				//Iterate over all reconstructions
				Iterator<Gene> iter=tree.toCollection().iterator();
				while(iter.hasNext()){
					Gene gene=iter.next();
					//For all assemblies downstream of this assembly in 10kB regions
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getEnd(), gene.getEnd()+constant);
					
					newGenes.add(gene);
					while(overlappers.hasNext()){
						Gene other = overlappers.next();
						if(isCandidate(gene,other)){
							if(pairedEndReadSpansTranscripts(gene, other)){ 
								//if(secondTranscriptIsSingleExon(gene,other)){
									logger.debug("Attempt to connect "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC());
									
									//Connect the genes
									Annotation connected = getConnectedTranscript(gene,other,medianInsertSize);
									if(connected!=null){
										somethingWasConnected = true;
										Gene newConnected = new Gene(connected);
										double[] scores = getScores(newConnected);
										double[] fields = new double[4];
										newConnected.setName(gene.getName()+"_"+other.getName());
										//[0] : sum
										fields[0] = scores[0];
										//[1] : p-value
										fields[1] = scores[1];
										//[2] : FPK
										fields[2] = (scores[0]*1000.0)/newConnected.getSize();
										//[3] : FPKM
										//Calculate FPKM
										fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
										//logger.debug("For isoform : "+isoform.getName()+"\tNum of exons: "+isoform.getSpliceConnections().size()+"\t"+fields[0]+"\t"+fields[1]);
										newConnected.setBedScore(fields[3]);
										newConnected.setExtraFields(fields);
										newGenes.remove(gene);
										newGenes.add(newConnected);
										conn.get(chr).remove(gene);
										conn.get(chr).remove(other);
										conn.get(chr).add(newConnected);
									}
								//}
							}
							else{
								//logger.info("The genes "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC()+" do not have paired reads");
							}
						}
					}				
				}
			}
			annotations = conn;
		}
		
		//FileWriter bw = new FileWriter(outName+".12connected.bed");
		FileWriter bw = new FileWriter(outName+".connected.bed");
		for(String name:conn.keySet()){
			Iterator<Gene> ter = conn.get(name).iterator();
			while(ter.hasNext()){
				Gene isoform = ter.next();
				Gene iso = new Gene(trimEnds(isoform,0.1));
				bw.write(iso.toBED()+"\n");
			}
		}
		bw.close();
	}
	
	/**
	 * Connect the two genes by fusing the last exon of the first and the first exon of the last gene
	 * @param gene
	 * @param other
	 * @return
	 * @throws IOException 
	 */
	private Annotation getConnectedTranscript(Gene gene,Gene other,double medianInsertSize) throws IOException{
		
		Pair<Gene> orderedGenes = ConnectDisconnectedTranscripts.getOrderedAssembly(gene, other);
		Annotation connected = null;
		/*
		 * If the distance between the transcripts is less than the insert size,
		 * paste through
		 */
		if(orderedGenes.getValue2().getStart()-orderedGenes.getValue1().getEnd() < medianInsertSize){
			Gene firstGene = orderedGenes.getValue1().copy();
			firstGene.setEnd(orderedGenes.getValue2().getStart());
			connected = firstGene.union(orderedGenes.getValue2());
		}
		/*
		 * If the distance between the transcripts is more than the insert size,
		 * 		if there is at least 1 splice read spanning the junction,
		 * 			union it
		 */
		else{
			boolean flag= false;
			Annotation junction = orderedGenes.getValue1().getLastExon().union(orderedGenes.getValue2().getFirstExon());
			CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(junction,false), new SplicedReadFilter());
			while(splicedIter.hasNext()){
				Alignment read = splicedIter.next();
				//if there is at least 1 splice read spanning the junction,
				if(BuildScriptureCoordinateSpace.compatible(read,junction)){
					flag=true;
					break;
				}
			}
			splicedIter.close();
			if(flag){
				connected = gene.union(other);				
			}
			/*
			 * Else NO CONNECTION
			 */
		}
		return connected;		
	}
	
	/**
	 * This function returns true if
	 * 		1. gene and other do not overlap
	 * 	    2. gene and other are in the same orientation
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean isCandidate(Gene gene,Gene other){ 
		
		if(gene.overlaps(other)){
			//logger.info(gene.getName()+" overlaps "+gene.toUCSC());
			return false;
		}
		else{
			//if transcript is in an intron of the other transcript
			if(gene.getEnd()>other.getStart() && other.getEnd()>gene.getStart()){
				//logger.info("One is inside the intron of the other");
				return false;
			}
			if(gene.getOrientation().equals(other.getOrientation())){
				//logger.debug("The orientation is same");
				return true;
			}
		}
		//logger.debug("Something else");
		return false;
	}
	
	/**
	 * This function will return true if gene and other have at least 1 paired end read in common
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean pairedEndReadSpansTranscripts(Gene gene,Gene other){
		
		boolean rtrn = false;
		//Get all overlapping paired end reads in same orientation as the gene
		CloseableFilterIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(new BasicAnnotation(gene.getChr(),gene.getStart(),other.getEnd()),false), new SameOrientationFilter(gene));
		while(iter.hasNext()){
			Alignment read = iter.next();
			List<Annotation> mates = (List<Annotation>) read.getReadAlignments(model.getCoordinateSpace());
			if(mates.get(0).getOrientation().equals(gene.getOrientation()) && mates.get(1).getOrientation().equals(gene.getOrientation())){
				if((gene.overlaps(mates.get(0)) && other.overlaps(mates.get(1)))
						||
						gene.overlaps(mates.get(1)) && other.overlaps(mates.get(0))){
					rtrn = true;
					break;
				}
			}
			
		}
		iter.close();
		return rtrn;
	}
	
	/**
	 * Returns true if the transcript at the oriented 3' end is a single exon transcript.
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean secondTranscriptIsSingleExon(Gene gene,Gene other){
		
		Pair<Gene> orderedGenes = ConnectDisconnectedTranscripts.getOrderedAssembly(gene, other);
		if(gene.isNegativeStrand()){
			if(orderedGenes.getValue1().getBlocks().size()==1){
				return true;
			}
		}
		else{
			if(orderedGenes.getValue2().getBlocks().size()==1){
				return true;
			}
		}
		return false;
	}
	
	public static void write(String save, Map<String, Collection<Gene>> rtrn) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String chr: rtrn.keySet()){
			Collection<Gene> genes=rtrn.get(chr);
			for(Gene gene: genes){
				writer.write(gene.toBED()+"\n");
			}
		}
		
		writer.close();
	}
	
	public static void memoryStats(){
		int mb = 1024*1024;
        
        //Getting the runtime reference from system
        Runtime runtime = Runtime.getRuntime();
         
        logger.info("##### Heap utilization statistics [MB] #####");
         
        //Print used memory
        logger.info("Used Memory:"
            + (runtime.totalMemory() - runtime.freeMemory()) / mb);
 
        //Print free memory
        logger.info("Free Memory:"
            + runtime.freeMemory() / mb);
         
        //Print total available memory
        logger.info("Total Memory:" + runtime.totalMemory() / mb);
 
        //Print Maximum available memory
        logger.info("Max Memory:" + runtime.maxMemory() / mb);
	}
	
	/**
	 * Makes assemblies using all reads and subjects all assemblies to a list of filters.
	 * @param chr
	 * @param strand
	 * @return
	 */
	private ChromosomeTranscriptGraph assembleDirectly(String chr,TranscriptionRead strand){
		
		model=new JCSAlignmentModel(bamfile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		model.addFilter(new ProperPairFilter());
		model.addFilter(new IndelFilter());
		model.addFilter(new GenomicSpanFilter(20000000));
		this.space=model.getCoordinateSpace();

		long S = System.currentTimeMillis();	
		logger.info("Assembling spliced reads");
		long start = System.currentTimeMillis();		
		//SPLICED READS
		CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new CanonicalSpliceFilter(genomeSeq));
		IntervalTree<Assembly> splicedAssemblies=assembleDirectly(splicedIter,strand);
		long end = System.currentTimeMillis();
		logger.debug("TIME: ASSEMBLE SPLICED: "+(end-start));
		//try{write(splicedAssemblies, outName+"."+chr+"."+"01splicedAssemblies.bed");}catch(IOException ex){}		
		logger.info("Size of spliced assemblies: "+splicedAssemblies.size());
		
		//NON_SPLICED READS
/*		CloseableIterator<Alignment> iter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new NoSpliceFilter());
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, splicedAssemblies,strand);
		end = System.currentTimeMillis();*/
		
		// ALL READS
		CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, splicedAssemblies,strand);
		end = System.currentTimeMillis();
		
		/*
		 * ONLY ALL READS
		 */ 
/*		long start = System.currentTimeMillis();
		logger.info("Assembling all reads");
		CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, strand);
		long end = System.currentTimeMillis();*/
		logger.debug("TIME: ASSEMBLE NON SPLICED: "+(end-start));
		logger.info("Size of direct assemblies: "+workingAssemblies.size());		
		//try{write(workingAssemblies, outName+"."+chr+"."+"02directAssemblies.bed");}catch(IOException ex){}
				
/*		model=new JCSAlignmentModel(bamFileName.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		model.addFilter(new ProperPairFilter());
		model.addFilter(new IndelFilter());
		model.addFilter(new GenomicSpanFilter(20000000));*/
		globalPairedLambda = model.getGlobalPairedFragments()/model.getGlobalLength();
		logger.info("Intron retention filter");
		//REMOVE SPURIOUS
		start = System.currentTimeMillis();
		IntervalTree<Assembly> unspuriousAssemblies = intronRetentionFilter(workingAssemblies,chr);
		end = System.currentTimeMillis();
		logger.debug("TIME: REMOVE SPURIOUS: "+(end-start));
		//try{write(unspuriousAssemblies, outName+"."+chr+"."+"03intronRetentionAssemblies.bed");}catch(IOException ex){}
									
		logger.info("Merge assemblies");
		start = System.currentTimeMillis();
		mergeAssembly(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: MERGE: "+(end-start));
		//try{write(unspuriousAssemblies, outName+"."+chr+"."+"04mergedAssemblies.bed");}catch(IOException ex){}

		logger.info("Extend compatible assemblies");
		start = System.currentTimeMillis();
		extendAssembly(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: BRANCHING: "+(end-start));
		//try{write(unspuriousAssemblies, outName+"."+chr+"."+"05extendedAssemblies.bed");}catch(IOException ex){}
				
		logger.info("Remove premature assemblies");
		start = System.currentTimeMillis();
		// Flag and remove premature
		IntervalTree<Assembly> filteredAssemblies = removePrematureAssemblies(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: FILTER PREMATURE: "+(end-start));
		//IntervalTree<Assembly> filteredAssemblies = removePrematureAssembliesUsingPairedEnds(unspuriousAssemblies);
		//try{write(filteredAssemblies, outName+"."+chr+"."+"06filteredAssemblies.bed");}catch(IOException ex){}

		start = System.currentTimeMillis();
		IntervalTree<Assembly> highSplicedAssemblies = removeLowSpliceJunctionAssemblies(filteredAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: REMOVE LOW SPLICED JUNCTIONS: "+(end-start));
		//try{write(highSplicedAssemblies, outName+"."+chr+"."+"07highSplicedAssemblies.bed");}catch(IOException ex){}
		
		//Lets split potential preprocessed and mature transcripts and test whether to include certain nodes
		
		start = System.currentTimeMillis();
		//Make graph
		ChromosomeTranscriptGraph graph=makeGraph(highSplicedAssemblies, chr); //TODO Should use the merged set
		end = System.currentTimeMillis();
		logger.debug("TIME: MAKE GRAPH: "+(end-start));
		long E = System.currentTimeMillis();	
		logger.debug("TIME: TOTAL "+(E-S));
		return graph;
	}

	private void extendAssembly(IntervalTree<Assembly> tree) {
		//We have a set of assemblies that are all incompatible
		//We want to link up parts
		logger.debug("Enter extend assembly");
		//iterate through and get overlapping assemblies
		Iterator<Assembly> iter=tree.toCollection().iterator();
		Collection<Assembly> considered = new HashSet<Assembly>();
		
		//IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();
		
		while(iter.hasNext()){
			//try and merge the non-overlapping portions
			Assembly assembly1=iter.next();
			//currentAssemblies.put(assembly1.getStart(), assembly1.getEnd(), assembly1);
			//boolean a1Changed = false;
			considered.add(assembly1);
			Iterator<Assembly> overlappers=tree.overlappingValueIterator(assembly1.getStart(), assembly1.getEnd());
			while(overlappers.hasNext()){
				Assembly assembly2=overlappers.next();
				//logger.debug("Assembly1: "+assembly1.getName()+" Assembly2: "+assembly2.getName());
				//logger.debug("overlaps "+assembly2.toUCSC());
				if(!considered.contains(assembly2)){
					//try and merge the non-overlapping portions
					Collection<Assembly> merged=branchAssemblies(assembly1, assembly2);
					if(!merged.isEmpty()){
						logger.debug("Branching "+assembly1.getName()+" and "+assembly2.getName());
						
						//remove assembly1
						tree.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//remove assembly2
						tree.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);						
						considered.remove(assembly1);
						considered.remove(assembly2);
						//currentAssemblies.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//currentAssemblies.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);
						//add merged
						//logger.debug("No. of assemblies "+ merged.size());
						for(Assembly merge: merged){
							tree.put(merge.getStart(), merge.getEnd(), merge);
							considered.add(merge);
							//currentAssemblies.put(merge.getStart(), merge.getEnd(), merge);
						}
					}
				}
			}
		}
		//return currentAssemblies;
	}

	private Collection<Assembly> branchAssemblies(Assembly assembly1, Assembly assembly2) {
		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		
		//order the two assemblies by which starts first
		Pair<Assembly> orderedAssembly=getOrderedAssembly(assembly1, assembly2);
		
		//go from first.getFirstExon() to second.getFirstExon
		Annotation portionToConsiderAdding=orderedAssembly.getValue1().intersect(new Alignments(assembly1.getChr(), orderedAssembly.getValue1().getStart(), orderedAssembly.getValue2().getBlocks().iterator().next().getEnd()));
//		Annotation portionToTest = portionToConsiderAdding.minus(orderedAssembly.getValue2());
		//ONLY IF THE REGION OF THE PORTIONTOCONSIDERADDING THAT DOES NOT OVERLAP THE SECOND ASSEMBLY IS NOT SPLICED
//		if(portionToTest.getSpliceConnections().isEmpty()){
			//if this is compatible with second assemebly
			if(compatible(orderedAssembly.getValue2(), portionToConsiderAdding)){
				Assembly merged=merge(portionToConsiderAdding, orderedAssembly.getValue2());
				merged.setName(orderedAssembly.getValue2().getName());
				rtrn.add(merged);
				rtrn.add(orderedAssembly.getValue1());
				//System.out.println(merged);
			}
//		}
		//get first exon of later start site exon
		//trim assembly 2 from start till this point
		//ask if compatible
		//if so, merge trim with rest of later assembly
		return rtrn;
	}

	private Pair<Assembly> getOrderedAssembly(Assembly assembly1, Assembly assembly2) {
		Pair<Assembly> rtrn=new Pair<Assembly>();
		//Order by CompareTo
		if(assembly1.compareTo(assembly2)<0){
			rtrn.setValue1(assembly1);
			rtrn.setValue2(assembly2);
		}
		else{
			rtrn.setValue1(assembly2);
			rtrn.setValue2(assembly1);
		}
		return rtrn;
	}

	private void write(IntervalTree<Assembly> tree, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Assembly> set=tree.toCollection();
		
		for(Assembly align: set){
			writer.write(align+"\n");
		}
		
		writer.close();
	}
		
	private IntervalTree<Assembly> assembleDirectly(CloseableIterator<Alignment> iter,TranscriptionRead strand) {
		IntervalTree<Assembly> workingAssemblies=new IntervalTree<Assembly>();
		return assembleDirectly(iter, workingAssemblies,strand);
	}


	/**
	 * Makes a graph using the specified assemblies for the specified chrosomome
	 * This function will go through each gene and only add those that pass the p-value threshold
	 * @param workingAssemblies
	 * @param chr
	 * @return
	 */
	private ChromosomeTranscriptGraph makeGraph(IntervalTree<Assembly> workingAssemblies, String chr) {
		ChromosomeTranscriptGraph graph=new ChromosomeTranscriptGraph(chr);
		Iterator<Node<Assembly>> iter=workingAssemblies.iterator();
		//For each assembly node
		while(iter.hasNext()){
			Collection<Assembly> genes=iter.next().getContainedValues();
			//For each gene(assembly)
			for(Annotation gene: genes){
				//Get all the exons of the gene
				List<? extends Annotation> blocks=gene.getBlocks();
				//if the gene has 1 exon
				if(blocks.size()==1){
					double pval = getScores(gene)[1];
					//System.err.println(gene.toUCSC()+" Count: "+score.getCount()+" pval:"+pval+"new count: "+s+"  new p-val "+pval2);
					if(pval<alpha){
						//System.err.println("passes");
						graph.connectVertexToGraph(blocks.get(0));
					}
					else{
						logger.debug(gene.getName()+" is not added to the graph because it does not pass the p-value threshold.");
					}
				}
				else{
			/*		double pval = getScores(gene)[1];
					if(pval<DEFAULT_ALPHA){*/
						graph.addAnnotationToGraph(gene);
			/*		}
					else{
						logger.debug("Gene "+gene.toUCSC()+" is filtered out because it does not meet the significance threshold : "+pval);
					}*/
				}
			}
		}
		return graph;
	}	

	/**
	 * Returns paired end counts and scan p-value for the specified gene
	 * @param gene
	 * @return
	 */
	private double[] getScores(Annotation gene){
		double[] scores = new double[2];
		scores[0] = 0.0;
		//Get all reads overlapping the transcript
		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(gene,true), new PairedAndProperFilter());
		//For each read,
		while(iter.hasNext()){					
			Alignment read = iter.next();
			boolean countRead = true;
			for(Annotation mate:read.getReadAlignments(space)){
				if(!compatible(gene,mate)){
					//logger.debug("Read "+mate.toUCSC()+" is not compatible with isoform with "+isoform.getExons().length);
					countRead=false;
					break;
				}
			}
			//For the assembly, we need to treat each read separately	
			if(countRead){
				scores[0] += read.getWeight();
			}
		}
		iter.close();
		//logger.debug("Count = "+scores[0]+" Int version "+new Double(scores[0]).intValue()+" global paired lambda = "+globalPairedLambda+" gene size = "+model.getCoordinateSpace().getSize(gene)+ " or "+gene.size()+" global length = "+model.getGlobalLength()+" global lambda = "+model.getGlobalLambda());
		scores[1] = AlignmentDataModelStats.calculatePVal(new Double(scores[0]).intValue(), globalPairedLambda,gene.size(), model.getGlobalLength());
		
//		scores[1] = AlignmentDataModelStats.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		
		return scores;
	}
	
	private void mergeAssembly(IntervalTree<Assembly> tree) {
		//Step 1: Iterate through the working assembly and find paths that overlap
		Iterator<Assembly> iter=tree.toCollection().iterator();
		//if overlaps but is incompatible, see if we can branch
		//iterate
		Collection<Assembly> considered = new HashSet<Assembly>();
 		logger.debug("Enter merging");
		while(iter.hasNext()){
			Assembly branch1=iter.next();
			considered.add(branch1);
			//get overlapping branches
			Iterator<Assembly> overlappers=tree.overlappingValueIterator(branch1.getStart(), branch1.getEnd());
			//Collection<Assembly> toRemove=new TreeSet<Assembly>();
			while(overlappers.hasNext()){
				Assembly branch2=overlappers.next();
				if(!considered.contains(branch2)){
					if(!branch1.equals(branch2) && compatible(branch1, branch2)){
						logger.debug("Merging: "+ branch1.getName()+ " and "+branch2.getName());
						Assembly merged=merge(branch1, branch2);
						merged.setName(branch1.getName());
						//remove annotation1 and annotation2
						tree.remove(branch1.getStart(), branch1.getEnd(), branch1);
						tree.remove(branch2.getStart(), branch2.getEnd(), branch2);
						//add merged
						tree.put(merged.getStart(), merged.getEnd(), merged);
					}
				}
			}
		}
	}
	
	/**
	 * This function removes any overlapping bad isoforms
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> intronRetentionFilter(IntervalTree<Assembly> assemblies,String chr) {
		//BAMFileWriter writer = new BAMFileWriter(new File(outName+"."+chr+".filtered.bam"));
		//writer.setSortOrder(SortOrder.coordinate, false);
		//writer.setHeader(header);
		//writer.setSortOrder(SortOrder.coordinate, false);
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
/*		double normalizationFactor = (double)spliceCount/(double)(spliceCount+nonSpliceCount);
		logger.error("R = "+normalizationFactor);*/
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();
		//For each assembly
		while(iter.hasNext()){
			Assembly assembly1=iter.next();
			//If assembly is not already marked as spurious
			if(!assembly1.isSpurious()){

				boolean	toRemoveAssembly1 = false;
				//For all overlapping assemblies
				Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly1.getStart(), assembly1.getEnd());
				while(overlappers.hasNext() && !toRemoveAssembly1){
					Assembly assembly2 = overlappers.next();
					if(assembly1.getOrientation().equals(assembly2.getOrientation())){
						//For all exons of this assembly
						for(Annotation intron1: assembly1.getSpliceConnections()){
							if(!toRemoveAssembly1){
								for(Annotation exon2:assembly2.getBlocks()){
									//If intron in assembly1 is contained in exon of assembly2
									if(exon2.contains(intron1)){
										//CHECK FOR COVERAGE
										double intronCount = model.getIntronCounts(intron1);
										double exonCount = (model.getCountStranded(new BasicAnnotation(exon2.getChr(), intron1.getStart(), intron1.getEnd(), exon2.getOrientation()), false)
																	/(double)(intron1.getEnd()-intron1.getStart()));
										double ratioI = intronCount/(intronCount+exonCount);
										double ratioE = exonCount/(intronCount+exonCount);
										//CASE 1 : Assembly with intron is mature and assembly with exon is immature
										if(ratioI>THRESHOLD_SPURIOUS){
											logger.debug("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.getName()+ " to intron "+intron1.toUCSC()+" in assembly "+assembly1.getName());
											assembly2.setSpurious(true);
											logger.debug("CASE 1: assembly "+ assembly2.getName()+" is spurious. Intron count "+intronCount+" > "+exonCount+" Ratio: "+ratioI);
										}
										else{
											//CASE 2: Assembly with exon is real and the other 
											if(ratioE>THRESHOLD_SPURIOUS){
												//logger.warn(assembly.toUCSC()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC());
												logger.debug("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.getName()+ " to intron "+intron1.toUCSC()+" in assembly "+assembly1.getName());
												assembly1.setSpurious(true);
												toRemoveAssembly1 = true;
												logger.debug("CASE 2: assembly "+ assembly1.getName()+" is spurious.Exon count "+exonCount+" > "+intronCount+" Ratio: "+(ratioE));
												break;
											}
											//else{
												//Both are real
											//}
										}
									}
								}
							}
							else{
								break;
							}
						}
					}
				}
				if(!toRemoveAssembly1){
					currentAssemblies.put(assembly1.getStart(), assembly1.getEnd(), assembly1);
					//for this assembly, get all the 
					//writer.addAlignment(alignment)
				}
			}
		}
		
		iter=currentAssemblies.toCollection().iterator();
		Set<Alignment> reads = new HashSet<Alignment>();
		while(iter.hasNext()){
			Assembly assembly = iter.next();
			CloseableIterator<Alignment> riter = model.getOverlappingReads(assembly, true);
			while(riter.hasNext()){
				Alignment read = riter.next();
				//System.out.println(read.getName());
				if(!reads.contains(read)){
					//writer.addAlignment(read.toSAMRecord());
					reads.add(read);
				}
			}
			riter.close();
		}
		//writer.close();
		return currentAssemblies;
	}
	
	/**
	 * This function removes any premature assemblies that have coverage conflicts in the exons
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removePrematureAssemblies(IntervalTree<Assembly> assemblies) {
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();

		while(iter.hasNext()){
			Assembly assembly=iter.next();
			boolean	toRemove = false;
			//For all overlapping assemblies
			Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly.getStart(), assembly.getEnd());
			while(overlappers.hasNext()){
				Assembly overlapper = overlappers.next();
				if(assembly.getOrientation().equals(overlapper.getOrientation())){
					
					//Check if the intron"ed" assembly is confident so it can be used
					if(!overlapper.isConfidentIsSet()){
						setConfidence(overlapper);
					}
					if(overlapper.isConfident()){
						//logger.warn("Compare " +assembly.getName()+" and "+overlapper.getName());
						//CHECK IF THE PARTS OF THE ASSEMBLY OTHER THAN THIS EXON OF THE ASSEMBLY) ARE COMPATIBLE
						//For all exons of this assembly
						for(Annotation exon1: assembly.getBlocks()){
							if(!toRemove){
								//If exon overlaps introns of overlapping assembly
								for(Annotation intron2:overlapper.getSpliceConnections()){
									//Flag as premature				
									if(exon1.overlaps(intron2)){
										Annotation p = assembly.minus(exon1);
										Annotation q = overlapper.minus(exon1);
										if(compatible(p, q) || p.length()==0){
											assembly.setPossiblePremature(true);
											//CHECK FOR COVERAGE
											if(coveragePassesCheck(exon1,intron2)){
												//toRemove = false;
											}
											else{
												//logger.warn(assembly.getName()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC()+" of "+overlapper.getName());
												toRemove = true;
												break;
											}
										}
									}
								}
							}
							else{
								break;
							}
						}
					}
				}
			}
			if(!toRemove){
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		
		return currentAssemblies;
	}
	
	private void setConfidence(Assembly overlapper){
		//for each intron
		for(Annotation intron2:overlapper.getSpliceConnections()){
			Annotation[] exons = overlapper.getFlankingBlocks(intron2);
			//exon[0] is left of intron
			//coverage to left of left exon
			//Since it is one base, we dont need coverage
			double leftScore = model.getCountStranded(new BasicAnnotation(exons[0].getChr(),exons[0].getEnd()-2,exons[0].getEnd()-1,overlapper.getOrientation()),false);
			//coverage to right of right exon
			double rightScore = model.getCountStranded(new BasicAnnotation(exons[1].getChr(),exons[1].getStart()+1,exons[1].getStart()+2,overlapper.getOrientation()),false);
			
			if(leftScore<rightScore){
				if(leftScore<rightScore*coveragePercentThreshold){ 
					logger.debug(overlapper.getName()+" does not pass confidence test because "+exons[0].toUCSC()+" has score "+leftScore+" compared to "+
								exons[1].toUCSC()+" which has "+rightScore);
					overlapper.setConfident(false);
				}
			}else{
				if(rightScore<leftScore*coveragePercentThreshold){
					logger.debug(overlapper.getName()+" does not pass confidence test because "+exons[1].toUCSC()+" has score "+rightScore+" compared to "+
							exons[0].toUCSC()+" which has "+leftScore);
					overlapper.setConfident(false);
				}
			}
		}
		if(!overlapper.isConfidentIsSet())
			overlapper.setConfident(true);
	}
	/**
	 * This function removes any premature assemblies that have coverage conflicts in the exons
	 * @param assemblies
	 * @return
	 */
/*	private IntervalTree<Assembly> removePrematureAssembliesUsingPairedEnds(IntervalTree<Assembly> assemblies) {
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();

		while(iter.hasNext()){
			Assembly assembly=iter.next();
			boolean	toRemove = false;
			//For all overlapping assemblies
			Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly.getStart(), assembly.getEnd());
			while(overlappers.hasNext()){
				Assembly overlapper = overlappers.next();
				//For all exons of this assembly
				for(Annotation exon1: assembly.getBlocks()){
					if(!toRemove){
						//If exon overlaps introns of overlapping assembly
						for(Annotation intron2:overlapper.getSpliceConnections()){
							//Flag as premature				
							if(exon1.overlaps(intron2)){
								assembly.setPossiblePremature(true);
								//CHECK FOR COVERAGE
								if(exonPassesPairedEndTest(exon1,overlapper,intron2)){
									//toRemove = false;
								}
								else{
									logger.warn(assembly.toUCSC()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC());
									toRemove = true;
									break;
								}
							}
						}
					}
					else{
						break;
					}
				}
			}
			if(!toRemove){
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		
		return currentAssemblies;
	}*/
	
	/**
	 * This function returns true if the number of paired end reads supporting exon1 and supporting the splice junction og intron2 are comparable, 
	 * that is, not below a certain threshold. HELPER function to removePrematureAssembliesUsingPairedEnds()
	 * @return
	 */
	private boolean exonPassesPairedEndTest(Annotation exon1,Assembly overlapper,Annotation intron2){
		
		Annotation overlap = exon1.intersect(intron2);
		Annotation nonoverlap = exon1.minus(intron2);
		
		//If exon is completely contained in the intron, pass
		if(overlap==null || overlap.equals(exon1)){
			return true;
		}
		boolean passes=true;
		//Get paired end reads for exon1
		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(exon1,true), new ProperPairFilter());
		double exonCount = 0.0;
		while(iter.hasNext()){
			Alignment read = iter.next();
			// if read is compatible with gene
			//If either read is not compatible flag will be set to false
			//If both mates of the pair are compatible
			if(boundaryIsCompatibleWithPairedEnd(overlap,nonoverlap,read)){
				exonCount = exonCount + 1.0;
			}
		}
		//exonCount = exonCount/(double)exon1.length();
		iter.close();
		
		//Get paired end reads for exons flanking intron2
		
		Annotation[] flankingExons = overlapper.getFlankingBlocks(intron2);
		Assembly sumAssembly = new Assembly(Arrays.asList(flankingExons),false);
		//int intronLength = flankingExons[0].length()+flankingExons[1].length();
		iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(sumAssembly,true), new ProperPairFilter());
		double intronCount = 0.0;
		while(iter.hasNext()){
			Alignment read = iter.next();
			// if read is compatible with gene
			//If either read is not compatible flag will be set to false
			//If both mates of the pair are compatible
			if(boundaryIsCompatibleWithPairedEnd(flankingExons[0],flankingExons[1],read)){
				intronCount = intronCount + 1.0;
			}
		}
		iter.close();
		//intronCount = intronCount / (double)intronLength;
		if(intronCount!=0.0){
			double ratio = exonCount/intronCount;
			if(ratio<0.8){
				logger.debug("Paired end count for exon " +exon1.toUCSC()+" : "+exonCount);
				logger.debug("Paired end count for intron " +intron2.toUCSC()+" : "+intronCount+ " ratio = "+ratio);
				return false;
			}
		}
		return passes;
	}
	
	/**
	 * Returns an array of size 2 where [0] is the first mate and [1] is the second mate
	 * @param read
	 * @return
	 */
	private Annotation[] getPairedEndReads(Alignment read){
		
		//We implement this specific to our purpose with paired ends
		Annotation[] annArr = new Annotation[2];
		int i=0;
		for(Annotation a: read.getReadAlignments(space)){
			if(i>1){
				logger.error("More than 2 mates are returned");
			}
			annArr[i] = a;
			i++;
		}
		return annArr;
	}
	
	/**
	 * This function returns true is the first mate is in the first exon and the second mate is in the second exon
	 * @param exon1
	 * @param exon2
	 * @param mates
	 * @return
	 */
	private boolean boundaryIsCompatibleWithPairedEnd(Annotation exon1,Annotation exon2,Alignment read){
		
		Annotation[] mates = getPairedEndReads(read);
		if((compatible(exon1,mates[0]) && 
				compatible(exon2,mates[1]))||
					(compatible(exon2,mates[0]) 
							&& compatible(exon1,mates[1])))
			return true;
		else
			return false;
	}
	/**
	 * This function will remove any assemblies that have introns with very low number of splice junctions
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removeLowSpliceJunctionAssemblies(IntervalTree<Assembly> assemblies){
	
		logger.info("Remove low splice junction assemblies");
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
				
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();//assemblies;

		//For each assembly
		while(iter.hasNext()){
			Assembly assembly=iter.next();
			
			boolean toRemove=false;
			//If the assembly does not pass the confidence test, remove it
			if(!assembly.isConfidentIsSet()){
				logger.debug("Confidence is set in the splice filter");
				setConfidence(assembly);
			}
			if(!assembly.isConfident()){
				logger.debug(assembly.getName()+" removed because does not pass the confidence test.");
			} else{		
				//IF THE GENE HAS ONE INTRON, check min #splice reads
				if(assembly.getSpliceConnections().size()==1){
					for(Annotation intron:assembly.getSpliceConnections()){
						double count=model.getIntronCounts(intron);
						if(count<minSpliceReads){
							toRemove = true;
						}
					}
					if(!toRemove)
						currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
				}
				else{
					Map<Annotation,Double> intronToSplicedCountMap = new TreeMap<Annotation,Double>();
					double avgCount = 0.0;
					//boolean toRemove = false;
					//for each intron
					for(Annotation intron:assembly.getSpliceConnections()){
						//find the number of supporting splice junctions
						double count=model.getIntronCounts(intron);											
						avgCount +=count;
						intronToSplicedCountMap.put(intron, count);
					}
					
					//Calculate average
					avgCount = (avgCount / (double)intronToSplicedCountMap.keySet().size());
					toRemove = false;
					//Check for all introns
					int cnt=0;
					//First check the first and last intron. If they do not pass the threshold
					//then trim and then go over this again.
					for(Annotation intron:intronToSplicedCountMap.keySet()){
						if(cnt==0 || cnt==intronToSplicedCountMap.keySet().size()-1){
							if(intronToSplicedCountMap.get(intron)<=avgCount*minSplicePercent){
								assembly = new Assembly(trimEnds(assembly,0.25));
							}
						}
						cnt++;
					}
					
					for(Annotation intron:assembly.getSpliceConnections()){
						if(intronToSplicedCountMap.get(intron)<=avgCount*minSplicePercent){
							//If  assembly is flagged to be removed because of an intron,
							logger.debug(assembly.getName()+" removed because intron "+intron.toUCSC()+" has coverage "+intronToSplicedCountMap.get(intron)+" compared to "+avgCount);
							toRemove = true; 
							//return true;
						}
					}
					if(!toRemove){
						//Removed by this filter but confidence was set to remove
						/*if(!assembly.isConfident())
							logger.debug(assembly.getName()+" is NOT confident but RETAINED.");*/
						currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
					}
					else{
						//Assembly flagged to remove
						logger.debug(assembly.getName()+" is confident but flagged to be REMOVED. Check trimming.");
						//trimEnds(assembly);
						
					}
				}
			}
		}
		return currentAssemblies;
	}

	/**
	 * Returns the assembly trimmed at both ends
	 * @param assembly
	 * @return
	 */
	private Annotation trimEnds(Annotation assembly,double pct){
		
		List<Double> counts = model.getCountsStrandedPerPosition(assembly);
		double[] cntArr = l2a(counts);

		Collections.sort(counts);		
		double cutoff = Math.max(2, Statistics.quantile(counts, pct));
/*		for(int i=0;i<assembly.size();i++){
			System.out.print(cntArr[i]+" ");
		}*/
		int trimStart = MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(cntArr, cutoff);
		int trimEnd   =  assembly.size() - MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(cntArr, cutoff);
		logger.debug(assembly.getName()+" trimStart " + trimStart + " trimEnd " + trimEnd + ", transcript length " + assembly.size()+ " with cutoff "+cutoff);

/*			int transcriptFirstExonEnd = assembly.getOrientation().equals(Strand.NEGATIVE) 
					? assembly.getPositionAtReferenceCoordinate(assembly.getStart())
					: assembly.getPositionAtReferenceCoordinate(assembly.getEnd()-1);

			int transcriptLastExonStart = assembly.getOrientation().equals(Strand.NEGATIVE)
						? assembly.getPositionAtReferenceCoordinate(assembly.getEnd()-1)
						: assembly.getPositionAtReferenceCoordinate(assembly.getStart());
				
			logger.trace("first exon end in transcript " + transcriptFirstExonEnd + " last exon start " + transcriptLastExonStart);

			if(trimStart > transcriptFirstExonEnd) {
				trimStart = Math.max(0, transcriptFirstExonEnd - 50);
			}	
			if(trimEnd < transcriptLastExonStart) {
				trimEnd = Math.min(align.getTranscriptLength(), transcriptLastExonStart + 50);
			}	
*/			
		if(trimStart>trimEnd){
			return assembly;
		}
		if(assembly.getOrientation().equals(Strand.NEGATIVE)){
			int temp = trimStart;
			trimStart = trimEnd;
			trimEnd = temp;
			logger.debug("Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
		}	
			
		Assembly newAssembly = new Assembly(assembly);
			
		if(((trimEnd+trimStart) < assembly.size()) && trimStart<trimEnd && trimStart<assembly.size() && trimEnd<assembly.size()){
			newAssembly.trim(trimStart, trimEnd);
			logger.debug("trimming ("+trimStart +" - "+ trimEnd+") gene was: " + assembly.toBED() + " and now is: " +newAssembly.toBED());
		}
		return newAssembly;
	}
	
	/**
	 * Helper function to removePrematureAssemblies
	 * @param exon1
	 * @param intron2
	 * @return
	 */
	private boolean coveragePassesCheck(Annotation exon1, Annotation intron2){
		Annotation overlap = exon1.intersect(intron2);
		//Annotation nonoverlap = exon1.minus(intron2);
		
		if(overlap==null){
			return true;
		}
		//If fully contained
		//TODO: A minimum coverage threshold?
		if(overlap.equals(exon1)){
			//logger.error("Fully Contained");
			return true;
		}
		
		if(overlap.getLengthOnReference()<=3){
			//logger.warn(exon1.toUCSC()+" overlaps "+intron2.toUCSC()+" with <=3");
			return false;
		}
		
		Annotation exonBoundary = null;
		Annotation intronBoundary = null;
		if(intron2.getEnd()>exon1.getEnd()){
			exonBoundary = new BasicAnnotation(exon1.getChr(),intron2.getStart()-2,intron2.getStart()-1,exon1.getOrientation());
			intronBoundary = new BasicAnnotation(intron2.getChr(),intron2.getStart()+1,intron2.getStart()+2,exon1.getOrientation());
		}
		else{
			exonBoundary = new BasicAnnotation(exon1.getChr(),intron2.getEnd()+1,intron2.getEnd()+2,exon1.getOrientation());
			intronBoundary = new BasicAnnotation(intron2.getChr(),intron2.getEnd()-2,intron2.getEnd()-1,exon1.getOrientation());
		}
		double overlapScore = model.getCountStranded(intronBoundary,false);
		double nonoverlapScore = model.getCountStranded(exonBoundary,false);
//		logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
		if(overlapScore>=nonoverlapScore*coveragePercentThreshold){
			//logger.error(exon1.toUCSC()+" passes coverage test with intron "+ intron2.toUCSC());
			//logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
			//logger.debug("Intron boundary score for "+intronBoundary.toUCSC()+ " : "+overlapScore+" Exon boundary score for "+exonBoundary.toUCSC()+ " : "+nonoverlapScore);
			return true;
		}
		else{
			//logger.error(exon1.toUCSC()+" does not pass coverage test with intron "+ intron2.toUCSC());
			logger.debug("REMOVED: Intron boundary score for "+intronBoundary.toUCSC()+ " : "+overlapScore+" Exon boundary score for "+exonBoundary.toUCSC()+ " : "+nonoverlapScore);
			return false;
		}		
	}

	/**
	 * 
	 * @param iter
	 * @param workingAssemblies
	 * @param strand
	 * @return
	 */
	private IntervalTree<Assembly> assembleDirectly(CloseableIterator<Alignment> iter, IntervalTree<Assembly> workingAssemblies,TranscriptionRead strand) {
		
		String linc="gene_v2_";
		long cnt=0;
		long start = System.currentTimeMillis();
		int c=0;
		boolean flagPremature=!workingAssemblies.isEmpty();
		while(iter.hasNext()){
			
			Alignment reads=iter.next();			
			//reads.setFragmentStrand(strand);				
			//For the assembly, we need to treat each read separately
			for(Annotation read: reads.getReadAlignments(space)){
				cnt++;
				if(cnt%100000==0){
					long end = System.currentTimeMillis();
					logger.debug("Processed "+cnt+" reads in "+(end-start)/1000+" seconds with "+workingAssemblies.size()+" "+c);
					start = System.currentTimeMillis();
				}
				
				//Find all compatible assemblies
				Collection<Assembly> compatibleAssemblies = new ArrayList<Assembly>();

/*				logger.debug("New read : "+read.toUCSC()+" "+!read.getSpliceConnections().isEmpty());
				for(Assembly as:workingAssemblies.toCollection()){
					logger.debug(as.toUCSC()+"\t"+as.toBED());
				
				}*/
				//EACH READ HAS THE FRAGMENT STRAND
				//for each read, get overlapping assemblies
				Iterator<Node<Assembly>> overlappers=workingAssemblies.overlappers(read.getStart(), read.getEnd());
				//if no overlappers add read as assembly
				if(!overlappers.hasNext()){
					//add the read as an annotation
					//Flag this as likely premature
					Assembly readAssembly=new Assembly(read, false);
					readAssembly.setName(linc+globalCounter);
					globalCounter++;
					if(flagPremature){
						readAssembly.setPossiblePremature(true);
					}
					workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
					c++;
				}
				else{
					boolean hasCompatible=false;
//					Collection<Assembly> overlapNotCompatible = new TreeSet<Assembly>();
//					long s1 = System.currentTimeMillis();
					while(overlappers.hasNext()){
						Collection<Assembly> assemblies=new TreeSet<Assembly>(overlappers.next().getContainedValues());
						//Find all compatible assemblies
						for(Assembly assembly: assemblies){
							if(compatible(assembly, read)){
								compatibleAssemblies.add(assembly);
								Assembly merged=merge(assembly, read);
								merged.setName(assembly.getName());
								//remove assembly
								workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
								//add merged
								workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
								hasCompatible=true;
							}
/*							else{
								if(assembly.overlaps(read)){
									overlapNotCompatible.add(assembly);	
								}
							}*/
						}
						
					}
/*					if(cnt>2900){
						logger.debug("compatible assemblies size "+compatibleAssemblies.size());
						logger.debug("overlap not compatible "+overlapNotCompatible.size());
					}*/
					
					//If some compatible assembly
/*					if(hasCompatible & compatibleAssemblies.size()>0){
						//MERGE NON-OVERLAPPING ASSEMBLIES
						if(compatibleAssemblies.size()>1){
							Collection<Assembly> alreadyMerged = new HashSet<Assembly>();
							Collection<String> names = new HashSet<String>();
							boolean mergedFlag = false;
							//Compare all compatible assemblies to each other
							for (Iterator<Assembly> iterator = compatibleAssemblies.iterator(); iterator.hasNext();) {
								Assembly assembly1 = iterator.next();
								for (Iterator<Assembly> iterator2 = compatibleAssemblies.iterator(); iterator2.hasNext();) {
									Assembly assembly2 = iterator2.next();
									if(!(assembly1.equals(assembly2)) && (!names.contains(assembly2.getName()+"_"+assembly1.getName()))){
										if(!assembly1.overlaps(assembly2) && assembly1.getOrientation().equals(assembly2.getOrientation())){
											alreadyMerged.add(assembly1);
											alreadyMerged.add(assembly2);
											mergedFlag=true;
											//Merge assembly1 and assembly2 and read
											Assembly merged=merge(assembly1, read);
											merged=merge(merged,assembly2);
											for(Assembly a:workingAssemblies.toCollection()){
												if(merged.equals(a)){
													logger.debug("Merged: "+merged.toBED());
													logger.debug("a = "+a.toBED());
												}
											}
											String name =assembly1.getName()+"_"+assembly2.getName();
											merged.setName(name);
											workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
											c++;
											names.add(name);
										}
									}
								}
								if(!mergedFlag){
									Assembly merged=merge(assembly1, read);
									merged.setName(assembly1.getName());
									//remove assembly
									workingAssemblies.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
									//add merged
									workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
								}
								iterator.remove();
							}
							
							for(Assembly assembly:alreadyMerged){
								c--;
								workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
							}
						}
						else{
							//Only 1 compatible assembly
							Iterator<Assembly> compIter=compatibleAssemblies.iterator();
							while(compIter.hasNext()){
								Assembly comp = compIter.next();
								//if compatible --> merge
								Assembly merged=merge(comp, read);
								merged.setName(comp.getName());
								//remove assembly
								workingAssemblies.remove(comp.getStart(), comp.getEnd(), comp);
								//add merged
								workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
							}
						}
					}*/
					if(!hasCompatible){
/*						boolean splittable =false;
						if(!overlapNotCompatible.isEmpty()){
//							s1 = System.currentTimeMillis();
							for(Assembly assembly:overlapNotCompatible){
								if(read.getSpliceConnections().isEmpty()){
									//OLD SPLIT
									Collection<Assembly> branches=oldSplitBranch(assembly, read);
									if(branches!=null && !branches.isEmpty()){
										splittable=true;
										for(Assembly branch: branches){
											branch.setName(linc+globalCounter);
											globalCounter++;
											//add merged assemblies to the tree
											workingAssemblies.put(branch.getStart(), branch.getEnd(), branch);
										}
									}
								}
								else{
								//NEW SPLIT
									Collection<Assembly> branches=newSplitBranch(assembly, read);
									if(branches!=null && !branches.isEmpty()){
										splittable=true;
										//remove assembly
										//workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
										for(Assembly branch: branches){
											//add merged assemblies to the tree
											branch.setName(linc+globalCounter);
											globalCounter++;
											workingAssemblies.put(branch.getStart(), branch.getEnd(), branch);
										}
									}
								}
								
								//SOPHISTICATED SPLITTING
//								Annotation portionToConsiderAdding=assembly.intersect(new Alignments(assembly.getChr(), assembly.getStart(), read.getBlocks().get(0).getEnd()));
								if(assembly.getOverlappingAssembly(read)==null){
									logger.debug("Assembly: "+assembly.toBED());
									logger.debug("Read: "+read.toUCSC()+" "+!read.getSpliceConnections().isEmpty()+" "+assembly.overlaps(read));
								}
								Annotation portionToConsiderAdding=assembly.intersect
										(new Alignments(assembly.getChr(), assembly.getStart(), 
												(assembly.getOverlappingAssembly(read)).getEnd()));
								if(compatible(portionToConsiderAdding,read)){
									Assembly branch = merge(portionToConsiderAdding,read);
									branch.setName(linc+globalCounter);
									globalCounter++;
									workingAssemblies.put(branch.getStart(), branch.getEnd(), branch);
									c++;
									splittable=true;
								}
							}
							if(!splittable){
								Assembly readAssembly=new Assembly(read, false);
								readAssembly.setName(linc+globalCounter);
								globalCounter++;
								if(flagPremature){
									readAssembly.setPossiblePremature(true);
								}
								workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
								c++;
							}
						}
						else{*/
							Assembly readAssembly=new Assembly(read, false);
							readAssembly.setName(linc+globalCounter);
							globalCounter++;
							if(flagPremature){
								readAssembly.setPossiblePremature(true);
							}
							workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
							c++;
//						}
					}
				}
				
			}
		}		
		iter.close(); //close the iterator
		
		return workingAssemblies;
	}

	/**
	 * See if the read is partially compatible with the assembly
	 * If so, branch and merge
	 * @param assembly
	 * @param read
	 * @return all branches
	 */
	private Collection<Assembly> newSplitBranch(Assembly assembly, Annotation read) {

		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		if(read.overlaps(assembly)){
			Assembly truncated=findMaxCompatibilityNew(assembly, read);
			if(truncated!=null){
				Assembly merged=merge(truncated, read);
				rtrn.add(merged);
				//Why add the assembly again? Just don't remove it above
				//rtrn.add(assembly);
			}
		}
		return rtrn;
	}

	/**
	 * See if the read is partially compatible with the assembly
	 * If so, branch and merge
	 * @param assembly
	 * @param read
	 * @return all branches
	 */
	private Collection<Assembly> oldSplitBranch(Assembly assembly, Annotation read) {

		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		if(read.overlaps(assembly)){
			Assembly truncated=findMaxCompatibility(assembly, read);
			if(truncated!=null){
				Assembly merged=merge(truncated, read);
				rtrn.add(merged);
				//Why add the assembly again? Just don't remove it above
				rtrn.add(assembly);
			}
		}
		return rtrn;
	}
	
	private Assembly findMaxCompatibilityNew(Assembly assembly, Annotation read) {
		/*Iterator<Assembly> iter=assembly.trimNodes(read);
		while(iter.hasNext()){
			Assembly truncated=iter.next();
			if(compatible(truncated, read)){
				return truncated;
			}
		}
		*/
		Assembly trimmedAssembly = assembly.trim(read);
		if(trimmedAssembly!=null){
			//TODO: doesnt the logic dictate this?
			if(compatible(trimmedAssembly.getLastBlock(), read)){
				return trimmedAssembly;
			}
		}
		return null;
	}
	
	private Assembly findMaxCompatibility(Assembly assembly, Annotation read) {
		Iterator<Assembly> iter=assembly.trimNodes();
		while(iter.hasNext()){
			Assembly truncated=iter.next();
			if(compatible(truncated, read)){
				return truncated;
			}
		}
		return null;
	}

	/*private boolean partiallyCompatible(Annotation assembly, Annotation read) {
		//The annotations are not compatible, so we will check if they are partially compatible
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//If they dont overlap they cant be partially compatible
		if(!assembly.overlaps(read)){return false;}
		
		//Otherwise, they are partially compatible if:
		//(i) some of the introns are compatible
		if(!assemblyIntrons.isEmpty() && !readIntrons.isEmpty()){
			//if there are some introns that
		}
		
		//(ii) or if exons of one overlap the other but dont overlap the intron of the other
		
	}*/

	private Collection<Annotation> branch(Annotation read, Annotation assembly) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//Define incompatible regions
		
		//sanity check to ensure they physically overlap
		if(read.overlaps(assembly)){
			//define compatible edges: find introns that overlap but arent equal
			//define compatible nodes: find exons that overlap and are compatible
			Collection<Annotation> incompatible=this.getIncompatibleRegions(read, assembly);
			try{write(incompatible, "incompatible.bed");}catch(IOException ex){}
			Collection<Annotation> compatible=getCompatibleRegions(read, assembly, incompatible); //simple subtraction of intervals
			for(Annotation incompatibleRegion: incompatible){
				Annotation branch=merge(compatible, incompatibleRegion);
				rtrn.add(branch);
			}
		}
		
		return rtrn;
		
	}

	
	private Collection<Annotation> getCompatibleRegions(Annotation read, Annotation assembly, Collection<Annotation> incompatible) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		//simply take read and assembly and subtract the incompatible regions
		if(!incompatible.isEmpty()){
			Annotation ann1=(read.minus(incompatible));
			Annotation ann2=(assembly.minus(incompatible));
			if(ann1.getSize()>0 && compatible(ann1, read) && compatible(ann1, assembly)){
				rtrn.add(ann1);
			}
			if(ann2.getSize()>0 && compatible(ann2, read) && compatible(ann2, assembly)){
				rtrn.add(ann2);
			}
			
		}
		else{
			rtrn.add(read);
			rtrn.add(assembly);
		}
		return rtrn;
	}

	

	private Annotation merge(Collection<Annotation> compatible, Annotation incompatibleRegion) {
		//do a simple block merge
		Collection<Annotation> blocks=new TreeSet<Annotation>();
		for(Annotation ann: compatible){
			blocks.addAll(ann.getBlocks());
		}
		blocks.addAll(incompatibleRegion.getBlocks());
		Annotation rtrn=new BasicAnnotation(blocks);
		return rtrn;
	}

	
	
	//MG: This was working well
	public static boolean compatible(Annotation assembly, Annotation read) {
		//Two alignments will be defined as compatible if:
		//(i) the intronic locations are exactly same
		//if both have introns, ensure that there are no non-overlapping introns
		//logger.debug(read.getName());
		if(!assembly.getSpliceConnections().isEmpty() 
				|| !read.getSpliceConnections().isEmpty()){
			//check if introns are compatible
			if(areIntronsCompatible(assembly, read)){
				//logger.debug("Introns are compatible");
				return true;
			}
			//logger.debug("Introns are NOT compatible");
			//TODO In the case of partial compatibility we should split and make new path
			return false;
		}
		//(ii) one is spliced, the other is not 
		/*if(assembly.getBlocks().size()>1 || !read.getSpliceConnections().isEmpty()){
			//and the non-spliced does not overlap the intron at all, but overlaps the exon
			boolean overlapsWithoutCrossingIntron=overlapsWithoutCrossingIntron(assembly, read);
			if(overlapsWithoutCrossingIntron){return true;}
			return false;
		}*/
		//(iii) both are unspliced and overlap
		boolean overlap=overlap(assembly, read);
		
		if(overlap){return true;}
		return false;
	}
	

	private static boolean areIntronsCompatible(Annotation assembly, Annotation read) {
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//introns are compatible if:
		//(i) all overlapping introns are identical
		for(Annotation intron1: assemblyIntrons){
			for(Annotation intron2: readIntrons){
				//if overlaps but not identical
				if(intron1.overlaps(intron2,forceStrandSpecificity)){
					if(!intron1.equals(intron2, forceStrandSpecificity)){
					//	logger.debug("Case1");
						return false;
					} //TODO: This should use strand info or not based on flag
				}
			}
		}
		
		//(ii) if none of the exons overlap an intron in the other
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation intron2: readIntrons){
				if(exon1.overlaps(intron2,forceStrandSpecificity)){
					//logger.debug("Case2a");
					return false;
				}
			}
		}
		
		for(Annotation exon2: read.getBlocks()){
			for(Annotation intron1: assemblyIntrons){
				if(exon2.overlaps(intron1,forceStrandSpecificity)){
					//logger.debug("Case2b");
					return false;
				}
			}
		}
		
		//(ii) the introns dont overlap but the exons do
		//just need to test that any exons overlap
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation exon2: read.getBlocks()){
				if(exon1.overlaps(exon2,forceStrandSpecificity)){
					return true;
				}
			}
		}
		//logger.debug("Case3");
		return false;
	}
	
	
	private Collection<Annotation> getIncompatibleRegions(Annotation assembly, Annotation read){
		
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		Collection<Annotation> rtrn1=new TreeSet<Annotation>();
		Collection<Annotation> rtrn2=new TreeSet<Annotation>();
		
		//introns are compatible if:
		//(i) all overlapping introns are identical
		for(Annotation intron1: assemblyIntrons){
			for(Annotation intron2: readIntrons){
				//if overlaps but not identical
				if(intron1.overlaps(intron2)){
					if(!intron1.equals(intron2, forceStrandSpecificity)){
						rtrn1.add(intron1);
						rtrn2.add(intron2);
						//intron1 and intron2 are not compatible
					} //This should use strand info or not based on flag
				}
			}
		}
		
		
		//now check if an exon overlaps an intron in the other
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation intron2: readIntrons){
				if(exon1.overlaps(intron2)){
					//exon1 is incompatible
					rtrn1.add(exon1);
					rtrn2.add(intron2);
				}
			}
		}
		
		for(Annotation exon2: read.getBlocks()){
			for(Annotation intron1: assemblyIntrons){
				if(exon2.overlaps(intron1)){
					rtrn2.add(exon2);
					rtrn1.add(intron1);
				}
			}
		}
		
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//merge rtrn1
		if(!rtrn1.isEmpty()){
			Annotation block1=new BasicAnnotation(rtrn1);
			//trim assembly
			Annotation tmp1=assembly.intersect(new Alignments(block1.getChr(), block1.getStart(), block1.getEnd()));
			if(tmp1!=null && tmp1.getSize()>0){rtrn.add(tmp1);}
		}
		
		//merge rtrn2
		if(!rtrn2.isEmpty()){
			Annotation block2=new BasicAnnotation(rtrn2);
			//trim read
			Annotation tmp2=read.intersect(new Alignments(block2.getChr(), block2.getStart(), block2.getEnd()));
			if(tmp2!=null && tmp2.getSize()>0){rtrn.add(tmp2);}
		}
		
		//return collection
		return rtrn;
	}

	private static boolean overlap(Annotation assembly, Annotation read) {
		return assembly.overlaps(read,forceStrandSpecificity);
	}

	private boolean overlapsWithoutCrossingIntron(Annotation assembly, Annotation read) {
		//one has introns, one does not
		//see if one without overlaps the exons of the other
		//AND does not overlap the intron
		if(!read.getSpliceConnections().isEmpty()){
			//read has introns
			//test if exons overlap the intron
			for(Annotation intron: read.getSpliceConnections()){
				if(assembly.overlaps(intron)){return false;}
			}
			//if not, test that exons overlap
			for(Annotation exon: read.getBlocks()){
				if(assembly.overlaps(exon)){return true;}
			}
			return false;
		}
		
		//assembly has intron
		Annotation complement=assembly.complement();
		for(Annotation intron: complement.getBlocks()){
			if(read.overlaps(intron)){return false;}
		}
		for(Annotation exon: assembly.getBlocks()){
			if(read.overlaps(exon)){return true;}
		}
		
		return false;
	}

	private boolean hasNonOverlapping(Annotation assembly, Alignment read) {
		//check if every intron in read is also an intron in assembly and vice versa
		
		Collection<? extends Annotation> assemblyIntrons=assembly.complement().getBlocks();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		for(Annotation intron: assemblyIntrons){
			if(!readIntrons.contains(intron)){return true;}
		}
		
		for(Annotation intron: readIntrons){
			if(!assemblyIntrons.contains(intron)){return true;}
		}
		
		return false;
	}

	private Annotation merge(Annotation assembly, Alignment read) {
		//since these are compatable, we will simply add the read to the assembly by merging exons
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		Collection<Annotation> readAlignments=read.getReadAlignments(this.space);
		for(Annotation align: readAlignments){
			rtrn.addAll(align.getBlocks());
		}
		rtrn.addAll(assembly.getBlocks());
		return new BasicAnnotation(rtrn);
	}

	private Assembly merge(Annotation assembly, Annotation read) {
		//since these are compatable, we will simply add the read to the assembly by merging exons
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		rtrn.addAll(read.getBlocks());
		rtrn.addAll(assembly.getBlocks());
		return new Assembly(rtrn);

	}
	

	private void write(Collection<? extends Annotation> blocks, String save) throws IOException {
		write(blocks, save, true);
	}
		

	private void write(Collection<? extends Annotation> blocks, String save, boolean append) throws IOException {
		FileWriter writer=new FileWriter(save, append);
		
		for(Annotation block: blocks){writer.write(block+"\n");}
		
		writer.close();
	}

	private Pair<Collection<Annotation>> splitExonsByJunctions(Annotation exonicRegion, Collection<Annotation> junctions) {
		Collection<Annotation> spliced=new TreeSet<Annotation>();
		Collection<Annotation> preSpliced=new TreeSet<Annotation>();
		
		//make junction tree
		IntervalTree<Annotation> junctionTree=makeTree(junctions);
		
		//go through exons and see if it overlaps a junction
		for(Annotation exon: exonicRegion.getBlocks()){
			//if overlaps a junction
			Iterator<Node<Annotation>> overlappers=junctionTree.overlappers(exon.getStart()-1, exon.getEnd()+1);
			Pair<Collection<Annotation>> exons=split(exon, overlappers);
			spliced.addAll(exons.getValue1());
			preSpliced.addAll(exons.getValue2());
		}
		
		Pair<Collection<Annotation>> rtrn=new Pair<Collection<Annotation>>(spliced, preSpliced);
		return rtrn;
	}

	private IntervalTree<Annotation> makeTree(Collection<Annotation> junctions) {
		IntervalTree<Annotation> rtrn=new IntervalTree<Annotation>();
		for(Annotation junction: junctions){
			rtrn.put(junction.getStart(), junction.getEnd(), junction);
		}
		return rtrn;
	}

	private Pair<Collection<Annotation>> split(Annotation exon, Iterator<Node<Annotation>> overlappers) {
		Collection<Annotation> spliced=new TreeSet<Annotation>();
		Collection<Annotation> preSpliced=new TreeSet<Annotation>();
		
		
		List<Annotation> introns=new ArrayList<Annotation>();
		
		while(overlappers.hasNext()){
			introns.add(overlappers.next().getValue());
		}
		
		if(introns.isEmpty()){
			spliced.add(exon);
		}
		else{
			boolean counted=false;
			//define overlapping introns
			List<List<Annotation>> overlappingIntrons=defineOverlappingIntrons(introns); //TODO Why would any list be empty?
			//create sets of introns
			Collection<Annotation> setOfIntrons=createSetsOfIntrons(overlappingIntrons);
			//subtract each set and add the resulting exons
			for(Annotation setOfIntron: setOfIntrons){
				Annotation minus=exon.minus(setOfIntron);
				if(!minus.getBlocks().isEmpty()){
					spliced.addAll(minus.getBlocks()); //These are legit
					counted=true;
				}
			}
			
			//We'll add all the complements as possible pre-mRNA background to test
			if(!spliced.isEmpty()){
				Annotation comp=exon.minus(spliced);
				preSpliced.addAll(comp.getBlocks()); //These are likely preprocessed mRNA
			}
						
			
			if(!counted && this.alignsWithEdge(exon, introns)){
				logger.debug("aligns with edge "+exon.toUCSC());
				spliced.add(exon); //These are likely legit
				counted=true;
			}
			
			//We'll add all exons fully contained in an intron as possible pre-mRNA background
			//TODO We may just want to add all exons not already counted
			if(!counted && fullyContained(exon, introns)){preSpliced.add(exon);} //These are likely preprocessed mRNA	
		}
			
		Pair<Collection<Annotation>> rtrn=new Pair<Collection<Annotation>>(spliced, preSpliced);
		return rtrn;
	}

	private boolean notOverlapping(Annotation exon, Collection<Annotation> spliced) {
		for(Annotation splic: spliced){
			if(exon.overlaps(splic)){
				logger.debug(exon.toUCSC()+" overlaps");
				return false;
			}
		}
		return true;
	}

	private Collection<Annotation> createSetsOfIntrons(List<List<Annotation>> overlappingIntrons) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		//go through the overlapping introns and create sets
		int[] sizes=new int[overlappingIntrons.size()];
		for(int i=0; i<overlappingIntrons.size(); i++){
			sizes[i]=overlappingIntrons.get(i).size();
		}
		
		CombinationGenerator combinations=new CombinationGenerator(sizes);
		int[][] possibleCombos=combinations.getPermutations();
		
		//Rows are the combos to use
		for(int i=0; i<possibleCombos.length; i++){
			int[] setPositions=possibleCombos[i];
			Annotation intronSet=makeIntronSet(overlappingIntrons, setPositions);
			rtrn.add(intronSet);
		}
		
		return rtrn;
	}

	private Annotation makeIntronSet(List<List<Annotation>> overlappingIntrons,	int[] setPositions) {
		Collection<Annotation> blocks=new TreeSet<Annotation>();
		for(int i=0; i<setPositions.length; i++){
			Annotation intron=overlappingIntrons.get(i).get(setPositions[i]);
			blocks.add(intron);
		}
		return new BasicAnnotation(blocks);
	}

	private void write(List<List<Annotation>> overlappingIntrons, String save, Annotation exon) throws IOException {
		FileWriter writer=new FileWriter(save, true);
		
		int i=0;
		for(List<Annotation> introns: overlappingIntrons){
			for(Annotation intron: introns){
				//intron.setName("set"+i);
				writer.write(intron.getChr()+"\t"+intron.getStart()+"\t"+intron.getEnd()+"\t"+("set:"+i+"-"+exon.toUCSC())+"\n");
			}
			i++;
		}
		
		writer.close();
	}

	private List<List<Annotation>> defineOverlappingIntrons(List<Annotation> introns) {
		//TODO This does not account for exon exclusion events
		//To account for this, we will want to include introns that overlap both sets to be included in both sets
		
		List<List<Annotation>> rtrn=new ArrayList<List<Annotation>>();
		IntervalTree<Annotation> intronTree=makeTree(introns);
		Set<Annotation> accountedFor=new TreeSet<Annotation>();
		
		for(Annotation intron: introns){
			if(!accountedFor.contains(intron)){
				List<Annotation> overlappingIntrons=new ArrayList<Annotation>();
				//get overlappers
				List<Annotation> overlappers=getList(intronTree.overlappers(intron.getStart(), intron.getEnd()));
				add(overlappingIntrons, overlappers, accountedFor);
				for(Annotation ann2: overlappers){
					List<Annotation> overlappers2=getList(intronTree.overlappers(ann2.getStart(), ann2.getEnd()));
					add(overlappingIntrons, overlappers2, accountedFor);
				}
				rtrn.add(overlappingIntrons);
			}
		}
		
		//Lets do a second split to take each list and divide up bins that overlap
		rtrn=splitOverlapping(rtrn);
		
		return rtrn;
	}

	private List<List<Annotation>> splitOverlapping(List<List<Annotation>> set) {
		List<List<Annotation>> rtrn=new ArrayList<List<Annotation>>();
		//Lets go through each list
		//If the reads within it have junctions that dont overlap each other then lets split into new bins
		for(List<Annotation> list: set){
			rtrn.addAll(splitOverlappingForList(list));
		}
		return rtrn;
	}


	private List<List<Annotation>> splitOverlappingForList(List<Annotation> list) {
		//take this list and make sure that every element overlaps every other element
		//else split into n bins of overlapping bins
		
		
		Collection<Annotation> newList=new TreeSet<Annotation>();
		
		//List<Annotation> originalList=new ArrayList<Annotation>();
		
		for(int i=0; i<list.size(); i++){
			for(int j=(i+1); j<list.size(); j++){
				Annotation intron1=list.get(i);
				Annotation intron2=list.get(j);
				//if intron 1 and intron 2 don't overlap then split
				if(!intron1.overlaps(intron2)){
					//logger.debug(intron1.toUCSC()+" "+intron2.toUCSC()+" dont overlap");
					newList.add(intron1);
					newList.add(intron2);
				}
			}
		}
		
		//partition the new list
		List<List<Annotation>> newLists=partition(newList, list);
		
		return newLists;
	}

	private List<List<Annotation>> partition(Collection<Annotation> newList, Collection<Annotation> originalList) {
		Collection<Annotation> accountedFor=new TreeSet<Annotation>();
		
		List<List<Annotation>> lists=new ArrayList<List<Annotation>>();
		Map<Annotation, Collection<Annotation>> rtrn=new TreeMap<Annotation, Collection<Annotation>>();
		for(Annotation intron1: newList){
			Collection<Annotation> overlaps=new TreeSet<Annotation>();
			for(Annotation intron2: originalList){
				if(intron1.overlaps(intron2)){
					overlaps.add(intron2);
				}
			}
			rtrn.put(intron1, overlaps);
		}
		
		for(Annotation ref: rtrn.keySet()){
			boolean fullyCompat=isFullyCompatible(rtrn.get(ref));
			if(fullyCompat){
				if(!rtrn.get(ref).isEmpty()){
					lists.add(new ArrayList(rtrn.get(ref))); //TODO We should figure out a way to exclude the exact same lists from being added multiple times
				}
				accountedFor.addAll(rtrn.get(ref));
				//logger.debug(ref.toUCSC()+" "+rtrn.get(ref).size());
				//for(Annotation temp: rtrn.get(ref)){logger.debug(temp.toUCSC());}
			}
		}
		//TODO Make sure every intron is accounted for
		for(Annotation newIntron: newList){
			if(!accountedFor.contains(newIntron)){
				logger.debug("MISSING: "+newIntron.toUCSC());
			}
		}
		
		
		//Add the remainder of the originalList as a separate list
		List<Annotation> remainderList=new ArrayList<Annotation>();
		for(Annotation original: originalList){
			if(!accountedFor.contains(original)){remainderList.add(original);}
		}
		
		if(!remainderList.isEmpty()){
			lists.add(remainderList);
		}
		
		return lists;
	}

	private boolean isFullyCompatible(Collection<Annotation> collection) {
		//go through each annotation and ensure that nothing doesnt overlap with something else
		for(Annotation an1: collection){
			for(Annotation an2: collection){
				if(!an1.overlaps(an2)){return false;}
			}
		}
		return true;
	}

	private void allOverlap(List<Annotation> newList, List<Annotation> originalList, Annotation reference) {
		//Check that all overlap in list
		//if not, exclude from newList
		//if so, exclude from original list
		
		for(int i=0; i<newList.size(); i++){
			for(int j=(i+1); j<newList.size(); j++){
				Annotation intron1=newList.get(i);
				Annotation intron2=newList.get(j);
				if(!intron1.overlaps(intron2)){
					throw new IllegalStateException();
				}
				else{
					originalList.remove(intron1);
					originalList.remove(intron2);
				}
			}
		}
		
	}

	private List<Annotation> getList(Iterator<Node<Annotation>> overlappers) {
		List<Annotation> rtrn=new ArrayList<Annotation>();
		
		while(overlappers.hasNext()){
			rtrn.add(overlappers.next().getValue());
		}
		
		return rtrn;
	}

	private void add(List<Annotation> overlappingIntrons, List<Annotation> overlappers, Set<Annotation> accountedFor) {
		for(Annotation intron: overlappers){
			if(!overlappingIntrons.contains(intron)){
				overlappingIntrons.add(intron);
			}
			if(!accountedFor.contains(intron)){
				accountedFor.add(intron);
			}
		}
		
	}

	private Collection<Annotation> getRemainingIntrons(List<Annotation> introns, Collection<Annotation> exons) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//TODO Make more efficient
		//go through the exon set and see what introns are fully explained
		for(Annotation intron: introns){
			if(!hasExons(intron, exons)){
				rtrn.add(intron);
			}
		}
		
		return rtrn;
	}

	private boolean hasExons(Annotation intron, Collection<Annotation> exons) {
		boolean left=false;
		boolean right=false;
		
		for(Annotation exon: exons){
			if(exon.getStart()==intron.getEnd()){left=true;}
			if(exon.getEnd()==intron.getStart()){right=true;}
		}
		
		return left && right;
	}

	private Collection<? extends Annotation> dissect(Annotation exon, List<Annotation> introns) {
		logger.debug("dissecting "+exon+" "+introns.size());
		return exon.intersect(introns);
	}

	private boolean spans1Intron(Annotation exon, List<Annotation> introns) {
		Annotation block=new BasicAnnotation(introns);
		
		int counter=0;
		for(Annotation intron: block.getBlocks()){
			if(exon.overlaps(intron)){counter++;}
		}
		
		if(counter==1){return true;}
		
		return false;
		
	}

	private boolean alignsWithEdge(Annotation exon, List<Annotation> introns) {
		for(Annotation intron: introns){
			if(exon.getEnd()== intron.getStart() || exon.getStart()==intron.getEnd()){return true;}
		}
		return false;
	}

	private boolean fullyContained(Annotation exon, List<Annotation> introns) {
		for(Annotation intron: introns){
			if(intron.fullyContains(exon)){
				return true;
			}
		}
		return false;
	}

	private Collection<Annotation> intersect(List<Annotation> introns) {
		Collection<Annotation> rtrn=CollapseByIntersection.collapseByIntersection(introns, true);
		return rtrn;
	}

	private Collection<? extends Annotation> splitByIntron(Annotation exon, Iterator<Node<Annotation>> overlappers) {
		List<Annotation> introns=new ArrayList<Annotation>();
		while(overlappers.hasNext()){introns.add(overlappers.next().getValue());}
		
		//TODO minus the consensus intronic region
		
		return exon.disect(introns);
		//return exon.minus(introns).getBlocks();
	}

	private Collection<? extends Annotation> getExons(Alignment read) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		Collection<Annotation> readPairs=read.getReadAlignments(space); // this should return a collection for left and right
		for(Annotation align: readPairs){
			rtrn.addAll(align.getBlocks());
		}
		return rtrn;
	}
	
	private Collection<? extends Annotation> getJunctions(Alignment read) {
		//TODO This should only get spliced reads NOT indels
			
		return read.getSpliceConnections();
	}

	public Map<String, Collection<Gene>> getPaths() {
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		
		//For each chromosome
		for(String chr: this.graphs.keySet()){
			List<GraphPath<Annotation, TranscriptGraphEdge>> paths=this.graphs.get(chr).getPaths();
			Collection<Gene> genes=new TreeSet<Gene>();
			for(GraphPath<Annotation, TranscriptGraphEdge> path: paths){
				//logger.debug(path.toString());
				Gene gene=this.graphs.get(chr).pathToGene(path);
				genes.add(gene);
			}
			rtrn.put(chr, genes);
		}
		//Add orphan genes
		for(String chr: this.graphs.keySet()){
			for(Gene g:this.graphs.get(chr).getOrphanGenes()){
				rtrn.get(chr).add(g);
			}			
		}
		return rtrn;
	}
	
	/**
	 * Converts a list to an array
	 * @param list
	 * @return
	 */
	private double[] l2a(List<Double> list){
		double[] rtrn=new double[list.size()];
	
		int i=0;
		for(Double val: list){rtrn[i++]=val;}
	
		return rtrn;
	}
	
	/**
	 * Converts a collection of paths to a collection of genes
	 * @param paths to be converted
	 * @return
	 */
	private Collection<Gene> convert(Collection<Path> paths) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		for(Path path: paths){
			rtrn.add(path.toGene());
		}
		
		return rtrn;
	}

	/**
	 * This function will go through each gene and calculate the number of paired end reads 
	 * fully contained within each isoform for the gene
	 * @return
	 * @throws IOException 
	 */
	private Map<String,Collection<Gene>> setFPKMScores(Map<String,Collection<Gene>> geneMap,FileWriter writer,FileWriter writer2){
		
		Map<String,Collection<Gene>> filteredGenes = new HashMap<String,Collection<Gene>>();
		String name = "gene.v2.1_";
		for(String chr:geneMap.keySet()){
			//logger.debug("For chromosome "+chr+" graph gave "+geneMap.get(chr).size()+" genes");
			Collection<Gene> filtered = new ArrayList<Gene>();
			try{
			//MAKE A MAP OF GENE TO ISOFORMS
			Map<Gene,Set<Gene>> isoformMap = getIsoformMap(geneMap.get(chr));
			Map<Gene,Double> geneToPairedCount = new HashMap<Gene,Double>();
			
			//For each gene
			for(Gene gene:isoformMap.keySet()){
				writer2.write("Gene:\n");
				//logger.debug("Starting new gene");
				//For each transcript
				for(Gene isoform:isoformMap.get(gene)){
					
					double[] scores = getScores(isoform);
					double[] fields = new double[4];
					logger.debug(gene.toUCSC()+" "+scores[0]+"\t"+scores[1]);
					if(scores[1]<alpha){
						isoform.setName(name+new Double(counter).toString()+"_"+isoform.getChr());
						//[0] : sum
						fields[0] = scores[0];
						//[1] : p-value
						fields[1] = scores[1];
						//[2] : FPK
						fields[2] = (scores[0]*1000.0)/isoform.getSize();
						//[3] : FPKM
						//Calculate FPKM
						fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
						//logger.debug("For isoform : "+isoform.getName()+"\tNum of exons: "+isoform.getSpliceConnections().size()+"\t"+fields[0]+"\t"+fields[1]);
						isoform.setBedScore(fields[3]);
						isoform.setExtraFields(fields);
						counter++;
						writer.write(isoform+"\n");
						writer2.write(isoform.toBED()+"\t"+fields[0]+"\n");
						geneToPairedCount.put(isoform, fields[0]);
						filtered.add(isoform);
					}
					else{
						logger.debug("Gene "+gene.toUCSC()+" is filtered out because it does not meet the significance threshold : "+scores[1]);
					}
				}
			}
			} catch (IOException e) {
				e.printStackTrace();
			}
			logger.info("After applyPairedEndFilter "+filtered.size()+" genes");
			filteredGenes.put(chr, filtered);
		}
		return filteredGenes;
		
	}

	/**
	 * 
	 * @param genes
	 * @return
	 */
	public static Map<Gene,Set<Gene>> getIsoformMap(Collection<Gene> genes){
				
		Map<Gene,Set<Gene>> isoformMap = new HashMap<Gene,Set<Gene>>();
		
		for(Gene g:genes){
			//If gene has not already been processed
			if(!isoformMapContains(isoformMap,g)){
				//Add gene to its own isoform map
				Set<Gene> set = new HashSet<Gene>();
				set.add(g);
				Collection<Gene> potentialIsoforms = getOverlappingGenes(g,genes,0.5);
				for(Gene p:potentialIsoforms){
					if(!isoformMapContains(isoformMap,p)){
						set.add(p);
					}
				}
				isoformMap.put(g, set);
			}			
		}
		return isoformMap;
	}
	
	public static boolean isoformMapContains(Map<Gene,Set<Gene>> isoformMap,Gene gene){
	
		for(Set<Gene> set:isoformMap.values()){
			for(Gene g:set)
				if(g.equals(gene))
					return true;
		}
		return false;
	}
	/**
	 * 
	 * @param gene
	 * @param allGenes
	 * @param minPctOverlap
	 * @return
	 */
	public static Collection<Gene> getOverlappingGenes(Gene gene,Collection<Gene> allGenes, double minPctOverlap){
		
		Collection<Gene> overlappers = new HashSet<Gene>();
		for(Gene g:allGenes){
			if(gene.overlaps(g, minPctOverlap) && g.getOrientation().equals(gene.getOrientation()))
				overlappers.add(g);
		}
		return overlappers;
	}
	
	public static void main(String[] args)throws IOException{
		
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"reconstruct");
		
		TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand").equalsIgnoreCase("first")){
			//System.out.println("First read");
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(argMap.get("strand").equalsIgnoreCase("second")){
			//System.out.println("Second read");
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}
		else
			logger.info("no strand");
		
		if(argMap.containsKey("chr")){
			new BuildScriptureCoordinateSpace(new File(argMap.getMandatory("alignment")),argMap.getMandatory("genome"),argMap.getOutput(),true, strand,argMap.getMandatory("chr"),argMap);
		}
		else{
			new BuildScriptureCoordinateSpace(new File(argMap.getMandatory("alignment")),argMap.getMandatory("genome"),argMap.getOutput(),true, strand,argMap);
		}
/*		if(args.length>1){
			File bamFile=new File(args[0]);
			String genomeSeqFile = null;
			double threshold = new Double(args[1]);
			genomeSeqFile = args[2];
			if(args.length==3)
				new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3]);
			if(args.length==3)
				new BuildScriptureCoordinateSpace(bamFile,threshold,null,args[2],true, TranscriptionRead.SECOND_OF_PAIR);
			if(args.length>=5){
				TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
				if(args[4].equalsIgnoreCase("first")){
					//System.out.println("First read");
					strand = TranscriptionRead.FIRST_OF_PAIR;
				}
				else if(args[4].equalsIgnoreCase("second")){
					//System.out.println("Second read");
					strand = TranscriptionRead.SECOND_OF_PAIR;
				}
				else
					System.out.println("no strand");
				if(args.length==5){
					new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3],true, strand);
				}
				else{
					new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3],true, strand,args[5]);
				}
			}
		}
		else{System.err.println(usage);}*/
	}
	
	static final String usage = "Usage: BuildScriptureCoordinateSpace -task reconstruct "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Alignment file to be used for reconstruction.> "+
			"\n\n\t\t-genome <Fasta file with Genome sequence to be used as reference.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+		
			"\n\t\t-strand <VALUES: first, second, unstranded. Specifies the mate that is in the direction of transcription DEFAULT: Unstranded> "+
			
			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\n\t\t-chr <If specified, Scripture will be run for this chromosome only.> "+
			"\n\t\t-coverage <Specifies the minimum percentage of drop in coverage allowed for an exon. DEFAULT: 0.2> "+
			"\n\t\t-minSpliceReads <The minimum number of splice reads allowed to support a single intron transcript. DEFAULT: 3> "+
			"\n\t\t-percentSpliceReads <The minimum percentage of the average splice counts for a transcript, that an intron can be supported by. DEFAULT: 0.05> "+
			"\n\t\t-alpha <The significance p-value threshold for reconstructions. DEFAULT: 0.01> "+
			"\n";
	
	//static String usage=" args[0]=bam file \n\t args[1]=minimum percentage threshold for coverage \n\t args[2]: Fasta file with Genome sequence"
	//					+"\n\t args[3]= outputName \n\t args[4] transcription strand \n\targs[5] If specified only this chromosome";
	
}
package nextgen.core.scripture;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.jgrapht.GraphPath;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.Pair;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CollapseByIntersection;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.Path;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.JCSAlignmentModel;
import nextgen.core.readFilters.CanonicalSpliceFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.IndelFilter;
import nextgen.core.readFilters.PairedAndProperFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.readFilters.SameOrientationFilter;
import nextgen.core.readFilters.SplicedReadFilter;
import nextgen.core.scripture.OrientedChromosomeTranscriptGraph.TranscriptGraphEdge;
import nextgen.core.scripture.statistics.ConnectDisconnectedTranscripts;

public class BuildScriptureCoordinateSpace {

	private static final TranscriptionRead DEFAULT_TXN_READ =  TranscriptionRead.UNSTRANDED;
	static Logger logger = Logger.getLogger(BuildScriptureCoordinateSpace.class.getName());
	private JCSAlignmentModel model;
	String genomeSeq = null;
	int windowSize=20000000;
	private CoordinateSpace space;
	private Map<String, ChromosomeTranscriptGraph> graphs;
	private static boolean forceStrandSpecificity=true; //TODO This should be passed or at least determined from data
	private static double DEFAULT_MIN_COV_THRESHOLD = 0.2;
	private static double MIN_SPLICE_PERCENT = 0.05;
	private double coveragePercentThreshold = DEFAULT_MIN_COV_THRESHOLD;
	//private static TranscriptionRead DEFAULT_TXN_READ = TranscriptionRead.UNSTRANDED;
	String outName = null;
	private double DEFAULT_ALPHA = 0.01;
	private double alpha=DEFAULT_ALPHA;
	private static double MIN_SPLICE_READS = 3.0;
	private double THRESHOLD_SPURIOUS = 0.95;
	private double minSpliceReads = MIN_SPLICE_READS;
	private double minSplicePercent = MIN_SPLICE_PERCENT;
	int counter = 1000;
	int globalCounter = 1000;
	File bamFileName;
	private int constant = 10000;
	private double globalPairedLambda=0.0;
	File bamfile;
	//double globalFragments;
	
/*	public BuildScriptureCoordinateSpace(File bamFile){
		this(bamFile,DEFAULT_MIN_COV_THRESHOLD,null,null,true,DEFAULT_TXN_READ);
		logger.info("Genome sequence has not been provided");
	}
	
	
	public BuildScriptureCoordinateSpace(File bamFile,double threshold,String genomeDir,String outputName){
		//By default first read is transcription read
		this(bamFile,threshold,genomeDir,outputName,true,DEFAULT_TXN_READ);
	}
	*/
	
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir){
		this(bamFile,genomeDir,bamFile.getName()+".reconstructions",true,DEFAULT_TXN_READ,null);
	}
	/**
	 * 
	 * @param bamFile
	 * @param threshold
	 * @param genomeDir
	 * @param outputName
	 * @param forceStrandedness
	 */
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead strand,ArgumentMap argMap){
			
		bamfile=bamFile;
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		bamFileName = bamFile;
		forceStrandSpecificity = forceStrandedness;
		outName = outputName;
		model=new JCSAlignmentModel(bamfile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		space=model.getCoordinateSpace();

		setThresholds(argMap);
//		globalFragments = calculateGlobalFragments();
		logger.info("Parameters used: " +
				"\nIntron intention Filter: "+THRESHOLD_SPURIOUS+
				"\nPremature Assembly Filter: "+coveragePercentThreshold+
				"\nSplice junction Filter : "+
				"\n\tNumber of spliced reads : "+minSpliceReads+
				"\n\tPercentage of total spliced reads: "+minSplicePercent+
				"\nAlpha for single exon assemblies : "+alpha
				);
		
		assemble(strand);
		
		Map<String, Collection<Gene>> rtrn=getPaths();
/*		try {
			FileWriter writer=new FileWriter(outName+".08graph.all.paths.bed");
			for(String chr: rtrn.keySet()){
				for(Gene g:rtrn.get(chr)){
					writer.write(g.toBED()+"\n");
				}
			}			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}*/
		
		try {
			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			Map<String, Collection<Gene>> genes = setFPKMScores(rtrn,writer,writer2);
			writer.close();
			writer2.close();

//			write(outName+".10final.paths.bed",genes);
			
			postProcess(genes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	} 
	
	/**
	 * 
	 * @param bamFile
	 * @param threshold
	 * @param genomeDir
	 * @param outputName
	 * @param forceStrandedness
	 * @param chr
	 */
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead strand,String chr,ArgumentMap argMap){
			
		bamfile=bamFile;
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		bamFileName = bamFile;
		forceStrandSpecificity = forceStrandedness;

		outName = outputName;
		
		setThresholds(argMap);
		//assemble(strand);
		logger.info("Parameters used: " +
				"\nIntron Retention Filter: "+THRESHOLD_SPURIOUS+
				"\nPremature Assembly Filter: "+coveragePercentThreshold+
				"\nSplice junction Filter : "+
				"\n\tNumber of spliced reads : "+minSpliceReads+
				"\n\tPercentage of total spliced reads: "+minSplicePercent+
				"\nAlpha for single exon assemblies : "+alpha
				);
		
		ChromosomeTranscriptGraph graph=assemble(chr,strand);
		graphs.put(chr, graph);
		
		Map<String, Collection<Gene>> rtrn=getPaths();
/*		try {
			FileWriter writer=new FileWriter(outName+".08graph.all.paths.bed");
			for(Gene g:rtrn.get(chr)){
				writer.write(g.toBED()+"\n");
			}			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
*/		
		try {
			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			Map<String, Collection<Gene>> genes = setFPKMScores(rtrn,writer,writer2);
			writer.close();
			writer2.close();

//			write(outName+".10final.paths.bed",genes);
			
			postProcess(genes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	} 
	
	private void setThresholds(ArgumentMap argMap){ 
		
		if(argMap!=null){
			coveragePercentThreshold = argMap.getDouble("coverage", DEFAULT_MIN_COV_THRESHOLD);
			alpha = argMap.getDouble("alpha", DEFAULT_ALPHA);
			minSpliceReads = argMap.getDouble("minSpliceReads", MIN_SPLICE_READS);
			minSplicePercent = argMap.getDouble("percentSpliceReads", MIN_SPLICE_PERCENT);
		}
	}

	private void assemble(TranscriptionRead strand) {
		//Iterate over all chromosomes
		for(String chr: space.getReferenceNames()){
			logger.info("Reference name: "+chr);
			if(model.getRefSequenceLambda(chr)==0.0){
				logger.info(chr+" is not expressed in the alignment file");
			}
			else{
				ChromosomeTranscriptGraph graph=assemble(chr,strand);
				this.graphs.put(chr, graph);
			}
		}
	}

	private void postProcess(Map<String, Collection<Gene>> oldGenes) throws IOException{
		
		Map<String, Collection<Gene>> newGenes = new HashMap<String, Collection<Gene>>();
		for(String chr:oldGenes.keySet()){
			//For all paths on this chromosome
			
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:oldGenes.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			/*
			 * MERGE THE GENES
			 */
			//Step 1: Iterate through the working assembly and find paths that overlap
			//if overlaps but is incompatible, see if we can branch
			//iterate
			Iterator<Gene> iter=tree.toCollection().iterator();			

			Collection<Gene> considered = new HashSet<Gene>();
	 		logger.info("Enter merging");
			while(iter.hasNext()){
				Gene branch1=iter.next();
				considered.add(branch1);
				//get overlapping branches
				Iterator<Gene> overlappers=tree.overlappingValueIterator(branch1.getStart(), branch1.getEnd());
				//Collection<Assembly> toRemove=new TreeSet<Assembly>();
				while(overlappers.hasNext()){
					Gene branch2=overlappers.next();
					if(!considered.contains(branch2)){
						if(!branch1.equals(branch2) && compatible(branch1, branch2)){
							logger.debug("Merging: "+ branch1.getName()+ " and "+branch2.getName());
							Collection<Annotation> rtrn=new TreeSet<Annotation>();
							rtrn.addAll(branch1.getBlocks());
							rtrn.addAll(branch2.getBlocks());
							Gene merged=new Gene(rtrn);
							merged.setName(branch1.getName());
							//SET THE SCORE TO THE MAX FPKM
							merged.setBedScore(Math.max(branch1.getBedScore(),branch2.getBedScore()));
							//remove annotation1 and annotation2
							tree.remove(branch1.getStart(), branch1.getEnd(), branch1);
							tree.remove(branch2.getStart(), branch2.getEnd(), branch2);
							//add merged
							tree.put(merged.getStart(), merged.getEnd(), merged);
						}
					}
				}
			}
			
			/*
			 * FPKM THRESHOLD
			 */
			//MAKE A MAP OF GENE TO ISOFORMS
			Map<Gene,Set<Gene>> isoformMap = getIsoformMap(tree.toCollection());		
			Set<Gene> finalSet = new HashSet<Gene>();
			//For each gene
			for(Gene gene:isoformMap.keySet()){
				//For each transcript
				double max = Double.MIN_VALUE;
				for(Gene isoform:isoformMap.get(gene)){					
					double score = isoform.getBedScore();
					if(score>max){
						max = score;
					}
				}
				//CHeck if the FPKM of each is at least 10% of the mac FPKM
				for(Gene isoform:isoformMap.get(gene)){		
					double score = isoform.getBedScore();
					if(score>0.1*max && score>0.0){
						//filtered.remove(isoform.getStart(), isoform.getEnd(), isoform);
						finalSet.add(isoform);
					}
				}
			}
			newGenes.put(chr, finalSet);
		}
	
		//write(outName+".11postprocessed.paths.bed",newGenes);
		write(outName+".scripture.paths.bed",newGenes);
		connectDisconnectedTranscripts(newGenes);
	}
	
	
	private ChromosomeTranscriptGraph assemble(String chr,TranscriptionRead strand) {
		
		//Option 1: Scan through the space and collapse into compatible edges and nodes
		ChromosomeTranscriptGraph graph=assembleDirectly(chr,strand);
		
		return graph;
	}

	/**
	 * This function connects disconnected reconstructions using paired end reads.
	 * @param annotations
	 * @throws IOException 
	 */
	private void connectDisconnectedTranscripts(Map<String, Collection<Gene>> annotations) throws IOException{
		
		model.addFilter(new PairedAndProperFilter());
		int loop=0;
		boolean somethingWasConnected = true;
		double medianInsertSize=0.0; 
		Map<String,Collection<Gene>> temp = new HashMap<String,Collection<Gene>>();
		for(String chr:annotations.keySet()){
			if(!annotations.get(chr).isEmpty()){
				temp.put(chr, annotations.get(chr));
			}
		}
		medianInsertSize += model.getReadSizeDistribution(new TranscriptomeSpace(temp), 800, 100).getMedianOfAllDataValues();
		//double medianInsertSize = 600;
		logger.info("Median size = "+medianInsertSize);

		Map<String,Collection<Gene>> conn = null;
		
		while(somethingWasConnected || loop<10){
			somethingWasConnected =false;
			loop++;
			
			logger.debug("Connected disconnected transcripts: Loop "+loop);
			conn = new HashMap<String,Collection<Gene>>();
			for(String chr:annotations.keySet()){
				conn.put(chr, new TreeSet<Gene>());			
			}
			
			for(String chr:annotations.keySet()){
				//For all genes on this chromosome
				logger.debug("Connecting, Processing "+chr);
				Collection<Gene> newGenes = new TreeSet<Gene>();
				//MAKE AN INTERVAL TREE OF THE GENES on this chr
				IntervalTree<Gene> tree = new IntervalTree<Gene>();
				for(Gene g:annotations.get(chr)){
					conn.get(chr).add(g);
					tree.put(g.getStart(), g.getEnd(), g);
				}
				//For each transcript
				//Iterate over all reconstructions
				Iterator<Gene> iter=tree.toCollection().iterator();
				while(iter.hasNext()){
					Gene gene=iter.next();
					//For all assemblies downstream of this assembly in 10kB regions
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getEnd(), gene.getEnd()+constant);
					
					newGenes.add(gene);
					while(overlappers.hasNext()){
						Gene other = overlappers.next();
						if(isCandidate(gene,other)){
							if(pairedEndReadSpansTranscripts(gene, other)){ 
								//if(secondTranscriptIsSingleExon(gene,other)){
									logger.debug("Attempt to connect "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC());
									
									//Connect the genes
									Annotation connected = getConnectedTranscript(gene,other,medianInsertSize);
									if(connected!=null){
										somethingWasConnected = true;
										Gene newConnected = new Gene(connected);
										double[] scores = getScores(newConnected);
										double[] fields = new double[4];
										newConnected.setName(gene.getName()+"_"+other.getName());
										//[0] : sum
										fields[0] = scores[0];
										//[1] : p-value
										fields[1] = scores[1];
										//[2] : FPK
										fields[2] = (scores[0]*1000.0)/newConnected.getSize();
										//[3] : FPKM
										//Calculate FPKM
										fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
										//logger.debug("For isoform : "+isoform.getName()+"\tNum of exons: "+isoform.getSpliceConnections().size()+"\t"+fields[0]+"\t"+fields[1]);
										newConnected.setBedScore(fields[3]);
										newConnected.setExtraFields(fields);
										newGenes.remove(gene);
										newGenes.add(newConnected);
										conn.get(chr).remove(gene);
										conn.get(chr).remove(other);
										conn.get(chr).add(newConnected);
									}
								//}
							}
							else{
								//logger.info("The genes "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC()+" do not have paired reads");
							}
						}
					}				
				}
			}
			annotations = conn;
		}
		
		//FileWriter bw = new FileWriter(outName+".12connected.bed");
		FileWriter bw = new FileWriter(outName+".connected.bed");
		for(String name:conn.keySet()){
			Iterator<Gene> ter = conn.get(name).iterator();
			while(ter.hasNext()){
				Gene isoform = ter.next();
				Gene iso = new Gene(trimEnds(isoform,0.1));
				bw.write(iso.toBED()+"\n");
			}
		}
		bw.close();
	}
	
	/**
	 * Connect the two genes by fusing the last exon of the first and the first exon of the last gene
	 * @param gene
	 * @param other
	 * @return
	 * @throws IOException 
	 */
	private Annotation getConnectedTranscript(Gene gene,Gene other,double medianInsertSize) throws IOException{
		
		Pair<Gene> orderedGenes = ConnectDisconnectedTranscripts.getOrderedAssembly(gene, other);
		Annotation connected = null;
		/*
		 * If the distance between the transcripts is less than the insert size,
		 * paste through
		 */
		if(orderedGenes.getValue2().getStart()-orderedGenes.getValue1().getEnd() < medianInsertSize){
			Gene firstGene = orderedGenes.getValue1().copy();
			firstGene.setEnd(orderedGenes.getValue2().getStart());
			connected = firstGene.union(orderedGenes.getValue2());
		}
		/*
		 * If the distance between the transcripts is more than the insert size,
		 * 		if there is at least 1 splice read spanning the junction,
		 * 			union it
		 */
		else{
			boolean flag= false;
			Annotation junction = orderedGenes.getValue1().getLastExon().union(orderedGenes.getValue2().getFirstExon());
			CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(junction,false), new SplicedReadFilter());
			while(splicedIter.hasNext()){
				Alignment read = splicedIter.next();
				//if there is at least 1 splice read spanning the junction,
				if(BuildScriptureCoordinateSpace.compatible(read,junction)){
					flag=true;
					break;
				}
			}
			splicedIter.close();
			if(flag){
				connected = gene.union(other);				
			}
			/*
			 * Else NO CONNECTION
			 */
		}
		return connected;		
	}
	
	/**
	 * This function returns true if
	 * 		1. gene and other do not overlap
	 * 	    2. gene and other are in the same orientation
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean isCandidate(Gene gene,Gene other){ 
		
		if(gene.overlaps(other)){
			//logger.info(gene.getName()+" overlaps "+gene.toUCSC());
			return false;
		}
		else{
			//if transcript is in an intron of the other transcript
			if(gene.getEnd()>other.getStart() && other.getEnd()>gene.getStart()){
				//logger.info("One is inside the intron of the other");
				return false;
			}
			if(gene.getOrientation().equals(other.getOrientation())){
				//logger.debug("The orientation is same");
				return true;
			}
		}
		//logger.debug("Something else");
		return false;
	}
	
	/**
	 * This function will return true if gene and other have at least 1 paired end read in common
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean pairedEndReadSpansTranscripts(Gene gene,Gene other){
		
		boolean rtrn = false;
		//Get all overlapping paired end reads in same orientation as the gene
		CloseableFilterIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(new BasicAnnotation(gene.getChr(),gene.getStart(),other.getEnd()),false), new SameOrientationFilter(gene));
		while(iter.hasNext()){
			Alignment read = iter.next();
			List<Annotation> mates = (List<Annotation>) read.getReadAlignments(model.getCoordinateSpace());
			if(mates.get(0).getOrientation().equals(gene.getOrientation()) && mates.get(1).getOrientation().equals(gene.getOrientation())){
				if((gene.overlaps(mates.get(0)) && other.overlaps(mates.get(1)))
						||
						gene.overlaps(mates.get(1)) && other.overlaps(mates.get(0))){
					rtrn = true;
					break;
				}
			}
			
		}
		iter.close();
		return rtrn;
	}
	
	/**
	 * Returns true if the transcript at the oriented 3' end is a single exon transcript.
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean secondTranscriptIsSingleExon(Gene gene,Gene other){
		
		Pair<Gene> orderedGenes = ConnectDisconnectedTranscripts.getOrderedAssembly(gene, other);
		if(gene.isNegativeStrand()){
			if(orderedGenes.getValue1().getBlocks().size()==1){
				return true;
			}
		}
		else{
			if(orderedGenes.getValue2().getBlocks().size()==1){
				return true;
			}
		}
		return false;
	}
	
	public static void write(String save, Map<String, Collection<Gene>> rtrn) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String chr: rtrn.keySet()){
			Collection<Gene> genes=rtrn.get(chr);
			for(Gene gene: genes){
				writer.write(gene.toBED()+"\n");
			}
		}
		
		writer.close();
	}
	
	public static void memoryStats(){
		int mb = 1024*1024;
        
        //Getting the runtime reference from system
        Runtime runtime = Runtime.getRuntime();
         
        logger.debug("##### Heap utilization statistics [MB] #####");
         
        //Print used memory
        logger.debug("Used Memory:"
            + (runtime.totalMemory() - runtime.freeMemory()) / mb);
 
        //Print free memory
        logger.debug("Free Memory:"
            + runtime.freeMemory() / mb);
         
        //Print total available memory
        logger.debug("Total Memory:" + runtime.totalMemory() / mb);
 
        //Print Maximum available memory
        logger.debug("Max Memory:" + runtime.maxMemory() / mb);
	}
	
	/**
	 * Makes assemblies using all reads and subjects all assemblies to a list of filters.
	 * @param chr
	 * @param strand
	 * @return
	 */
	private ChromosomeTranscriptGraph assembleDirectly(String chr,TranscriptionRead strand){
		
		model=new JCSAlignmentModel(bamfile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		model.addFilter(new ProperPairFilter());
		model.addFilter(new IndelFilter());
		model.addFilter(new GenomicSpanFilter(20000000));
		this.space=model.getCoordinateSpace();

		long S = System.currentTimeMillis();	
		logger.info("Assembling spliced reads");
		long start = System.currentTimeMillis();		
		//SPLICED READS
		CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new CanonicalSpliceFilter(genomeSeq));
		IntervalTree<Assembly> splicedAssemblies=assembleDirectly(splicedIter,strand);
		long end = System.currentTimeMillis();
		logger.debug("TIME: ASSEMBLE SPLICED: "+(end-start));
		//try{write(splicedAssemblies, outName+"."+chr+"."+"01splicedAssemblies.bed");}catch(IOException ex){}		
		logger.info("Size of spliced assemblies: "+splicedAssemblies.size());
		
		//NON_SPLICED READS
/*		CloseableIterator<Alignment> iter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new NoSpliceFilter());
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, splicedAssemblies,strand);
		end = System.currentTimeMillis();*/
		
		// ALL READS
		CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, splicedAssemblies,strand);
		end = System.currentTimeMillis();
		
		/*
		 * ONLY ALL READS
		 */ 
/*		long start = System.currentTimeMillis();
		logger.info("Assembling all reads");
		CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, strand);
		long end = System.currentTimeMillis();*/
		logger.debug("TIME: ASSEMBLE NON SPLICED: "+(end-start));
		logger.info("Size of direct assemblies: "+workingAssemblies.size());		
		//try{write(workingAssemblies, outName+"."+chr+"."+"02directAssemblies.bed");}catch(IOException ex){}
				
/*		model=new JCSAlignmentModel(bamFileName.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		model.addFilter(new ProperPairFilter());
		model.addFilter(new IndelFilter());
		model.addFilter(new GenomicSpanFilter(20000000));*/
		globalPairedLambda = model.getGlobalPairedFragments()/model.getGlobalLength();
		logger.info("Intron retention filter");
		//REMOVE SPURIOUS
		start = System.currentTimeMillis();
		IntervalTree<Assembly> unspuriousAssemblies = intronRetentionFilter(workingAssemblies,chr);
		end = System.currentTimeMillis();
		logger.debug("TIME: REMOVE SPURIOUS: "+(end-start));
		//try{write(unspuriousAssemblies, outName+"."+chr+"."+"03intronRetentionAssemblies.bed");}catch(IOException ex){}
									
		logger.info("Merge assemblies");
		start = System.currentTimeMillis();
		mergeAssembly(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: MERGE: "+(end-start));
		//try{write(unspuriousAssemblies, outName+"."+chr+"."+"04mergedAssemblies.bed");}catch(IOException ex){}

		logger.info("Extend compatible assemblies");
		start = System.currentTimeMillis();
		extendAssembly(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: BRANCHING: "+(end-start));
		//try{write(unspuriousAssemblies, outName+"."+chr+"."+"05extendedAssemblies.bed");}catch(IOException ex){}
				
		logger.info("Remove premature assemblies");
		start = System.currentTimeMillis();
		// Flag and remove premature
		IntervalTree<Assembly> filteredAssemblies = removePrematureAssemblies(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: FILTER PREMATURE: "+(end-start));
		//IntervalTree<Assembly> filteredAssemblies = removePrematureAssembliesUsingPairedEnds(unspuriousAssemblies);
		//try{write(filteredAssemblies, outName+"."+chr+"."+"06filteredAssemblies.bed");}catch(IOException ex){}

		start = System.currentTimeMillis();
		IntervalTree<Assembly> highSplicedAssemblies = removeLowSpliceJunctionAssemblies(filteredAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: REMOVE LOW SPLICED JUNCTIONS: "+(end-start));
		//try{write(highSplicedAssemblies, outName+"."+chr+"."+"07highSplicedAssemblies.bed");}catch(IOException ex){}
		
		//Lets split potential preprocessed and mature transcripts and test whether to include certain nodes
		
		start = System.currentTimeMillis();
		//Make graph
		ChromosomeTranscriptGraph graph=makeGraph(highSplicedAssemblies, chr); //TODO Should use the merged set
		end = System.currentTimeMillis();
		logger.debug("TIME: MAKE GRAPH: "+(end-start));
		long E = System.currentTimeMillis();	
		logger.debug("TIME: TOTAL "+(E-S));
		return graph;
	}

	private void extendAssembly(IntervalTree<Assembly> tree) {
		//We have a set of assemblies that are all incompatible
		//We want to link up parts
		logger.debug("Enter extend assembly");
		//iterate through and get overlapping assemblies
		Iterator<Assembly> iter=tree.toCollection().iterator();
		Collection<Assembly> considered = new HashSet<Assembly>();
		
		//IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();
		
		while(iter.hasNext()){
			//try and merge the non-overlapping portions
			Assembly assembly1=iter.next();
			//currentAssemblies.put(assembly1.getStart(), assembly1.getEnd(), assembly1);
			//boolean a1Changed = false;
			considered.add(assembly1);
			Iterator<Assembly> overlappers=tree.overlappingValueIterator(assembly1.getStart(), assembly1.getEnd());
			while(overlappers.hasNext()){
				Assembly assembly2=overlappers.next();
				//logger.debug("Assembly1: "+assembly1.getName()+" Assembly2: "+assembly2.getName());
				//logger.debug("overlaps "+assembly2.toUCSC());
				if(!considered.contains(assembly2)){
					//try and merge the non-overlapping portions
					Collection<Assembly> merged=branchAssemblies(assembly1, assembly2);
					if(!merged.isEmpty()){
						logger.debug("Branching "+assembly1.getName()+" and "+assembly2.getName());
						
						//remove assembly1
						tree.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//remove assembly2
						tree.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);						
						considered.remove(assembly1);
						considered.remove(assembly2);
						//currentAssemblies.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//currentAssemblies.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);
						//add merged
						//logger.debug("No. of assemblies "+ merged.size());
						for(Assembly merge: merged){
							tree.put(merge.getStart(), merge.getEnd(), merge);
							considered.add(merge);
							//currentAssemblies.put(merge.getStart(), merge.getEnd(), merge);
						}
					}
				}
			}
		}
		//return currentAssemblies;
	}

	private Collection<Assembly> branchAssemblies(Assembly assembly1, Assembly assembly2) {
		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		
		//order the two assemblies by which starts first
		Pair<Assembly> orderedAssembly=getOrderedAssembly(assembly1, assembly2);
		
		//go from first.getFirstExon() to second.getFirstExon
		Annotation portionToConsiderAdding=orderedAssembly.getValue1().intersect(new Alignments(assembly1.getChr(), orderedAssembly.getValue1().getStart(), orderedAssembly.getValue2().getBlocks().iterator().next().getEnd()));
//		Annotation portionToTest = portionToConsiderAdding.minus(orderedAssembly.getValue2());
		//ONLY IF THE REGION OF THE PORTIONTOCONSIDERADDING THAT DOES NOT OVERLAP THE SECOND ASSEMBLY IS NOT SPLICED
//		if(portionToTest.getSpliceConnections().isEmpty()){
			//if this is compatible with second assemebly
			if(compatible(orderedAssembly.getValue2(), portionToConsiderAdding)){
				Assembly merged=merge(portionToConsiderAdding, orderedAssembly.getValue2());
				merged.setName(orderedAssembly.getValue2().getName());
				rtrn.add(merged);
				rtrn.add(orderedAssembly.getValue1());
				//System.out.println(merged);
			}
//		}
		//get first exon of later start site exon
		//trim assembly 2 from start till this point
		//ask if compatible
		//if so, merge trim with rest of later assembly
		return rtrn;
	}

	private Pair<Assembly> getOrderedAssembly(Assembly assembly1, Assembly assembly2) {
		Pair<Assembly> rtrn=new Pair<Assembly>();
		//Order by CompareTo
		if(assembly1.compareTo(assembly2)<0){
			rtrn.setValue1(assembly1);
			rtrn.setValue2(assembly2);
		}
		else{
			rtrn.setValue1(assembly2);
			rtrn.setValue2(assembly1);
		}
		return rtrn;
	}

	private void write(IntervalTree<Assembly> tree, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Assembly> set=tree.toCollection();
		
		for(Assembly align: set){
			writer.write(align+"\n");
		}
		
		writer.close();
	}
		
	private IntervalTree<Assembly> assembleDirectly(CloseableIterator<Alignment> iter,TranscriptionRead strand) {
		IntervalTree<Assembly> workingAssemblies=new IntervalTree<Assembly>();
		return assembleDirectly(iter, workingAssemblies,strand);
	}


	/**
	 * Makes a graph using the specified assemblies for the specified chrosomome
	 * This function will go through each gene and only add those that pass the p-value threshold
	 * @param workingAssemblies
	 * @param chr
	 * @return
	 */
	private ChromosomeTranscriptGraph makeGraph(IntervalTree<Assembly> workingAssemblies, String chr) {
		ChromosomeTranscriptGraph graph=new ChromosomeTranscriptGraph(chr);
		Iterator<Node<Assembly>> iter=workingAssemblies.iterator();
		//For each assembly node
		while(iter.hasNext()){
			Collection<Assembly> genes=iter.next().getContainedValues();
			//For each gene(assembly)
			for(Annotation gene: genes){
				//Get all the exons of the gene
				List<? extends Annotation> blocks=gene.getBlocks();
				//if the gene has 1 exon
				if(blocks.size()==1){
					double pval = getScores(gene)[1];
					//System.err.println(gene.toUCSC()+" Count: "+score.getCount()+" pval:"+pval+"new count: "+s+"  new p-val "+pval2);
					if(pval<alpha){
						//System.err.println("passes");
						graph.connectVertexToGraph(blocks.get(0));
					}
					else{
						logger.debug(gene.getName()+" is not added to the graph because it does not pass the p-value threshold.");
					}
				}
				else{
			/*		double pval = getScores(gene)[1];
					if(pval<DEFAULT_ALPHA){*/
						graph.addAnnotationToGraph(gene);
			/*		}
					else{
						logger.debug("Gene "+gene.toUCSC()+" is filtered out because it does not meet the significance threshold : "+pval);
					}*/
				}
			}
		}
		return graph;
	}	

	/**
	 * Returns paired end counts and scan p-value for the specified gene
	 * @param gene
	 * @return
	 */
	private double[] getScores(Annotation gene){
		double[] scores = new double[2];
		scores[0] = 0.0;
		//Get all reads overlapping the transcript
		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(gene,true), new PairedAndProperFilter());
		//For each read,
		while(iter.hasNext()){					
			Alignment read = iter.next();
			boolean countRead = true;
			for(Annotation mate:read.getReadAlignments(space)){
				if(!compatible(gene,mate)){
					//logger.debug("Read "+mate.toUCSC()+" is not compatible with isoform with "+isoform.getExons().length);
					countRead=false;
					break;
				}
			}
			//For the assembly, we need to treat each read separately	
			if(countRead){
				scores[0] += read.getWeight();
			}
		}
		iter.close();
		//logger.debug("Count = "+scores[0]+" Int version "+new Double(scores[0]).intValue()+" global paired lambda = "+globalPairedLambda+" gene size = "+model.getCoordinateSpace().getSize(gene)+ " or "+gene.size()+" global length = "+model.getGlobalLength()+" global lambda = "+model.getGlobalLambda());
		scores[1] = AlignmentDataModelStats.calculatePVal(new Double(scores[0]).intValue(), globalPairedLambda,gene.size(), model.getGlobalLength());
		
//		scores[1] = AlignmentDataModelStats.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		
		return scores;
	}
	
	private void mergeAssembly(IntervalTree<Assembly> tree) {
		//Step 1: Iterate through the working assembly and find paths that overlap
		Iterator<Assembly> iter=tree.toCollection().iterator();
		//if overlaps but is incompatible, see if we can branch
		//iterate
		Collection<Assembly> considered = new HashSet<Assembly>();
 		logger.debug("Enter merging");
		while(iter.hasNext()){
			Assembly branch1=iter.next();
			considered.add(branch1);
			//get overlapping branches
			Iterator<Assembly> overlappers=tree.overlappingValueIterator(branch1.getStart(), branch1.getEnd());
			//Collection<Assembly> toRemove=new TreeSet<Assembly>();
			while(overlappers.hasNext()){
				Assembly branch2=overlappers.next();
				if(!considered.contains(branch2)){
					if(!branch1.equals(branch2) && compatible(branch1, branch2)){
						logger.debug("Merging: "+ branch1.getName()+ " and "+branch2.getName());
						Assembly merged=merge(branch1, branch2);
						merged.setName(branch1.getName());
						//remove annotation1 and annotation2
						tree.remove(branch1.getStart(), branch1.getEnd(), branch1);
						tree.remove(branch2.getStart(), branch2.getEnd(), branch2);
						//add merged
						tree.put(merged.getStart(), merged.getEnd(), merged);
					}
				}
			}
		}
	}
	
	/**
	 * This function removes any overlapping bad isoforms
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> intronRetentionFilter(IntervalTree<Assembly> assemblies,String chr) {
		//BAMFileWriter writer = new BAMFileWriter(new File(outName+"."+chr+".filtered.bam"));
		//writer.setSortOrder(SortOrder.coordinate, false);
		//writer.setHeader(header);
		//writer.setSortOrder(SortOrder.coordinate, false);
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
/*		double normalizationFactor = (double)spliceCount/(double)(spliceCount+nonSpliceCount);
		logger.error("R = "+normalizationFactor);*/
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();
		//For each assembly
		while(iter.hasNext()){
			Assembly assembly1=iter.next();
			//If assembly is not already marked as spurious
			if(!assembly1.isSpurious()){

				boolean	toRemoveAssembly1 = false;
				//For all overlapping assemblies
				Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly1.getStart(), assembly1.getEnd());
				while(overlappers.hasNext() && !toRemoveAssembly1){
					Assembly assembly2 = overlappers.next();
					if(assembly1.getOrientation().equals(assembly2.getOrientation())){
						//For all exons of this assembly
						for(Annotation intron1: assembly1.getSpliceConnections()){
							if(!toRemoveAssembly1){
								for(Annotation exon2:assembly2.getBlocks()){
									//If intron in assembly1 is contained in exon of assembly2
									if(exon2.contains(intron1)){
										//CHECK FOR COVERAGE
										double intronCount = model.getIntronCounts(intron1);
										double exonCount = (model.getCountStranded(new BasicAnnotation(exon2.getChr(), intron1.getStart(), intron1.getEnd(), exon2.getOrientation()), false)
																	/(double)(intron1.getEnd()-intron1.getStart()));
										double ratioI = intronCount/(intronCount+exonCount);
										double ratioE = exonCount/(intronCount+exonCount);
										//CASE 1 : Assembly with intron is mature and assembly with exon is immature
										if(ratioI>THRESHOLD_SPURIOUS){
											logger.debug("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.getName()+ " to intron "+intron1.toUCSC()+" in assembly "+assembly1.getName());
											assembly2.setSpurious(true);
											logger.debug("CASE 1: assembly "+ assembly2.getName()+" is spurious. Intron count "+intronCount+" > "+exonCount+" Ratio: "+ratioI);
										}
										else{
											//CASE 2: Assembly with exon is real and the other 
											if(ratioE>THRESHOLD_SPURIOUS){
												//logger.warn(assembly.toUCSC()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC());
												logger.debug("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.getName()+ " to intron "+intron1.toUCSC()+" in assembly "+assembly1.getName());
												assembly1.setSpurious(true);
												toRemoveAssembly1 = true;
												logger.debug("CASE 2: assembly "+ assembly1.getName()+" is spurious.Exon count "+exonCount+" > "+intronCount+" Ratio: "+(ratioE));
												break;
											}
											//else{
												//Both are real
											//}
										}
									}
								}
							}
							else{
								break;
							}
						}
					}
				}
				if(!toRemoveAssembly1){
					currentAssemblies.put(assembly1.getStart(), assembly1.getEnd(), assembly1);
					//for this assembly, get all the 
					//writer.addAlignment(alignment)
				}
			}
		}
		
		iter=currentAssemblies.toCollection().iterator();
		Set<Alignment> reads = new HashSet<Alignment>();
		while(iter.hasNext()){
			Assembly assembly = iter.next();
			CloseableIterator<Alignment> riter = model.getOverlappingReads(assembly, true);
			while(riter.hasNext()){
				Alignment read = riter.next();
				//System.out.println(read.getName());
				if(!reads.contains(read)){
					//writer.addAlignment(read.toSAMRecord());
					reads.add(read);
				}
			}
			riter.close();
		}
		//writer.close();
		return currentAssemblies;
	}
	
	/**
	 * This function removes any premature assemblies that have coverage conflicts in the exons
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removePrematureAssemblies(IntervalTree<Assembly> assemblies) {
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();

		while(iter.hasNext()){
			Assembly assembly=iter.next();
			boolean	toRemove = false;
			//For all overlapping assemblies
			Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly.getStart(), assembly.getEnd());
			while(overlappers.hasNext()){
				Assembly overlapper = overlappers.next();
				if(assembly.getOrientation().equals(overlapper.getOrientation())){
					
					//Check if the intron"ed" assembly is confident so it can be used
					if(!overlapper.isConfidentIsSet()){
						setConfidence(overlapper);
					}
					if(overlapper.isConfident()){
						//logger.warn("Compare " +assembly.getName()+" and "+overlapper.getName());
						//CHECK IF THE PARTS OF THE ASSEMBLY OTHER THAN THIS EXON OF THE ASSEMBLY) ARE COMPATIBLE
						//For all exons of this assembly
						for(Annotation exon1: assembly.getBlocks()){
							if(!toRemove){
								//If exon overlaps introns of overlapping assembly
								for(Annotation intron2:overlapper.getSpliceConnections()){
									//Flag as premature				
									if(exon1.overlaps(intron2)){
										Annotation p = assembly.minus(exon1);
										Annotation q = overlapper.minus(exon1);
										if(compatible(p, q) || p.length()==0){
											assembly.setPossiblePremature(true);
											//CHECK FOR COVERAGE
											if(coveragePassesCheck(exon1,intron2)){
												//toRemove = false;
											}
											else{
												//logger.warn(assembly.getName()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC()+" of "+overlapper.getName());
												toRemove = true;
												break;
											}
										}
									}
								}
							}
							else{
								break;
							}
						}
					}
				}
			}
			if(!toRemove){
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		
		return currentAssemblies;
	}
	
	private void setConfidence(Assembly overlapper){
		//for each intron
		for(Annotation intron2:overlapper.getSpliceConnections()){
			Annotation[] exons = overlapper.getFlankingBlocks(intron2);
			//exon[0] is left of intron
			//coverage to left of left exon
			//Since it is one base, we dont need coverage
			double leftScore = model.getCountStranded(new BasicAnnotation(exons[0].getChr(),exons[0].getEnd()-2,exons[0].getEnd()-1,overlapper.getOrientation()),false);
			//coverage to right of right exon
			double rightScore = model.getCountStranded(new BasicAnnotation(exons[1].getChr(),exons[1].getStart()+1,exons[1].getStart()+2,overlapper.getOrientation()),false);
			
			if(leftScore<rightScore){
				if(leftScore<rightScore*coveragePercentThreshold){ 
					logger.debug(overlapper.getName()+" does not pass confidence test because "+exons[0].toUCSC()+" has score "+leftScore+" compared to "+
								exons[1].toUCSC()+" which has "+rightScore);
					overlapper.setConfident(false);
				}
			}else{
				if(rightScore<leftScore*coveragePercentThreshold){
					logger.debug(overlapper.getName()+" does not pass confidence test because "+exons[1].toUCSC()+" has score "+rightScore+" compared to "+
							exons[0].toUCSC()+" which has "+leftScore);
					overlapper.setConfident(false);
				}
			}
		}
		if(!overlapper.isConfidentIsSet())
			overlapper.setConfident(true);
	}
	/**
	 * This function removes any premature assemblies that have coverage conflicts in the exons
	 * @param assemblies
	 * @return
	 */
/*	private IntervalTree<Assembly> removePrematureAssembliesUsingPairedEnds(IntervalTree<Assembly> assemblies) {
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();

		while(iter.hasNext()){
			Assembly assembly=iter.next();
			boolean	toRemove = false;
			//For all overlapping assemblies
			Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly.getStart(), assembly.getEnd());
			while(overlappers.hasNext()){
				Assembly overlapper = overlappers.next();
				//For all exons of this assembly
				for(Annotation exon1: assembly.getBlocks()){
					if(!toRemove){
						//If exon overlaps introns of overlapping assembly
						for(Annotation intron2:overlapper.getSpliceConnections()){
							//Flag as premature				
							if(exon1.overlaps(intron2)){
								assembly.setPossiblePremature(true);
								//CHECK FOR COVERAGE
								if(exonPassesPairedEndTest(exon1,overlapper,intron2)){
									//toRemove = false;
								}
								else{
									logger.warn(assembly.toUCSC()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC());
									toRemove = true;
									break;
								}
							}
						}
					}
					else{
						break;
					}
				}
			}
			if(!toRemove){
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		
		return currentAssemblies;
	}*/
	
	/**
	 * This function returns true if the number of paired end reads supporting exon1 and supporting the splice junction og intron2 are comparable, 
	 * that is, not below a certain threshold. HELPER function to removePrematureAssembliesUsingPairedEnds()
	 * @return
	 */
	private boolean exonPassesPairedEndTest(Annotation exon1,Assembly overlapper,Annotation intron2){
		
		Annotation overlap = exon1.intersect(intron2);
		Annotation nonoverlap = exon1.minus(intron2);
		
		//If exon is completely contained in the intron, pass
		if(overlap==null || overlap.equals(exon1)){
			return true;
		}
		boolean passes=true;
		//Get paired end reads for exon1
		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(exon1,true), new ProperPairFilter());
		double exonCount = 0.0;
		while(iter.hasNext()){
			Alignment read = iter.next();
			// if read is compatible with gene
			//If either read is not compatible flag will be set to false
			//If both mates of the pair are compatible
			if(boundaryIsCompatibleWithPairedEnd(overlap,nonoverlap,read)){
				exonCount = exonCount + 1.0;
			}
		}
		//exonCount = exonCount/(double)exon1.length();
		iter.close();
		
		//Get paired end reads for exons flanking intron2
		
		Annotation[] flankingExons = overlapper.getFlankingBlocks(intron2);
		Assembly sumAssembly = new Assembly(Arrays.asList(flankingExons),false);
		//int intronLength = flankingExons[0].length()+flankingExons[1].length();
		iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(sumAssembly,true), new ProperPairFilter());
		double intronCount = 0.0;
		while(iter.hasNext()){
			Alignment read = iter.next();
			// if read is compatible with gene
			//If either read is not compatible flag will be set to false
			//If both mates of the pair are compatible
			if(boundaryIsCompatibleWithPairedEnd(flankingExons[0],flankingExons[1],read)){
				intronCount = intronCount + 1.0;
			}
		}
		iter.close();
		//intronCount = intronCount / (double)intronLength;
		if(intronCount!=0.0){
			double ratio = exonCount/intronCount;
			if(ratio<0.8){
				logger.debug("Paired end count for exon " +exon1.toUCSC()+" : "+exonCount);
				logger.debug("Paired end count for intron " +intron2.toUCSC()+" : "+intronCount+ " ratio = "+ratio);
				return false;
			}
		}
		return passes;
	}
	
	/**
	 * Returns an array of size 2 where [0] is the first mate and [1] is the second mate
	 * @param read
	 * @return
	 */
	private Annotation[] getPairedEndReads(Alignment read){
		
		//We implement this specific to our purpose with paired ends
		Annotation[] annArr = new Annotation[2];
		int i=0;
		for(Annotation a: read.getReadAlignments(space)){
			if(i>1){
				logger.error("More than 2 mates are returned");
			}
			annArr[i] = a;
			i++;
		}
		return annArr;
	}
	
	/**
	 * This function returns true is the first mate is in the first exon and the second mate is in the second exon
	 * @param exon1
	 * @param exon2
	 * @param mates
	 * @return
	 */
	private boolean boundaryIsCompatibleWithPairedEnd(Annotation exon1,Annotation exon2,Alignment read){
		
		Annotation[] mates = getPairedEndReads(read);
		if((compatible(exon1,mates[0]) && 
				compatible(exon2,mates[1]))||
					(compatible(exon2,mates[0]) 
							&& compatible(exon1,mates[1])))
			return true;
		else
			return false;
	}
	/**
	 * This function will remove any assemblies that have introns with very low number of splice junctions
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removeLowSpliceJunctionAssemblies(IntervalTree<Assembly> assemblies){
	
		logger.info("Remove low splice junction assemblies");
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
				
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();//assemblies;

		//For each assembly
		while(iter.hasNext()){
			Assembly assembly=iter.next();
			
			boolean toRemove=false;
			//If the assembly does not pass the confidence test, remove it
			if(!assembly.isConfidentIsSet()){
				logger.debug("Confidence is set in the splice filter");
				setConfidence(assembly);
			}
			if(!assembly.isConfident()){
				logger.debug(assembly.getName()+" removed because does not pass the confidence test.");
			} else{		
				//IF THE GENE HAS ONE INTRON, check min #splice reads
				if(assembly.getSpliceConnections().size()==1){
					for(Annotation intron:assembly.getSpliceConnections()){
						double count=model.getIntronCounts(intron);
						if(count<minSpliceReads){
							toRemove = true;
						}
					}
					if(!toRemove)
						currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
				}
				else{
					Map<Annotation,Double> intronToSplicedCountMap = new TreeMap<Annotation,Double>();
					double avgCount = 0.0;
					//boolean toRemove = false;
					//for each intron
					for(Annotation intron:assembly.getSpliceConnections()){
						//find the number of supporting splice junctions
						double count=model.getIntronCounts(intron);											
						avgCount +=count;
						intronToSplicedCountMap.put(intron, count);
					}
					
					//Calculate average
					avgCount = (avgCount / (double)intronToSplicedCountMap.keySet().size());
					toRemove = false;
					//Check for all introns
					int cnt=0;
					//First check the first and last intron. If they do not pass the threshold
					//then trim and then go over this again.
					for(Annotation intron:intronToSplicedCountMap.keySet()){
						if(cnt==0 || cnt==intronToSplicedCountMap.keySet().size()-1){
							if(intronToSplicedCountMap.get(intron)<=avgCount*minSplicePercent){
								assembly = new Assembly(trimEnds(assembly,0.25));
							}
						}
						cnt++;
					}
					
					for(Annotation intron:assembly.getSpliceConnections()){
						if(intronToSplicedCountMap.get(intron)<=avgCount*minSplicePercent){
							//If  assembly is flagged to be removed because of an intron,
							logger.debug(assembly.getName()+" removed because intron "+intron.toUCSC()+" has coverage "+intronToSplicedCountMap.get(intron)+" compared to "+avgCount);
							toRemove = true; 
							//return true;
						}
					}
					if(!toRemove){
						//Removed by this filter but confidence was set to remove
						/*if(!assembly.isConfident())
							logger.debug(assembly.getName()+" is NOT confident but RETAINED.");*/
						currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
					}
					else{
						//Assembly flagged to remove
						logger.debug(assembly.getName()+" is confident but flagged to be REMOVED. Check trimming.");
						//trimEnds(assembly);
						
					}
				}
			}
		}
		return currentAssemblies;
	}

	/**
	 * Returns the assembly trimmed at both ends
	 * @param assembly
	 * @return
	 */
	private Annotation trimEnds(Annotation assembly,double pct){
		
		List<Double> counts = model.getCountsStrandedPerPosition(assembly);
		double[] cntArr = l2a(counts);

		Collections.sort(counts);		
		double cutoff = Math.max(2, Statistics.quantile(counts, pct));
/*		for(int i=0;i<assembly.size();i++){
			System.out.print(cntArr[i]+" ");
		}*/
		int trimStart = MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(cntArr, cutoff);
		int trimEnd   =  assembly.size() - MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(cntArr, cutoff);
		logger.debug(assembly.getName()+" trimStart " + trimStart + " trimEnd " + trimEnd + ", transcript length " + assembly.size()+ " with cutoff "+cutoff);

/*			int transcriptFirstExonEnd = assembly.getOrientation().equals(Strand.NEGATIVE) 
					? assembly.getPositionAtReferenceCoordinate(assembly.getStart())
					: assembly.getPositionAtReferenceCoordinate(assembly.getEnd()-1);

			int transcriptLastExonStart = assembly.getOrientation().equals(Strand.NEGATIVE)
						? assembly.getPositionAtReferenceCoordinate(assembly.getEnd()-1)
						: assembly.getPositionAtReferenceCoordinate(assembly.getStart());
				
			logger.trace("first exon end in transcript " + transcriptFirstExonEnd + " last exon start " + transcriptLastExonStart);

			if(trimStart > transcriptFirstExonEnd) {
				trimStart = Math.max(0, transcriptFirstExonEnd - 50);
			}	
			if(trimEnd < transcriptLastExonStart) {
				trimEnd = Math.min(align.getTranscriptLength(), transcriptLastExonStart + 50);
			}	
*/			
		if(trimStart>trimEnd){
			return assembly;
		}
		if(assembly.getOrientation().equals(Strand.NEGATIVE)){
			int temp = trimStart;
			trimStart = trimEnd;
			trimEnd = temp;
			logger.debug("Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
		}	
			
		Assembly newAssembly = new Assembly(assembly);
			
		if(((trimEnd+trimStart) < assembly.size()) && trimStart<trimEnd && trimStart<assembly.size() && trimEnd<assembly.size()){
			newAssembly.trim(trimStart, trimEnd);
			logger.debug("trimming ("+trimStart +" - "+ trimEnd+") gene was: " + assembly.toBED() + " and now is: " +newAssembly.toBED());
		}
		return newAssembly;
	}
	
	/**
	 * Helper function to removePrematureAssemblies
	 * @param exon1
	 * @param intron2
	 * @return
	 */
	private boolean coveragePassesCheck(Annotation exon1, Annotation intron2){
		Annotation overlap = exon1.intersect(intron2);
		//Annotation nonoverlap = exon1.minus(intron2);
		
		if(overlap==null){
			return true;
		}
		//If fully contained
		//TODO: A minimum coverage threshold?
		if(overlap.equals(exon1)){
			//logger.error("Fully Contained");
			return true;
		}
		
		if(overlap.getLengthOnReference()<=3){
			//logger.warn(exon1.toUCSC()+" overlaps "+intron2.toUCSC()+" with <=3");
			return false;
		}
		
		Annotation exonBoundary = null;
		Annotation intronBoundary = null;
		if(intron2.getEnd()>exon1.getEnd()){
			exonBoundary = new BasicAnnotation(exon1.getChr(),intron2.getStart()-2,intron2.getStart()-1,exon1.getOrientation());
			intronBoundary = new BasicAnnotation(intron2.getChr(),intron2.getStart()+1,intron2.getStart()+2,exon1.getOrientation());
		}
		else{
			exonBoundary = new BasicAnnotation(exon1.getChr(),intron2.getEnd()+1,intron2.getEnd()+2,exon1.getOrientation());
			intronBoundary = new BasicAnnotation(intron2.getChr(),intron2.getEnd()-2,intron2.getEnd()-1,exon1.getOrientation());
		}
		double overlapScore = model.getCountStranded(intronBoundary,false);
		double nonoverlapScore = model.getCountStranded(exonBoundary,false);
//		logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
		if(overlapScore>=nonoverlapScore*coveragePercentThreshold){
			//logger.error(exon1.toUCSC()+" passes coverage test with intron "+ intron2.toUCSC());
			//logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
			//logger.debug("Intron boundary score for "+intronBoundary.toUCSC()+ " : "+overlapScore+" Exon boundary score for "+exonBoundary.toUCSC()+ " : "+nonoverlapScore);
			return true;
		}
		else{
			//logger.error(exon1.toUCSC()+" does not pass coverage test with intron "+ intron2.toUCSC());
			logger.debug("REMOVED: Intron boundary score for "+intronBoundary.toUCSC()+ " : "+overlapScore+" Exon boundary score for "+exonBoundary.toUCSC()+ " : "+nonoverlapScore);
			return false;
		}		
	}

	/**
	 * 
	 * @param iter
	 * @param workingAssemblies
	 * @param strand
	 * @return
	 */
	private IntervalTree<Assembly> assembleDirectly(CloseableIterator<Alignment> iter, IntervalTree<Assembly> workingAssemblies,TranscriptionRead strand) {
		
		String linc="gene_v2_";
		long cnt=0;
		long start = System.currentTimeMillis();
		int c=0;
		boolean flagPremature=!workingAssemblies.isEmpty();
		while(iter.hasNext()){
			
			Alignment reads=iter.next();			
			//reads.setFragmentStrand(strand);				
			//For the assembly, we need to treat each read separately
			for(Annotation read: reads.getReadAlignments(space)){
				cnt++;
				if(cnt%100000==0){
					long end = System.currentTimeMillis();
					logger.debug("Processed "+cnt+" reads in "+(end-start)/1000+" seconds with "+workingAssemblies.size()+" "+c);
					start = System.currentTimeMillis();
				}
				
				//Find all compatible assemblies
				Collection<Assembly> compatibleAssemblies = new ArrayList<Assembly>();

/*				logger.debug("New read : "+read.toUCSC()+" "+!read.getSpliceConnections().isEmpty());
				for(Assembly as:workingAssemblies.toCollection()){
					logger.debug(as.toUCSC()+"\t"+as.toBED());
				
				}*/
				//EACH READ HAS THE FRAGMENT STRAND
				//for each read, get overlapping assemblies
				Iterator<Node<Assembly>> overlappers=workingAssemblies.overlappers(read.getStart(), read.getEnd());
				//if no overlappers add read as assembly
				if(!overlappers.hasNext()){
					//add the read as an annotation
					//Flag this as likely premature
					Assembly readAssembly=new Assembly(read, false);
					readAssembly.setName(linc+globalCounter);
					globalCounter++;
					if(flagPremature){
						readAssembly.setPossiblePremature(true);
					}
					workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
					c++;
				}
				else{
					boolean hasCompatible=false;
//					Collection<Assembly> overlapNotCompatible = new TreeSet<Assembly>();
//					long s1 = System.currentTimeMillis();
					while(overlappers.hasNext()){
						Collection<Assembly> assemblies=new TreeSet<Assembly>(overlappers.next().getContainedValues());
						//Find all compatible assemblies
						for(Assembly assembly: assemblies){
							if(compatible(assembly, read)){
								compatibleAssemblies.add(assembly);
								Assembly merged=merge(assembly, read);
								merged.setName(assembly.getName());
								//remove assembly
								workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
								//add merged
								workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
								hasCompatible=true;
							}
/*							else{
								if(assembly.overlaps(read)){
									overlapNotCompatible.add(assembly);	
								}
							}*/
						}
						
					}
/*					if(cnt>2900){
						logger.debug("compatible assemblies size "+compatibleAssemblies.size());
						logger.debug("overlap not compatible "+overlapNotCompatible.size());
					}*/
					
					//If some compatible assembly
/*					if(hasCompatible & compatibleAssemblies.size()>0){
						//MERGE NON-OVERLAPPING ASSEMBLIES
						if(compatibleAssemblies.size()>1){
							Collection<Assembly> alreadyMerged = new HashSet<Assembly>();
							Collection<String> names = new HashSet<String>();
							boolean mergedFlag = false;
							//Compare all compatible assemblies to each other
							for (Iterator<Assembly> iterator = compatibleAssemblies.iterator(); iterator.hasNext();) {
								Assembly assembly1 = iterator.next();
								for (Iterator<Assembly> iterator2 = compatibleAssemblies.iterator(); iterator2.hasNext();) {
									Assembly assembly2 = iterator2.next();
									if(!(assembly1.equals(assembly2)) && (!names.contains(assembly2.getName()+"_"+assembly1.getName()))){
										if(!assembly1.overlaps(assembly2) && assembly1.getOrientation().equals(assembly2.getOrientation())){
											alreadyMerged.add(assembly1);
											alreadyMerged.add(assembly2);
											mergedFlag=true;
											//Merge assembly1 and assembly2 and read
											Assembly merged=merge(assembly1, read);
											merged=merge(merged,assembly2);
											for(Assembly a:workingAssemblies.toCollection()){
												if(merged.equals(a)){
													logger.debug("Merged: "+merged.toBED());
													logger.debug("a = "+a.toBED());
												}
											}
											String name =assembly1.getName()+"_"+assembly2.getName();
											merged.setName(name);
											workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
											c++;
											names.add(name);
										}
									}
								}
								if(!mergedFlag){
									Assembly merged=merge(assembly1, read);
									merged.setName(assembly1.getName());
									//remove assembly
									workingAssemblies.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
									//add merged
									workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
								}
								iterator.remove();
							}
							
							for(Assembly assembly:alreadyMerged){
								c--;
								workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
							}
						}
						else{
							//Only 1 compatible assembly
							Iterator<Assembly> compIter=compatibleAssemblies.iterator();
							while(compIter.hasNext()){
								Assembly comp = compIter.next();
								//if compatible --> merge
								Assembly merged=merge(comp, read);
								merged.setName(comp.getName());
								//remove assembly
								workingAssemblies.remove(comp.getStart(), comp.getEnd(), comp);
								//add merged
								workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
							}
						}
					}*/
					if(!hasCompatible){
/*						boolean splittable =false;
						if(!overlapNotCompatible.isEmpty()){
//							s1 = System.currentTimeMillis();
							for(Assembly assembly:overlapNotCompatible){
								if(read.getSpliceConnections().isEmpty()){
									//OLD SPLIT
									Collection<Assembly> branches=oldSplitBranch(assembly, read);
									if(branches!=null && !branches.isEmpty()){
										splittable=true;
										for(Assembly branch: branches){
											branch.setName(linc+globalCounter);
											globalCounter++;
											//add merged assemblies to the tree
											workingAssemblies.put(branch.getStart(), branch.getEnd(), branch);
										}
									}
								}
								else{
								//NEW SPLIT
									Collection<Assembly> branches=newSplitBranch(assembly, read);
									if(branches!=null && !branches.isEmpty()){
										splittable=true;
										//remove assembly
										//workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
										for(Assembly branch: branches){
											//add merged assemblies to the tree
											branch.setName(linc+globalCounter);
											globalCounter++;
											workingAssemblies.put(branch.getStart(), branch.getEnd(), branch);
										}
									}
								}
								
								//SOPHISTICATED SPLITTING
//								Annotation portionToConsiderAdding=assembly.intersect(new Alignments(assembly.getChr(), assembly.getStart(), read.getBlocks().get(0).getEnd()));
								if(assembly.getOverlappingAssembly(read)==null){
									logger.debug("Assembly: "+assembly.toBED());
									logger.debug("Read: "+read.toUCSC()+" "+!read.getSpliceConnections().isEmpty()+" "+assembly.overlaps(read));
								}
								Annotation portionToConsiderAdding=assembly.intersect
										(new Alignments(assembly.getChr(), assembly.getStart(), 
												(assembly.getOverlappingAssembly(read)).getEnd()));
								if(compatible(portionToConsiderAdding,read)){
									Assembly branch = merge(portionToConsiderAdding,read);
									branch.setName(linc+globalCounter);
									globalCounter++;
									workingAssemblies.put(branch.getStart(), branch.getEnd(), branch);
									c++;
									splittable=true;
								}
							}
							if(!splittable){
								Assembly readAssembly=new Assembly(read, false);
								readAssembly.setName(linc+globalCounter);
								globalCounter++;
								if(flagPremature){
									readAssembly.setPossiblePremature(true);
								}
								workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
								c++;
							}
						}
						else{*/
							Assembly readAssembly=new Assembly(read, false);
							readAssembly.setName(linc+globalCounter);
							globalCounter++;
							if(flagPremature){
								readAssembly.setPossiblePremature(true);
							}
							workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
							c++;
//						}
					}
				}
				
			}
		}		
		iter.close(); //close the iterator
		
		return workingAssemblies;
	}

	/**
	 * See if the read is partially compatible with the assembly
	 * If so, branch and merge
	 * @param assembly
	 * @param read
	 * @return all branches
	 */
	private Collection<Assembly> newSplitBranch(Assembly assembly, Annotation read) {

		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		if(read.overlaps(assembly)){
			Assembly truncated=findMaxCompatibilityNew(assembly, read);
			if(truncated!=null){
				Assembly merged=merge(truncated, read);
				rtrn.add(merged);
				//Why add the assembly again? Just don't remove it above
				//rtrn.add(assembly);
			}
		}
		return rtrn;
	}

	/**
	 * See if the read is partially compatible with the assembly
	 * If so, branch and merge
	 * @param assembly
	 * @param read
	 * @return all branches
	 */
	private Collection<Assembly> oldSplitBranch(Assembly assembly, Annotation read) {

		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		if(read.overlaps(assembly)){
			Assembly truncated=findMaxCompatibility(assembly, read);
			if(truncated!=null){
				Assembly merged=merge(truncated, read);
				rtrn.add(merged);
				//Why add the assembly again? Just don't remove it above
				rtrn.add(assembly);
			}
		}
		return rtrn;
	}
	
	private Assembly findMaxCompatibilityNew(Assembly assembly, Annotation read) {
		/*Iterator<Assembly> iter=assembly.trimNodes(read);
		while(iter.hasNext()){
			Assembly truncated=iter.next();
			if(compatible(truncated, read)){
				return truncated;
			}
		}
		*/
		Assembly trimmedAssembly = assembly.trim(read);
		if(trimmedAssembly!=null){
			//TODO: doesnt the logic dictate this?
			if(compatible(trimmedAssembly.getLastBlock(), read)){
				return trimmedAssembly;
			}
		}
		return null;
	}
	
	private Assembly findMaxCompatibility(Assembly assembly, Annotation read) {
		Iterator<Assembly> iter=assembly.trimNodes();
		while(iter.hasNext()){
			Assembly truncated=iter.next();
			if(compatible(truncated, read)){
				return truncated;
			}
		}
		return null;
	}

	/*private boolean partiallyCompatible(Annotation assembly, Annotation read) {
		//The annotations are not compatible, so we will check if they are partially compatible
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//If they dont overlap they cant be partially compatible
		if(!assembly.overlaps(read)){return false;}
		
		//Otherwise, they are partially compatible if:
		//(i) some of the introns are compatible
		if(!assemblyIntrons.isEmpty() && !readIntrons.isEmpty()){
			//if there are some introns that
		}
		
		//(ii) or if exons of one overlap the other but dont overlap the intron of the other
		
	}*/

	private Collection<Annotation> branch(Annotation read, Annotation assembly) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//Define incompatible regions
		
		//sanity check to ensure they physically overlap
		if(read.overlaps(assembly)){
			//define compatible edges: find introns that overlap but arent equal
			//define compatible nodes: find exons that overlap and are compatible
			Collection<Annotation> incompatible=this.getIncompatibleRegions(read, assembly);
			try{write(incompatible, "incompatible.bed");}catch(IOException ex){}
			Collection<Annotation> compatible=getCompatibleRegions(read, assembly, incompatible); //simple subtraction of intervals
			for(Annotation incompatibleRegion: incompatible){
				Annotation branch=merge(compatible, incompatibleRegion);
				rtrn.add(branch);
			}
		}
		
		return rtrn;
		
	}

	
	private Collection<Annotation> getCompatibleRegions(Annotation read, Annotation assembly, Collection<Annotation> incompatible) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		//simply take read and assembly and subtract the incompatible regions
		if(!incompatible.isEmpty()){
			Annotation ann1=(read.minus(incompatible));
			Annotation ann2=(assembly.minus(incompatible));
			if(ann1.getSize()>0 && compatible(ann1, read) && compatible(ann1, assembly)){
				rtrn.add(ann1);
			}
			if(ann2.getSize()>0 && compatible(ann2, read) && compatible(ann2, assembly)){
				rtrn.add(ann2);
			}
			
		}
		else{
			rtrn.add(read);
			rtrn.add(assembly);
		}
		return rtrn;
	}

	

	private Annotation merge(Collection<Annotation> compatible, Annotation incompatibleRegion) {
		//do a simple block merge
		Collection<Annotation> blocks=new TreeSet<Annotation>();
		for(Annotation ann: compatible){
			blocks.addAll(ann.getBlocks());
		}
		blocks.addAll(incompatibleRegion.getBlocks());
		Annotation rtrn=new BasicAnnotation(blocks);
		return rtrn;
	}

	
	
	//MG: This was working well
	public static boolean compatible(Annotation assembly, Annotation read) {
		//Two alignments will be defined as compatible if:
		//(i) the intronic locations are exactly same
		//if both have introns, ensure that there are no non-overlapping introns
		//logger.debug(read.getName());
		if(!assembly.getSpliceConnections().isEmpty() 
				|| !read.getSpliceConnections().isEmpty()){
			//check if introns are compatible
			if(areIntronsCompatible(assembly, read)){
				//logger.debug("Introns are compatible");
				return true;
			}
			//logger.debug("Introns are NOT compatible");
			//TODO In the case of partial compatibility we should split and make new path
			return false;
		}
		//(ii) one is spliced, the other is not 
		/*if(assembly.getBlocks().size()>1 || !read.getSpliceConnections().isEmpty()){
			//and the non-spliced does not overlap the intron at all, but overlaps the exon
			boolean overlapsWithoutCrossingIntron=overlapsWithoutCrossingIntron(assembly, read);
			if(overlapsWithoutCrossingIntron){return true;}
			return false;
		}*/
		//(iii) both are unspliced and overlap
		boolean overlap=overlap(assembly, read);
		
		if(overlap){return true;}
		return false;
	}
	

	private static boolean areIntronsCompatible(Annotation assembly, Annotation read) {
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//introns are compatible if:
		//(i) all overlapping introns are identical
		for(Annotation intron1: assemblyIntrons){
			for(Annotation intron2: readIntrons){
				//if overlaps but not identical
				if(intron1.overlaps(intron2,forceStrandSpecificity)){
					if(!intron1.equals(intron2, forceStrandSpecificity)){
					//	logger.debug("Case1");
						return false;
					} //TODO: This should use strand info or not based on flag
				}
			}
		}
		
		//(ii) if none of the exons overlap an intron in the other
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation intron2: readIntrons){
				if(exon1.overlaps(intron2,forceStrandSpecificity)){
					//logger.debug("Case2a");
					return false;
				}
			}
		}
		
		for(Annotation exon2: read.getBlocks()){
			for(Annotation intron1: assemblyIntrons){
				if(exon2.overlaps(intron1,forceStrandSpecificity)){
					//logger.debug("Case2b");
					return false;
				}
			}
		}
		
		//(ii) the introns dont overlap but the exons do
		//just need to test that any exons overlap
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation exon2: read.getBlocks()){
				if(exon1.overlaps(exon2,forceStrandSpecificity)){
					return true;
				}
			}
		}
		//logger.debug("Case3");
		return false;
	}
	
	
	private Collection<Annotation> getIncompatibleRegions(Annotation assembly, Annotation read){
		
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		Collection<Annotation> rtrn1=new TreeSet<Annotation>();
		Collection<Annotation> rtrn2=new TreeSet<Annotation>();
		
		//introns are compatible if:
		//(i) all overlapping introns are identical
		for(Annotation intron1: assemblyIntrons){
			for(Annotation intron2: readIntrons){
				//if overlaps but not identical
				if(intron1.overlaps(intron2)){
					if(!intron1.equals(intron2, forceStrandSpecificity)){
						rtrn1.add(intron1);
						rtrn2.add(intron2);
						//intron1 and intron2 are not compatible
					} //This should use strand info or not based on flag
				}
			}
		}
		
		
		//now check if an exon overlaps an intron in the other
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation intron2: readIntrons){
				if(exon1.overlaps(intron2)){
					//exon1 is incompatible
					rtrn1.add(exon1);
					rtrn2.add(intron2);
				}
			}
		}
		
		for(Annotation exon2: read.getBlocks()){
			for(Annotation intron1: assemblyIntrons){
				if(exon2.overlaps(intron1)){
					rtrn2.add(exon2);
					rtrn1.add(intron1);
				}
			}
		}
		
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//merge rtrn1
		if(!rtrn1.isEmpty()){
			Annotation block1=new BasicAnnotation(rtrn1);
			//trim assembly
			Annotation tmp1=assembly.intersect(new Alignments(block1.getChr(), block1.getStart(), block1.getEnd()));
			if(tmp1!=null && tmp1.getSize()>0){rtrn.add(tmp1);}
		}
		
		//merge rtrn2
		if(!rtrn2.isEmpty()){
			Annotation block2=new BasicAnnotation(rtrn2);
			//trim read
			Annotation tmp2=read.intersect(new Alignments(block2.getChr(), block2.getStart(), block2.getEnd()));
			if(tmp2!=null && tmp2.getSize()>0){rtrn.add(tmp2);}
		}
		
		//return collection
		return rtrn;
	}

	private static boolean overlap(Annotation assembly, Annotation read) {
		return assembly.overlaps(read,forceStrandSpecificity);
	}

	private boolean overlapsWithoutCrossingIntron(Annotation assembly, Annotation read) {
		//one has introns, one does not
		//see if one without overlaps the exons of the other
		//AND does not overlap the intron
		if(!read.getSpliceConnections().isEmpty()){
			//read has introns
			//test if exons overlap the intron
			for(Annotation intron: read.getSpliceConnections()){
				if(assembly.overlaps(intron)){return false;}
			}
			//if not, test that exons overlap
			for(Annotation exon: read.getBlocks()){
				if(assembly.overlaps(exon)){return true;}
			}
			return false;
		}
		
		//assembly has intron
		Annotation complement=assembly.complement();
		for(Annotation intron: complement.getBlocks()){
			if(read.overlaps(intron)){return false;}
		}
		for(Annotation exon: assembly.getBlocks()){
			if(read.overlaps(exon)){return true;}
		}
		
		return false;
	}

	private boolean hasNonOverlapping(Annotation assembly, Alignment read) {
		//check if every intron in read is also an intron in assembly and vice versa
		
		Collection<? extends Annotation> assemblyIntrons=assembly.complement().getBlocks();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		for(Annotation intron: assemblyIntrons){
			if(!readIntrons.contains(intron)){return true;}
		}
		
		for(Annotation intron: readIntrons){
			if(!assemblyIntrons.contains(intron)){return true;}
		}
		
		return false;
	}

	private Annotation merge(Annotation assembly, Alignment read) {
		//since these are compatable, we will simply add the read to the assembly by merging exons
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		Collection<Annotation> readAlignments=read.getReadAlignments(this.space);
		for(Annotation align: readAlignments){
			rtrn.addAll(align.getBlocks());
		}
		rtrn.addAll(assembly.getBlocks());
		return new BasicAnnotation(rtrn);
	}

	private Assembly merge(Annotation assembly, Annotation read) {
		//since these are compatable, we will simply add the read to the assembly by merging exons
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		rtrn.addAll(read.getBlocks());
		rtrn.addAll(assembly.getBlocks());
		return new Assembly(rtrn);

	}
	

	private void write(Collection<? extends Annotation> blocks, String save) throws IOException {
		write(blocks, save, true);
	}
		

	private void write(Collection<? extends Annotation> blocks, String save, boolean append) throws IOException {
		FileWriter writer=new FileWriter(save, append);
		
		for(Annotation block: blocks){writer.write(block+"\n");}
		
		writer.close();
	}

	private Pair<Collection<Annotation>> splitExonsByJunctions(Annotation exonicRegion, Collection<Annotation> junctions) {
		Collection<Annotation> spliced=new TreeSet<Annotation>();
		Collection<Annotation> preSpliced=new TreeSet<Annotation>();
		
		//make junction tree
		IntervalTree<Annotation> junctionTree=makeTree(junctions);
		
		//go through exons and see if it overlaps a junction
		for(Annotation exon: exonicRegion.getBlocks()){
			//if overlaps a junction
			Iterator<Node<Annotation>> overlappers=junctionTree.overlappers(exon.getStart()-1, exon.getEnd()+1);
			Pair<Collection<Annotation>> exons=split(exon, overlappers);
			spliced.addAll(exons.getValue1());
			preSpliced.addAll(exons.getValue2());
		}
		
		Pair<Collection<Annotation>> rtrn=new Pair<Collection<Annotation>>(spliced, preSpliced);
		return rtrn;
	}

	private IntervalTree<Annotation> makeTree(Collection<Annotation> junctions) {
		IntervalTree<Annotation> rtrn=new IntervalTree<Annotation>();
		for(Annotation junction: junctions){
			rtrn.put(junction.getStart(), junction.getEnd(), junction);
		}
		return rtrn;
	}

	private Pair<Collection<Annotation>> split(Annotation exon, Iterator<Node<Annotation>> overlappers) {
		Collection<Annotation> spliced=new TreeSet<Annotation>();
		Collection<Annotation> preSpliced=new TreeSet<Annotation>();
		
		
		List<Annotation> introns=new ArrayList<Annotation>();
		
		while(overlappers.hasNext()){
			introns.add(overlappers.next().getValue());
		}
		
		if(introns.isEmpty()){
			spliced.add(exon);
		}
		else{
			boolean counted=false;
			//define overlapping introns
			List<List<Annotation>> overlappingIntrons=defineOverlappingIntrons(introns); //TODO Why would any list be empty?
			//create sets of introns
			Collection<Annotation> setOfIntrons=createSetsOfIntrons(overlappingIntrons);
			//subtract each set and add the resulting exons
			for(Annotation setOfIntron: setOfIntrons){
				Annotation minus=exon.minus(setOfIntron);
				if(!minus.getBlocks().isEmpty()){
					spliced.addAll(minus.getBlocks()); //These are legit
					counted=true;
				}
			}
			
			//We'll add all the complements as possible pre-mRNA background to test
			if(!spliced.isEmpty()){
				Annotation comp=exon.minus(spliced);
				preSpliced.addAll(comp.getBlocks()); //These are likely preprocessed mRNA
			}
						
			
			if(!counted && this.alignsWithEdge(exon, introns)){
				logger.debug("aligns with edge "+exon.toUCSC());
				spliced.add(exon); //These are likely legit
				counted=true;
			}
			
			//We'll add all exons fully contained in an intron as possible pre-mRNA background
			//TODO We may just want to add all exons not already counted
			if(!counted && fullyContained(exon, introns)){preSpliced.add(exon);} //These are likely preprocessed mRNA	
		}
			
		Pair<Collection<Annotation>> rtrn=new Pair<Collection<Annotation>>(spliced, preSpliced);
		return rtrn;
	}

	private boolean notOverlapping(Annotation exon, Collection<Annotation> spliced) {
		for(Annotation splic: spliced){
			if(exon.overlaps(splic)){
				logger.debug(exon.toUCSC()+" overlaps");
				return false;
			}
		}
		return true;
	}

	private Collection<Annotation> createSetsOfIntrons(List<List<Annotation>> overlappingIntrons) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		//go through the overlapping introns and create sets
		int[] sizes=new int[overlappingIntrons.size()];
		for(int i=0; i<overlappingIntrons.size(); i++){
			sizes[i]=overlappingIntrons.get(i).size();
		}
		
		CombinationGenerator combinations=new CombinationGenerator(sizes);
		int[][] possibleCombos=combinations.getPermutations();
		
		//Rows are the combos to use
		for(int i=0; i<possibleCombos.length; i++){
			int[] setPositions=possibleCombos[i];
			Annotation intronSet=makeIntronSet(overlappingIntrons, setPositions);
			rtrn.add(intronSet);
		}
		
		return rtrn;
	}

	private Annotation makeIntronSet(List<List<Annotation>> overlappingIntrons,	int[] setPositions) {
		Collection<Annotation> blocks=new TreeSet<Annotation>();
		for(int i=0; i<setPositions.length; i++){
			Annotation intron=overlappingIntrons.get(i).get(setPositions[i]);
			blocks.add(intron);
		}
		return new BasicAnnotation(blocks);
	}

	private void write(List<List<Annotation>> overlappingIntrons, String save, Annotation exon) throws IOException {
		FileWriter writer=new FileWriter(save, true);
		
		int i=0;
		for(List<Annotation> introns: overlappingIntrons){
			for(Annotation intron: introns){
				//intron.setName("set"+i);
				writer.write(intron.getChr()+"\t"+intron.getStart()+"\t"+intron.getEnd()+"\t"+("set:"+i+"-"+exon.toUCSC())+"\n");
			}
			i++;
		}
		
		writer.close();
	}

	private List<List<Annotation>> defineOverlappingIntrons(List<Annotation> introns) {
		//TODO This does not account for exon exclusion events
		//To account for this, we will want to include introns that overlap both sets to be included in both sets
		
		List<List<Annotation>> rtrn=new ArrayList<List<Annotation>>();
		IntervalTree<Annotation> intronTree=makeTree(introns);
		Set<Annotation> accountedFor=new TreeSet<Annotation>();
		
		for(Annotation intron: introns){
			if(!accountedFor.contains(intron)){
				List<Annotation> overlappingIntrons=new ArrayList<Annotation>();
				//get overlappers
				List<Annotation> overlappers=getList(intronTree.overlappers(intron.getStart(), intron.getEnd()));
				add(overlappingIntrons, overlappers, accountedFor);
				for(Annotation ann2: overlappers){
					List<Annotation> overlappers2=getList(intronTree.overlappers(ann2.getStart(), ann2.getEnd()));
					add(overlappingIntrons, overlappers2, accountedFor);
				}
				rtrn.add(overlappingIntrons);
			}
		}
		
		//Lets do a second split to take each list and divide up bins that overlap
		rtrn=splitOverlapping(rtrn);
		
		return rtrn;
	}

	private List<List<Annotation>> splitOverlapping(List<List<Annotation>> set) {
		List<List<Annotation>> rtrn=new ArrayList<List<Annotation>>();
		//Lets go through each list
		//If the reads within it have junctions that dont overlap each other then lets split into new bins
		for(List<Annotation> list: set){
			rtrn.addAll(splitOverlappingForList(list));
		}
		return rtrn;
	}


	private List<List<Annotation>> splitOverlappingForList(List<Annotation> list) {
		//take this list and make sure that every element overlaps every other element
		//else split into n bins of overlapping bins
		
		
		Collection<Annotation> newList=new TreeSet<Annotation>();
		
		//List<Annotation> originalList=new ArrayList<Annotation>();
		
		for(int i=0; i<list.size(); i++){
			for(int j=(i+1); j<list.size(); j++){
				Annotation intron1=list.get(i);
				Annotation intron2=list.get(j);
				//if intron 1 and intron 2 don't overlap then split
				if(!intron1.overlaps(intron2)){
					//logger.debug(intron1.toUCSC()+" "+intron2.toUCSC()+" dont overlap");
					newList.add(intron1);
					newList.add(intron2);
				}
			}
		}
		
		//partition the new list
		List<List<Annotation>> newLists=partition(newList, list);
		
		return newLists;
	}

	private List<List<Annotation>> partition(Collection<Annotation> newList, Collection<Annotation> originalList) {
		Collection<Annotation> accountedFor=new TreeSet<Annotation>();
		
		List<List<Annotation>> lists=new ArrayList<List<Annotation>>();
		Map<Annotation, Collection<Annotation>> rtrn=new TreeMap<Annotation, Collection<Annotation>>();
		for(Annotation intron1: newList){
			Collection<Annotation> overlaps=new TreeSet<Annotation>();
			for(Annotation intron2: originalList){
				if(intron1.overlaps(intron2)){
					overlaps.add(intron2);
				}
			}
			rtrn.put(intron1, overlaps);
		}
		
		for(Annotation ref: rtrn.keySet()){
			boolean fullyCompat=isFullyCompatible(rtrn.get(ref));
			if(fullyCompat){
				if(!rtrn.get(ref).isEmpty()){
					lists.add(new ArrayList(rtrn.get(ref))); //TODO We should figure out a way to exclude the exact same lists from being added multiple times
				}
				accountedFor.addAll(rtrn.get(ref));
				//logger.debug(ref.toUCSC()+" "+rtrn.get(ref).size());
				//for(Annotation temp: rtrn.get(ref)){logger.debug(temp.toUCSC());}
			}
		}
		//TODO Make sure every intron is accounted for
		for(Annotation newIntron: newList){
			if(!accountedFor.contains(newIntron)){
				logger.debug("MISSING: "+newIntron.toUCSC());
			}
		}
		
		
		//Add the remainder of the originalList as a separate list
		List<Annotation> remainderList=new ArrayList<Annotation>();
		for(Annotation original: originalList){
			if(!accountedFor.contains(original)){remainderList.add(original);}
		}
		
		if(!remainderList.isEmpty()){
			lists.add(remainderList);
		}
		
		return lists;
	}

	private boolean isFullyCompatible(Collection<Annotation> collection) {
		//go through each annotation and ensure that nothing doesnt overlap with something else
		for(Annotation an1: collection){
			for(Annotation an2: collection){
				if(!an1.overlaps(an2)){return false;}
			}
		}
		return true;
	}

	private void allOverlap(List<Annotation> newList, List<Annotation> originalList, Annotation reference) {
		//Check that all overlap in list
		//if not, exclude from newList
		//if so, exclude from original list
		
		for(int i=0; i<newList.size(); i++){
			for(int j=(i+1); j<newList.size(); j++){
				Annotation intron1=newList.get(i);
				Annotation intron2=newList.get(j);
				if(!intron1.overlaps(intron2)){
					throw new IllegalStateException();
				}
				else{
					originalList.remove(intron1);
					originalList.remove(intron2);
				}
			}
		}
		
	}

	private List<Annotation> getList(Iterator<Node<Annotation>> overlappers) {
		List<Annotation> rtrn=new ArrayList<Annotation>();
		
		while(overlappers.hasNext()){
			rtrn.add(overlappers.next().getValue());
		}
		
		return rtrn;
	}

	private void add(List<Annotation> overlappingIntrons, List<Annotation> overlappers, Set<Annotation> accountedFor) {
		for(Annotation intron: overlappers){
			if(!overlappingIntrons.contains(intron)){
				overlappingIntrons.add(intron);
			}
			if(!accountedFor.contains(intron)){
				accountedFor.add(intron);
			}
		}
		
	}

	private Collection<Annotation> getRemainingIntrons(List<Annotation> introns, Collection<Annotation> exons) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//TODO Make more efficient
		//go through the exon set and see what introns are fully explained
		for(Annotation intron: introns){
			if(!hasExons(intron, exons)){
				rtrn.add(intron);
			}
		}
		
		return rtrn;
	}

	private boolean hasExons(Annotation intron, Collection<Annotation> exons) {
		boolean left=false;
		boolean right=false;
		
		for(Annotation exon: exons){
			if(exon.getStart()==intron.getEnd()){left=true;}
			if(exon.getEnd()==intron.getStart()){right=true;}
		}
		
		return left && right;
	}

	private Collection<? extends Annotation> dissect(Annotation exon, List<Annotation> introns) {
		logger.debug("dissecting "+exon+" "+introns.size());
		return exon.intersect(introns);
	}

	private boolean spans1Intron(Annotation exon, List<Annotation> introns) {
		Annotation block=new BasicAnnotation(introns);
		
		int counter=0;
		for(Annotation intron: block.getBlocks()){
			if(exon.overlaps(intron)){counter++;}
		}
		
		if(counter==1){return true;}
		
		return false;
		
	}

	private boolean alignsWithEdge(Annotation exon, List<Annotation> introns) {
		for(Annotation intron: introns){
			if(exon.getEnd()== intron.getStart() || exon.getStart()==intron.getEnd()){return true;}
		}
		return false;
	}

	private boolean fullyContained(Annotation exon, List<Annotation> introns) {
		for(Annotation intron: introns){
			if(intron.fullyContains(exon)){
				return true;
			}
		}
		return false;
	}

	private Collection<Annotation> intersect(List<Annotation> introns) {
		Collection<Annotation> rtrn=CollapseByIntersection.collapseByIntersection(introns, true);
		return rtrn;
	}

	private Collection<? extends Annotation> splitByIntron(Annotation exon, Iterator<Node<Annotation>> overlappers) {
		List<Annotation> introns=new ArrayList<Annotation>();
		while(overlappers.hasNext()){introns.add(overlappers.next().getValue());}
		
		//TODO minus the consensus intronic region
		
		return exon.disect(introns);
		//return exon.minus(introns).getBlocks();
	}

	private Collection<? extends Annotation> getExons(Alignment read) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		Collection<Annotation> readPairs=read.getReadAlignments(space); // this should return a collection for left and right
		for(Annotation align: readPairs){
			rtrn.addAll(align.getBlocks());
		}
		return rtrn;
	}
	
	private Collection<? extends Annotation> getJunctions(Alignment read) {
		//TODO This should only get spliced reads NOT indels
			
		return read.getSpliceConnections();
	}

	public Map<String, Collection<Gene>> getPaths() {
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		
		//For each chromosome
		for(String chr: this.graphs.keySet()){
			List<GraphPath<Annotation, TranscriptGraphEdge>> paths=this.graphs.get(chr).getPaths();
			Collection<Gene> genes=new TreeSet<Gene>();
			for(GraphPath<Annotation, TranscriptGraphEdge> path: paths){
				//logger.debug(path.toString());
				Gene gene=this.graphs.get(chr).pathToGene(path);
				genes.add(gene);
			}
			rtrn.put(chr, genes);
		}
		//Add orphan genes
		for(String chr: this.graphs.keySet()){
			for(Gene g:this.graphs.get(chr).getOrphanGenes()){
				rtrn.get(chr).add(g);
			}			
		}
		return rtrn;
	}
	
	/**
	 * Converts a list to an array
	 * @param list
	 * @return
	 */
	private double[] l2a(List<Double> list){
		double[] rtrn=new double[list.size()];
	
		int i=0;
		for(Double val: list){rtrn[i++]=val;}
	
		return rtrn;
	}
	
	/**
	 * Converts a collection of paths to a collection of genes
	 * @param paths to be converted
	 * @return
	 */
	private Collection<Gene> convert(Collection<Path> paths) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		for(Path path: paths){
			rtrn.add(path.toGene());
		}
		
		return rtrn;
	}

	/**
	 * This function will go through each gene and calculate the number of paired end reads 
	 * fully contained within each isoform for the gene
	 * @return
	 * @throws IOException 
	 */
	private Map<String,Collection<Gene>> setFPKMScores(Map<String,Collection<Gene>> geneMap,FileWriter writer,FileWriter writer2){
		
		Map<String,Collection<Gene>> filteredGenes = new HashMap<String,Collection<Gene>>();
		String name = "gene.v2.1_";
		for(String chr:geneMap.keySet()){
			//logger.debug("For chromosome "+chr+" graph gave "+geneMap.get(chr).size()+" genes");
			Collection<Gene> filtered = new ArrayList<Gene>();
			try{
			//MAKE A MAP OF GENE TO ISOFORMS
			Map<Gene,Set<Gene>> isoformMap = getIsoformMap(geneMap.get(chr));
			Map<Gene,Double> geneToPairedCount = new HashMap<Gene,Double>();
			
			//For each gene
			for(Gene gene:isoformMap.keySet()){
				writer2.write("Gene:\n");
				//logger.debug("Starting new gene");
				//For each transcript
				for(Gene isoform:isoformMap.get(gene)){
					
					double[] scores = getScores(isoform);
					double[] fields = new double[4];
					logger.debug(gene.toUCSC()+" "+scores[0]+"\t"+scores[1]);
					if(scores[1]<alpha){
						isoform.setName(name+new Double(counter).toString()+"_"+isoform.getChr());
						//[0] : sum
						fields[0] = scores[0];
						//[1] : p-value
						fields[1] = scores[1];
						//[2] : FPK
						fields[2] = (scores[0]*1000.0)/isoform.getSize();
						//[3] : FPKM
						//Calculate FPKM
						fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
						//logger.debug("For isoform : "+isoform.getName()+"\tNum of exons: "+isoform.getSpliceConnections().size()+"\t"+fields[0]+"\t"+fields[1]);
						isoform.setBedScore(fields[3]);
						isoform.setExtraFields(fields);
						counter++;
						writer.write(isoform+"\n");
						writer2.write(isoform.toBED()+"\t"+fields[0]+"\n");
						geneToPairedCount.put(isoform, fields[0]);
						filtered.add(isoform);
					}
					else{
						logger.debug("Gene "+gene.toUCSC()+" is filtered out because it does not meet the significance threshold : "+scores[1]);
					}
				}
			}
			} catch (IOException e) {
				e.printStackTrace();
			}
			logger.info("After applyPairedEndFilter "+filtered.size()+" genes");
			filteredGenes.put(chr, filtered);
		}
		return filteredGenes;
		
	}

	/**
	 * 
	 * @param genes
	 * @return
	 */
	public static Map<Gene,Set<Gene>> getIsoformMap(Collection<Gene> genes){
				
		Map<Gene,Set<Gene>> isoformMap = new HashMap<Gene,Set<Gene>>();
		
		for(Gene g:genes){
			//If gene has not already been processed
			if(!isoformMapContains(isoformMap,g)){
				//Add gene to its own isoform map
				Set<Gene> set = new HashSet<Gene>();
				set.add(g);
				Collection<Gene> potentialIsoforms = getOverlappingGenes(g,genes,0.5);
				for(Gene p:potentialIsoforms){
					if(!isoformMapContains(isoformMap,p)){
						set.add(p);
					}
				}
				isoformMap.put(g, set);
			}			
		}
		return isoformMap;
	}
	
	public static boolean isoformMapContains(Map<Gene,Set<Gene>> isoformMap,Gene gene){
	
		for(Set<Gene> set:isoformMap.values()){
			for(Gene g:set)
				if(g.equals(gene))
					return true;
		}
		return false;
	}
	/**
	 * 
	 * @param gene
	 * @param allGenes
	 * @param minPctOverlap
	 * @return
	 */
	public static Collection<Gene> getOverlappingGenes(Gene gene,Collection<Gene> allGenes, double minPctOverlap){
		
		Collection<Gene> overlappers = new HashSet<Gene>();
		for(Gene g:allGenes){
			if(gene.overlaps(g, minPctOverlap) && g.getOrientation().equals(gene.getOrientation()))
				overlappers.add(g);
		}
		return overlappers;
	}
	
	public static void main(String[] args)throws IOException{
		
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"reconstruct");
		
		TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand").equalsIgnoreCase("first")){
			//System.out.println("First read");
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(argMap.get("strand").equalsIgnoreCase("second")){
			//System.out.println("Second read");
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}
		else
			logger.info("no strand");
		
		if(argMap.containsKey("chr")){
			new BuildScriptureCoordinateSpace(new File(argMap.getMandatory("alignment")),argMap.getMandatory("genome"),argMap.getOutput(),true, strand,argMap.getMandatory("chr"),argMap);
		}
		else{
			new BuildScriptureCoordinateSpace(new File(argMap.getMandatory("alignment")),argMap.getMandatory("genome"),argMap.getOutput(),true, strand,argMap);
		}
/*		if(args.length>1){
			File bamFile=new File(args[0]);
			String genomeSeqFile = null;
			double threshold = new Double(args[1]);
			genomeSeqFile = args[2];
			if(args.length==3)
				new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3]);
			if(args.length==3)
				new BuildScriptureCoordinateSpace(bamFile,threshold,null,args[2],true, TranscriptionRead.SECOND_OF_PAIR);
			if(args.length>=5){
				TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
				if(args[4].equalsIgnoreCase("first")){
					//System.out.println("First read");
					strand = TranscriptionRead.FIRST_OF_PAIR;
				}
				else if(args[4].equalsIgnoreCase("second")){
					//System.out.println("Second read");
					strand = TranscriptionRead.SECOND_OF_PAIR;
				}
				else
					System.out.println("no strand");
				if(args.length==5){
					new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3],true, strand);
				}
				else{
					new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3],true, strand,args[5]);
				}
			}
		}
		else{System.err.println(usage);}*/
	}
	
	static final String usage = "Usage: BuildScriptureCoordinateSpace -task reconstruct "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Alignment file to be used for reconstruction.> "+
			"\n\n\t\t-genome <Fasta file with Genome sequence to be used as reference.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+		
			"\n\t\t-strand <VALUES: first, second, unstranded. Specifies the mate that is in the direction of transcription DEFAULT: Unstranded> "+
			
			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\n\t\t-chr <If specified, Scripture will be run for this chromosome only.> "+
			"\n\t\t-coverage <Specifies the minimum percentage of drop in coverage allowed for an exon. DEFAULT: 0.2> "+
			"\n\t\t-minSpliceReads <The minimum number of splice reads allowed to support a single intron transcript. DEFAULT: 3> "+
			"\n\t\t-percentSpliceReads <The minimum percentage of the average splice counts for a transcript, that an intron can be supported by. DEFAULT: 0.05> "+
			"\n\t\t-alpha <The significance p-value threshold for reconstructions. DEFAULT: 0.01> "+
			"\n";
	
	//static String usage=" args[0]=bam file \n\t args[1]=minimum percentage threshold for coverage \n\t args[2]: Fasta file with Genome sequence"
	//					+"\n\t args[3]= outputName \n\t args[4] transcription strand \n\targs[5] If specified only this chromosome";
	
}

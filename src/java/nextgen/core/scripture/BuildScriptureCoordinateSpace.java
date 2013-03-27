package nextgen.core.scripture;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
import org.jgrapht.GraphPath;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.Pair;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.Path;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.PairedEndAlignment.TranscriptionRead;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readFilters.CanonicalSpliceFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.IndelFilter;
import nextgen.core.readFilters.NoSpliceFilter;
import nextgen.core.readFilters.PairedAndProperFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.scripture.OrientedChromosomeTranscriptGraph.TranscriptGraphEdge;

public class BuildScriptureCoordinateSpace {

	static Logger logger = Logger.getLogger(BuildScriptureCoordinateSpace.class.getName());
	private AlignmentModel model;
	String genomeSeq = null;
	int windowSize=20000000;
	private CoordinateSpace space;
	private Map<String, ChromosomeTranscriptGraph> graphs;
	private boolean forceStrandSpecificity=true; //TODO This should be passed or at least determined from data
	private static double DEFAULT_MIN_COV_THRESHOLD = 0.5;
	private static double MIN_SPLICE_PERCENT = 0.40;
	private double coveragePercentThreshold = DEFAULT_MIN_COV_THRESHOLD;
	private static TranscriptionRead DEFAULT_TXN_READ = TranscriptionRead.UNSTRANDED;
	String outName = null;
	private SAMFileHeader header;
	private long spliceCount;
	private long nonSpliceCount;
	private double DEFAULT_ALPHA = 0.05;
	private static double MIN_SPLICE_READS = 2.0;
	private double THRESHOLD_SPURIOUS = 95.0;
	
	//TODO: Handle genome sequence file
	public BuildScriptureCoordinateSpace(File bamFile){
		this(bamFile,DEFAULT_MIN_COV_THRESHOLD,null,null,true,DEFAULT_TXN_READ);
		logger.info("Genome sequence has not been provided");
	}
	
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir){
		this(bamFile,DEFAULT_MIN_COV_THRESHOLD,genomeDir,null,true,DEFAULT_TXN_READ);
	}
	
	public BuildScriptureCoordinateSpace(File bamFile,double threshold,String genomeDir,String outputName){
		//By default first read is transcription read
		this(bamFile,threshold,genomeDir,outputName,true,DEFAULT_TXN_READ);
	}
	
	/**
	 * 
	 * @param bamFile
	 * @param threshold
	 * @param genomeDir
	 * @param outputName
	 * @param forceStrandedness
	 */
	public BuildScriptureCoordinateSpace(File bamFile,double threshold,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead strand){
			
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		forceStrandSpecificity = forceStrandedness;
		this.model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand);
		header = new SAMFileReader(bamFile).getFileHeader();
		model.addFilter(new ProperPairFilter());
		model.addFilter(new IndelFilter());
		model.addFilter(new GenomicSpanFilter(20000000));
		this.space=model.getCoordinateSpace();
		spliceCount = 0;
		nonSpliceCount = 0;
		//logger.info("Count: "+model.getCount(new BasicAnnotation("chr6",52062677,52064112,Strand.NEGATIVE),true));
		outName = outputName;
		coveragePercentThreshold = threshold;
		
		logger.info("Parameters used: " +
				"\nSpurious Filter: "+THRESHOLD_SPURIOUS+
				"\nPremature Assembly Filter: "+coveragePercentThreshold+
				"\nSplice junction Filter : "+
				"\n\tNumber of spliced reads : "+MIN_SPLICE_READS+
				"\n\tPercentage of total spliced reads: "+MIN_SPLICE_PERCENT+
				"\nAlpha for single exon assemblies : "+this.DEFAULT_ALPHA
				);
		
		assemble(strand);
	//	ChromosomeTranscriptGraph graph=assemble("chr19",strand);
	//	graphs.put("chr19", graph);
	//	assemble();
		
		Map<String, Collection<Gene>> rtrn=this.getPaths();
		try {
			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			for(String chr: rtrn.keySet()){
				applyPairedEndFilter(rtrn.get(chr),writer,writer2);
			}
			writer.close();
			writer2.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		try {
			write(outName+".paths.bed",rtrn);
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
	public BuildScriptureCoordinateSpace(File bamFile,double threshold,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead strand,String chr){
			
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		forceStrandSpecificity = forceStrandedness;
		this.model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand);
		header = new SAMFileReader(bamFile).getFileHeader();
		model.addFilter(new ProperPairFilter());
		model.addFilter(new IndelFilter());
		model.addFilter(new GenomicSpanFilter(20000000));
		this.space=model.getCoordinateSpace();
		spliceCount = 0;
		nonSpliceCount = 0;
		//logger.info("Count: "+model.getCount(new BasicAnnotation("chr6",52062677,52064112,Strand.NEGATIVE),true));

		outName = outputName;
		coveragePercentThreshold = threshold;
		//assemble(strand);
		ChromosomeTranscriptGraph graph=assemble(chr,strand);
		graphs.put(chr, graph);
		
		Map<String, Collection<Gene>> rtrn=this.getPaths();
		try {
			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			for(String ch: rtrn.keySet()){
				applyPairedEndFilter(rtrn.get(ch),writer,writer2);
			}
			writer.close();
			writer2.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		try {
			write(outName+".paths.bed",rtrn);
		} catch (IOException e) {
			e.printStackTrace();
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

	private ChromosomeTranscriptGraph assemble(String chr,TranscriptionRead strand) {
		
		//Option 1: Scan through the space and collapse into compatible edges and nodes
		ChromosomeTranscriptGraph graph=assembleDirectly(chr,strand);
		
		return graph;
	}


	private void write(String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<String, Collection<Gene>> rtrn=this.getPaths();
		
		for(String chr: rtrn.keySet()){
			Collection<Gene> genes=rtrn.get(chr);
			for(Gene gene: genes){
				writer.write(gene+"\n");
			}
		}
		
		writer.close();
	}

	private void write(String save, Map<String, Collection<Gene>> rtrn) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String chr: rtrn.keySet()){
			Collection<Gene> genes=rtrn.get(chr);
			for(Gene gene: genes){
				writer.write(gene+"\n");
			}
		}
		
		writer.close();
	}
	
	private ChromosomeTranscriptGraph assembleDirectly(String chr,TranscriptionRead strand){
		
		logger.info("Assembling spliced reads");
		CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new CanonicalSpliceFilter(genomeSeq));

		IntervalTree<Assembly> splicedAssemblies=assembleDirectly(splicedIter,strand);
		//logger.info(model.getCount());
		try{
			write(splicedAssemblies, outName+"."+chr+"."+"01splicedAssemblies.bed");
			}catch(IOException ex){}
		
		logger.info("Size of spliced assemblies: "+splicedAssemblies.size());
		
		logger.info("Remove low splice read-supported assemblies");
		//REMOVE SPURIOUS
		IntervalTree<Assembly> notLowSplicedAssemblies = removeLowSplicedAssemblies(splicedAssemblies,chr);
		try{
			write(notLowSplicedAssemblies, outName+"."+chr+"."+"01BsplicedAssemblies.bed");
			}catch(IOException ex){}
		
		logger.info("Assembling all reads : working assemblies");
		CloseableIterator<Alignment> iter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new NoSpliceFilter());
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, notLowSplicedAssemblies,strand);
		
		logger.info("Size of direct assemblies: "+workingAssemblies.size());
		
		
		try{
		write(workingAssemblies, outName+"."+chr+"."+"02directAssemblies.bed");
		}catch(IOException ex){}
		
		logger.info("Remove spurious assemblies");
		//REMOVE SPURIOUS
		IntervalTree<Assembly> unspuriousAssemblies = removeSpuriousAssemblies(workingAssemblies,chr);
		try{
			write(unspuriousAssemblies, outName+"."+chr+"."+"03unspuriousAssemblies.bed");
			}catch(IOException ex){}
									
		logger.info("Merge assemblies");
		mergeAssembly(unspuriousAssemblies);
		try{write(unspuriousAssemblies, outName+"."+chr+"."+"04mergedAssemblies.bed");}catch(IOException ex){}
		
		logger.info("Branch compatible assemblies");
		//TODO Branch incompatible assemblies
		branchAssembly(unspuriousAssemblies);
		try{write(unspuriousAssemblies, outName+"."+chr+"."+"05branchedAssemblies.bed");}catch(IOException ex){}
				
		logger.info("Remove premature assemblies");
		// Flag and remove premature
		IntervalTree<Assembly> filteredAssemblies = removePrematureAssemblies(unspuriousAssemblies);
		//IntervalTree<Assembly> filteredAssemblies = removePrematureAssembliesUsingPairedEnds(unspuriousAssemblies);
		try{
			write(filteredAssemblies, outName+"."+chr+"."+"06filteredAssemblies.bed");
			}catch(IOException ex){}

		
		IntervalTree<Assembly> highSplicedAssemblies = removeLowSpliceJunctionAssemblies(filteredAssemblies);
		try{write(highSplicedAssemblies, outName+"."+chr+"."+"07highSplicedAssemblies.bed");}catch(IOException ex){}
		
		//Lets split potential preprocessed and mature transcripts and test whether to include certain nodes
		
		//Make graph
		ChromosomeTranscriptGraph graph=makeGraph(highSplicedAssemblies, chr); //TODO Should use the merged set
		return graph;
	}

	private void branchAssembly(IntervalTree<Assembly> tree) {
		//We have a set of assemblies that are all incompatible
		//We want to link up parts
		logger.info("Enter branch assembly");
		//iterate through and get overlapping assemblies
		Iterator<Assembly> iter=tree.toCollection().iterator();
		Collection<Assembly> considered = new ArrayList<Assembly>();
		
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
				//logger.info("overlaps "+assembly2.toUCSC());
				if(!considered.contains(assembly2)){
					//try and merge the non-overlapping portions
					Collection<Assembly> merged=branchAssemblies(assembly1, assembly2);
					if(!merged.isEmpty()){
						logger.info("Branching "+assembly1.toUCSC()+" and "+assembly2.toUCSC());
						
						//remove assembly1
						tree.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//remove assembly2
						tree.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);
						
						considered.remove(assembly1);
						considered.remove(assembly2);
						//currentAssemblies.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//currentAssemblies.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);
						//add merged
						//logger.info("No. of assemblies "+ merged.size());
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
		
		//if this is compatible with second assemebly
		if(compatible(orderedAssembly.getValue2(), portionToConsiderAdding)){
			Assembly merged=merge(portionToConsiderAdding, orderedAssembly.getValue2());
			rtrn.add(merged);
			rtrn.add(orderedAssembly.getValue1());
			//System.out.println(merged);
		}
		
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
					//Add the exon as a vertex to the graph
					//logger.error("Assembly "+gene.toUCSC()+" has one exon. Attempting to add the one exon to the graph");
					//CHECK IF THE VERTEX MEETS THE SIGNIFICANCE THRESHOLD
					double pval = new ScanStatisticScore(model,blocks.get(0)).getScanPvalue();
					if(pval<DEFAULT_ALPHA){
						graph.connectVertexToGraph(blocks.get(0));
					}
					else{
						logger.info(blocks.get(0).toUCSC()+" is not added to the graph because it does not pass the p-value threshold.");
					}
				}
				else{
					graph.addAnnotationToGraph(gene);
				}
			}
		}
		return graph;
	}

	private void mergeAssembly(IntervalTree<Assembly> tree) {
		//Step 1: Iterate through the working assembly and find paths that overlap
		Iterator<Assembly> iter=tree.toCollection().iterator();
		//if overlaps but is incompatible, see if we can branch
		//iterate
		Collection<Assembly> considered = new ArrayList<Assembly>();
 		logger.info("Enter merging");
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
						logger.info("Merging: "+ branch1.toUCSC()+ " and "+branch2.toUCSC());
						Assembly merged=merge(branch1, branch2);
						//remove annotation1 and annotation2
						tree.remove(branch1.getStart(), branch1.getEnd(), branch1);
						tree.remove(branch2.getStart(), branch2.getEnd(), branch2);
						//add merged
						tree.put(merged.getStart(), merged.getEnd(), merged);
					}
				}
			}
			considered.add(branch1);
		}
	}
	
	/**
	 * Remove all assemblies that are supported by 
	 * @param assemblies
	 * @param chr
	 * @return
	 */
	private IntervalTree<Assembly> removeLowSplicedAssemblies(IntervalTree<Assembly> assemblies,String chr) {
		
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		//For each assembly
		while(iter.hasNext()){
			Assembly assembly=iter.next();

			//IF THE GENE HAS ONE INTRON, CHECK IF IT HAS >1 SPLICED READS
			if(assembly.getSpliceConnections().size()==1){
				for(Annotation intron:assembly.getSpliceConnections()){
					double count=model.getIntronCounts(intron);
					if(count>MIN_SPLICE_READS){
						currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
					}
					else{
						logger.info(assembly.toUCSC()+" is removed because intron "+intron.toUCSC()+" has "+count+" supporting read.");
					}
				}
			}
			else{
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		return currentAssemblies;
	}
	
	/**
	 * This function removes any overlapping bad isoforms
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removeSpuriousAssemblies(IntervalTree<Assembly> assemblies,String chr) {
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
				while(overlappers.hasNext()){
					Assembly assembly2 = overlappers.next();
					//For all exons of this assembly
					for(Annotation intron1: assembly1.getSpliceConnections()){
						if(!toRemoveAssembly1){
							for(Annotation exon2:assembly2.getBlocks()){
								//If intron in assembly1 is contained in exon of assembly2
								if(exon2.contains(intron1)){
									//CHECK FOR COVERAGE
									double intronCount = model.getIntronCounts(intron1);
									double exonCount = (model.getCount(new BasicAnnotation(exon2.getChr(), intron1.getStart(), intron1.getEnd(), exon2.getOrientation()), false)
																/(double)(intron1.getEnd()-intron1.getStart()))*1000;
									double ratio = intronCount/exonCount;
									//CASE 1 : Assembly with intron is mature and assembly with exon is immature
									if(ratio>THRESHOLD_SPURIOUS){
										logger.info("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.toUCSC()+ "to intron "+intron1.toUCSC()+" in assembly "+assembly1.toUCSC());
										assembly2.setSpurious(true);
										logger.info("CASE 1: assembly "+ assembly2.toUCSC()+" is spurious. Intron count "+intronCount+" > "+exonCount+" Ratio: "+ratio);
									}
									else{
										//CASE 2: Assembly with exon is real and the other 
										if(ratio<(1.0/THRESHOLD_SPURIOUS)){
											//logger.warn(assembly.toUCSC()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC());
											logger.info("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.toUCSC()+ "to intron "+intron1.toUCSC()+" in assembly "+assembly1.toUCSC());
											assembly1.setSpurious(true);
											toRemoveAssembly1 = true;
											logger.info("CASE 2: assembly "+ assembly1.toUCSC()+" is spurious.Exon count"+(model.getCount(new BasicAnnotation(exon2.getChr(), intron1.getStart(), intron1.getEnd(), exon2.getOrientation()), false))+" > "+intronCount+" Ratio: "+(ratio));
											break;
										}
										else{
											//Both are real
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
				//For all exons of this assembly
				for(Annotation exon1: assembly.getBlocks()){
					if(!toRemove){
						//If exon overlaps introns of overlapping assembly
						for(Annotation intron2:overlapper.getSpliceConnections()){
							//Flag as premature				
							if(exon1.overlaps(intron2)){
								assembly.setPossiblePremature(true);
								//CHECK FOR COVERAGE
								if(coveragePassesCheck(exon1,intron2)){
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
	}
	
	/**
	 * This function removes any premature assemblies that have coverage conflicts in the exons
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removePrematureAssembliesUsingPairedEnds(IntervalTree<Assembly> assemblies) {
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
	}
	
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
		int intronLength = flankingExons[0].length()+flankingExons[1].length();
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
				logger.info("Paired end count for exon " +exon1.toUCSC()+" : "+exonCount);
				logger.info("Paired end count for intron " +intron2.toUCSC()+" : "+intronCount+ " ratio = "+ratio);
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
			
			Map<Annotation,Double> intronToSplicedCountMap = new HashMap<Annotation,Double>();
			double avgCount = 0.0;
			
			boolean toRemove = false;
			//IF THE GENE HAS ONE INTRON, NOTHING
			if(assembly.getSpliceConnections().size()==1){
			/*	for(Annotation intron:assembly.getSpliceConnections()){
					double count=model.getIntronCounts(intron);
					if(count<=MIN_SPLICE_READS){
						toRemove = true;
					}
				}*/
			}
			else{
				//for each intron
				for(Annotation intron:assembly.getSpliceConnections()){
					//find the number of supporting splice junctions
					double count=model.getIntronCounts(intron);
												
					avgCount +=count;
					intronToSplicedCountMap.put(intron, count);
				}
				
				avgCount = (avgCount / (double)intronToSplicedCountMap.keySet().size());
				
				//logger.error(assembly.toUCSC()+" has average count "+avgCount);
				
				toRemove = false;
				//Check for all introns
				for(Annotation intron:intronToSplicedCountMap.keySet()){
					if(intronToSplicedCountMap.get(intron)<=avgCount*MIN_SPLICE_PERCENT){
						logger.error(assembly.toUCSC()+" removed because intron "+intron.toUCSC()+" has coverage "+intronToSplicedCountMap.get(intron)+" compared to "+avgCount);
						toRemove = true; 
						break;
					}
					else{
						//logger.error(intron.toUCSC()+" has coverage "+intronToSplicedCountMap.get(intron)+" passes");
					}
				}
			}
			if(!toRemove){
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		return currentAssemblies;
	}

	/**
	 * Helper function to removePrematureAssemblies
	 * @param exon1
	 * @param intron2
	 * @return
	 */
	private boolean coveragePassesCheck(Annotation exon1, Annotation intron2){
		Annotation overlap = exon1.intersect(intron2);
		Annotation nonoverlap = exon1.minus(intron2);
		
		if(overlap==null){
			return true;
		}
		//If fully contained
		//TODO: A minimum coverage threshold?
		if(overlap.equals(exon1)){
			//logger.error("Fully Contained");
			return true;
		}
		
		double overlapScore = new ScanStatisticScore(model, overlap).getAverageCoverage();
		double nonoverlapScore = new ScanStatisticScore(model, nonoverlap).getAverageCoverage();
		
		//logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
		if(overlapScore>=nonoverlapScore*coveragePercentThreshold){
			//logger.error(exon1.toUCSC()+" passes coverage test with intron "+ intron2.toUCSC());
			//logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
			return true;
		}
		else{
			//logger.error(exon1.toUCSC()+" does not pass coverage test with intron "+ intron2.toUCSC());
			return false;
		}
	}

	private void remove(IntervalTree<Assembly> tree, Collection<Assembly> toRemove) {
		for(Assembly assembly: toRemove){
			tree.remove(assembly.getStart(), assembly.getEnd(), assembly);
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
		boolean flagPremature=!workingAssemblies.isEmpty();
		while(iter.hasNext()){
			
			Alignment reads=iter.next();
			reads.setFragmentStrand(strand);	
			//boolean isfragmentspliced = false;
			
			//For the assembly, we need to treat each read separately
			for(Annotation read: reads.getReadAlignments(space)){
				
				//EACH READ HAS THE FRAGMENT STRAND
				//for each read, get overlapping assemblies
				Iterator<Node<Assembly>> overlappers=workingAssemblies.overlappers(read.getStart(), read.getEnd());

				/*if(read.getSpliceConnections().isEmpty())
					isfragmentspliced=true;*/
				//if no overlappers add read as assembly
				if(!overlappers.hasNext()){
					//add the read as an annotation
					//Flag this as likely premature
					Assembly readAssembly=new Assembly(read, false);
					if(flagPremature){
						readAssembly.setPossiblePremature(true);
						//logger.info("premature");
					}
					workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
				}
				else{
					boolean hasCompatible=false;
					//store the overlappers
					//else, test compatability between read and overlapping assembly
					while(overlappers.hasNext()){
						Collection<Assembly> assemblies=new TreeSet<Assembly>(overlappers.next().getContainedValues());
						for(Assembly assembly: assemblies){
							boolean isCompatible=compatible(assembly, read);
							if(isCompatible){
								//if compatable --> merge
								Assembly merged=merge(assembly, read);
								//remove assembly
								workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
								//add merged
								workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
								//set has compatable to true
								hasCompatible=true;
							}
						}
					}
					if(!hasCompatible){
						//add the read as an annotation
						//Flag this as likely premature
						Assembly readAssembly=new Assembly(read, false);
						if(flagPremature){
							readAssembly.setPossiblePremature(true);
						}
						workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
					}
				}
			}
/*			if(isfragmentspliced){
				spliceCount++;
			}
			else{
				nonSpliceCount++;
			}*/
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
	private Collection<Assembly> splitBranch(Assembly assembly, Annotation read) {
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		
		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		if(read.overlaps(assembly)){
			Assembly truncated=findMaxCompatibility(assembly, read);
			if(truncated!=null){
				Assembly merged=merge(truncated, read);
				rtrn.add(merged);
				rtrn.add(assembly);
			}
		}
		return rtrn;
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
	private boolean compatible(Annotation assembly, Annotation read) {
		//Two alignments will be defined as compatible if:
		//(i) the intronic locations are exactly same
		//if both have introns, ensure that there are no non-overlapping introns
		//logger.info(read.getName());
		if(!assembly.getSpliceConnections().isEmpty() 
				|| !read.getSpliceConnections().isEmpty()){
			//check if introns are compatible
			if(areIntronsCompatible(assembly, read)){
				if(!assembly.getOrientation().equals(read.getOrientation()))
					logger.error("Non compatible units are returning TRUE");

				return true;
			}
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
	

	private boolean areIntronsCompatible(Annotation assembly, Annotation read) {
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//introns are compatible if:
		//(i) all overlapping introns are identical
		for(Annotation intron1: assemblyIntrons){
			for(Annotation intron2: readIntrons){
				//if overlaps but not identical
				if(intron1.overlaps(intron2,forceStrandSpecificity)){
					if(!intron1.equals(intron2, forceStrandSpecificity)){return false;} //This should use strand info or not based on flag
				}
			}
		}
		
		//(ii) if none of the exons overlap an intron in the other
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation intron2: readIntrons){
				if(exon1.overlaps(intron2,forceStrandSpecificity)){return false;}
			}
		}
		
		for(Annotation exon2: read.getBlocks()){
			for(Annotation intron1: assemblyIntrons){
				if(exon2.overlaps(intron1,forceStrandSpecificity)){return false;}
			}
		}
		
		//(ii) the introns dont overlap but the exons do
		//just need to test that any exons overlap
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation exon2: read.getBlocks()){
				if(exon1.overlaps(exon2,forceStrandSpecificity)){return true;}
			}
		}
		
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
					if(!intron1.equals(intron2, this.forceStrandSpecificity)){
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

	private boolean overlap(Annotation assembly, Annotation read) {
		return assembly.overlaps(read);
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

	private Assembly merge(Assembly assembly, Annotation read) {
		//since these are compatable, we will simply add the read to the assembly by merging exons
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		rtrn.addAll(read.getBlocks());
		rtrn.addAll(assembly.getBlocks());
		
		return new Assembly(rtrn, assembly.getPossiblePremature());
	}
	
	private Assembly merge(Annotation assembly, Annotation read) {
		//since these are compatable, we will simply add the read to the assembly by merging exons
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		rtrn.addAll(read.getBlocks());
		rtrn.addAll(assembly.getBlocks());
		
		return new Assembly(rtrn);
	}
	
	private OrientedChromosomeTranscriptGraph assemble(CloseableIterator<Alignment> reads, String chr) {
		//Convert reads into collection of reads
		Collection<Alignment> readCollection=new TreeSet<Alignment>();
		
		//Also keep collection of spans from paired ends
		Collection<Annotation> spans=new TreeSet<Annotation>();
		
		//Also keep collection of junctions
		Collection<Annotation> junctions=new TreeSet<Annotation>();
		
		//Also keep just the exonic regions
		Collection<Annotation> exons=new TreeSet<Annotation>();
		
		logger.info("Iterating through reads");
		int counter=0;
		while(reads.hasNext()){
			Alignment read=reads.next();
			
			//the raw read
			readCollection.add(read);
			
			//get span from paired end
			spans.addAll(read.getFragment(space));
			
			//get junctions from reads
			junctions.addAll(getJunctions(read));
			
			//get exons from reads
			exons.addAll(getExons(read));
			counter++;
			if(counter%1000 ==0){logger.info(counter);}
		}
		
		reads.close();
		
		//First, lets collapse all exons
		Annotation exonicRegion=new BasicAnnotation(exons);
		
		logger.info("Splitting exons by intron locations");
		//Second, lets split exons by junctions
		Pair<Collection<Annotation>> exonAssembly=splitExonsByJunctions(exonicRegion, junctions);
				
		//Third, lets filter pre-spliced mRNA
		exons=filter(exonAssembly); //TODO Finish implementing this
		
		logger.info(exons.size()+" "+junctions.size());
		
		//Fourth, build a graph from the exons and introns
		/*ChromosomeTranscriptGraph graph=createGraph(chr, exons, junctions);
		return graph;*/
		
		//Fifth, lets partition seperate gene graphs by paired end spans
		//if has pairs then partition the graph by spans
		//else simply return the contiguous spans
		
		//TEMP
		try {
			write(exonicRegion.getBlocks(), "allExons.bed");
			write(exons, "decollapsedExons.bed");
			write(junctions, "junctions.bed");
			write(spans, "spans.bed");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		//TODO Temporary return
		return null;
	}

	

	private Collection<Annotation> filter( Pair<Collection<Annotation>> exonAssembly) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		//We will compute the coverage for the spliced regions and for each exon in prespliced
		Collection<Annotation> spliced=exonAssembly.getValue1();
	
		//TODO If the prespliced connects to a existing exon we will keep current and make a merged version
		
		//Add all from spliced
		rtrn.addAll(spliced);
		
		//return the collection
		return rtrn;
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
				logger.info("aligns with edge "+exon.toUCSC());
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
				logger.info(exon.toUCSC()+" overlaps");
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
					//logger.info(intron1.toUCSC()+" "+intron2.toUCSC()+" dont overlap");
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
				//logger.info(ref.toUCSC()+" "+rtrn.get(ref).size());
				//for(Annotation temp: rtrn.get(ref)){logger.info(temp.toUCSC());}
			}
		}
		//TODO Make sure every intron is accounted for
		for(Annotation newIntron: newList){
			if(!accountedFor.contains(newIntron)){
				logger.error("MISSING: "+newIntron.toUCSC());
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
		logger.info("dissecting "+exon+" "+introns.size());
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

	private Collection<? extends Annotation> split(Annotation exon,	Annotation intron) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		rtrn.add(exon.minus(intron));
		//rtrn.add(exon);
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
				logger.info(path.toString());
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
	 * This function will go through each gene
	 * @return
	 * @throws IOException 
	 */
	private void applyPairedEndFilter(Collection<Gene> genes,FileWriter writer,FileWriter writer2){
		
		try{
		//MAKE A MAP OF GENE TO ISOFORMS
		Map<Gene,Set<Gene>> isoformMap = getIsoformMap(genes);
		Map<Gene,Double> geneToPairedCount = new HashMap<Gene,Double>();
		
		//For each gene
		for(Gene gene:isoformMap.keySet()){
			writer2.write("Gene:\n");
			//For each transcript
			for(Gene isoform:isoformMap.get(gene)){
				
				logger.info("For isoform: "+isoform.toUCSC());
				double count = 0.0;
				//Get all reads overlapping the transcript
				CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(isoform,true), new PairedAndProperFilter());
				//For each read,
				while(iter.hasNext()){					
					Alignment read = iter.next();
					boolean countRead = true;
					for(Annotation mate:read.getReadAlignments(space)){
						if(!compatible(isoform,mate)){
							logger.info("Read "+mate.toUCSC()+" is not compatible with isoform with "+isoform.getExons().length);
							countRead=false;
						}
					}
					//For the assembly, we need to treat each read separately	
					if(countRead){
						count++;
					}
				}
				iter.close();
				logger.info(count);
				isoform.setName(new Double(count).toString());
				writer.write(isoform+"\n");
				writer2.write(isoform.toUCSC()+"\t"+count+"\n");
				geneToPairedCount.put(isoform, count);
			}
		}
		writer.close();
		writer2.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	/**
	 * 
	 * @param genes
	 * @return
	 */
	private Map<Gene,Set<Gene>> getIsoformMap(Collection<Gene> genes){
		
		Map<Gene,Set<Gene>> isoformMap = new HashMap<Gene,Set<Gene>>();
		
		for(Gene g:genes){
			//If gene has not already been processed
			if(!isoformMapContains(isoformMap,g)){
				//Add gene to its own isoform map
				Set<Gene> set = new HashSet<Gene>();
				set.add(g);
				Collection<Gene> potentialIsoforms = getOverlappingGenes(g,genes,0.5);
				for(Gene p:potentialIsoforms){
					if(!isoformMap.values().contains(p)){
						set.add(p);
					}
				}
				isoformMap.put(g, set);
			}
		}
		return isoformMap;
	}
	
	private boolean isoformMapContains(Map<Gene,Set<Gene>> isoformMap,Gene gene){
	
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
	private Collection<Gene> getOverlappingGenes(Gene gene,Collection<Gene> allGenes, double minPctOverlap){
		
		Collection<Gene> overlappers = new HashSet<Gene>();
		for(Gene g:allGenes){
			if(gene.overlaps(g, minPctOverlap))
				overlappers.add(g);
		}
		return overlappers;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File bamFile=new File(args[0]);
			String genomeSeqFile = null;
			double threshold = new Double(args[1]);
			genomeSeqFile = args[2];
			if(args.length==3)
				new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3]);
			if(args.length>=5){
				TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
				if(args[4].equalsIgnoreCase("first")){
					System.out.println("First read");
					strand = TranscriptionRead.FIRST_OF_PAIR;
				}
				else if(args[4].equalsIgnoreCase("second")){
					System.out.println("Second read");
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
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n\t args[1]=minimum percentage threshold for coverage \n\t args[2]: Fasta file with Genome sequence"
						+"\n\t args[3]= outputName \n\t args[4] transcription strand \n\targs[5] If specified only this chromosome";
	
}

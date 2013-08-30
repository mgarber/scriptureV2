package nextgen.core.scripture.statistics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.PairedAndProperFilter;
import nextgen.core.readFilters.SameOrientationFilter;
import nextgen.core.readFilters.SplicedReadFilter;
import nextgen.core.scripture.BuildScriptureCoordinateSpace;
import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.Pair;
import broad.core.error.ParseException;
import broad.core.math.ScanStatistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

public class CalculateSignificance {
	
	static Logger logger = Logger.getLogger(CalculateSignificance.class.getName());
	Map<String,Collection<Gene>> annotations;
	private AlignmentModel model;
	private static double DEFAULT_ALPHA = 0.05;
	private CoordinateSpace space;
	int counter = 1000;
	
	static final String usage = "Usage: CalculateSignific -task doWork "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-in <Reconstruction bed file. [BED by default]> "+
			"\n\t\t-alignment <Alignment bam file with data on which reconstructions were calculated> "+
			"\n\t\t-strand <first/second : mate in the direction of transcription> "+
			"\n\t\t-maxSize <Maximum insert size. Default=500bp> "+
			"\n";

	public CalculateSignificance(String annotationFile,String outputName,File bamFile,TranscriptionRead strand) throws IOException{
		
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		model.addFilter(new PairedAndProperFilter());
		
		FileWriter writer=new FileWriter(outputName+".09pairedGenes.bed");
		Map<String, Collection<Gene>> genes = setFPKMScores(annotations,writer);
		writer.close();
		space=model.getCoordinateSpace();

		BuildScriptureCoordinateSpace.write(outputName+".10final.paths.bed",genes);
		
		postProcess(genes,outputName);
		
	}
	
	/**
	 * This function will go through each gene and calculate the number of paired end reads 
	 * fully contained within each isoform for the gene
	 * @return
	 * @throws IOException 
	 */
	private Map<String,Collection<Gene>> setFPKMScores(Map<String,Collection<Gene>> geneMap,FileWriter writer){
		
		Map<String,Collection<Gene>> filteredGenes = new HashMap<String,Collection<Gene>>();
		String name = "gene.v1.3_";
		for(String chr:geneMap.keySet()){
			//logger.info("For chromosome "+chr+" graph gave "+geneMap.get(chr).size()+" genes");
			Collection<Gene> filtered = new ArrayList<Gene>();
			try{
			//MAKE A MAP OF GENE TO ISOFORMS
//			Map<Gene,Set<Gene>> isoformMap = BuildScriptureCoordinateSpace.getIsoformMap(geneMap.get(chr));			
			//For each gene
			for(Gene gene:geneMap.get(chr)){					
					double[] scores = getScores(gene);
					double[] fields = new double[4];
					if(scores[1]<DEFAULT_ALPHA){
						//logger.info(count);
						gene.setName(name+new Double(counter).toString());
						//[0] : sum
						fields[0] = scores[0];
						//[1] : p-value
						fields[1] = scores[1];
						//[2] : FPK
						fields[2] = (scores[0]*1000.0)/gene.getSize();
						//[3] : FPKM
						//Calculate FPKM
						fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
						//logger.info("For isoform : "+isoform.getName()+"\tNum of exons: "+isoform.getSpliceConnections().size()+"\t"+fields[0]+"\t"+fields[1]);
						gene.setBedScore(fields[3]);
						gene.setExtraFields(fields);
						counter++;
						writer.write(gene+"\n");
						filtered.add(gene);
					}
					else{
						logger.info("Gene "+gene.toUCSC()+" is filtered out because it does not meet the significance threshold : "+scores[1]);
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
				if(!BuildScriptureCoordinateSpace.compatible(gene,mate)){
					//logger.info("Read "+mate.toUCSC()+" is not compatible with isoform with "+isoform.getExons().length);
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
		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		
		return scores;
	}
	
	private void postProcess(Map<String, Collection<Gene>> oldGenes,String outputName) throws IOException{
		
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
						if(!branch1.equals(branch2) && BuildScriptureCoordinateSpace.compatible(branch1, branch2)){
							logger.info("Merging: "+ branch1.getName()+ " and "+branch2.getName());
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
			Map<Gene,Set<Gene>> isoformMap = BuildScriptureCoordinateSpace.getIsoformMap(tree.toCollection());		
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
	
		BuildScriptureCoordinateSpace.write(outputName+".11postprocessed.paths.bed",newGenes);
	}
	
	
	public static void main (String [] args) throws ParseException, IOException {
		
		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"doWork");
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
			System.out.println("no strand");
		
		
		new CalculateSignificance(argMap.getInput(),argMap.getOutput(),new File(argMap.getMandatory("alignment")),strand);
		
	}
}

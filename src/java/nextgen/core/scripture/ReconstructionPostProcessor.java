package nextgen.core.scripture;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;


public class ReconstructionPostProcessor {
	
	static Logger logger = Logger.getLogger(ReconstructionPostProcessor.class.getName());
	Map<String,Collection<Gene>> annotations;
	
	static final String usage = "Usage: ReconstructionPostProcessor -task doWork "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-input <Reconstruction bed file. [BED by default]> "+
			"\n";
	
	public ReconstructionPostProcessor(String annotationFile,String outputName) throws IOException{
		
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		doWork(outputName);
	}
	
	
	public void doWork(String outputName) throws IOException{
		
		Map<String, Collection<Gene>> newGenes = new HashMap<String, Collection<Gene>>();
		for(String chr:annotations.keySet()){
			//For all paths on this chromosome
			
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			/*
			 * REMOVE ANY GENES WITH 1-3BP OVERHANGS
			 */
			//Iterate over all assemblies
			Iterator<Gene> iter=tree.toCollection().iterator();			
			//Check for coverage and add only those that pass the coverage mark
			IntervalTree<Gene> filtered = new IntervalTree<Gene>();

			while(iter.hasNext()){
				Gene assembly=iter.next();
				boolean	toRemove = false;
				//For all overlapping assemblies
				Iterator<Gene> overlappers=tree.overlappingValueIterator(assembly.getStart(), assembly.getEnd());
				while(overlappers.hasNext()){
					Gene overlapper = overlappers.next();
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
									
									if(compatible(p, q)){
										Annotation overlap = exon1.intersect(intron2);
										if(overlap.getLengthOnReference()<=3){
											logger.warn(assembly.getName()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC()+" of "+overlapper.getName());
											toRemove = true;
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
				if(!toRemove){
					filtered.put(assembly.getStart(), assembly.getEnd(), assembly);
				}
			}
			
			
			/*
			 * MERGE THE GENES
			 */
			//Step 1: Iterate through the working assembly and find paths that overlap
			iter=filtered.toCollection().iterator();
			//if overlaps but is incompatible, see if we can branch
			//iterate
			Collection<Gene> considered = new HashSet<Gene>();
	 		logger.info("Enter merging");
			while(iter.hasNext()){
				Gene branch1=iter.next();
				considered.add(branch1);
				//get overlapping branches
				Iterator<Gene> overlappers=filtered.overlappingValueIterator(branch1.getStart(), branch1.getEnd());
				//Collection<Assembly> toRemove=new TreeSet<Assembly>();
				while(overlappers.hasNext()){
					Gene branch2=overlappers.next();
					if(!considered.contains(branch2)){
						if(!branch1.equals(branch2) && compatible(branch1, branch2)){
							logger.info("Merging: "+ branch1.getName()+ " and "+branch2.getName());
							Collection<Annotation> rtrn=new TreeSet<Annotation>();
							rtrn.addAll(branch1.getBlocks());
							rtrn.addAll(branch2.getBlocks());
							Gene merged=new Gene(rtrn);
							merged.setName(branch1.getName());
							//SET THE SCORE TO THE MAX FPKM
							merged.setBedScore(Math.max(branch1.getBedScore(),branch2.getBedScore()));
							//remove annotation1 and annotation2
							filtered.remove(branch1.getStart(), branch1.getEnd(), branch1);
							filtered.remove(branch2.getStart(), branch2.getEnd(), branch2);
							//add merged
							filtered.put(merged.getStart(), merged.getEnd(), merged);
						}
					}
				}
			}
			
			/*
			 * FPKM THRESHOLD
			 */
			//MAKE A MAP OF GENE TO ISOFORMS
			Map<Gene,Set<Gene>> isoformMap = BuildScriptureCoordinateSpace.getIsoformMap(filtered.toCollection());		
			Set<Gene> finalSet = new HashSet<Gene>();
			//For each gene
			for(Gene gene:isoformMap.keySet()){
				//For each transcript
				double max = Double.MIN_VALUE;
				String maxName = null;
				
				for(Gene isoform:isoformMap.get(gene)){					
					double score = isoform.getBedScore();
					if(score>max){
						max = score;
						maxName = isoform.getName();
					}
				}
				int counter=1;
				//CHeck if the FPKM of each is at least 10% of the mac FPKM
				for(Gene isoform:isoformMap.get(gene)){		
					double score = isoform.getBedScore();
					if(score>0.1*max && score>0.0){
						//filtered.remove(isoform.getStart(), isoform.getEnd(), isoform);
						isoform.setName(maxName+"."+counter);
						counter++;
						finalSet.add(isoform);
					}
					else{
						filtered.remove(isoform.getStart(), isoform.getEnd(), isoform);
						logger.info(isoform.getName()+" is removed because FPKM = "+score);
					}
				}
			}
			newGenes.put(chr, finalSet);		
		}
		
	
		write(outputName,newGenes);
	}
	
	private void write(String save, Map<String, Collection<Gene>> rtrn) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String chr: rtrn.keySet()){
			Collection<Gene> genes=rtrn.get(chr);
			for(Gene gene: genes){
				writer.write(gene.toBED()+"\n");
			}
		}
		
		writer.close();
	}
	
	private boolean compatible(Annotation assembly, Annotation read) {
		//Two alignments will be defined as compatible if:
		//(i) the intronic locations are exactly same
		//if both have introns, ensure that there are no non-overlapping introns
		//logger.info(read.getName());
		if(!assembly.getSpliceConnections().isEmpty() 
				|| !read.getSpliceConnections().isEmpty()){
			//check if introns are compatible
			if(areIntronsCompatible(assembly, read)){
				//logger.info("Introns are compatible");
				return true;
			}
			//logger.info("Introns are NOT compatible");
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
	
	private boolean overlap(Annotation assembly, Annotation read) {
		return assembly.overlaps(read,true);
	}
	
	private boolean areIntronsCompatible(Annotation assembly, Annotation read) {
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//introns are compatible if:
		//(i) all overlapping introns are identical
		for(Annotation intron1: assemblyIntrons){
			for(Annotation intron2: readIntrons){
				//if overlaps but not identical
				if(intron1.overlaps(intron2,true)){
					if(!intron1.equals(intron2, true)){
					//	logger.info("Case1");
						return false;
					} //TODO: This should use strand info or not based on flag
				}
			}
		}
		
		//(ii) if none of the exons overlap an intron in the other
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation intron2: readIntrons){
				if(exon1.overlaps(intron2,true)){
					//logger.info("Case2a");
					return false;
				}
			}
		}
		
		for(Annotation exon2: read.getBlocks()){
			for(Annotation intron1: assemblyIntrons){
				if(exon2.overlaps(intron1,true)){
					//logger.info("Case2b");
					return false;
				}
			}
		}
		
		//(ii) the introns dont overlap but the exons do
		//just need to test that any exons overlap
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation exon2: read.getBlocks()){
				if(exon1.overlaps(exon2,true)){
					return true;
				}
			}
		}
		//logger.info("Case3");
		return false;
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
		new ReconstructionPostProcessor(argMap.getMandatory("annotations"),argMap.getOutput());
		
	}
	
}

package nextgen.core.scripture.statistics;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

public class FilterDuplicateGenes {

	static Logger logger = Logger.getLogger(FilterDuplicateGenes.class.getName());
	Map<String,Collection<Gene>> annotations;
	String output;
	
	public FilterDuplicateGenes(String task,String inputFile,String outputFile) throws IOException{
		
		annotations= BEDFileParser.loadDataByChr(new File(inputFile));
		output = outputFile;
		
		if(task.equalsIgnoreCase("fullLength")){
			logger.info("Filtering for full length");
			filterFullLength();
		}
		else if(task.equalsIgnoreCase("3p")){
			logger.info("Filtering for 3p");
			filter3p();
		}
		else if(task.equalsIgnoreCase("5p")){
			logger.info("Filtering for 5p");
			filter5p();
		}
		else if(task.equalsIgnoreCase("overlapEnds")){
			logger.info("Removing all annotations where genes of opposite orientation overlap");
			removeOverlap();
		}
		else if(task.equalsIgnoreCase("singleIsoform")){
			logger.info("Removing all annotations with more than 1 isoform");
			filterSingleIsoforms();
		}
		
	}
	
	/**
	 * Returns on only those annotations without any overlapping isoforms
	 */
	private void filterSingleIsoforms() throws IOException{
		
		FileWriter writer=new FileWriter(output);
		for(String chr:annotations.keySet()){
			
			Collection<Gene> fl= new TreeSet<Gene>();
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
				fl.add(g);
			}
			
			Iterator<Gene> iter=tree.toCollection().iterator();
			
			Collection<Gene> considered = new HashSet<Gene>();
			
			while(iter.hasNext()){
				Gene gene=iter.next();
				
				boolean overlaps=false;
				Collection<Gene> listToWrite = new HashSet<Gene>();
				if(!considered.contains(gene)){
					considered.add(gene);
					//get overlapping branches
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getStart(), gene.getEnd());
					while(overlappers.hasNext()){
						Gene gene2=overlappers.next();
						//If gene not already considered
						if(!considered.contains(gene2)){
							//If both genes are not same and OPPOSITE orientation
							if(!gene.equals(gene2) && gene.getOrientation().equals(gene2.getOrientation()) && gene.overlaps(gene2, 0.3)){
								overlaps=true;
								listToWrite.add(gene2);
								considered.add(gene2);
							}
						}
					}
					
					//Now compare all overlapping genes
					for(Gene g:listToWrite){
						fl.remove(g);
					}										
				}
				if(overlaps){
					fl.remove(gene);
				}
			}
			for(Gene gene:fl){
				writer.write(gene.toBED()+"\n");
			}
		}
		writer.close();
	}
	/**
	 * filter out overlapping genes, keeping the more abundant one.
	 * @throws IOException 
	 */
	private void filterFullLength() throws IOException{
		
		FileWriter writer=new FileWriter(output);
		
		for(String chr:annotations.keySet()){
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			
			Iterator<Gene> iter=tree.toCollection().iterator();			

			Collection<Gene> considered = new HashSet<Gene>();
			
			while(iter.hasNext()){
				Gene gene=iter.next();
				Collection<Gene> listToWrite = new HashSet<Gene>();
				listToWrite.add(gene);
				if(!considered.contains(gene)){
					considered.add(gene);
					//get overlapping branches
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getStart(), gene.getEnd());
					while(overlappers.hasNext()){
						Gene gene2=overlappers.next();
						//If gene not already considered
						if(!considered.contains(gene2)){
							//If both genes are not same and same orientation
							if(!gene.equals(gene2) && gene.getOrientation().equals(gene2.getOrientation())){
								listToWrite.add(gene2);
								considered.add(gene2);
							}
						}
					}
					
					//Now compare all overlapping genes
					Gene best = gene;					
					for(Gene g:listToWrite){
						if(g.getBedScore()>best.getBedScore()){
							best=g;
						}
					}
					//Write the best gene
					writer.write(best.toBED()+"\n");					
				}
			}

		}
		writer.close();
	}
	
	/**
	 * filter out overlapping genes, keeping the more abundant one.
	 * @throws IOException 
	 */
	private void removeOverlap() throws IOException{
		
		FileWriter writer=new FileWriter(output);
		for(String chr:annotations.keySet()){
			
			Collection<Gene> fl= new TreeSet<Gene>();
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
				fl.add(g);
			}
			
			Iterator<Gene> iter=tree.toCollection().iterator();
			
			Collection<Gene> considered = new HashSet<Gene>();
			
			while(iter.hasNext()){
				Gene gene=iter.next();
				
				boolean overlaps=false;
				Collection<Gene> listToWrite = new HashSet<Gene>();
				if(!considered.contains(gene)){
					considered.add(gene);
					//get overlapping branches
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getStart(), gene.getEnd());
					while(overlappers.hasNext()){
						Gene gene2=overlappers.next();
						//If gene not already considered
						if(!considered.contains(gene2)){
							//If both genes are not same and OPPOSITE orientation
							if(!gene.equals(gene2) && !gene.getOrientation().equals(gene2.getOrientation())){
								overlaps=true;
								fl.remove(gene2);
								considered.add(gene2);
							}
						}
					}
					
					for(Gene g:listToWrite){
						fl.remove(g);
					}
				}
				if(overlaps){
					fl.remove(gene);
				}
			}
			for(Gene gene:fl){
				writer.write(gene.toBED()+"\n");
			}
		}
		writer.close();
	}
	
	/**
	 * filter out overlapping genes with same 3' end, keeping the more abundant one.
	 * @throws IOException 
	 */
	private void filter3p() throws IOException{
		
		FileWriter writer=new FileWriter(output);
		
		for(String chr:annotations.keySet()){
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			
			Iterator<Gene> iter=tree.toCollection().iterator();			

			Collection<Gene> considered = new HashSet<Gene>();
			
			while(iter.hasNext()){
				Gene gene=iter.next();
				Collection<Gene> listToWrite = new HashSet<Gene>();
				listToWrite.add(gene);
				if(!considered.contains(gene)){
					considered.add(gene);
					//get overlapping branches
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getStart(), gene.getEnd());
					while(overlappers.hasNext()){
						Gene gene2=overlappers.next();
						//If gene not already considered
						if(!considered.contains(gene2)){
							//If both genes are not same and same orientation and their 3' ends are the same
							if(!gene.equals(gene2) && gene.getOrientation().equals(gene2.getOrientation()) && gene.getOrientedEnd()==gene2.getOrientedEnd()){
								listToWrite.add(gene2);
								considered.add(gene2);
							}
						}
					}
					
					//Now compare all overlapping genes
					Gene best = gene;					
					for(Gene g:listToWrite){
						if(g.getBedScore()>best.getBedScore()){
							best=g;
						}
					}
					//Write the best gene
					writer.write(best.toBED()+"\n");
				}
			}

		}
		writer.close();
	}
	
	/**
	 * filter out overlapping genes, keeping the more abundant one.
	 * @throws IOException 
	 */
	private void filter5p() throws IOException{
		
		FileWriter writer=new FileWriter(output);
		
		for(String chr:annotations.keySet()){
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			
			Iterator<Gene> iter=tree.toCollection().iterator();			

			Collection<Gene> considered = new HashSet<Gene>();
			
			while(iter.hasNext()){
				Gene gene=iter.next();
				Collection<Gene> listToWrite = new HashSet<Gene>();
				listToWrite.add(gene);
				if(!considered.contains(gene)){
					considered.add(gene);
					//get overlapping branches
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getStart(), gene.getEnd());
					while(overlappers.hasNext()){
						Gene gene2=overlappers.next();
						//If gene not already considered
						if(!considered.contains(gene2)){
							//If both genes are not same and same orientation
							if(!gene.equals(gene2) && gene.getOrientation().equals(gene2.getOrientation()) && gene.getOrientedStart()==gene2.getOrientedStart()){
								listToWrite.add(gene2);
								considered.add(gene2);
							}
						}
					}
					
					//Now compare all overlapping genes
					Gene best = gene;					
					for(Gene g:listToWrite){
						if(g.getBedScore()>best.getBedScore()){
							best=g;
						}
					}
					//Write the best gene
					writer.write(best.toBED()+"\n");
				}
			}

		}
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		 
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"fullLength");
		
		new FilterDuplicateGenes(argMap.getTask(),argMap.getInput(),argMap.getOutput());
	}
	static final String usage = "Usage: FilterDuplicateGenes " +
			"\n\t-task fullLength "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\t\t-annotations <Annotations bed file. [BED by default]> "+
			"\n\t\t-out <Output file [Defaults to stdout]> ";		

}

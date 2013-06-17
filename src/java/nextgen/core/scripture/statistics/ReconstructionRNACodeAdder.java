package nextgen.core.scripture.statistics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;
import nextgen.core.scripture.BuildScriptureCoordinateSpace;

import broad.core.datastructures.IntervalTree;
import broad.pda.annotation.BEDFileParser;

public class ReconstructionRNACodeAdder {

	private Map<String, String> classFiles;
	private Map<String, Boolean> stranded;
	private Map<String, Boolean> overlap;
	Map<String,Collection<Gene>> reconstructions;
	private String outputName;
	private double THRESHOLD=25.0;
	
	static Logger logger = Logger.getLogger(ClassifyReconstructions.class.getName());
	
	public ReconstructionRNACodeAdder(String reconstructFile,String fileName,String outName) throws IOException{
		
		reconstructions = BEDFileParser.loadDataByChr(reconstructFile);
		
		outputName = outName;
		
		calculate(fileName);
	}
	
	/**
	 * Calculates and separates the reconstructions that map to each annotation class
	 * @throws IOException 
	 */
	private void calculate(String fileName) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".summary.txt"));
		bw.write("AnnotationClass : ReconstructionsThatOverlap\tRemainingReconstructions+\n");
		
		//Map<String,Collection<Gene>> newGenes = new HashMap<String,Collection<Gene>>();
		//Count initial reconstructions
		int recon1 = 0;
		for(String chr:reconstructions.keySet()){
			for(Gene gene:reconstructions.get(chr)){
				recon1++;
			}
		}
		logger.info("Starting with "+recon1+" reconstructions");
						
		int counter = 0;
		
		Map<String,Collection<Gene>> annotations = BEDFileParser.loadDataByChr(fileName); 
		Map<String,Collection<Gene>> overlappingGenes = new HashMap<String,Collection<Gene>>();

		String outFile = outputName+".overlapping.bed";
		
		// Find the overlapping reconstructions
		for(String chr:reconstructions.keySet()){
			
			overlappingGenes.put(chr, new TreeSet<Gene>());
			//newGenes.put(chr, new TreeSet<Gene>());
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			
	/*		if(!annotations.containsKey(chr)){
				for(Gene gene:reconstructions.get(chr)){
					double[] extras = new double[1];
					extras[0] = 0.0;
					gene.setExtraFields(extras);
					//newGenes.get(chr).add(gene);
				}
			}
			else{*/
				//MAKE AN INTERVAL TREE OF THE GENES on this chr in the annotations
				for(Gene g:annotations.get(chr)){
					tree.put(g.getStart(), g.getEnd(), g);
				}
			
				for(Gene gene:reconstructions.get(chr)){
					
					boolean overlaps=false;
					double score = 0.0;
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getStart(), gene.getEnd());
					
					//Collection<Assembly> toRemove=new TreeSet<Assembly>();
					while(overlappers.hasNext()){
						//For each overlapping gene
						Gene other=overlappers.next();
						
						score+=(other.getBedScore());
	
					}
					logger.info(gene.getName()+" score: "+score);
					if(score>THRESHOLD){
						overlaps=true;
					}
					
					if(overlaps){
						overlappingGenes.get(chr).add(gene);
						counter++;
					}
					//now set extra fields
					/*double[] extras = new double[1];
					extras[0] = score;
					gene.setExtraFields(extras);
					newGenes.get(chr).add(gene);*/	
				}
				
				logger.info(chr);
				// Get the remaining reconstructions
				for(Gene g:overlappingGenes.get(chr)){
					reconstructions.get(chr).remove(g);
				}
				
			}
		//Count remainging reconstructions
		int recon = 0;
		for(String chr:reconstructions.keySet()){
			for(Gene gene:reconstructions.get(chr)){
				recon++;
			}
		}
		// Output to a separate file - write out the number
		BuildScriptureCoordinateSpace.write(outFile, overlappingGenes);			
		bw.write(outputName+" : "+counter+"\t"+recon+"\n");
		
		bw.close();
//		}
		
		//The final reconstructions could be lincRNAs
		BuildScriptureCoordinateSpace.write(outputName, reconstructions);
//		BuildScriptureCoordinateSpace.write(outputName, newGenes);	
	}
	
	public static void main(String[] args) throws IOException{
		
		if(args.length<3)
			System.err.println(usage);
		else
			new ReconstructionRNACodeAdder(args[0],args[1],args[2]);
	}
	
	static String usage=" args[0]=reconstructions bed file \n\t args[1]=Bed file with RNACode for genome in bed format \n\t args[2]: output file name";

	public static String whitespaceDelimiter = "\\s++";
}


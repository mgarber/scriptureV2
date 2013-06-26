
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

public class RemoveIntronicReconstructions {

	private Map<String, String> classFiles;
	private Map<String, Boolean> stranded;
	private Map<String, Boolean> overlap;
	Map<String,Collection<Gene>> reconstructions;
	private String outputName;
	
	static Logger logger = Logger.getLogger(ClassifyReconstructions.class.getName());
	
	public RemoveIntronicReconstructions(String reconstructFile,String fileName,String outName) throws IOException{
		
		reconstructions = BEDFileParser.loadDataByChr(reconstructFile);
		
		outputName = outName;
		getClassFileNames(fileName);
		
		calculate();
	}
	
	/**
	 * Calculates and separates the reconstructions that map to each annotation class
	 * @throws IOException 
	 */
	private void calculate() throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".summary.txt"));
		bw.write("AnnotationClass : ReconstructionsThatOverlap\tRemainingReconstructions+\n");
		//Count initial reconstructions
		int recon1 = 0;
		for(String chr:reconstructions.keySet()){
			for(Gene gene:reconstructions.get(chr)){
				recon1++;
			}
		}
		logger.info("Starting with "+recon1+"reconstructions");
		
		
		for(String annClass:classFiles.keySet()){
			String outFile = outputName+"."+annClass+".overlapping.bed";
			
			Map<String,Collection<Gene>> annotations = BEDFileParser.loadDataByChr(classFiles.get(annClass)); 
			
			Map<String,Collection<Gene>> overlappingGenes = new HashMap<String,Collection<Gene>>();
			
			logger.info("For "+annClass+" checking overlap "+overlap.get(annClass)+" and stranded: "+stranded.get(annClass));
			int counter = 0;
			// Find the overlapping reconstructions
			for(String chr:reconstructions.keySet()){
				//MAKE AN INTERVAL TREE OF THE GENES on this chr in the annotations
				
				overlappingGenes.put(chr, new TreeSet<Gene>());
				IntervalTree<Gene> tree = new IntervalTree<Gene>();
				
				if(!annotations.containsKey(chr))
					continue;
				for(Gene g:annotations.get(chr)){
					tree.put(g.getStart(), g.getEnd(), g);
				}
			
				for(Gene gene:reconstructions.get(chr)){
					
					boolean overlaps = false;
					
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getStart(), gene.getEnd());
					
					//Collection<Assembly> toRemove=new TreeSet<Assembly>();
					while(overlappers.hasNext()){
						//For each overlapping gene
						Gene other=overlappers.next();
						
						Gene intron = other.getIntrons();
						
						if(intron!=null){
							if(intron.contains(gene) && gene.getOrientation().equals(other.getOrientation())){
								overlaps = true;
							}
						}
					}
					if(overlaps){
						overlappingGenes.get(chr).add(gene);
						counter++;
					}
				}
				
				// Get the remaining reconstructions
				for(Gene gene:overlappingGenes.get(chr)){
					reconstructions.get(chr).remove(gene);
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
			bw.write(annClass+" : "+counter+"\t"+recon+"\n");
			
		}
		
		//The final reconstructions could be lincRNAs
		BuildScriptureCoordinateSpace.write(outputName, reconstructions);
		
		bw.write("\n\n");
		//Write the lincs on each chromosomes
		for(String chr:reconstructions.keySet()){
			bw.write(chr+"\t"+reconstructions.get(chr).size()+"\n");
		}
		bw.close();
	}
	
	
	/**
	 * Reads a 2-column file with name of the annotation class mapped to its path
	 * @param fileName
	 * @throws IOException
	 */
	private void getClassFileNames(String fileName) throws IOException{
		
		classFiles = new HashMap<String,String>();
		stranded = new HashMap<String,Boolean>();
		overlap = new HashMap<String,Boolean>();
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String s;
		while((s = br.readLine())!= null){
			boolean value = false;
			classFiles.put(s.split(whitespaceDelimiter)[0],s.split(whitespaceDelimiter)[1]);	
			
			logger.info(s.split(whitespaceDelimiter)[0]+":");
			
			if(s.split(whitespaceDelimiter)[2].equals("stranded")){
				value = true;
				logger.info("is stranded " );
			}	
			else if(s.split(whitespaceDelimiter)[2].equals("unstranded")){
				logger.info("is UNstranded " );
				value = false;
			}
				else
					logger.info("Problem in value of strand in "+s);
			stranded.put(s.split(whitespaceDelimiter)[0],value);
			
			if(s.split(whitespaceDelimiter)[3].equals("overlap")){
				logger.info("is OVERLAP " );
				value = true;
			}	
			else
				if(s.split(whitespaceDelimiter)[3].equals("contain")){
					logger.info("is CONTAIN " );
					value = false;
				}
				else
					logger.info("Problem in value of overlap in "+s);
			overlap.put(s.split(whitespaceDelimiter)[0],value);
		}		
	}
	
	public static void main(String[] args) throws IOException{
		
		if(args.length<3)
			System.err.println(usage);
		else
			new RemoveIntronicReconstructions(args[0],args[1],args[2]);
	}
	
	static String usage=" args[0]=reconstructions bed file \n\t args[1]=list of bed files to filter \n\t args[2]: output file name";

	public static String whitespaceDelimiter = "\\s++";
}


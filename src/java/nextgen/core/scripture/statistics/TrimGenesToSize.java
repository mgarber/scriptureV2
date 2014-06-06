package nextgen.core.scripture.statistics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;

import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

public class TrimGenesToSize {

	static final String usage = "Usage: TrimGenesToSize -task trim3p/trim5p "+
			"\n\t\t-in <Bed file with annotations> "+
			"\n\t\t-length <Size to which gene end must be trimmed>"+
			"\n\t\t-out <Output file [Defaults to stdout]> "
			;

	Map<String,Collection<Gene>> annotations;
	static Logger logger = Logger.getLogger(TrimGenesToSize.class.getName());
	int length=0;
	
	public TrimGenesToSize(File input,String outputFile,int o,String task) throws IOException{
		
		annotations= BEDFileParser.loadDataByChr(input);
		length = o;
		if(task.equalsIgnoreCase("trim3p"))
			trimGenes(outputFile);
		else
			trim5pGenes(outputFile);
	}
	
	private void trimGenes(String outputFile) throws IOException{
		
		logger.info("Trimming genes to length "+length);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		for(String chr:annotations.keySet()){
			logger.info("Processing "+chr);
			
			for(Gene gene:annotations.get(chr)){
				//if gene is shorter than trim length
				if(gene.length()<length){
					//Write as it is
					bw.write(gene.toBED()+"\n");
				}
				else{
					if(gene.isNegativeStrand()){
						Gene gs = gene.trimGene(0, length);
						gs.setName(gene.getName());
						bw.write(gs.toBED()+"\n");
					}
					else{
						//Trim gene is strand insensitive
						Gene gs = gene.trimGene(gene.length()-length, gene.length());
						gs.setName(gene.getName());
						bw.write(gs.toBED()+"\n");
					}
				}
			}
		}
		bw.close();
	}
	
	private void trim5pGenes(String outputFile) throws IOException{
		
		logger.info("Trimming genes to 5p length "+length);
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		for(String chr:annotations.keySet()){
			logger.info("Processing "+chr);
			
			for(Gene gene:annotations.get(chr)){
				//if gene is shorter than trim length
				if(gene.length()<length){
					//Write as it is
					bw.write(gene.toBED()+"\n");
				}
				else{
					if(!gene.isNegativeStrand()){
						Gene gs = gene.trimGene(0, length);
						gs.setName(gene.getName());
						bw.write(gs.toBED()+"\n");
					}
					else{
						//trimGene is strand insensitive
						Gene gs = gene.trimGene(gene.length()-length, gene.length());
						gs.setName(gene.getName());
						bw.write(gs.toBED()+"\n");
					}
				}
			}
		}
		bw.close();
	}
	
public static void main (String [] args) throws ParseException, IOException {
		
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"trim");

		//Test trim genes
/*		Gene g = new Gene("chr1",100,5000,"gene",Strand.POSITIVE);
		System.out.println(g.toBED());
		Gene gs = g.trimGene(g.length()-1000, g.length());
		System.out.println(gs.toBED());
		
		g = new Gene("chr1",100,5000,"gene",Strand.NEGATIVE);
		System.out.println(g.toBED());
		gs = g.trimGene(0, 1000);
		System.out.println(gs.toBED());
*/		
		new TrimGenesToSize(new File(argMap.getInput()),argMap.getOutput(),argMap.getInteger("length"),argMap.getTask());
		
		
	}
}

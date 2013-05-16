package nextgen.core.scripture.statistics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.scripture.AddEndRNASeqToScripture;
import nextgen.core.scripture.AddEndRNASeqToScripture.IsoformMap;

public class OldReconstructionCompare {
	
	static final String usage = "Usage: OldReconstructionCompare -task <task name> "+
			"\n\tcompare: Compares the reference with another annotation file and outputs a table with results" +
			"\n\t\t-old <Bed file with old reconstructions> "+
			"\n\t\t-new <Bed file with new reconstructions> "+
			"\n\t\t-minOverlap MinimumOverlap required for comparison. [Default is 50%] "+
			"\n\t\t-out <Output file [Defaults to stdout]> "
			;


	Map<String,Collection<Gene>> oldAnnotations;
	Map<String,Collection<Gene>> newAnnotations;
	static Logger logger = Logger.getLogger(OldReconstructionCompare.class.getName());
	double minOverlap;
	private static double DEFAULT_OVERLAP=0.5;

	public OldReconstructionCompare(File oldR,File newR,String outputFile,double o) throws IOException{
		
		oldAnnotations= BEDFileParser.loadDataByChr(oldR);
		newAnnotations= BEDFileParser.loadDataByChr(newR);
		minOverlap = o;
		calculateGlobalStatsOnIsoformNumber(outputFile+".isoforms.stats");
	}
	
	/**
	 * Calculates the global statistics on differences in isoform number between new and old reconstructions
	 * @param outputName
	 * @throws IOException 
	 */
	private void calculateGlobalStatsOnIsoformNumber(String outputName) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName));
		BufferedWriter bwSum = new BufferedWriter(new FileWriter(outputName+".summary"));
		bw.write("LineName\t#IsoformsInV2\t#IsoformsInNew\n");
		//MAKE ISOFORM MAPS FOR BOTH ANNOTATION SETS
		IsoformMap oldMap = AddEndRNASeqToScripture.buildIsoformMap(oldAnnotations);
		IsoformMap newMap = AddEndRNASeqToScripture.buildIsoformMap(newAnnotations);
		
		int lincsWithIncreasedNumIsoforms= 0;
		int lincsWithDecreasedNumIsoforms=0;
		int lincsWithSameNumIsoforms=0;
		//For each gene in the old reconstructions
		for(Gene gene:oldMap.getAllGenes()){
			//Find the gene in new reconstructions that it overlaps
			for(Gene gene2:newMap.getGenesForChromosome(gene.getChr())){
				if(gene.overlapsStranded(gene2,minOverlap)){
					int oldNum = oldMap.getNumOfIsoformsForGene(gene);
					int newNum = newMap.getNumOfIsoformsForGene(gene2);
					bw.write(gene.getName()+"\t"+oldNum+"\t"+newNum+"\n"); 
					
					if(oldNum>newNum){
						lincsWithDecreasedNumIsoforms++;
					}
					else{
						if(oldNum<newNum){
							lincsWithIncreasedNumIsoforms++;
						}
						else{
							lincsWithSameNumIsoforms++;
						}
					}
				}
			}
		}
		//Summarize
		bwSum.write("LincRNAs with an increase in the number of Isoforms reported: "+lincsWithIncreasedNumIsoforms+"\n");
		bwSum.write("LincRNAs with an decrease in the number of Isoforms reported: "+lincsWithDecreasedNumIsoforms+"\n");
		bwSum.write("LincRNAs with same number of Isoforms reported: "+lincsWithSameNumIsoforms+"\n");
		
		bw.close();
		bwSum.close();
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
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"compare");

		new OldReconstructionCompare(new File(argMap.getMandatory("old")),new File(argMap.getMandatory("new")),argMap.getOutput(),argMap.getDouble("minOverlap", DEFAULT_OVERLAP));
		
	}

}

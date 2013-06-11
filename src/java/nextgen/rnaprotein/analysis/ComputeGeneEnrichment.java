package nextgen.rnaprotein.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.RatioScore;

public class ComputeGeneEnrichment {

	boolean useFragments=true;
	static Logger logger = Logger.getLogger(ComputeGeneEnrichment.class.getName());
	
	public ComputeGeneEnrichment(String bam1, String bam2, Map<String, Collection<Gene>> genesByChr, String save) throws IOException{
		CoordinateSpace cs=new TranscriptomeSpace(genesByChr);
		
		AlignmentModel model1=new AlignmentModel(bam1, cs, useFragments);
		AlignmentModel model2=new AlignmentModel(bam2, cs, useFragments);
		
		logger.info("Made alignment models");
		
		Map<Gene, RatioScore> scores=scoreGenes(genesByChr, model1, model2);
		
		write(save, scores);
	}
		
	private void write(String save, Map<Gene, RatioScore> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: scores.keySet()){
			RatioScore score=scores.get(gene);
			writer.write(gene.getChr()+"\t"+gene.getStart()+"\t"+gene.getEnd()+"\t"+score.getRatio()+"\n");
		}
		
		writer.close();
	}

	private Map<Gene, RatioScore> scoreGenes(Map<String, Collection<Gene>> genesByChr, AlignmentModel model1, AlignmentModel model2) {
		Map<Gene, RatioScore> rtrn=new TreeMap<Gene, RatioScore>();
		RatioScore.Processor processor=new RatioScore.Processor(model1, model2);
		for(String chr: genesByChr.keySet()){
			logger.info(chr);
			Collection<Gene> genes=genesByChr.get(chr);
			for(Gene gene: genes){
				RatioScore score=model1.getScore(gene, processor);
				rtrn.put(gene, score);
			}
		}
		return rtrn;
	}



	public static void main(String[] args) throws IOException{
		if(args.length>3){
			String bam1=args[0];
			String bam2=args[1];
			Map<String, Collection<Gene>> genesByChr=BEDFileParser.loadDataByChr(args[2]);
			String save=args[3];
			new ComputeGeneEnrichment(bam1, bam2, genesByChr, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam1 \n args[1]=bam2 \n args[2]=genes \n args[3]=save";
	
}

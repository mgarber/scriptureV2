package nextgen.core.normalize;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.rnaprotein.analysis.ComputeGeneEnrichment;

/**
 * 
 * @author mguttman
 * The purpose of this class is to take a BAM file and normalize the counts to show enrichment over the transcript average
 * Will write a big wig
 */
public class NormalizeCountsToTranscriptAverage {

	static Logger logger = Logger.getLogger(NormalizeCountsToTranscriptAverage.class.getName());
	private boolean useFragments=true;
	
	public NormalizeCountsToTranscriptAverage(String bamFile, Map<String, Collection<Gene>> genes, String sizes, String save) throws IOException, InterruptedException{
		CoordinateSpace cs=new TranscriptomeSpace(genes);
		AlignmentModel model=new AlignmentModel(bamFile, cs, useFragments);
		model.addFilter(new GenomicSpanFilter(50000));
		
		FileWriter writer=new FileWriter(save);
		
		//normalize genes
		normalizeGenes(model, genes.get("chr19"), writer);
		
		//close file
		writer.close();
		
		//convert to bigwig
		convertToBW(save, sizes);
	}
	
	private void convertToBW(String save, String sizes) throws IOException, InterruptedException {
		String cmd="/seq/lincRNA/scripts/UCSC/wigToBigWig save "+sizes+" "+save+".bw";
		Runtime run=Runtime.getRuntime();
		Process p=run.exec(cmd);
		p.waitFor();
	}

	private void normalizeGenes(AlignmentModel model, Map<String, Collection<Gene>> genes, FileWriter writer) throws IOException {
		for(String chr: genes.keySet()){
			logger.info(chr);
			normalizeGenes(model, genes.get(chr), writer);
		}
	}
	
	

	private void normalizeGenes(AlignmentModel model, Collection<Gene> genes, FileWriter writer) throws IOException {
		for(Gene gene: genes){
			
			WindowScoreIterator<CountScore> scoreIter = model.scan(gene, 1, 0);
			double average=model.getCount(gene)/gene.size();
			//logger.info(gene.toUCSC()+"\t"+average);
			
			writeNorm(scoreIter, average, writer);
		}
		
	}

	private void writeNorm(WindowScoreIterator<CountScore> iter, double average, FileWriter writer) throws IOException {
		double actualAverage=0;
		int num=0;
		ArrayList<CountScore> list=new ArrayList<CountScore>();
		while(iter.hasNext()){
			CountScore score=iter.next();
			list.add(score);
			actualAverage+=score.getCount();
			num++;
		}
		actualAverage=actualAverage/(double)num;
		logger.info(average+" "+actualAverage);
		
		for(CountScore score: list){
			double norm=score.getCount()/actualAverage;
			writer.write(score.getAnnotation().getChr()+"\t"+score.getAnnotation().getStart()+"\t"+score.getAnnotation().getEnd()+"\t"+norm+"\n");
		}
		
		/*
		//fixedStep  chrom=chrN  start=position  step=stepInterval
		writer.write("fixedStep chrom="+gene.getChr()+"\tstart="+gene.getStart()+"\tstep=1\n");
		
		for(int i=0; i<countsPerPosition.length; i++){
			if(countsPerPosition[i]>0){
				logger.info(gene.toUCSC()+" "+countsPerPosition[i]);
			}
			
			double normScore=countsPerPosition[i]/average;
			writer.write(normScore+"\n");
		}*/
	}

	public static void main(String[] args) throws IOException, InterruptedException{
		if (args.length>2){
			String bam=args[0];
			Map<String, Collection<Gene>> genes=BEDFileParser.loadDataByChr(args[1]);
			String sizes=args[2];
			String save=args[3];
			new NormalizeCountsToTranscriptAverage(bam, genes, sizes, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=bam file \n args[1]=genes \n args[2]=sizes \n args[3]=save";
}

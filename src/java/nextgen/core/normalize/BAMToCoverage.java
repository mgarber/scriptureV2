package nextgen.core.normalize;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GeneScore;

import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.Window;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;

//Convert a BAM File to a coverage profile
public class BAMToCoverage {

	static Logger logger = Logger.getLogger(BAMToCoverage.class.getName());
	int windowSize;
	
	public BAMToCoverage(ScanStatisticDataAlignmentModel data, int windowSize, String sizes, String save, String chr) throws IOException, InterruptedException{
		this.windowSize=windowSize;
		
		logger.info("Started the window slider");
		//For each window compute median at center point
		windowSlider(data, save+".wig", chr);
				
		logger.info("Converting to BigWig");
		//make bigWig
		makeBigWig(save, sizes);
	}
	
	//TODO This is not memory efficient, seems that the jump to new gene is screwing up
	private void windowSlider(ScanStatisticDataAlignmentModel data, String save, String chr) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Runtime run=Runtime.getRuntime();
		
		logger.info("Created the iterator");
		//Iterator<Window> iter=data.getCoordinateSpace().getWindowIterator(windowSize);
		Iterator<ScanStatisticScore> iter=data.scan(chr, this.windowSize, this.windowSize-1);
		
		logger.info("starting the iterations and scoring");
		int counter=0;
		while(iter.hasNext()){
			ScanStatisticScore score=iter.next();
			writer.write(score.getAnnotation().getChr()+"\t"+score.getAnnotation().getStart()+"\t"+score.getAnnotation().getEnd()+"\t"+score.getCount()/windowSize+"\n");
			//if(counter%1000 ==0){logger.info("At iteration "+counter+" "+score.getAlignment().toUCSC()+" "+score.getNumberOfReads());}
			counter++;
		}
		
		logger.info("done writing");
		writer.close();
	}

	private void makeBigWig(String save, String sizes) throws IOException, InterruptedException {
		String cmd="/seq/lincRNA/scripts/UCSC/bedGraphToBigWig "+save+".wig "+sizes+" "+save+".bw";
		Runtime run=Runtime.getRuntime();
		Process p=run.exec(cmd);
		p.waitFor();
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>5){
			ScanStatisticDataAlignmentModel data=new ScanStatisticDataAlignmentModel(args[0], new GenomicSpace(args[3]));
			int windowSize=new Integer(args[2]);
			String chrSizes=args[3];
			String save=args[4];
			String chr=args[5];
			new BAMToCoverage(data, windowSize, chrSizes, save, chr);
		}
		else{System.err.println(usage);}
	}
	
	static String usage= " args[0]=BAM file \n args[1]=genes \n args[2]=window size \n args[3]=chromosome sizes file \n args[4]=save \n args[5]=chr";
}

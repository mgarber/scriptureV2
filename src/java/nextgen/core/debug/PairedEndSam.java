package nextgen.core.debug;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.ReadFilter;

import nextgen.core.annotation.*;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.feature.Window;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.readFilters.ProperPairFilter;

public class PairedEndSam {

	
	public PairedEndSam(String bam, Alignments region, CoordinateSpace space) throws IOException {
		ScanStatisticDataAlignmentModel data=new ScanStatisticDataAlignmentModel(bam, space);
		data.addFilter(new ProperPairFilter());
		
		Window w=new GeneWindow(region);
		double count=data.getCount(w);
		System.err.println(region.toUCSC()+" "+count);
	}

	public PairedEndSam(Map<String, Collection<Gene>> genes, Alignments region) throws IOException {
		FileWriter writer=new FileWriter("w.bed");
		
		TranscriptomeSpace space=new TranscriptomeSpace(genes);
		
		Iterator<Window> iter=space.getWindowIterator(30, region.getChr(), region.getStart(), region.getEnd(), 29);
		
		//Iterator<? extends Window> iter=space.getOverlappingRegion(region.getChr(), region.getStart(), region.getEnd()).iterator();
		
		while(iter.hasNext()){
			Window w=iter.next();
			writer.write(w+"\n");
		}
		writer.close();
	}

	public static void main(String[] args)throws IOException{
		if(args.length>1){
			//String bam=args[0];
			Alignments region=new Alignments(args[1]);
			//CoordinateSpace cs=new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(args[2])));
			new PairedEndSam(BEDFileParser.loadDataByChr(new File(args[0])), region);			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes \n args[1]=ucsc string";
	
}

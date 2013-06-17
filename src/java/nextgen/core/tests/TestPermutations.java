package nextgen.core.tests;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;

import broad.pda.annotation.BEDFileParser;

public class TestPermutations {

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			String bamFile=(args[0]);
			Map<String, Collection<Gene>> genes=BEDFileParser.loadDataByChr(new File(args[1]));
			String save=args[2];
			CoordinateSpace space=new TranscriptomeSpace(genes);
			AlignmentModel model=new AlignmentModel(bamFile, space);
			
			Gene gene=genes.get(genes.keySet().iterator().next()).iterator().next();
			
			CloseableIterator<Alignment> shuffle=model.getPermutedAnnotations(gene);
			FileWriter writer=new FileWriter(save);
			while(shuffle.hasNext()){
				Alignment next=shuffle.next();
				writer.write(next+"\n");
			}
			writer.close();
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=genes \n args[2]=save";
}

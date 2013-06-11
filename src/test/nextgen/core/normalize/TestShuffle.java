package nextgen.core.normalize;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;

public class TestShuffle {

	public TestShuffle(AlignmentModel model, Gene gene, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		CloseableIterator<Alignment> iter=model.getPermutedAnnotations(gene);
		while(iter.hasNext()){
			writer.write(iter.next()+"\n");
		}
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Map<String, Collection<Gene>> genes=BEDFileParser.loadDataByChr(args[1]);
			CoordinateSpace cs=new TranscriptomeSpace(genes);
			AlignmentModel model=new AlignmentModel(args[0], cs);
			Gene g=BEDFileParser.getGene(args[1], "NM_013633");
			new TestShuffle(model, g, args[2]);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam \n args[1]=genes \n args[2]=save";
}

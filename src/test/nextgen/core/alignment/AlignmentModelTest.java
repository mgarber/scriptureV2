package nextgen.core.alignment;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.model.AlignmentModel;

import org.apache.log4j.Logger;

public class AlignmentModelTest {
	
	private static Logger logger = Logger.getLogger(AlignmentModelTest.class.getName());
	
	public static void main(String[] args) throws IOException {
		
		//String bamFile = "ProtectSeq_March_Rnase1K_merged.bam";
		//String bedFile = "RefSeq.nonrandom.bed";
		String bamFile = args[0];
		String bedFile = args[1];
		//String chromosome = args[2];
		Map<String, Collection<Gene>> genesByChromosome = BEDFileParser.loadDataByChr(new File (bedFile));
		TranscriptomeSpace transcriptomeSpace = new TranscriptomeSpace(genesByChromosome);
		AlignmentModel data = new AlignmentModel(bamFile, transcriptomeSpace);
		for(String name: genesByChromosome.keySet()) {
		//long startTime=System.currentTimeMillis();	
			for (Gene gene: genesByChromosome.get(name)) {
				data.scan(gene, 10, 1);
				
			//gene.getWindows(10);
			} 
		//long endTime=System.currentTimeMillis();
		//logger.info((endTime-startTime));
		logger.info("");
		}
		
		
		/*Gene gene=genesByChromosome.values().iterator().next().iterator().next();
		Collection<GeneWindow> windows=gene.getWindows(10);
		FileWriter writer=new FileWriter("windowTest.bed");
		for(GeneWindow window: windows) {
			writer.write(window+"\n");
		}
		writer.close();
		*/
	}
}
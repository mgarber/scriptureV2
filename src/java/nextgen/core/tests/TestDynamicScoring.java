package nextgen.core.tests;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.readFilters.SplicedReadFilter;
import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class TestDynamicScoring {

	public TestDynamicScoring(){}
	
	
	public static void main(String[] args) throws IOException{
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Bed gene annotation", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-chr", "Chromosome name", true);
		p.addIntegerArg("-s", "Start position", true);
		p.addIntegerArg("-e", "End position", true);
		p.addBooleanArg("-pe", "Use fragments", true);
		p.addStringArg("-save", "save file", true);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-c");
		String chrName = p.getStringArg("-chr");
		int start = p.getIntArg("-s");
		int end = p.getIntArg("-e");
		boolean useFragments = p.getBooleanArg("-pe");
		FileWriter writer = new FileWriter(p.getStringArg("-save"));
		
		TranscriptomeSpace t = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(bedFile)));
		GenomicSpace g = new GenomicSpace(chrSizeFile);
		
		AlignmentModel transcriptomeAlignments = new AlignmentModel(bamFile, t, useFragments);
		AlignmentModel genomeAlignments = new AlignmentModel(bamFile, g, useFragments);
		
		WindowScoreIterator<CountScore> scores=genomeAlignments.scan(new Alignments(chrName, start,end), 90, 89);
		
		
		genomeAlignments.addFilter(new SplicedReadFilter());
		System.err.println(genomeAlignments.getCount(new Alignments(chrName, start,end)));
		
		/*while(scores.hasNext()){
			CountScore score=scores.next();
			Annotation region=score.getAnnotation();
			//region.setScore(score.getCount());
			writer.write(region.getChr()+"\t"+region.getStart()+"\t"+region.getEnd()+"\t"+score.getCount()+"\n");
		}*/
		
		
		
		writer.close();
		
	}
	
}

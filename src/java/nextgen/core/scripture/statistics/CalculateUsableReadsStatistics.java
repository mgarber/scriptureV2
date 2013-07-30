package nextgen.core.scripture.statistics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.alignment.SingleEndAlignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.PairedAndProperFilter;
import nextgen.core.readFilters.ProperPairFilter;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.pda.annotation.BEDFileParser;

public class CalculateUsableReadsStatistics {
	
	Map<String,Collection<Gene>> annotations;
	private AlignmentModel model;
	private String outputFileName;
	
	static Logger logger = Logger.getLogger(CalculateUsableReadsStatistics.class.getName());

	public CalculateUsableReadsStatistics(String bedFile, String bamFile,String outputName,TranscriptionRead strand) throws IOException {
		 outputFileName = outputName;
		 annotations = BEDFileParser.loadDataByChr(bedFile);
		 model=new AlignmentModel(new File(bamFile).getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		 model.addFilter(new GenomicSpanFilter(20000000));
		 calculate();
	}
	
	private void calculate() throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
		bw.write("Total_Proper_Fragments\tTotal_PairedReads\tTotalSingleReads\t" +
				"Usable_PairedReads\tUsable_SingleReads\tTotal_UsableFragments" +
				"\tBothContained_PairedReads\tTotal_SingleMappedFragments\tTotal_MultiMappedFragments" +
				"\tUsable_SingleMappedFragments\tUsable_MultiMappedFragments\n");

		long total_paired_reads=0;
		long total_single_reads=0;
		long total_single_mapped_fragments=0;
		long total_multi_mapped_fragments=0;
		
		long usable_paired_reads=0;
		long usable_single_reads=0;
		long contained_paired_end=0;
		long usable_single_mapped_fragments=0;
		long usable_multimapped_fragments=0;
		
		//For each chromosome
		for(String chr:annotations.keySet()){
			logger.info("Processing "+chr);
			//MAKE AN INTERVAL TREE OF THE GENES on this chr in the annotations
			if(!model.containsReference(chr))
				continue;
			//For each annotation
			for(Gene gene:annotations.get(chr)){
				//Get all reads overlapping the transcript
				CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(gene,false), new ProperPairFilter());
				while(iter.hasNext()){	
					Alignment read = iter.next();
					//USABLE SINGLE OR PAIRED FRAGMENTS
					if(SingleEndAlignment.class.isInstance(read)){
						usable_single_reads++;
					}
					else{
						usable_paired_reads++;
						boolean contained=true;
						boolean singlyMapped=true;
						for(Annotation mate:read.getReadAlignments(model.getCoordinateSpace())){
							if(!gene.contains(mate)){
								contained=false;
							}
						}
						
						if(contained){
							contained_paired_end++;
						}
					}
					//SINGLE AND MULTI MAPPERS
					if(read.getWeight()==1.0){
						usable_single_mapped_fragments++;
					}
					else{
						usable_multimapped_fragments++;
					}
					
				}
				iter.close();
			}
			
			//NOW TOTAL
			CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
			while(iter.hasNext()){
				Alignment read=iter.next();
				if(SingleEndAlignment.class.isInstance(read)){
					total_single_reads++;
				}
				else{
					total_paired_reads++;
				}
				//SINGLE AND MULTI MAPPERS
				if(read.getWeight()==1.0){
					total_single_mapped_fragments++;
				}
				else{
					total_multi_mapped_fragments++;
				}
			}
			iter.close();
			
			logger.info((total_single_reads+total_paired_reads)+"\t"+total_paired_reads+"\t"+total_single_reads+"\t" +
					usable_paired_reads+"\t"+usable_single_reads+"\t"+(usable_paired_reads+usable_single_reads)+"\t"+
					contained_paired_end+"\t"+total_single_mapped_fragments+"\t"+total_multi_mapped_fragments+"\t"+
					usable_single_mapped_fragments+"\t"+usable_multimapped_fragments);
			
		}
		bw.write((total_single_reads+total_paired_reads)+"\t"+total_paired_reads+"\t"+total_single_reads+"\t" +
				usable_paired_reads+"\t"+usable_single_reads+"\t"+(usable_paired_reads+usable_single_reads)+"\t"+
				contained_paired_end+"\t"+total_single_mapped_fragments+"\t"+total_multi_mapped_fragments+"\t"+
				usable_single_mapped_fragments+"\t"+usable_multimapped_fragments+"\n");
		bw.close();
	}

	public static void main(String[] args) throws IOException{
		
		if(args.length<4)
			System.err.println(usage);
		else{
			TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
			if(args[4].equalsIgnoreCase("first")){
				//System.out.println("First read");
				strand = TranscriptionRead.FIRST_OF_PAIR;
			}
			else if(args[4].equalsIgnoreCase("second")){
				//System.out.println("Second read");
				strand = TranscriptionRead.SECOND_OF_PAIR;
			}
			else
				System.out.println("no strand");
			new CalculateUsableReadsStatistics(args[0],args[1],args[2],strand);
		}
	}
	
	static String usage=" args[0]=annotation bed file \n\t args[1]=bam file \n\t args[2]: output file name \n\t args[3] = strand";

}

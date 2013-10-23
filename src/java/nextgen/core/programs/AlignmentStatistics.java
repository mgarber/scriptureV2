package nextgen.core.programs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections15.Predicate;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.MatrixWithHeaders;
import broad.pda.annotation.BEDFileParser;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;

import nextgen.core.alignment.AbstractPairedEndAlignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.SingleEndAlignment;
import nextgen.core.annotation.Gene;
import nextgen.core.model.AlignmentModel;

public class AlignmentStatistics extends CommandLineProgram {
    private static final Log log = Log.getInstance(AlignmentStatistics.class);
    
    @Option(doc="Input BAM file (default format)", shortName="I")
    String INPUT;
    
    @Option(doc="Output File name", shortName="O")
    String OUTPUT;
    
    @Option(doc="Strand", shortName="S")
    String STRAND;
    
    @Option(doc="Annotations", shortName="A")
    String ANNOTATIONS;
    
    @Option(doc="Paired", shortName="P")
    String PAIRED;
    
    TranscriptionRead strand;
    Map<String,Collection<Gene>> annotations;
    boolean pairedData;
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new AlignmentStatistics().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			
			//Transcription strand
			strand = TranscriptionRead.UNSTRANDED;
			if(STRAND.equalsIgnoreCase("first")){
				strand = TranscriptionRead.FIRST_OF_PAIR;
			}
			else if(STRAND.equalsIgnoreCase("second")){
				strand = TranscriptionRead.SECOND_OF_PAIR;
			}
			else
				System.out.println("no strand");
			
			if(PAIRED.equals("true"))
				pairedData = true;
			else
				pairedData = false;
			
			/*
			 * Read the names of the alignment file into an array
			 */
			BufferedReader br = new BufferedReader(new FileReader(INPUT));
			List<String> alignmentFiles = new ArrayList<String>();
			String s;
			while((s = br.readLine())!= null){
				alignmentFiles.add(s);
			}
			
			/*
			 * The columns of the matrix will be each metrics
			 */
			List<String> cols = new ArrayList<String>();
			cols.add("TotalPairedReadsAligned");
			cols.add("TotalFirstMateAligned");
			cols.add("TotalSecondMateAligned");
			cols.add("TotalAligned");
			cols.add("TotalPairedCounts");
			cols.add("TotalFirstMateCounts");
			cols.add("TotalSecondMateCounts");
			cols.add("TotalCounts");
			
			MatrixWithHeaders countsMatrix = new MatrixWithHeaders(alignmentFiles,cols);
			
			for(String alignmentFile: alignmentFiles){
				AlignmentModel model=new AlignmentModel(new File(alignmentFile).getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),pairedData,strand,true);
							
				long pairedT = 0;
				long firstSingleT = 0; 
				long secondSingleT = 0;
				
				double weightedPairedCounterT=0.0;
				double weightedFirstSingleCounterT=0.0;
				double weightedSecondSingleCounterT=0.0;
				
				FileWriter writer2=new FileWriter(alignmentFile+".statistics");
				writer2.write("Chromosome\tPairedReads\tFirstMate\tSecondMate\tPairedCounts\tFirstMateCounts\tSecondMateCounts\n");
				for(String chr:model.getCoordinateSpace().getChromosomeNames()){
					long paired = 0;
					long firstSingle = 0; 
					long secondSingle = 0;
					
					double weightedPairedCounter=0.0;
					double weightedFirstSingleCounter=0.0;
					double weightedSecondSingleCounter=0.0;
				
					CloseableIterator<Alignment> iter= model.getOverlappingReads(chr);
					while(iter.hasNext()){
						Alignment read = iter.next();
						
						if(read instanceof AbstractPairedEndAlignment){
							paired++;
							weightedPairedCounter+=read.getWeight();
						}else{
							if(read instanceof SingleEndAlignment){
								if(((SingleEndAlignment) read).getIsFirstMate()){
									firstSingle++;
									weightedFirstSingleCounter+=read.getWeight();
								}
								else{
									secondSingle++;
									weightedSecondSingleCounter+=read.getWeight();
								}
							}
							else{
								System.out.println("Error! Neither paired non single alignment");
							}
						}
					}
					iter.close();
					pairedT +=paired;
					firstSingleT+=firstSingle;
					secondSingleT+=secondSingle;
					weightedPairedCounterT+=weightedPairedCounter;
					weightedFirstSingleCounterT += weightedFirstSingleCounter;
					weightedSecondSingleCounterT += weightedSecondSingleCounter;
					
					writer2.write(chr+"\t"+paired+"\t"+firstSingle+"\t"+secondSingle+"\t"+weightedPairedCounter+"\t"+weightedFirstSingleCounter+"\t"+weightedSecondSingleCounter+"\n");
				}		
				countsMatrix.set(alignmentFile, "TotalPairedReadsAligned", weightedPairedCounterT);
				countsMatrix.set(alignmentFile, "TotalFirstMateAligned", weightedFirstSingleCounterT);
				countsMatrix.set(alignmentFile, "TotalSecondMateAligned", weightedSecondSingleCounterT);
				countsMatrix.set(alignmentFile, "TotalAligned", (weightedPairedCounterT+weightedFirstSingleCounterT+weightedSecondSingleCounterT));
				
				writer2.write("\nTotalPairedReads\tTotalFirstMate\tTotalSecondMate\tTotalPairedCounts\tTotalFirstMateCounts\tTotalSecondMateCounts\n");
				writer2.write(pairedT+"\t"+firstSingleT+"\t"+secondSingleT+"\t"+weightedPairedCounterT+"\t"+weightedFirstSingleCounterT+"\t"+weightedSecondSingleCounterT+"\n");

				
				//ANNOTATIONS
				annotations= BEDFileParser.loadDataByChr(new File(ANNOTATIONS));
				pairedT = 0;
				firstSingleT = 0; 
				secondSingleT = 0;
				
				weightedPairedCounterT=0.0;
				weightedFirstSingleCounterT=0.0;
				weightedSecondSingleCounterT=0.0; 
				
				writer2.write("Chromosome\tGenes-PairedReads\tFirstMate\tSecondMate\tPairedCounts\tFirstMateCounts\tSecondMateCounts\n");

				for(String chr:annotations.keySet()){
					long paired = 0;
					long firstSingle = 0; 
					long secondSingle = 0;
					
					double weightedPairedCounter=0.0;
					double weightedFirstSingleCounter=0.0;
					double weightedSecondSingleCounter=0.0;
								
					for(Gene g:annotations.get(chr)){
						//ASSUME THE GENES ARE ALREADY COLLAPSED
						CloseableIterator<Alignment> iter = model.getOverlappingReads(g, false);
						while(iter.hasNext()){
							Alignment read = iter.next();
							if(read instanceof AbstractPairedEndAlignment){
								paired++;
								weightedPairedCounter+=read.getWeight();
							}else{
								if(read instanceof SingleEndAlignment){
									if(((SingleEndAlignment) read).getIsFirstMate()){
										firstSingle++;
										weightedFirstSingleCounter+=read.getWeight();
									}
									else{
										secondSingle++;
										weightedSecondSingleCounter+=read.getWeight();
									}
								}
								else{
									System.out.println("Error! Neither paired non single alignment");
								}
							}
						}
						iter.close();
					} 
					writer2.write(chr+"\t"+paired+"\t"+firstSingle+"\t"+secondSingle+"\t"+weightedPairedCounter+"\t"+weightedFirstSingleCounter+"\t"+weightedSecondSingleCounter+"\n");
					pairedT +=paired;
					firstSingleT+=firstSingle;
					secondSingleT+=secondSingle;
					weightedPairedCounterT+=weightedPairedCounter;
					weightedFirstSingleCounterT += weightedFirstSingleCounter;
					weightedSecondSingleCounterT += weightedSecondSingleCounter;
				}
				writer2.write("\nTotalPairedReads\tTotalFirstMate\tTotalSecondMate\tTotalPairedCounts\tTotalFirstMateCounts\tTotalSecondMateCounts\n");
				writer2.write(pairedT+"\t"+firstSingleT+"\t"+secondSingleT+"\t"+weightedPairedCounterT+"\t"+weightedFirstSingleCounterT+"\t"+weightedSecondSingleCounterT+"\n");
				countsMatrix.set(alignmentFile, "TotalPairedCounts", weightedPairedCounterT);
				countsMatrix.set(alignmentFile, "TotalFirstMateCounts", weightedFirstSingleCounterT);
				countsMatrix.set(alignmentFile, "TotalSecondMateCounts", weightedSecondSingleCounterT);
				countsMatrix.set(alignmentFile, "TotalCounts", (weightedPairedCounterT+weightedFirstSingleCounterT+weightedSecondSingleCounterT));
				
				writer2.close();
				
			}
			BufferedWriter outBw = new BufferedWriter(new FileWriter(OUTPUT));
			countsMatrix.write(outBw);
			outBw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return 0;
	}

	
	
}

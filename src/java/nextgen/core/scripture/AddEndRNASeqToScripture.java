package nextgen.core.scripture;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.alignment.PairedReadAlignment;
import nextgen.core.alignment.SingleEndAlignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.Window;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;

public class AddEndRNASeqToScripture {
	
	static Logger logger = Logger.getLogger(AddEndRNASeqToScripture.class.getName());
	
	private AlignmentModel model5p;
	private AlignmentModel model3p;
	private AlignmentModel model;
	//private BEDFileParser annotationParser;
	private int windowSize;
	private int extension;
	Map<String,Collection<Gene>> annotations;
	private TranscriptionRead strand;
	private static int DEFAULT_EXTENSION = 0;
	private static int DEFAULT_WINDOW_SIZE = 2;
	
	static final String usage = "Usage: AddEndRNASeqToScripture -task <task name> "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-3p <3P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\n\t\t-5p <5P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-annotations <Reconstruction bed file. [BED by default]> "+

			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\\nn\t\t-window <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends. Default = 5bp> "+
			"\n\t\t-extension <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends> "+

			"\n";
	
	public AddEndRNASeqToScripture(File bamFile5p,File bamFile3p,TranscriptionRead str,String annotationFile,String outputName,int windowS,File fullBam,int ext) throws IOException{
		
		model5p=new AlignmentModel(bamFile5p.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,str,false); 
		
		TranscriptionRead oppStrand = TranscriptionRead.UNSTRANDED;
		if(str == (TranscriptionRead.FIRST_OF_PAIR))
			oppStrand = TranscriptionRead.SECOND_OF_PAIR;
		else 
			if(str == (TranscriptionRead.SECOND_OF_PAIR))
				oppStrand = TranscriptionRead.FIRST_OF_PAIR;
		strand = str;
		model3p=new AlignmentModel(bamFile3p.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,oppStrand,false);
		model=new AlignmentModel(fullBam.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand);
		//Read annotation file
		//annotationParser = new BEDFileParser(annotationFile);	
				
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		windowSize = windowS;
		extension = ext;
		
		numberOfIsoformsPerGene(outputName);
		findCompleteTranscripts(outputName);
	}
	
	public void numberOfIsoformsPerGene(String outputName) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".isoforms.info"));
		BufferedWriter bwCov = new BufferedWriter(new FileWriter(outputName+".isoforms.coverage.info"));
		int singleExonGenes = 0;
		int multiExonGenes = 0;
		Map<Integer,Integer> singleIsoforms = new HashMap<Integer,Integer>();
		Map<Integer,Integer> multiIsoforms = new HashMap<Integer,Integer>();
		//For each chromosome in the annotation set
		/*Iterator<String> iter = annotationParser.getChromosomeIterator();
		while(iter.hasNext()){
			String chr = iter.next();*/
		for(String chr:annotations.keySet()){
			// Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
			//Iterator<GeneWithIsoforms> annotation_iter = annotationParser.getChrTree(chr).valueIterator();
			
			logger.info("Processing : "+chr);
			Map<Gene,Set<Gene>> isoformMap = BuildScriptureCoordinateSpace.getIsoformMap(annotations.get(chr));
			logger.info("Built isoform map");
			for(Gene gene:isoformMap.keySet()){
				
				if(gene.getBlocks().size()==1){
					//Get number of isoforms
					int cnt = 0;
					if(singleIsoforms.containsKey(isoformMap.get(gene).size())){
						cnt = singleIsoforms.get(isoformMap.get(gene).size());
					}	
					cnt++;
					singleIsoforms.put(isoformMap.get(gene).size(), cnt);
				}
				else{
					 multiExonGenes++;
					 int cnt = 0;
					 if(multiIsoforms.containsKey(isoformMap.get(gene).size())){
							cnt = multiIsoforms.get(isoformMap.get(gene).size());
					}	
					cnt++;
					multiIsoforms.put(isoformMap.get(gene).size(), cnt);
				}
				//For each isoform calculate coverage
				double avg=0.0;
				for(Gene isoform:isoformMap.get(gene)){
					avg+=new ScanStatisticScore(model,isoform).getAverageCoverage(model);
				}
				avg = avg/(double)isoformMap.get(gene).size();
				bwCov.write(isoformMap.get(gene).size()+"\t"+avg+"\n");
				logger.info(isoformMap.get(gene).size()+"\t"+avg+"\n");
			}
			
			logger.info("For chromosome "+chr);
			logger.info("Single Exon genes: "+singleExonGenes+"\n");
			logger.info("Multiple Exon genes: "+multiExonGenes+"\n");
			logger.info("\nSingle Exon genes: \n");
			for(int num:singleIsoforms.keySet()){
				logger.info(num+"\t"+singleIsoforms.get(num)+"\n");
			}
			logger.info("\nMultiple Exon genes: \n");
			for(int num:multiIsoforms.keySet()){
				logger.info(num+"\t"+multiIsoforms.get(num)+"\n");
			}
		}
		
			//While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
			//while(annotation_iter.hasNext()){
/*			for(Annotation gene: annotations.get(chr)){
				Gene gene = annotation_iter.next();
				
				if(gene.numBlocks()==1){
					singleExonGenes++;
					//Get number of isoforms
					int cnt = 0;
					if(singleIsoforms.containsKey(gene.getIsoforms().size())){
						cnt = singleIsoforms.get(gene.getIsoforms().size());
					}	
					cnt++;
					singleIsoforms.put(gene.getIsoforms().size(), cnt);
				}
				else{
					 multiExonGenes++;
					 int cnt = 0;
					 if(multiIsoforms.containsKey(gene.getIsoforms().size())){
							cnt = multiIsoforms.get(gene.getIsoforms().size());
					}	
					cnt++;
					multiIsoforms.put(gene.getIsoforms().size(), cnt);
				}
			}
		}*/
		bw.write("Single Exon genes: "+singleExonGenes+"\n");
		bw.write("Multiple Exon genes: "+multiExonGenes+"\n");
		bw.write("\nSingle Exon genes: \n");
		for(int num:singleIsoforms.keySet()){
			bw.write(num+"\t"+singleIsoforms.get(num)+"\n");
		}
		bw.write("\nMultiple Exon genes: \n");
		for(int num:multiIsoforms.keySet()){
			bw.write(num+"\t"+multiIsoforms.get(num)+"\n");
		}
		bw.close();
		bwCov.close();
	}
	
	public void findCompleteTranscripts(String outputName) throws IOException{
		
		BufferedWriter bw5p = new BufferedWriter(new FileWriter(outputName+".peaks.5p"));
		BufferedWriter bw5pBed = new BufferedWriter(new FileWriter(outputName+".peaks.5p.bed"));
		BufferedWriter bw3p = new BufferedWriter(new FileWriter(outputName+".peaks.3p"));
		BufferedWriter bw3pBed = new BufferedWriter(new FileWriter(outputName+".peaks.3p.bed"));
		
		int count = 0;
		//For each chromosome in the annotation set
		for(String chr:annotations.keySet()){
			
			int numFullSing = 0;
			int numPartialSing = 0;
			int num5pPartialSing = 0;
			int num3pPartialSing = 0;
			int numFullMult = 0;
			int numPartialMult = 0;
			int num5pPartialMult = 0;
			int num3pPartialMult = 0;
			//If 5' or 3' end RNA-seq does not have data for it, dont run
			if(model5p.getRefSequenceLambda(chr)==0.0 || model3p.getRefSequenceLambda(chr)==0.0){
				logger.warn(chr +" is not expressed in end RNA-seq");
			}
			else{
				logger.info("Processing " + chr);
				
				/*
				 * FOR THIS GENE
				 */
				for(Gene gene:annotations.get(chr)){
					
					//Make a coordinate space with the gene
					Map<String,Collection<Gene>> chrToGenesMap = new HashMap<String,Collection<Gene>>();
					List<Gene> g = new ArrayList<Gene>();
					g.add(gene);
					chrToGenesMap.put(gene.getChr(), g);
					CoordinateSpace space = new TranscriptomeSpace(chrToGenesMap);
					
					count++;
					/*
					 * IS THERE AT LEAST 1 5P END
					 */
					boolean has5pPeak = false;
					boolean has3pPeak = false;
					double[] nulls5p = get5pNullDistribution(gene,space);
					double[] nulls3p = get3pNullDistribution(gene,space);
/*					Collection<double[]> nulls = new ArrayList<double[]>();
					//Compute nulls for each isoform transcript
					for(Gene isoform:gene.getIsoforms()){
						double[] rtrn = get5pNullDistribution(isoform);
						if(rtrn[0]>0.0 || rtrn[1]>0.0)
							nulls.add(rtrn);
					}
					
					//For each isoform,
					for(Gene isoform:gene.getIsoforms()){
*/
						Collection<Annotation> peaks = new ArrayList<Annotation>();
						//Compute the z-score and take the minimum of the z-scores
//						Iterator<? extends Window> witer = model5p.getCoordinateSpace().getWindowIterator(windowSize, chr, isoform.getStart(), isoform.getEnd(), 0);
						//Iterator<? extends Window> witer = model5p.getCoordinateSpace().getWindowIterator(windowSize, chr, gene.getStart(), gene.getEnd(), 0);
						
						/*
						 * ITERATE IN THE TRANSCRIPTOME SPACE ONLY
						 */
						int start = 0;
						int end =0;
						if(gene.isNegativeStrand()){
							end = extension;
							//logger.info("Start: "+start+" End: "+end);
						}
						else{
							start = extension;
						}
						Annotation ge = gene.copy();
						ge.expand(start, end);

						Iterator<? extends Window> giter = space.getWindowIterator(ge, windowSize, 0);
						
						boolean flag5p = false;
						Annotation prev5p = null;
						Annotation peak5p =null;
						boolean flag3p = false;
						Annotation prev3p = null;
						Annotation peak3p =null;
						
						List<Annotation> this5pPeaks = new ArrayList<Annotation>();
						List<Annotation> this3pPeaks = new ArrayList<Annotation>();
						//For every window in the transcript
 						while(giter.hasNext()){
							Window window = giter.next();
							
							/*
							 * 5P 
							 */
							double windowCount5p = get5pWindowCount(window,gene.getOrientation());
							//Get the z-score of each window
//							double zscore = getMinimumZScore(windowCount,nulls,window.getSize());
							double zscore5p = Statistics.zScore(windowCount5p, nulls5p[0],nulls5p[1],window.getSize());
							//Associate the high z-scores with gene
							//If window is significant
							if(zscore5p>=7){
								//if flag=false, that is, no peak found before this(?)
								if(!flag5p){
									has5pPeak=true;
									//Set flag to true, start a new peak
									flag5p = true;
									peak5p = new BasicAnnotation(window);
									peak5p.setName(gene.getName());
									
									//IF THERE WAS A PREVIOUS PEAK AND THE DISTANCE BET THE TWO IS LESS THAN 25bp, merge
									if(!(prev5p==null)){
										if(((peak5p.getStart()-prev5p.getEnd())<=25 && (peak5p.getStart()-prev5p.getEnd())>=0) 
												|| ((prev5p.getStart()-peak5p.getEnd())<=25 && (prev5p.getStart()-peak5p.getEnd())>=0) ){
											peak5p.setStart(Math.min(peak5p.getStart(),prev5p.getStart()));
											peak5p.setEnd(Math.max(peak5p.getEnd(), prev5p.getEnd()));
											this5pPeaks.remove(prev5p);
											prev5p=null;
										}
									}
								}
								//peak was already started
								else{
									//extend it
									if(window.getStart()<peak5p.getStart())
										peak5p.setStart(window.getStart());
									if(window.getEnd()>peak5p.getEnd())
										peak5p.setEnd(window.getEnd());
								}
								peaks.add(window);
							}
							//Either end of window or no peak found yet
							else{
								//if flag=true
								if(flag5p){
									//Set flag = false, end peak and report to bed file
									flag5p = false;
									prev5p = peak5p;
									this5pPeaks.add(peak5p);
									peak5p = null;
								}
								else{
									//nothing
								}
							}
							
							/*
							 * 3p
							 */
							double windowCount3p = get3pWindowCount(window,gene.getOrientation());
							//Get the z-score of each window
//							double zscore = getMinimumZScore(windowCount,nulls,window.getSize());
							double zscore3p = Statistics.zScore(windowCount3p, nulls3p[0],nulls3p[1],window.getSize());
							//If window is significant
							if(zscore3p>=7){
								//if flag=false, that is, no peak found before this(?)
								if(!flag3p){
									has3pPeak=true;
									//Set flag to true, start a new peak
									flag3p = true;
									peak3p = new BasicAnnotation(window);
									peak3p.setName(gene.getName());
									
									//IF THERE WAS A PREVIOUS PEAK AND THE DISTANCE BET THE TWO IS LESS THAN 25bp, merge
									if(!(prev3p==null)){
										if(((peak3p.getStart()-prev3p.getEnd())<=25 && (peak3p.getStart()-prev3p.getEnd())>=0) 
												|| ((prev3p.getStart()-peak3p.getEnd())<=25 && (prev3p.getStart()-peak3p.getEnd())>=0) ){
											peak3p.setStart(Math.min(peak3p.getStart(),prev3p.getStart()));
											peak3p.setEnd(Math.max(peak3p.getEnd(), prev3p.getEnd()));
											this3pPeaks.remove(prev3p);
											prev3p=null;
										}
									}
								}
								//peak was already started
								else{
									//extend it
									if(window.getStart()<peak3p.getStart())
										peak3p.setStart(window.getStart());
									if(window.getEnd()>peak3p.getEnd())
										peak3p.setEnd(window.getEnd());
								}
								peaks.add(window);
							}
							//Either end of window or no peak found yet
							else{
								//if flag=true
								if(flag3p){
									//Set flag = false, end peak and report to bed file
									flag3p = false;
									prev3p = peak3p;
									this3pPeaks.add(peak3p);
									peak3p = null;
								}
								else{
									//nothing
								}
							}
						}
 						//Last peak
 						if(flag5p){
							this5pPeaks.add(peak5p);
						}
 						for(Annotation p:this5pPeaks){
							bw5pBed.write(p.toBED()+"\n");
							bw5p.write(gene.getName()+"\t"+p.toUCSC()+"\n");
						}
 						//Last peak
 						if(flag3p){
							this3pPeaks.add(peak3p);
						}
 						for(Annotation p:this3pPeaks){
							bw3pBed.write(p.toBED()+"\n");
							bw3p.write(gene.getName()+"\t"+p.toUCSC()+"\n");
						}
 						
 						if(gene.getBlocks().size()==1){
	 						if(has5pPeak){
	 							if(has3pPeak)
	 								numFullSing++;
	 							else
	 								//5p but no 3p
	 								num5pPartialSing++;
	 						}
	 						else{
	 							if(has3pPeak)
	 								num3pPartialSing++;
	 							else
	 								numPartialSing++;
	 						}
 						}
 						else{
 							if(has5pPeak){
	 							if(has3pPeak)
	 								numFullMult++;
	 							else
	 								//5p but no 3p
	 								num5pPartialMult++;
	 						}
	 						else{
	 							if(has3pPeak)
	 								num3pPartialMult++;
	 							else
	 								numPartialMult++;
	 						}
 						}
 						if(count%1000.0==0.0){
 							logger.info("Single: Number of Full= "+numFullSing+" Incomplete= "+numPartialSing+
 									" 5p No 3p= "+num5pPartialSing+" 3p No 5p= "+num3pPartialSing);
 							logger.info("Multiple: Number of Full= "+numFullMult+" Incomplete= "+numPartialMult+
 									" 5p No 3p= "+num5pPartialMult+" 3p No 5p= "+num3pPartialMult);
 						}
					}					
				}
			}
		bw5p.close();
		bw3p.close();
		bw5pBed.close();
		bw3pBed.close();
	}
	
	private double get5pWindowCount(Window window,Strand orientation){
		double windowCount = 0.0;
		//Get the reads in the window
		
		Iterator<Alignment> readiter = model5p.getOverlappingReads(window,false);
		//for all reads in the window
		while(readiter.hasNext()){
			Alignment read = readiter.next();
			if(passesChecks(read,window,orientation)){
				windowCount +=1.0;
			}
		}
		return windowCount;
	}
	
	private double get3pWindowCount(Window window,Strand orientation){
		double windowCount = 0.0;
		//Get the reads in the window
		
		Iterator<Alignment> readiter = model3p.getOverlappingReads(window,false);
		//for all reads in the window
		while(readiter.hasNext()){
			Alignment read = readiter.next();
			if(passesChecks(read,window,orientation)){
				windowCount +=1.0;
			}
		}
		return windowCount;
	}
	
	private double getMinimumZScore(double windowCount,Collection<double[]> nulls,int windowS){
		double zscore = Double.MAX_VALUE;
		for(double[] nullDist:nulls){
			zscore = Math.min(zscore, Statistics.zScore(windowCount, nullDist[0],nullDist[1],windowS));
		}
		return zscore;
	}

	/**
	 * CALCULATE THE NULL OVER THE ISOFORM
	 * @param annotation
	 * @return 	[0]: mean
	 * 			[1]: variance
	 */
	private double[] get5pNullDistribution(Gene annotation,CoordinateSpace space){
		
		List<Double> values = new ArrayList<Double>();
				
		Iterator<? extends Window> giter = space.getWindowIterator(annotation, windowSize, 0);
		//For each window
		while(giter.hasNext()){
			Window window = giter.next();
			double windowCount = 0.0;
			//Get the reads in the window
			//For each block in the window
			for(Annotation block: window.getBlocks()){
				
				Iterator<Alignment> readiter = model5p.getOverlappingReads(block,false);
				//for all reads in the window
				while(readiter.hasNext()){
					Alignment read = readiter.next();					
					if(passesChecks(read,window,annotation.getOrientation())){
						windowCount += 1.0;
					}				
				}
			}
			if(windowCount>0.0){
				values.add(windowCount);
			}
		}
		double[] rtrn = new double[2];
		if(values.size()>0){
			rtrn[0] = Statistics.mean(values);
			rtrn[1] = Statistics.variance(values);
		}
		else{ 
			rtrn[0] = 0.0;
			rtrn[1] = 0.0;
		}
		return rtrn;
	}
	
	/**
	 * CALCULATE THE NULL OVER THE ISOFORM
	 * @param annotation
	 * @return 	[0]: mean
	 * 			[1]: variance
	 */
	private double[] get3pNullDistribution(Gene annotation,CoordinateSpace space){
		
		List<Double> values = new ArrayList<Double>();
				
		Iterator<? extends Window> giter = space.getWindowIterator(annotation, windowSize, 0);
		//For each window
		while(giter.hasNext()){
			Window window = giter.next();
			double windowCount = 0.0;
			//Get the reads in the window
			//For each block in the window
			for(Annotation block: window.getBlocks()){
				
				Iterator<Alignment> readiter = model3p.getOverlappingReads(block,false);
				//for all reads in the window
				while(readiter.hasNext()){
					Alignment read = readiter.next();					
					if(passesChecks(read,window,annotation.getOrientation())){
						windowCount += 1.0;
					}				
				}
			}
			if(windowCount>0.0){
				values.add(windowCount);
			}
		}
		double[] rtrn = new double[2];
		if(values.size()>0){
			rtrn[0] = Statistics.mean(values);
			rtrn[1] = Statistics.variance(values);
		}
		else{ 
			rtrn[0] = 0.0;
			rtrn[1] = 0.0;
		}
		return rtrn;
	}
	
	/**
	 * Returns true is the specified alignment starts in window and 
	 * if Single ended, matches the mate of transcription
	 * if paired ended, the mate in the direction of transcription starts in the window
	 * @param read
	 * @param window
	 * @return
	 */
	private boolean passesChecks(Alignment read,Window window,Strand orientation){
		
		if(SingleEndAlignment.class.isInstance(read)){
			SingleEndAlignment align = (SingleEndAlignment) read;
			//Check if read is the correct read
			//if read starts in window
			if(((strand==(TranscriptionRead.FIRST_OF_PAIR) && align.getIsFirstMate()) || 
					(strand==(TranscriptionRead.SECOND_OF_PAIR) && !align.getIsFirstMate())) 
						&& (readStartFallsInWindow(read,window))
							&& (read.getOrientation().equals(orientation))){
				return true;
			}
		}
		//ELSE PAIRED
		else{
			PairedReadAlignment align = (PairedReadAlignment) read;
			Annotation mate;
			if(strand==TranscriptionRead.FIRST_OF_PAIR){
				mate = align.getFirstMate();
			}
			else{
				mate = align.getSecondMate();
			}
			if(readStartFallsInWindow(mate,window) && (read.getOrientation().equals(orientation))){
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns true if the oriented start of the read falls in the window
	 * @param align
	 * @param window
	 * @return
	 */
	private boolean readStartFallsInWindow(Annotation align,Window window){
		
		int start;
		if(align.isNegativeStrand()){
			start = align.getEnd();
		}
		else{
			start = align.getStart();
		}
		if(start>=window.getStart() && start<=window.getEnd())
			return true;
		else
			return false;
	}

	
	public static void main (String [] args) throws ParseException, IOException {
		
		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"5p");
		TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand").equalsIgnoreCase("first")){
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else{ 
			if(argMap.get("strand").equalsIgnoreCase("second")){
				strand = TranscriptionRead.SECOND_OF_PAIR;
			}
			else
				logger.warn("Strand is not first or second");
		}
		
		logger.info("Checking strand equality:");

		new AddEndRNASeqToScripture(new File(argMap.getMandatory("5p")),new File(argMap.getMandatory("3p")),strand,argMap.getMandatory("annotations"),argMap.getOutput(),argMap.getInteger("window", DEFAULT_WINDOW_SIZE),new File(argMap.getMandatory("full")),argMap.getInteger("extension", DEFAULT_EXTENSION));
		
	}
	
	/*static String usage=" args[0]=bam file \n\t args[1]=annotation file \n\t args[2]: Output name"
			+"\n\t args[3]= window size \n\t args[4] transcription strand";*/
	
}

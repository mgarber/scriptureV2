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

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.GeneWithIsoforms;

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

public class EndPileUpScanner {
	
	static Logger logger = Logger.getLogger(EndPileUpScanner.class.getName());
	
	private AlignmentModel model;
	private BEDFileParser annotationParser;
	private int window;
	private int extension;
	private TranscriptionRead strand;
	private static int DEFAULT_EXTENSION = 0;
	private static int DEFAULT_WINDOW_SIZE = 2;
	
	static final String usage = "Usage: EndPileUpScanner -task <task name> "+
			"\n\tTASK 1: 5p: Identifies the significant 5' gene ends in the end RNA-Seq data for a given annotation set." + 
			"\n\n\tTASK 2 3p: Identifies the significant 3' gene ends in the end RNA-Seq data for a given annotation set." + 
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-annotations <FOR ANNOTATE TASK ONLYAnnotation file for which to calculate expression. [BED by default]> "+

			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\\nn\t\t-window <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends. Default = 5bp> "+
			"\n\t\t-extension <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends> "+

			"\n";
	
	public EndPileUpScanner(File bamFile,TranscriptionRead str,String annotationFile,String outputName,int windowSize,int ext,String task) throws IOException{
		
		model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,str);
		
		//Read annotation file
		annotationParser = new BEDFileParser(annotationFile);	
		extension = ext;
		window = windowSize;
		strand = str;
		if(task.equalsIgnoreCase("5p"))
			assign5pPeaksToGenes(outputName);
		else{
			if(task.equalsIgnoreCase("3p"))
				assign3pPeaksToGenes(outputName);
			else
				logger.error(usage);
		}
	}
	
	public void assign5pPeaksToGenes(String outputName) throws IOException{
		//Output files
		BufferedWriter bwP = new BufferedWriter(new FileWriter(outputName+".plus.wig"));
		BufferedWriter bwM = new BufferedWriter(new FileWriter(outputName+".minus.wig"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName));
		BufferedWriter bwWindows = new BufferedWriter(new FileWriter(outputName+".windows"));
		BufferedWriter bwPeaks = new BufferedWriter(new FileWriter(outputName+".peaks.bed"));
		BufferedWriter bwPeaksTxt = new BufferedWriter(new FileWriter(outputName+".peaks.txt"));
		
		Map<Gene,Collection<Annotation>> geneToPeakMap = new HashMap<Gene,Collection<Annotation>>();
		 
		//For each chromosome in the annotation set
		Iterator<String> iter = annotationParser.getChromosomeIterator();
		while(iter.hasNext()){
			String chr = iter.next();
			
			//If the chromosome is not expressed in the sample
			double lambda = model.getRefSequenceLambda(chr);
			if(lambda==0.0){
				logger.warn(chr +" is not expressed in sample");
			}
			//If the chromosome is expressed in the sample
			else{
				logger.info("Processing " + chr);
				//Start of new chromosome
				bwP.write("variableStep chrom="+chr+" span="+window+"\n");
				bwM.write("variableStep chrom="+chr+" span="+window+"\n");
				
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 */
				Iterator<GeneWithIsoforms> annotation_iter = annotationParser.getChrTree(chr).valueIterator();
				/*
				 * Parse annotation tree
				 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
				 * Thus, for each gene
				 */
				while(annotation_iter.hasNext()){
					
					Gene gene = annotation_iter.next();
					//For each gene
					for(Gene annotation : gene.getIsoforms()){	
					Collection<Annotation> peaks = new ArrayList<Annotation>();
					List<Double> values = new ArrayList<Double>();
					
					//Get a gene window iterator over the entire gene
					Collection<Gene> genes = new ArrayList<Gene>();
					genes.add(annotation);
					
					Strand orientation = annotation.getOrientation();
					/*
					 * CALCULATE THE NULL OVER THE ISOFORM
					 */
					Map<String,Collection<Gene>> chrToGenesMap = new HashMap<String,Collection<Gene>>();
					List<Gene> g = new ArrayList<Gene>();
					g.add(annotation);
					chrToGenesMap.put(annotation.getChr(), g);
					//Make a coordinate space with the gene
					CoordinateSpace space = new TranscriptomeSpace(chrToGenesMap);
					
					Iterator<? extends Window> giter = space.getWindowIterator(annotation, window, 0);
					//For each window
					while(giter.hasNext()){
						Window window = giter.next();
						double windowCount = 0.0;
						//Get the reads in the window
						//For each block in the window
						for(Annotation block: window.getBlocks()){
							
							Iterator<Alignment> readiter = model.getOverlappingReads(block,false);
							//for all reads in the window
							while(readiter.hasNext()){
								Alignment read = readiter.next();
								
								if(passesChecks(read,window,orientation)){
									windowCount += 1.0;
								}
							
							}
						}
						if(windowCount>0.0){
							values.add(windowCount);
						}
					}
					
					if(values.size()>0){
						//Get the median and standard deviation
						double median = Statistics.median(values);
						double variance = Statistics.variance(values);
						double mean=Statistics.mean(values);
						logger.info("For "+annotation.getName()+" Median = "+median+" Variance = "+variance+ " Mean = "+mean);

						//Get a genome window iterator over the genomic region of gene and upstream region
						int start = annotation.getStart();
						int end  = annotation.getEnd();
						if(annotation.isNegativeStrand()){
							end = annotation.getEnd() + getUpstreamNonOverlappingDistance(annotationParser, annotation, extension);
							//logger.info("Start: "+start+" End: "+end);
						}
						else{
							start = annotation.getStart() - getUpstreamNonOverlappingDistance(annotationParser, annotation, extension);
						}
						
						/*
						 * SCAN THE GENOMIC (model's) COORDINATE SPACE
						 */
						Iterator<? extends Window> witer = model.getCoordinateSpace().getWindowIterator(window, chr, start, end, 0);
						boolean flag=false;
						Annotation prev = null;
						Annotation peak =null;
						List<Annotation> thisPeaks = new ArrayList<Annotation>();
 						while(witer.hasNext()){
							Window window = witer.next();
							double windowCount = 0.0;
							//Get the reads in the window
							
							Iterator<Alignment> readiter = model.getOverlappingReads(window,false);
							//for all reads in the window
							while(readiter.hasNext()){
								Alignment read = readiter.next();
								if(passesChecks(read,window,orientation)){
									windowCount +=1.0;
								}
							}

							//Get the z-score of each window
							double zscore = Statistics.zScore(windowCount, mean, variance,window.getSize());
							double pval = Statistics.zscoreToPvalue(zscore);
							
							//Associate the high z-scores with gene
							//If window is significant
							if(zscore>=7){
								//if flag=false, that is, no peak found before this(?)
								if(!flag){
									//Set flag to true, start a new peak
									flag = true;
									peak = new BasicAnnotation(window);
									peak.setName(annotation.getName());
									
									//IF THERE WAS A PREVIOUS PEAK AND THE DISTANCE BET THE TWO IS LESS THAN 25bp, merge
									if(!(prev==null)){
										if(((peak.getStart()-prev.getEnd())<=25 && (peak.getStart()-prev.getEnd())>=0) 
												|| ((prev.getStart()-peak.getEnd())<=25 && (prev.getStart()-peak.getEnd())>=0) ){
											peak.setStart(Math.min(peak.getStart(),prev.getStart()));
											peak.setEnd(Math.max(peak.getEnd(), prev.getEnd()));
											thisPeaks.remove(prev);
											prev=null;
										}
									}
								}
								//peak was already started
								else{
									//extend it
									if(window.getStart()<peak.getStart())
										peak.setStart(window.getStart());
									if(window.getEnd()>peak.getEnd())
										peak.setEnd(window.getEnd());
								}
								peaks.add(window);
							}
							//Either end of window or no peak found yet
							else{
								//if flag=true
								if(flag){
									//Set flag = false, end peak and report to bed file
									flag = false;
									prev = peak;
									thisPeaks.add(peak);
									peak = null;
								}
								else{
									//nothing
								}
							}
							if(annotation.isNegativeStrand()){
								if(windowCount>0.0 && !Double.isInfinite(zscore))
									bwM.write(window.getStart()+"\t"+zscore+"\n");							
							}
							else{
								if(windowCount>0.0 && !Double.isInfinite(zscore))
									bwP.write(window.getStart()+"\t"+zscore+"\n");
							}
							if(windowCount>0.0 && (!Double.isInfinite(zscore)||!Double.isNaN(zscore))){
								bw.write(annotation.getName()+"\t"+window.toUCSC()+"\t"+windowCount+"\t"+zscore+"\t"+pval+"\n");
							}
						}
						if(flag){
							thisPeaks.add(peak);
						}
						if(peaks.size()>0)
							geneToPeakMap.put(annotation, peaks);
						for(Annotation p:thisPeaks){
							bwPeaks.write(p.toBED()+"\n");
							bwPeaksTxt.write(annotation.getName()+"\t"+p.toUCSC());
						}
					}
				}
				}
			}
		}
		
		
		bwP.close();
		bwM.close();
		BufferedWriter bwbed = new BufferedWriter(new FileWriter(outputName+".bed"));
		for(Gene g:geneToPeakMap.keySet()){
			for(Annotation p:geneToPeakMap.get(g)){
				bwWindows.write(g.getName()+"\t"+p.toUCSC()+"\n");
				bwbed.write(new Gene(p).toBED()+"\n");
			}
		}
		bw.close();
		bwbed.close();
		bwWindows.close();
					
	}
	
	public void assign3pPeaksToGenes(String outputName) throws IOException{
		//Output files
		BufferedWriter bwP = new BufferedWriter(new FileWriter(outputName+".plus.wig"));
		BufferedWriter bwM = new BufferedWriter(new FileWriter(outputName+".minus.wig"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName));
		BufferedWriter bwPeaks = new BufferedWriter(new FileWriter(outputName+".peaks"));
		
		Map<Gene,Collection<Annotation>> geneToPeakMap = new HashMap<Gene,Collection<Annotation>>();
		 
		//For each chromosome in the annotation set
		Iterator<String> iter = annotationParser.getChromosomeIterator();
		while(iter.hasNext()){
			String chr = iter.next();
			
			//If the chromosome is not expressed in the sample
			double lambda = model.getRefSequenceLambda(chr);
			if(lambda==0.0){
				logger.warn(chr +" is not expressed in sample");
			}
			//If the chromosome is expressed in the sample
			else{
				logger.info("Processing " + chr);
				//Start of new chromosome
				bwP.write("variableStep chrom="+chr+" span="+window+"\n");
				bwM.write("variableStep chrom="+chr+" span="+window+"\n");
				
				/*
				 * Obtain an iterator over the interval tree values built from the annotations (genes) on that particular chromosome.
				 */
				Iterator<GeneWithIsoforms> annotation_iter = annotationParser.getChrTree(chr).valueIterator();
				/*
				 * Parse annotation tree
				 * While there is an annotated RefSeqGeneWithIsoforms in Interval tree to analyze
				 * Thus, for each gene
				 */
				while(annotation_iter.hasNext()){
					
					Gene gene = annotation_iter.next();
					//For each gene
					for(Gene annotation : gene.getIsoforms()){	
					Collection<Annotation> peaks = new ArrayList<Annotation>();
					List<Double> values = new ArrayList<Double>();
					
					//Get a gene window iterator over the entire gene
					Collection<Gene> genes = new ArrayList<Gene>();
					genes.add(annotation);
					
					Strand orientation = annotation.getOrientation();
					/*
					 * CALCULATE THE NULL OVER THE ISOFORM
					 */
					Map<String,Collection<Gene>> chrToGenesMap = new HashMap<String,Collection<Gene>>();
					List<Gene> g = new ArrayList<Gene>();
					g.add(annotation);
					chrToGenesMap.put(annotation.getChr(), g);
					//Make a coordinate space with the gene
					CoordinateSpace space = new TranscriptomeSpace(chrToGenesMap);
					
					Iterator<? extends Window> giter = space.getWindowIterator(annotation, window, 0);
					//For each window
					while(giter.hasNext()){
						Window window = giter.next();
						double windowCount = 0.0;
						//Get the reads in the window
						//For each block in the window
						for(Annotation block: window.getBlocks()){
							
							Iterator<Alignment> readiter = model.getOverlappingReads(block,false);
							//for all reads in the window
							while(readiter.hasNext()){
								Alignment read = readiter.next();
								
								if(passesChecks(read,window,orientation)){
									windowCount += 1.0;
								}
							}
						}
						if(windowCount>0.0){
							values.add(windowCount);
						}
					}
					
					if(values.size()>0){
						//Get the median and standard deviation
						double median = Statistics.median(values);
						double variance = Statistics.variance(values);
						double mean=Statistics.mean(values);
						logger.info("For "+annotation.getName()+" Median = "+median+" Variance = "+variance+ " Mean = "+mean);

						//Get a genome window iterator over the genomic region of gene and upstream region
						int start = annotation.getStart();
						int end  = annotation.getEnd();
						if(annotation.isNegativeStrand()){
							end = annotation.getStart() - getDownstreamNonOverlappingDistance(annotationParser, annotation, extension);
							//logger.info("Start: "+start+" End: "+end);
						}
						else{
							start = annotation.getEnd() + getDownstreamNonOverlappingDistance(annotationParser, annotation, extension);
						}
						
						/*
						 * SCAN THE GENOMIC (model's) COORDINATE SPACE
						 */
						Iterator<? extends Window> witer = model.getCoordinateSpace().getWindowIterator(window, chr, start, end, 0);
						while(witer.hasNext()){
							Window window = witer.next();
							double windowCount = 0.0;
							//Get the reads in the window
							
							Iterator<Alignment> readiter = model.getOverlappingReads(window,false);
							//for all reads in the window
							while(readiter.hasNext()){
								Alignment read = readiter.next();
								if(passesChecks(read,window,orientation)){
									windowCount +=1.0;
								}
							}

							//Get the z-score of each window
							double zscore = Statistics.zScore(windowCount, mean, variance,window.getSize());
							double pval = Statistics.zscoreToPvalue(zscore);
							
							//Associate the high z-scores with gene
							if(zscore>6){
								peaks.add(window);
							}
							if(annotation.isNegativeStrand()){
								if(windowCount>0.0 && !Double.isInfinite(zscore))
									bwM.write(window.getStart()+"\t"+zscore+"\n");							
							}
							else{
								if(windowCount>0.0 && !Double.isInfinite(zscore))
									bwP.write(window.getStart()+"\t"+zscore+"\n");
							}
							if(windowCount>0.0 && (!Double.isInfinite(zscore)||!Double.isNaN(zscore)))
								bw.write(annotation.getName()+"\t"+window.toUCSC()+"\t"+windowCount+"\t"+zscore+"\t"+pval+"\n");
						}
						
						if(peaks.size()>0)
							geneToPeakMap.put(annotation, peaks);
					}
				}
				}
			}
		}
		
		
		bwP.close();
		bwM.close();
		BufferedWriter bwbed = new BufferedWriter(new FileWriter(outputName+".bed"));
		for(Gene g:geneToPeakMap.keySet()){
			for(Annotation p:geneToPeakMap.get(g)){
				bwPeaks.write(g.getName()+"\t"+p.toUCSC()+"\n");
				bwbed.write(new Gene(p).toBED()+"\n");
			}
		}
		bw.close();
		bwbed.close();
		bwPeaks.close();
					
	}

	
	private int getUpstreamNonOverlappingDistance(BEDFileParser annotationCollection,Gene thisGene,int maxRegion){
		
		if(thisGene.isNegativeStrand()){
			Gene closestUpstreamGene = annotationCollection.getClosestDownstream(thisGene);
			if(closestUpstreamGene==null){
				return maxRegion;
			}
			else{
				int distanceFromThisGene = closestUpstreamGene.getStart()-thisGene.getEnd();
				if(distanceFromThisGene>maxRegion){
					return maxRegion;
				}
				else{
					return distanceFromThisGene;
				}
			}
		}
		else{
			Gene closestUpstreamGene = annotationCollection.getClosestUpstream(thisGene);
			if(closestUpstreamGene==null){
				return maxRegion;
			}
			else{
				int distanceFromThisGene = thisGene.getStart() - closestUpstreamGene.getEnd();
				if(distanceFromThisGene>maxRegion){
					return maxRegion;
				}
				else{
					return distanceFromThisGene;
				}
			}
		}
	}
	
	private int getDownstreamNonOverlappingDistance(BEDFileParser annotationCollection,Gene thisGene,int maxRegion){
		
		if(!thisGene.isNegativeStrand()){
			Gene closestDownstreamGene = annotationCollection.getClosestDownstream(thisGene);
			if(closestDownstreamGene==null){
				return maxRegion;
			}
			else{
				int distanceFromThisGene = closestDownstreamGene.getStart()-thisGene.getEnd();
				if(distanceFromThisGene>maxRegion){
					return maxRegion;
				}
				else{
					return distanceFromThisGene;
				}
			}
		}
		else{
			Gene closestDownstreamGene = annotationCollection.getClosestUpstream(thisGene);
			if(closestDownstreamGene==null){
				return maxRegion;
			}
			else{
				int distanceFromThisGene = thisGene.getStart() - closestDownstreamGene.getEnd();
				if(distanceFromThisGene>maxRegion){
					return maxRegion;
				}
				else{
					return distanceFromThisGene;
				}
			}
		}
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
			if(((strand==TranscriptionRead.FIRST_OF_PAIR && align.getIsFirstMate()) || 
					(strand==TranscriptionRead.SECOND_OF_PAIR && !align.getIsFirstMate())) 
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
			if(readStartFallsInWindow(mate,window)){
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns true is the specified alignment starts in window and 
	 * if Single ended, matches the opposite mate of transcription
	 * if paired ended, the mate not in the direction of transcription starts in the window
	 * @param read
	 * @param window
	 * @return
	 */
	private boolean passesOppositeChecks(Alignment read,Window window,Strand orientation){
		
		if(SingleEndAlignment.class.isInstance(read)){
			SingleEndAlignment align = (SingleEndAlignment) read;
			//Check if read is the correct read
			//if read starts in window
			if(((strand==TranscriptionRead.FIRST_OF_PAIR && !align.getIsFirstMate()) || 
					(strand==TranscriptionRead.SECOND_OF_PAIR && align.getIsFirstMate())) 
						&& (readStartFallsInWindow(read,window))
							&& (!read.getOrientation().equals(orientation))){
				return true;
			}
		}
		//ELSE PAIRED
		else{
			PairedReadAlignment align = (PairedReadAlignment) read;
			Annotation mate;
			if(strand==TranscriptionRead.FIRST_OF_PAIR){
				mate = align.getSecondMate();
			}
			else{
				mate = align.getFirstMate();
			}
			if(readStartFallsInWindow(mate,window)){
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
		
		/*File bamFile=new File(args[0]);
		TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
		if(args[4].equalsIgnoreCase("first")){
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(args[4].equalsIgnoreCase("second")){
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}
		new EndPileUpScanner(bamFile,strand,args[1],args[2],new Integer(args[3]));*/
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
		
		new EndPileUpScanner(new File(argMap.getMandatory("alignment")),strand,argMap.getMandatory("annotations"),argMap.getOutput(),argMap.getInteger("window", DEFAULT_WINDOW_SIZE),argMap.getInteger("extension", DEFAULT_EXTENSION),argMap.getTask());
		
	}
	
	/*static String usage=" args[0]=bam file \n\t args[1]=annotation file \n\t args[2]: Output name"
			+"\n\t args[3]= window size \n\t args[4] transcription strand";*/
	
}

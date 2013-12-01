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
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

import net.sf.samtools.util.CloseableIterator;
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

/**
 * This class uses end RNASeq data to trim an annotation set
 * 1. To find complete transcripts using 3' and 5' data
 * 2. Only 3' ends using RNaseH data
 * 3. Only 5' ends using ExoCAGE data.
 * 
 * This class uses a window to scan the transcript and measures "pile ups" of reads starting/ending in each window
 * This is used as the null distribution and a zscore is calculated for each window to determine significant pileups
 * These significant pileups are used to trim annotations.
 * 
 * @author skadri
 *
 */
public class AddEndRNASeqToScripture {
	
	static Logger logger = Logger.getLogger(AddEndRNASeqToScripture.class.getName());
	
	private AlignmentModel model5p;
	private AlignmentModel model3p;
	//private AlignmentModel model;
	//private BEDFileParser annotationParser;
	private int windowSize;
	private int extension;
	Map<String,Collection<Gene>> annotations;
	Map<String,Collection<Gene>> complete;
	private TranscriptionRead strand5p;
	private TranscriptionRead strand3p;
	private static int DEFAULT_EXTENSION = 0;
	private static int DEFAULT_WINDOW_SIZE = 2;
	Map<String, IntervalTree<Gene>> intervalTrees;
	private double THRESHOLD = 7.0;
	private static final double DEFAULT_THRESHOLD = 7.0;
	private boolean singleEnd3p;
	
	static final String usage = "Usage: AddEndRNASeqToScripture -task <task name> "+
			"\n\tcompleteTranscripts"+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-3p <3P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\n\t\t-5p <5P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-annotations <Reconstruction bed file. [BED by default]> "+
			
			"\n\t trim3p"+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-3p <3P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-annotations <Reconstruction bed file. [BED by default]> "+

			"\n\t trim5p"+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-5p <5P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-annotations <Reconstruction bed file. [BED by default]> "+

//			"\n\n\t\t-full <full length Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+

			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\n\t\t-window <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends. Default = 5bp> "+
			"\n\t\t-extension <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends> "+
			"\n\t\t-strand3p <Specifies the mate that is in the direction of trnascription for 3' data. Options: first,second,unstranded. Default: Unstranded> "+
			"\n\t\t-strand5p <Specifies the mate that is in the direction of trnascription for 5' data. Options: first,second,unstranded. Default: Unstranded> "+
			"\n\t\t-singleEnd <if specified, means that 3' data is single end"+
			"\n";
	
	/**
	 * Instantiates the class and calls the function that finds complete transcripts and trims the reconstructions
	 * @throws IOException
	 */
	public AddEndRNASeqToScripture(ArgumentMap argMap) throws IOException{
					
		//Read the annotations
		annotations= BEDFileParser.loadDataByChr(new File(argMap.getMandatory("annotations")));
		initiateIntervalTrees(annotations);
		windowSize = argMap.getInteger("window", DEFAULT_WINDOW_SIZE);
		extension = argMap.getInteger("extension", DEFAULT_EXTENSION);

		//1. if task is complete transcripts
		if(argMap.getTask().equalsIgnoreCase("completeTranscripts")){
			
			model3p=load3pData(argMap);
			model5p=load5pData(argMap);
			findCompleteTranscripts(argMap.getOutput());
		}
		else{
			if(argMap.getTask().equalsIgnoreCase("trim3p")){
				
				model3p=load3pData(argMap);
				trim3pEnds(argMap.getOutput());
			}
			else{
				if(argMap.getTask().equalsIgnoreCase("trim5p")){
					model5p=load5pData(argMap);
					trim5pEnds(argMap.getOutput());
				}
				else{
					logger.info("Incorrect task. Please choose from:\n"+usage);
				}
			}
		}
		
		//numberOfIsoformsPerGene(outputName,fullBam);
	}
	
	/**
	 * Load the 3p model
	 * @param argMap
	 * @return
	 */
	private AlignmentModel load3pData(ArgumentMap argMap){
		
		//Strand for 3' data
		strand3p = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand3p").equalsIgnoreCase("first")){
			strand3p = TranscriptionRead.FIRST_OF_PAIR;
		}
		else{ 
			if(argMap.get("strand3p").equalsIgnoreCase("second")){
				strand3p = TranscriptionRead.SECOND_OF_PAIR;
			}
			else
				logger.warn("3' data is unstranded");
		}	

		//Determine if 3' end is single ended
		if(argMap.isPresent("singleEnd")){
			logger.info("3' data is single ended");
			singleEnd3p=true;
		}
		else{
			logger.info("3' data is paired ended");
			singleEnd3p =false;
		}
		
		AlignmentModel model=new AlignmentModel(new File(argMap.getMandatory("3p")).getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),!singleEnd3p,strand3p,false);
		
		return model;
	}
	
	
	/**
	 * Load the 5p model
	 * @param argMap
	 * @return
	 */
	private AlignmentModel load5pData(ArgumentMap argMap){
		
		//Strand for 5' data
		strand5p = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand5p").equalsIgnoreCase("first")){
			strand5p = TranscriptionRead.FIRST_OF_PAIR;
		}
		else{ 
			if(argMap.get("strand5p").equalsIgnoreCase("second")){
				strand5p = TranscriptionRead.SECOND_OF_PAIR;
			}
			else
				logger.warn("5' data is unstranded");
		}	

		AlignmentModel model=new AlignmentModel(new File(argMap.getMandatory("5p")).getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand5p,false); 
		
		return model;
	}
	
	/**
	 * Calls the function that calculates the number of isoforms per gene and outputs #isoforms to coverage map
	 * @param str
	 * @param annotationFile
	 * @param outputName
	 * @param windowS
	 * @param fullBam
	 * @param ext
	 * @throws IOException
	 */
	public AddEndRNASeqToScripture(TranscriptionRead str,String annotationFile,String outputName,int windowS,File fullBam,int ext) throws IOException{
				
		TranscriptionRead oppStrand = TranscriptionRead.UNSTRANDED;
		if(str == (TranscriptionRead.FIRST_OF_PAIR))
			oppStrand = TranscriptionRead.SECOND_OF_PAIR;
		else 
			if(str == (TranscriptionRead.SECOND_OF_PAIR))
				oppStrand = TranscriptionRead.FIRST_OF_PAIR;
		//strand = str;
		//Read annotation file
		//annotationParser = new BEDFileParser(annotationFile);	
				
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		windowSize = windowS;
		extension = ext;
						
		numberOfIsoformsPerGene(outputName,fullBam);
		//findCompleteTranscripts(outputName);
	}

	/**
	 * Builds the isoform map for each chromosome
	 * @param ann
	 * @return
	 */
	public static IsoformMap buildIsoformMap(Map<String,Collection<Gene>> ann){
	
		IsoformMap map = new IsoformMap();
		for(String chr:ann.keySet()){
			logger.info("Processing : "+chr);
			map.addChrMap(chr, BuildScriptureCoordinateSpace.getIsoformMap(ann.get(chr)));
			logger.info("Built isoform map for "+ chr);
		}
		
		return map;
	}
	
	/**
	 * Calculates the number of isoforms per gene and a map of #isforms to coverage
	 * @param outputName
	 * @param fullBam
	 * @throws IOException
	 */
	public void numberOfIsoformsPerGene(String outputName,File fullBam) throws IOException{
		
		IsoformMap isoformMap = buildIsoformMap(annotations);
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".isoforms.info"));
		BufferedWriter bwCov = new BufferedWriter(new FileWriter(outputName+".isoforms.coverage.info"));
		//model=new AlignmentModel(fullBam.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand5p);
		
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
			for(Gene gene:isoformMap.getGenesForChromosome(chr)){
				//gene.overlapsStranded(other);
				if(gene.getBlocks().size()==1){
					//Get number of isoforms
					int cnt = 0;
					if(singleIsoforms.containsKey(isoformMap.getNumOfIsoformsForGene(gene))){
						cnt = singleIsoforms.get(isoformMap.getNumOfIsoformsForGene(gene));
					}	
					cnt++;
					singleIsoforms.put(isoformMap.getNumOfIsoformsForGene(gene), cnt);
				}
				else{
					 multiExonGenes++;
					 int cnt = 0;
					 if(multiIsoforms.containsKey(isoformMap.getNumOfIsoformsForGene(gene))){
							cnt = multiIsoforms.get(isoformMap.getNumOfIsoformsForGene(gene));
					}	
					cnt++;
					multiIsoforms.put(isoformMap.getNumOfIsoformsForGene(gene), cnt);
				}
				//For each isoform calculate coverage
				double avg=0.0;
				for(Gene isoform:isoformMap.getIsoformsForGene(gene)){
//					avg+=new ScanStatisticScore(model,isoform).getAverageCoverage(model);
				}
				avg = avg/(double)isoformMap.getNumOfIsoformsForGene(gene);
				bwCov.write(isoformMap.getNumOfIsoformsForGene(gene)+"\t"+avg+"\n");
				//logger.info(isoformMap.getNumOfIsoformsForGene(gene)+"\t"+avg);
			}
			
			logger.info("For chromosome "+chr);
			logger.info("Single Exon genes: "+singleExonGenes);
			logger.info("Multiple Exon genes: "+multiExonGenes);
			logger.info("\nSingle Exon genes: \n");
			for(int num:singleIsoforms.keySet()){
				logger.info(num+"\t"+singleIsoforms.get(num));
			}
			logger.info("\nMultiple Exon genes: \n");
			for(int num:multiIsoforms.keySet()){
				logger.info(num+"\t"+multiIsoforms.get(num));
			}
		}
		
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
	
	/**
	 * Calculates number of complete, incomplete and partial transcripts based on end RNA-Seq data
	 * @param outputName
	 * @throws IOException
	 */
	public void findCompleteTranscripts(String outputName) throws IOException{
		
		BufferedWriter bw5p = new BufferedWriter(new FileWriter(outputName+".peaks.5p"));
		BufferedWriter bw5pBed = new BufferedWriter(new FileWriter(outputName+".peaks.5p.bed"));
		BufferedWriter bw3p = new BufferedWriter(new FileWriter(outputName+".peaks.3p"));
		BufferedWriter bw3pBed = new BufferedWriter(new FileWriter(outputName+".peaks.3p.bed"));
		BufferedWriter bwComplete = new BufferedWriter(new FileWriter(outputName+".complete.bed"));
		
		Map<Gene,List<Annotation>> geneTo5pPeakMap = new HashMap<Gene,List<Annotation>>();
		Map<Gene,List<Annotation>> geneTo3pPeakMap = new HashMap<Gene,List<Annotation>>();
		int count = 0;
		complete = new TreeMap<String,Collection<Gene>>();
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
			logger.info("Processing "+chr);
			complete.put(chr, new TreeSet<Gene>());
			//If 5' or 3' end RNA-seq does not have data for it, dont run
			if(!model5p.containsReference(chr)
					&& !model3p.containsReference(chr)){
				logger.warn(chr +" is not expressed in end RNA-seq");
			}
			else{				
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
					
					int start = 0;
					int end =0;
					int bestExtension = getDistanceToClosestSameOrientation5pGene(gene);
					//logger.info(gene.getName()+" getDistanceToClosestSameOrientation5pGene "+bestExtension);
					if(bestExtension>extension || bestExtension<0){
						bestExtension = extension;
					}
					if(gene.isNegativeStrand()){
						if((gene.getEnd()+bestExtension)>model3p.getCoordinateSpace().getReferenceAnnotation(chr).getEnd()){
							bestExtension=model3p.getCoordinateSpace().getReferenceAnnotation(chr).getEnd()-gene.getEnd();
						}
						end = bestExtension;
					}
					else{
						if((gene.getStart()-bestExtension)<0){
							bestExtension=gene.getStart();
						}
						start = bestExtension;
					}
					Gene ge = gene.copy();
					//EXPAND IS A STRAND-INDEPENDENT FUNCTION
					ge.expand(start, end);
					//logger.info(ge.toBED());
					chrToGenesMap = new HashMap<String,Collection<Gene>>();
					g = new ArrayList<Gene>();
					g.add(ge);
					chrToGenesMap.put(gene.getChr(), g);
					space = new TranscriptomeSpace(chrToGenesMap);
					
					double[] nulls5p = get5pNullDistribution(ge,space);
					//logger.info("5p null: "+nulls5p[0]+" "+nulls5p[1]+" "+nulls5p[2]+" ");
					if(nulls5p[2]<11){
						THRESHOLD=15;
						
					}
					else{
						THRESHOLD = DEFAULT_THRESHOLD;
					}
					start = 0;
					end =0;
					bestExtension = getDistanceToClosestSameOrientation3pGene(gene);
					//logger.info(gene.getName()+" getDistanceToClosestSameOrientation3pGene "+bestExtension);
					if(bestExtension>extension || bestExtension<0){
						bestExtension = extension;
					}
					if(gene.isNegativeStrand()){
						if((gene.getStart()-bestExtension)<0){
							bestExtension=gene.getStart();
						}
						start = bestExtension;
						//
					}
					else{
						if((gene.getEnd()+bestExtension)>model3p.getCoordinateSpace().getReferenceAnnotation(chr).getEnd()){
							bestExtension=model3p.getCoordinateSpace().getReferenceAnnotation(chr).getEnd()-gene.getEnd();
						}
						end = bestExtension;
					}
					//logger.info("Start: "+start+" End: "+end);
					Gene gs = gene.copy();
					//EXPAND IS A STRAND-INDEPENDENT FUNCTION
					gs.expand(start, end);
					chrToGenesMap = new HashMap<String,Collection<Gene>>();
					g = new ArrayList<Gene>();
					g.add(gs);
					chrToGenesMap.put(gene.getChr(), g);
					space = new TranscriptomeSpace(chrToGenesMap);
					
					double[] nulls3p = get3pNullDistribution(gs,space);
						/*
						 * RE-DEFINE COORDINATE SPACE
						 */
						chrToGenesMap = new HashMap<String,Collection<Gene>>();
						g = new ArrayList<Gene>();
						g.add(ge);
						chrToGenesMap.put(gene.getChr(), g);
						space = new TranscriptomeSpace(chrToGenesMap);
						Iterator<? extends Window> giter = space.getWindowIterator(ge, windowSize, 0);
						
						boolean flag5p = false;
						Annotation prev5p = null;
						Annotation peak5p =null;
						
						List<Annotation> this5pPeaks = new ArrayList<Annotation>();
						//For every window in the transcript
 						while(giter.hasNext()){
							Window window = giter.next();
							
							/*
							 * 5P 
							 */
							double windowCount5p = get5pWindowCount(window,gene.getOrientation());
							double zscore5p = Statistics.zScore(windowCount5p, nulls5p[0],nulls5p[1],window.getSize());
							//Associate the high z-scores with gene
							//If window is significant
							if(zscore5p>=THRESHOLD){
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
								//this5pPeaks.add(window);
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
						}
 						//Last peak
 						if(flag5p){
							this5pPeaks.add(peak5p);
						}
 						for(Annotation p:this5pPeaks){
 							double windowCount = get5pWindowCount(p,gene.getOrientation());
							double zscore = Statistics.zScore(windowCount, nulls5p[0],nulls5p[1],p.getSize());
							p.setScore(zscore);
							bw5pBed.write(p.toBED()+"\n");
							bw5p.write(gene.getName()+"\t"+p.toUCSC()+"\t"+windowCount+"\t"+zscore+"\t"+calculate5pDistance(gene,p)+"\n");
						}
 						
 						/**
 						 * 3P 
 						 */
 						if(nulls3p[2]<11){
 							THRESHOLD=15;
 						}
 						else{
 							THRESHOLD = DEFAULT_THRESHOLD;
 						}
 						chrToGenesMap = new HashMap<String,Collection<Gene>>();
						g = new ArrayList<Gene>();
						g.add(gs);
						chrToGenesMap.put(gene.getChr(), g);
						space = new TranscriptomeSpace(chrToGenesMap);
						
						Iterator<? extends Window> iter = space.getWindowIterator(gs, windowSize, 0);
						boolean flag3p = false;
						Annotation prev3p = null;
						Annotation peak3p =null;
						
						List<Annotation> this3pPeaks = new ArrayList<Annotation>();
						//For every window in the transcript
 						while(iter.hasNext()){
							Window window = iter.next();							
							/*
							 * 3p
							 */
							double windowCount3p = get3pWindowCount(window,gene.getOrientation());
//							if(windowCount3p>0)
//								logger.info(window.toUCSC()+" "+windowCount3p);
							//Get the z-score of each window
//							double zscore = getMinimumZScore(windowCount,nulls,window.getSize());
							double zscore3p = Statistics.zScore(windowCount3p, nulls3p[0],nulls3p[1],window.getSize());
//							if(windowCount3p>0)
//								logger.info(window.toUCSC()+" Count = "+windowCount3p+" zscore = "+zscore3p);
							//If window is significant
							if(zscore3p>=THRESHOLD){
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
								//this3pPeaks.add(window);
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
 						if(flag3p){
							this3pPeaks.add(peak3p);
						}
 						for(Annotation p:this3pPeaks){
 							double windowCount = get3pWindowCount(p,gene.getOrientation());
							double zscore = Statistics.zScore(windowCount, nulls3p[0],nulls3p[1],p.getSize());
							p.setScore(zscore);							
							bw3pBed.write(p.toBED()+"\n");
							bw3p.write(gene.getName()+"\t"+p.toUCSC()+"\t"+windowCount+"\t"+zscore+"\t"+calculate3pDistance(gene,p)+"\n");
						}
 						
 						if(gene.getBlocks().size()==1){
	 						if(has5pPeak){
	 							if(has3pPeak){
	 								numFullSing++;
	 								bwComplete.write(gene.toBED()+"\n");
	 								logger.info(gene.getName()+" is complete");
	 								complete.get(chr).add(gene);
	 							}
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
	 							if(has3pPeak){
	 								numFullMult++;
	 								bwComplete.write(gene.toBED()+"\n");
	 								logger.info(gene.getName()+" is complete");
	 								complete.get(chr).add(gene);
	 							}
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
 						
 						geneTo5pPeakMap.put(gene, this5pPeaks);
 						geneTo3pPeakMap.put(gene, this3pPeaks);
					}					
				}
			logger.info("Single: Number of Full= "+numFullSing+" Incomplete= "+numPartialSing+
						" 5p No 3p= "+num5pPartialSing+" 3p No 5p= "+num3pPartialSing);
			logger.info("Multiple: Number of Full= "+numFullMult+" Incomplete= "+numPartialMult+
						" 5p No 3p= "+num5pPartialMult+" 3p No 5p= "+num3pPartialMult);
		}
		bw5p.close();
		bw3p.close();
		bw5pBed.close();
		bw3pBed.close();
		bwComplete.close();
		
		Map<Gene,List<Double>> mapp = trimAndExtendBestIsoform(geneTo5pPeakMap,geneTo3pPeakMap,outputName);
		trimAndExtendAllIsoforms(geneTo5pPeakMap,geneTo3pPeakMap,outputName,mapp);
	}
	
	/**
	 * Trims or extends transcripts at 5' and/or 3' end to the best peak.
	 * @param geneTo5pPeakMap
	 * @param geneTo3pPeakMap
	 * @param outputName
	 * @return returns a map of each gene to a list of double values, where list[0] = max z score in 5' list[1] = max z score in 3'
	 * @throws IOException
	 */
	private Map<Gene,List<Double>> trimAndExtendBestIsoform(Map<Gene,List<Annotation>> geneTo5pPeakMap,Map<Gene,List<Annotation>> geneTo3pPeakMap,String outputName) throws IOException{
		
		Map<Gene,List<Double>> geneToMaxZScoreMap = new TreeMap<Gene,List<Double>>();
		
		BufferedWriter bwBed = new BufferedWriter(new FileWriter(outputName+".trimmed.best.bed"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".trimmed.summary.txt"));
		
		int trim5p =0;
		int extend5p = 0;
		int trim3p = 0;
		int extend3p=0;
		//ONLY TRIM COMPLETE TRANSCRIPTS
		//annotations = complete;
		
		//For each chromosome in the annotation set
		for(String chr:annotations.keySet()){
			//For each gene
			for(Gene gene:annotations.get(chr)){
				Gene edited = gene.copy();
				
				List<Double> list = new ArrayList<Double>(2);
				//logger.info(edited.toBED());
				//if there is a 5' peak for this gene
				if(geneTo5pPeakMap.containsKey(gene)){
					//Get the max z-score peak
					Annotation maxPeak = getMaxZScorePeak(geneTo5pPeakMap.get(gene));
					if(maxPeak !=null){
						list.add(0, maxPeak.getScore());
						if(gene.isNegativeStrand()){ 
							if(maxPeak.getEnd()>gene.getEnd()){
								extend5p++;
								edited.setEnd(maxPeak.getEnd());
								edited.setName(edited.getName()+".adj5p");
							}
							else{
								if(maxPeak.getEnd()<gene.getEnd()){
									if(maxPeak.getEnd()>gene.getStart()){
										trim5p++;
										edited.setEnd(maxPeak.getEnd());
										edited.setName(edited.getName()+".adj5p");
									}
									else{
										logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
									}
								}
							}
						}
						else{
							if(maxPeak.getStart()<gene.getStart()){
								extend5p++;
								edited.setStart(maxPeak.getStart());
								edited.setName(edited.getName()+".adj5p");
							}
							else{
								if(maxPeak.getStart()>gene.getStart()){
									if(maxPeak.getStart()<gene.getEnd()){
										trim5p++; 
										edited.setStart(maxPeak.getStart());
										edited.setName(edited.getName()+".adj5p");
									}
									else{
										logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
									}
								}
							}
						}
					}
					else{
						list.add(0, Double.NaN);
					}
				}
				else{
					list.add(0, Double.NaN);
				}
				//if there is a 3' peak for this gene
				if(geneTo3pPeakMap.containsKey(gene)){
					//Get the max z-score peak
					Annotation maxPeak = getMaxZScorePeak(geneTo3pPeakMap.get(gene));
					if(maxPeak !=null){
						list.add(1, maxPeak.getScore());
						if(gene.isNegativeStrand()){ 
							if(maxPeak.getStart()<edited.getStart()){
								extend3p++;
								edited.setStart(maxPeak.getStart());
								edited.setName(edited.getName()+".adj3p");
							}
							else{
								if(maxPeak.getStart()>edited.getStart()){
									if(maxPeak.getStart()<edited.getEnd()){
										trim3p++;
										edited.setStart(maxPeak.getStart());
										edited.setName(edited.getName()+".adj3p");
									}
									else{
										logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
									}
								}
							}
						}
						else{
							if(maxPeak.getEnd()>edited.getEnd()){
								extend3p++;
								edited.setEnd(maxPeak.getEnd());
								edited.setName(edited.getName()+".adj3p");
							}
							else{
								if(maxPeak.getEnd()<edited.getEnd()){
									if(maxPeak.getEnd()>edited.getStart()){
										trim3p++;
										edited.setEnd(maxPeak.getEnd());
										edited.setName(edited.getName()+".adj3p");
									}
									else{
										logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
									}
								}
							}
							
						}
					}
					else{
						list.add(1, Double.NaN);
					}
					//edited.setOrientedEnd(maxPeak.getOrientedEnd());
				}
				else{
					list.add(1, Double.NaN);
				}
				geneToMaxZScoreMap.put(gene, list);
				//logger.info(edited.toBED());
				bwBed.write(edited.toBED()+"\n");
			}
		}
		bw.write("Genes with trimmed 5p ends: "+trim5p+"\n");
		bw.write("Genes with extended 5p ends: "+extend5p+"\n");
		bw.write("Genes with trimmed 3p ends: "+trim3p+"\n");
		bw.write("Genes with extended 3p ends: "+extend3p+"\n");
		bwBed.close();
		bw.close();
		
		return geneToMaxZScoreMap;
	}
	
	/**
	 * Removes all peaks with z-score less than 10% of the max z-score
	 * @param peaks
	 * @param prevZscore
	 */
	private List<Annotation> cleanPeaks(List<Annotation> peaks,double prevZscore){
		List<Annotation> rtrn = new ArrayList<Annotation>();
		for(Annotation p:peaks){
			if(p.getScore()>0.1*prevZscore){
				rtrn.add(p);
			}
		}
		return rtrn;
	}
	
	private Annotation pop(List<Annotation> peaks,double prevZscore){
		
		double max = Double.MIN_VALUE;
		Annotation rtrn = null;
		//Find the peak with the next highest score in peaks
		for(Annotation p:peaks){
			if(p.getScore()>max){
				rtrn = p;
				max = p.getScore();
			}
		}
		return rtrn;
	}
	/**
	 * Trims or extends transcripts at 5' and/or 3' end to the best peak.
	 * @param geneTo5pPeakMap
	 * @param geneTo3pPeakMap
	 * @param outputName
	 * @throws IOException
	 */
	private void trimAndExtendAllIsoforms(Map<Gene,List<Annotation>> geneTo5pPeakMap,Map<Gene,List<Annotation>> geneTo3pPeakMap,String outputName,Map<Gene,List<Double>> mapp) throws IOException{
		
		BufferedWriter bwBed = new BufferedWriter(new FileWriter(outputName+".trimmed.all.bed"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".trimmed.all.summary.txt"));
		
		int trim5p =0;
		int extend5p = 0;
		int trim3p = 0;
		int extend3p=0;
		//ONLY TRIM COMPLETE TRANSCRIPTS
		//annotations = complete;
		//For each chromosome in the annotation set
		for(String chr:annotations.keySet()){
			//For each gene
			for(Gene gene:annotations.get(chr)){
				
				boolean geneTrimmed = false;
				boolean peakTrimmed = false;
				
				Collection<Gene> editedGenes = new TreeSet<Gene>();
				//if there is a 5' peak for this gene
				if(geneTo5pPeakMap.containsKey(gene)){
					
					double prevZScore = mapp.get(gene).get(0);
					List<Annotation> peaks = geneTo5pPeakMap.get(gene);
					//Removes all peaks with z-score less than 10% of the max z-score
					peaks = cleanPeaks(peaks,prevZScore);
					//While peaks is not empty
					while(!peaks.isEmpty()){
						peakTrimmed = false;
						//Get the next highest peak
						Annotation peak5 = pop(peaks,prevZScore);
						peaks.remove(peak5);
						//A new gene copy
						Gene edited = gene.copy();
						if(gene.isNegativeStrand()){ 
							if(peak5.getEnd()>gene.getEnd()){
								extend5p++;
								edited.setEnd(peak5.getEnd());
								geneTrimmed=true;
								peakTrimmed=true;
								edited.setName(edited.getName()+".adj5p");
							}
							else{
								if(peak5.getEnd()<gene.getEnd()){
									if(peak5.getEnd()>gene.getStart()){
										trim5p++;
										edited.setEnd(peak5.getEnd());
										geneTrimmed=true;
										peakTrimmed=true;
										edited.setName(edited.getName()+".adj5p");
									}
									else{
										logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
									}
								}
							}
						}
						else{
							if(peak5.getStart()<gene.getStart()){
								extend5p++;
								edited.setStart(peak5.getStart());
								geneTrimmed=true;
								peakTrimmed=true;
								edited.setName(edited.getName()+".adj5p");
							}
							else{
								if(peak5.getStart()>gene.getStart()){
									if(peak5.getStart()<gene.getEnd()){
										trim5p++; 
										edited.setStart(peak5.getStart());
										geneTrimmed=true;
										peakTrimmed=true;
										edited.setName(edited.getName()+".adj5p");
									}
									else{
										logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
									}
								}
							}
						}
						if(peakTrimmed){
							editedGenes.add(edited);
						}
					}

					//If none of the peaks were used for trimming
					if(!geneTrimmed)
						editedGenes.add(gene.copy());

				}
				//if there is a 3' peak for this gene
				if(geneTo3pPeakMap.containsKey(gene)){
					
					geneTrimmed = false;
					peakTrimmed = false;
					
					double prevZScore = mapp.get(gene).get(1);
					List<Annotation> peaks = geneTo3pPeakMap.get(gene);
					//Removes all peaks with z-score less than 10% of the max z-score
					peaks = cleanPeaks(peaks,prevZScore);
					//While peaks is not empty
					while(!peaks.isEmpty()){
						peakTrimmed = false;
						//Get the next highest peak
						Annotation peak3 = pop(peaks,prevZScore);
						peaks.remove(peak3);
						
						//For each edited gene
						for(Gene edited:editedGenes){
							
							Gene editedAgain = edited.copy();
							if(gene.isNegativeStrand()){ 
								if(peak3.getStart()<edited.getStart() && (edited.getEnd()-peak3.getStart())>200){
									extend3p++;
									editedAgain.setStart(peak3.getStart());
									editedAgain.setName(edited.getName()+".adj3p");
								}
								else{
									if(peak3.getStart()>edited.getStart()){
										if(peak3.getStart()<edited.getEnd() && (edited.getEnd()-peak3.getStart())>200){
											trim3p++;
											editedAgain.setStart(peak3.getStart());
											editedAgain.setName(edited.getName()+".adj3p");
										}
										else{
											logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
										}
									}
								}
							}
							else{
								if(peak3.getEnd()>edited.getEnd() && (peak3.getEnd()-edited.getStart())>200){
									extend3p++;
									editedAgain.setEnd(peak3.getEnd());
									editedAgain.setName(edited.getName()+".adj3p");
								}
								else{
									if(peak3.getEnd()<edited.getEnd()){
										if(peak3.getEnd()>edited.getStart() && (peak3.getEnd()-edited.getStart())>200){
											trim3p++;
											editedAgain.setEnd(peak3.getEnd());
											editedAgain.setName(edited.getName()+".adj3p");
										}
										else{
											logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
										}
									}
								}							
							}
							bwBed.write(editedAgain.toBED()+"\n");
						}

					}
				}
				//logger.info(edited.toBED());
			}
		}
		bw.write("Genes with trimmed 5p ends: "+trim5p+"\n");
		bw.write("Genes with extended 5p ends: "+extend5p+"\n");
		bw.write("Genes with trimmed 3p ends: "+trim3p+"\n");
		bw.write("Genes with extended 3p ends: "+extend3p+"\n");
		bwBed.close();
		bw.close();
	}
	
	/**
	 * Returns the peak with the maximum z-score
	 * @param peaks
	 * @return
	 */
	private Annotation getMaxZScorePeak(List<Annotation> peaks){
		
		if(peaks.size()==1){
			return peaks.get(0);
		}
		double max = Double.MIN_VALUE;
		Annotation maxA=null;
		for(Annotation p:peaks){
			if(p.getScore()>max){
				max = p.getScore();
				maxA = p;
			}
		}
		return maxA;
	}
	
	private double get5pWindowCount(Annotation window,Strand orientation){
		double windowCount = 0.0;
		//Get the reads in the window
		window.setOrientation(orientation);
		//Iterator<AlignmentCount> readiter = model5p.getOverlappingReadCountsStranded(window, false);
		CloseableIterator<Alignment> readiter = model5p.getOverlappingReads(window, false);
		//for all reads in the window
		while(readiter.hasNext()){
			//AlignmentCount read = readiter.next();
			Alignment read = readiter.next();
			if(passesStartChecks(read,window,orientation)){
				windowCount += read.getWeight();
			}
		}
		readiter.close();
		return windowCount;
	}
	
	private double get3pWindowCount(Annotation window,Strand orientation){
		double windowCount = 0.0;
		if(singleEnd3p){
			if(orientation.equals(Strand.POSITIVE))
				window.setOrientation(Strand.NEGATIVE);
			else
				if(orientation.equals(Strand.NEGATIVE))
					window.setOrientation(Strand.POSITIVE);
		}
		else{
			//Get the reads in the window
			window.setOrientation(orientation);
		}
		//Iterator<AlignmentCount> readiter = model3p.getOverlappingReadCountsStranded(window, false);
		CloseableIterator<Alignment> readiter = model3p.getOverlappingReads(window, false);
		//for all reads in the window
		while(readiter.hasNext()){
			//AlignmentCount read = readiter.next();
			Alignment read = readiter.next();
			if(passesEndChecks(read,window,orientation)){
				windowCount += read.getWeight();
			}
		}
		readiter.close();
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
	 * 			[2]: highest pileup
	 */
	private double[] get5pNullDistribution(Gene annotation,CoordinateSpace space){
		
		List<Double> values = new ArrayList<Double>();
		double max = Double.MIN_VALUE;
		Iterator<? extends Window> giter = space.getWindowIterator(annotation, windowSize, 0);
		//For each window
		while(giter.hasNext()){
			Window window = giter.next();
			double windowCount = 0.0;
			//Get the reads in the window
			//For each block in the window
			for(Annotation block: window.getBlocks()){
				
				CloseableIterator<Alignment> readiter = model5p.getOverlappingReads(block,false);
				//for all reads in the window
				while(readiter.hasNext()){
					Alignment read = readiter.next();					
					if(passesStartChecks(read,window,annotation.getOrientation())){
						windowCount += read.getWeight();
					}				
				}
				readiter.close();
			}
//			if(windowCount>0.0){
				values.add(windowCount);
				if(windowCount>max){
					max = windowCount;
				}
//			}
		}
		double[] rtrn = new double[3];
		if(values.size()>0){
			rtrn[0] = Statistics.mean(values);
			//rtrn[0] = Statistics.median(values);
			//rtrn[0] = Statistics.geometricMean(l2a(values));
			rtrn[1] = Statistics.variance(values);
			rtrn[2] = max;
		}
		else{ 
			rtrn[0] = 0.0;
			rtrn[1] = 0.0;
			rtrn[2] = 0.0;
		}
		return rtrn;
	}
	
	/**
	 * CALCULATE THE NULL OVER THE ISOFORM
	 * @param annotation
	 * @return 	[0]: mean
	 * 			[1]: variance
	 * 			[2]: highest pileup
	 */
	private double[] get3pNullDistribution(Gene annotation,CoordinateSpace space){
		
		List<Double> values = new ArrayList<Double>();
		double max = Double.MIN_VALUE;
		Iterator<? extends Window> giter = space.getWindowIterator(annotation, windowSize, 0);
		//For each window
		while(giter.hasNext()){
			Window window = giter.next();
			double windowCount = 0.0;
			//Get the reads in the window
			//For each block in the window
			for(Annotation block: window.getBlocks()){
				
				if(singleEnd3p){
					Strand s = window.getOrientation();
					if(s.equals(Strand.POSITIVE))
						block.setOrientation(Strand.NEGATIVE);
					else
						if(s.equals(Strand.NEGATIVE))
							block.setOrientation(Strand.POSITIVE);
				}
				CloseableIterator<Alignment> readiter = model3p.getOverlappingReads(block,false);
				//for all reads in the window
				while(readiter.hasNext()){
					Alignment read = readiter.next();
					//Orientation for the 3p models will be set already on opposite mate of 5p models OR if SINGLE END, orientation NOT changed.
					if(passesEndChecks(read,window,annotation.getOrientation())){
						windowCount += read.getWeight();
					}				
				}
				readiter.close();
			}
//			if(windowCount>0.0){
				values.add(windowCount);
				if(windowCount>max){
					max = windowCount;
				}
//			}
		}
		double[] rtrn = new double[3];
		if(values.size()>0){
			//rtrn[0] = Statistics.mean(values);
			//rtrn[0] = Statistics.median(values);
			rtrn[0] = Statistics.geometricMean(l2a(values));
			rtrn[1] = Statistics.variance(values);
			rtrn[2] = max;
		}
		else{ 
			rtrn[0] = 0.0;
			rtrn[1] = 0.0;
			rtrn[2] = 0.0;
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
	private boolean passesStartChecks(Alignment read,Annotation window,Strand orientation){
		
		if(SingleEndAlignment.class.isInstance(read)){
			SingleEndAlignment align = (SingleEndAlignment) read;
			//Check if read is the correct read
			//if read starts in window
			if(((strand5p==(TranscriptionRead.FIRST_OF_PAIR) && align.getIsFirstMate()) || 
					(strand5p==(TranscriptionRead.SECOND_OF_PAIR) && !align.getIsFirstMate())) 
						&& (readStartFallsInWindow(read,window))
							&& (read.getOrientation().equals(orientation))){
				return true;
			}
		}
		//ELSE PAIRED
		else{
			PairedReadAlignment align = (PairedReadAlignment) read;
			Annotation mate;
			if(strand5p==TranscriptionRead.FIRST_OF_PAIR){
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
	 * Returns true is the specified alignment starts in window and 
	 * same mate as is in direction of transcription for model 5p is in opp direction of transcription
	 * thus,
	 * opposite orientation as transcript but other checks are same
	 * @param read
	 * @param window
	 * @return
	 */
	private boolean passesEndChecks(Alignment read,Annotation window,Strand orientation){
		
		if(SingleEndAlignment.class.isInstance(read)){
			SingleEndAlignment align = (SingleEndAlignment) read;			
			//if data is single end
			if(singleEnd3p){
				if((readStartFallsInWindow(read,window))
						//Single end alignment. Orientation is already corrected.
						&& (!read.getOrientation().equals(orientation)))
					return true;
			}
			//Single end but data is paired. Thus single mates have same orientation as the transcript
			else{
				//Check if read is the correct read
				//if read starts in window
				if(((strand3p==(TranscriptionRead.FIRST_OF_PAIR) && !align.getIsFirstMate()) || 
						(strand3p==(TranscriptionRead.SECOND_OF_PAIR) && align.getIsFirstMate())) 
						//We used readEndFallsInWindow if the read would be in the same orientation as the gene
						//Here is it oppostite so
						//	&& (readEndFallsInWindow(read,window))
							&& (readEndFallsInWindow(read,window))
							//Single end alignment. Orientation is already corrected.
								&& (read.getOrientation().equals(orientation))){
					return true;
				}
			}
		}
		//ELSE PAIRED
		else{
			PairedReadAlignment align = (PairedReadAlignment) read;
			Annotation mate;
			if(strand3p==TranscriptionRead.SECOND_OF_PAIR){
				mate = align.getFirstMate();
			}
			else{
				mate = align.getSecondMate();
			}
			if(readStartFallsInWindow(mate,window) && (!read.getOrientation().equals(orientation))){
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
	private boolean readStartFallsInWindow(Annotation align,Annotation window){
		
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
	
	/**
	 * Returns true if the oriented end of the read falls in the window
	 * @param align
	 * @param window
	 * @return
	 */
	private boolean readEndFallsInWindow(Annotation align,Annotation window){
		
		int end;
		if(align.isNegativeStrand()){
			end = align.getStart();
		}
		else{
			end = align.getEnd();
		}
		if(end>=window.getStart() && end<=window.getEnd())
			return true;
		else
			return false;
	}
	
	/**
	 * Calculates the distance of the Annotation from the 5' end of the gene
	 * @param gene
	 * @param p
	 * @return
	 */
	private int calculate5pDistance(Gene gene,Annotation p){
		
		if(gene.isNegativeStrand()){
			return p.getEnd()-gene.getEnd();
		}
		else{
			return gene.getStart()-p.getStart();
		}
	}

	/**
	 * Calculates the distance of the Annotation from the 3' end of the gene
	 * @param gene
	 * @param p
	 * @return
	 */
	private int calculate3pDistance(Gene gene,Annotation p){
		
		if(!gene.isNegativeStrand()){
			return p.getEnd()-gene.getEnd();
		}
		else{
			return gene.getStart()-p.getStart();
		}
	}
	
	/**
	 * Returns a map of chromosome name to interval tree for that chromosome
	 * @param annotations
	 * @return
	 */
	public void initiateIntervalTrees(Map<String,Collection<Gene>> annotations){
		intervalTrees = new HashMap<String,IntervalTree<Gene>>();
		for(String chr:annotations.keySet()){
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			intervalTrees.put(chr, tree);
		}
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
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"completeTranscripts");
		
		new AddEndRNASeqToScripture(argMap);

	}
	
	public static class IsoformMap{
		
		//Map of chromosome to genes
		Map<String,Collection<Gene>> chrToGeneMap;
		//Map of gene to isoforms
		Map<Gene,Set<Gene>> geneToIsoformMap;
		
		public IsoformMap(){
			chrToGeneMap = new HashMap<String,Collection<Gene>>();
			geneToIsoformMap = new HashMap<Gene,Set<Gene>>();
		}
		
		/**
		 * ADDS MAPS FOR A CHROMOSOME
		 * If maps already exist, they will be replaced
		 * @param chr
		 * @param map
		 */
		public void addChrMap(String chr,Map<Gene,Set<Gene>> map){
			chrToGeneMap.put(chr, map.keySet());
			geneToIsoformMap.putAll(map);
		}
		
		public Collection<String> getChromosomes(){
			return chrToGeneMap.keySet();
		}

		public Collection<Gene> getAllGenes(){
			return geneToIsoformMap.keySet();
		}
		
		public Collection<Gene> getGenesForChromosome(String chr){
			return chrToGeneMap.get(chr);
		}
		
		public Collection<Gene> getIsoformsForGene(Gene gene){
			return geneToIsoformMap.get(gene);
		}
		
		public int getNumOfIsoformsForGene(Gene gene){
			return geneToIsoformMap.get(gene).size();
		}
	}
	
	public int getDistanceToClosestSameOrientation5pGene(Gene gene) {
		IntervalTree<Gene> chrTree = intervalTrees.get(gene.getChr());
		int closestUpstream = Integer.MAX_VALUE;
		if(chrTree != null) {
			if(gene.isNegativeStrand()){
				Node<Gene> max  = chrTree.max(gene.getEnd()+1, gene.getEnd()+2 );
				if(max != null&& !max.getValue().overlaps(gene) && max.getValue().getOrientation()==gene.getOrientation()) {
					//logger.info("Gene = "+gene.toUCSC()+" max is "+max.getValue().toUCSC());
					closestUpstream = max.getValue().getStart() - gene.getEnd();
				}
			}
			else{
				Node<Gene> max  = chrTree.max(gene.getStart()-1, gene.getStart() );
				if(max != null && !max.getValue().overlaps(gene) && max.getValue().getOrientation()==gene.getOrientation()) {
					closestUpstream = gene.getStart() - max.getValue().getEnd();
				}
			}
			
		}
		return closestUpstream;
		
	}
	
	/**
	 * Returns the distance to the next closest gene downstream of the current gene's 3' end
	 * @param gene
	 * @return
	 */
	public int getDistanceToClosestSameOrientation3pGene(Gene gene) {
		IntervalTree<Gene> chrTree = intervalTrees.get(gene.getChr());
		int closestUpstream = Integer.MAX_VALUE;
		if(chrTree != null) {
			if(!gene.isNegativeStrand()){
				Node<Gene> max  = chrTree.max(gene.getEnd()+1, gene.getEnd()+2 );
				if(max != null && !max.getValue().overlaps(gene) && max.getValue().getOrientation()==gene.getOrientation()) {
					closestUpstream = max.getValue().getStart() - gene.getEnd();
					//logger.info(gene.toBED());
					//logger.info(max.getValue().toBED());
				}
			}
			else{
				Node<Gene> max  = chrTree.max(gene.getStart()-1, gene.getStart() );
				if(max != null && !max.getValue().overlaps(gene) && max.getValue().getOrientation()==gene.getOrientation()) {
					closestUpstream = gene.getStart() - max.getValue().getEnd();
					//logger.info(gene.toBED());
					//logger.info(max.getValue().toBED());
				}
			}
			
		}
		return closestUpstream;
		
	}
	
	private double[] l2a(List<Double> list){
		double[] rtrn=new double[list.size()];
	
		int i=0;
		for(Double val: list){rtrn[i++]=val;}
	
		return rtrn;
	}
	
	private class NullDistribution{
		
		double[] distribution;
		int highestPileup;
		
		NullDistribution(double[] d,int h){
			distribution = d;
			highestPileup = h;
		}
		
		public double[] getDistribution(){
			return distribution;
		}
		
		public int getHighestPileup(){
			return highestPileup;
		}
	}
	
	
	/**
	 * Calculates number of complete, incomplete and partial transcripts based on end RNA-Seq data
	 * @param outputName
	 * @throws IOException
	 */
	public void trim3pEnds(String outputName) throws IOException{
		
		BufferedWriter bw3p = new BufferedWriter(new FileWriter(outputName+".peaks.3p"));
		BufferedWriter bw3pBed = new BufferedWriter(new FileWriter(outputName+".peaks.3p.bed"));
		
		Map<Gene,List<Annotation>> geneTo3pPeakMap = new HashMap<Gene,List<Annotation>>();
		int count = 0;
		complete = new TreeMap<String,Collection<Gene>>();

		//For each chromosome in the annotation set
		for(String chr:annotations.keySet()){
			
			int num3pPartialSing = 0;
			int num3pPartialMult = 0;

			logger.info("Processing "+chr);
			complete.put(chr, new TreeSet<Gene>());
			//If 5' or 3' end RNA-seq does not have data for it, dont run
			if(!model3p.containsReference(chr)){
				logger.warn("Warning: "+chr +" not expressed in 3' data.");
			}
			else{				
				/*
				 * FOR THIS GENE
				 */
				for(Gene gene:annotations.get(chr)){
					//logger.info("Processing "+gene.getName());
					//Make a coordinate space with the gene
					Map<String,Collection<Gene>> chrToGenesMap = new HashMap<String,Collection<Gene>>();
					List<Gene> g = new ArrayList<Gene>();
					g.add(gene);
					chrToGenesMap.put(gene.getChr(), g);
					CoordinateSpace space = new TranscriptomeSpace(chrToGenesMap);
					
					count++;
					/*
					 * IS THERE AT LEAST 1 3P END
					 */
					boolean has3pPeak = false;
										
					int start = 0;
					int end =0;
					int bestExtension = getDistanceToClosestSameOrientation3pGene(gene);
					//logger.info(gene.getName()+" getDistanceToClosestSameOrientation3pGene "+bestExtension);
					if(bestExtension>extension || bestExtension<0){
						bestExtension = extension;
					}
					if(gene.isNegativeStrand()){
						start = bestExtension;
						//
					}
					else{
						end = bestExtension;
					}
					//logger.info("Start: "+start+" End: "+end);
					Gene gs = gene.copy();
					//EXPAND IS A STRAND-INDEPENDENT FUNCTION
					gs.expand(start, end);
					chrToGenesMap = new HashMap<String,Collection<Gene>>();
					g = new ArrayList<Gene>();
					g.add(gs);
					chrToGenesMap.put(gene.getChr(), g);
					space = new TranscriptomeSpace(chrToGenesMap);
					
					double[] nulls3p = get3pNullDistribution(gs,space);
						
					/**
					 * 3P: For pileup less than 5bp, dont run
					 */
					if(nulls3p[2]<=5){
						
					}
					else{ if(nulls3p[2]<=10){
						THRESHOLD=15;
					}
					else{
						THRESHOLD = DEFAULT_THRESHOLD;
//					}		
 					
					chrToGenesMap = new HashMap<String,Collection<Gene>>();
					g = new ArrayList<Gene>();
					g.add(gs);
					chrToGenesMap.put(gene.getChr(), g);
					space = new TranscriptomeSpace(chrToGenesMap);
					
					Iterator<? extends Window> iter = space.getWindowIterator(gs, windowSize, 0);
					boolean flag3p = false;
					Annotation prev3p = null;
					Annotation peak3p =null;
					
					List<Annotation> this3pPeaks = new ArrayList<Annotation>();
					//For every window in the transcript
					while(iter.hasNext()){
						Window window = iter.next();							
						/*
						 * 3p
						 */
						double windowCount3p = get3pWindowCount(window,gene.getOrientation());
//						if(windowCount3p>0)
//							logger.info(window.toUCSC()+" "+windowCount3p);
						//Get the z-score of each window
//						double zscore = getMinimumZScore(windowCount,nulls,window.getSize());
						double zscore3p = Statistics.zScore(windowCount3p, nulls3p[0],nulls3p[1],window.getSize());
//						if(windowCount3p>0)
//							logger.info(window.toUCSC()+" Count = "+windowCount3p+" zscore = "+zscore3p);
						//If window is significant
						if(zscore3p>=THRESHOLD){
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
							//this3pPeaks.add(window);
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
					if(flag3p){
						this3pPeaks.add(peak3p);
					}
					for(Annotation p:this3pPeaks){
						double windowCount = get3pWindowCount(p,gene.getOrientation());
						double zscore = Statistics.zScore(windowCount, nulls3p[0],nulls3p[1],p.getSize());
						p.setScore(zscore);							
						bw3pBed.write(p.toBED()+"\n");
						bw3p.write(gene.getName()+"\t"+p.toUCSC()+"\t"+windowCount+"\t"+zscore+"\t"+calculate3pDistance(gene,p)+"\n");
					}
				
					if(gene.getBlocks().size()==1){
 						if(has3pPeak){
 							num3pPartialSing++;
 						}
					}
					else{
 						if(has3pPeak){
 							num3pPartialMult++;
 						}
					}
					geneTo3pPeakMap.put(gene, this3pPeaks);
				}}
					if(count%1000.0==0.0){
						logger.info("Single: 3p = "+num3pPartialSing);
						logger.info("Multiple: 3p = "+num3pPartialMult);
					}
					if(count%1000.0==0.0){
						logger.info("Processed "+count+" genes");
					}
					
				}					
			}
			logger.info("Single: Number of 3p = "+num3pPartialSing);
			logger.info("Multiple: Number of 3p = "+num3pPartialMult);
		}
		bw3p.close();
		bw3pBed.close();
		
		Map<Gene,List<Double>> mapp = trimAndExtendBestIsoform(geneTo3pPeakMap,outputName,"3p");
		trimAndExtendAllIsoforms(geneTo3pPeakMap,outputName,mapp,"3p");
	}
	
	/**
	 * Trims the 5' ends of transcripts based on end RNA-Seq data
	 * @param outputName
	 * @throws IOException
	 */
	public void trim5pEnds(String outputName) throws IOException{
		
		BufferedWriter bw5p = new BufferedWriter(new FileWriter(outputName+".peaks.5p"));
		BufferedWriter bw5pBed = new BufferedWriter(new FileWriter(outputName+".peaks.5p.bed"));
		
		Map<Gene,List<Annotation>> geneTo5pPeakMap = new HashMap<Gene,List<Annotation>>();
		int count = 0;
		complete = new TreeMap<String,Collection<Gene>>();
		//For each chromosome in the annotation set
		for(String chr:annotations.keySet()){
			
			int num5pPartialSing = 0;
			int num5pPartialMult = 0;
			logger.info("Processing "+chr);
			complete.put(chr, new TreeSet<Gene>());
			//If 5' or 3' end RNA-seq does not have data for it, dont run
			if(!model5p.containsReference(chr)){
				logger.warn("Warning:"+chr +" not present in 5' end RNA-seq");
			}
			else{				
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
					
					int start = 0;
					int end =0;
					int bestExtension = getDistanceToClosestSameOrientation5pGene(gene);
					//logger.info(gene.getName()+" getDistanceToClosestSameOrientation5pGene "+bestExtension);
					if(bestExtension>extension || bestExtension<0){
						bestExtension = extension;
					}
					if(gene.isNegativeStrand()){
						end = bestExtension;
					}
					else{
						start = bestExtension;
					}
					Gene ge = gene.copy();
					//EXPAND IS A STRAND-INDEPENDENT FUNCTION
					ge.expand(start, end);
					//logger.info(ge.toBED());
					chrToGenesMap = new HashMap<String,Collection<Gene>>();
					g = new ArrayList<Gene>();
					g.add(ge);
					chrToGenesMap.put(gene.getChr(), g);
					space = new TranscriptomeSpace(chrToGenesMap);
					
					double[] nulls5p = get5pNullDistribution(ge,space);
					
					if(nulls5p[2]<11){
						THRESHOLD=15;						
					}
					else{
						THRESHOLD = DEFAULT_THRESHOLD;
					}
						/*
						 * RE-DEFINE COORDINATE SPACE
						 */
						chrToGenesMap = new HashMap<String,Collection<Gene>>();
						g = new ArrayList<Gene>();
						g.add(ge);
						chrToGenesMap.put(gene.getChr(), g);
						space = new TranscriptomeSpace(chrToGenesMap);
						Iterator<? extends Window> giter = space.getWindowIterator(ge, windowSize, 0);
						
						boolean flag5p = false;
						Annotation prev5p = null;
						Annotation peak5p =null;
						
						List<Annotation> this5pPeaks = new ArrayList<Annotation>();
						//For every window in the transcript
 						while(giter.hasNext()){
							Window window = giter.next();
							
							/*
							 * 5P 
							 */
							double windowCount5p = get5pWindowCount(window,gene.getOrientation());
							double zscore5p = Statistics.zScore(windowCount5p, nulls5p[0],nulls5p[1],window.getSize());
							//Associate the high z-scores with gene
							//If window is significant
							if(zscore5p>=THRESHOLD){
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
								//this5pPeaks.add(window);
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
						}
 						//Last peak
 						if(flag5p){
							this5pPeaks.add(peak5p);
						}
 						for(Annotation p:this5pPeaks){
 							double windowCount = get5pWindowCount(p,gene.getOrientation());
							double zscore = Statistics.zScore(windowCount, nulls5p[0],nulls5p[1],p.getSize());
							p.setScore(zscore);
							bw5pBed.write(p.toBED()+"\n");
							bw5p.write(gene.getName()+"\t"+p.toUCSC()+"\t"+windowCount+"\t"+zscore+"\t"+calculate5pDistance(gene,p)+"\n");
						}
 						
  						if(gene.getBlocks().size()==1){
	 						if(has5pPeak){
	 								num5pPartialSing++;
	 						}
 						}
 						else{
 							if(has5pPeak){
	 							num5pPartialMult++;
	 						}
 						}
 						if(count%1000.0==0.0){
 							logger.info("Single: Number of 5p = "+num5pPartialSing);
 							logger.info("Multiple: Number of 5p = "+num5pPartialMult);
 						}
 						
 						geneTo5pPeakMap.put(gene, this5pPeaks);
					}					
				}
			logger.info("Single: Number of 5p = "+num5pPartialSing);
			logger.info("Multiple: Number of 5p = "+num5pPartialMult);
		}
		bw5p.close();
		bw5pBed.close();
		
		Map<Gene,List<Double>> mapp = trimAndExtendBestIsoform(geneTo5pPeakMap,outputName,"5p");
		trimAndExtendAllIsoforms(geneTo5pPeakMap,outputName,mapp,"5p");
	}

	/**
	 * Trims or extends transcripts at 5' and/or 3' end to the best peak.
	 * @param geneTo5pPeakMap
	 * @param geneTo3pPeakMap
	 * @param outputName
	 * @return returns a map of each gene to a list of double values, where list[0] = max z score in 5' list[1] = max z score in 3'
	 * @throws IOException
	 */
	private Map<Gene,List<Double>> trimAndExtendBestIsoform(Map<Gene,List<Annotation>> geneToPeakMap,String outputName,String whichEnd) throws IOException{
		
		Map<Gene,List<Double>> geneToMaxZScoreMap = new TreeMap<Gene,List<Double>>();
		
		BufferedWriter bwBed = new BufferedWriter(new FileWriter(outputName+".trimmed.best.bed"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".trimmed.summary.txt"));
		
		int count=0;
		if(whichEnd.equalsIgnoreCase("3p")){
			logger.info("Trimming Best isoform with 3p ends");
			int trim = 0;
			int extend=0;
			//For each chromosome in the annotation set
			for(String chr:annotations.keySet()){
				//For each gene
				for(Gene gene:annotations.get(chr)){
					
					count++;
					Gene edited = gene.copy();
					
					List<Double> list = new ArrayList<Double>(1);
					//logger.info(edited.toBED());
					//if there is a peak for this gene
					if(geneToPeakMap.containsKey(gene)){
						//Get the max z-score peak
						Annotation maxPeak = getMaxZScorePeak(geneToPeakMap.get(gene));
						if(maxPeak !=null){
							list.add(0, maxPeak.getScore());
							if(gene.isNegativeStrand()){ 
								if(maxPeak.getStart()<edited.getStart()){
									extend++;
									edited.setStart(maxPeak.getStart());
									edited.setName(edited.getName()+".adj3p");
								}
								else{
									if(maxPeak.getStart()>edited.getStart()){
										if(maxPeak.getStart()<edited.getEnd()){
											trim++;
											edited.setStart(maxPeak.getStart());
											edited.setName(edited.getName()+".adj3p");
										}
										else{
											logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
										}
									}
								}
							}
							else{
								if(maxPeak.getEnd()>edited.getEnd()){
									extend++;
									edited.setEnd(maxPeak.getEnd());
									edited.setName(edited.getName()+".adj3p");
								}
								else{
									if(maxPeak.getEnd()<edited.getEnd()){
										if(maxPeak.getEnd()>edited.getStart()){
											trim++;
											edited.setEnd(maxPeak.getEnd());
											edited.setName(edited.getName()+".adj3p");
										}
										else{
											logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
										}
									}
								}
								
							}
						}
						else{
							list.add(0, Double.NaN);
						}
						//edited.setOrientedEnd(maxPeak.getOrientedEnd());
					}
					else{
						list.add(0, Double.NaN);
					}
					geneToMaxZScoreMap.put(gene, list);
					//logger.info(edited.toBED());
					bwBed.write(edited.toBED()+"\n");
					
					if(count%500==0){
						logger.info("Trimming best isoform done for "+count+" genes");
					}
				}
			}
			bw.write("Genes with trimmed 3p ends: "+trim+"\n");
			bw.write("Genes with extended 3p ends: "+extend+"\n");
			bwBed.close();
			bw.close();
		}
		else{
			if(whichEnd.equalsIgnoreCase("5p")){
				logger.info("Trimming best isoforms with 5p ends");
				int trim = 0;
				int extend=0;
				//For each chromosome in the annotation set
				for(String chr:annotations.keySet()){
					//For each gene
					for(Gene gene:annotations.get(chr)){
						Gene edited = gene.copy();
						
						List<Double> list = new ArrayList<Double>(1);
						//logger.info(edited.toBED());
						//if there is a peak for this gene
						if(geneToPeakMap.containsKey(gene)){
							
							//Get the max z-score peak
							Annotation maxPeak = getMaxZScorePeak(geneToPeakMap.get(gene));
							if(maxPeak !=null){
								list.add(0, maxPeak.getScore());
								if(gene.isNegativeStrand()){ 
									if(maxPeak.getEnd()>gene.getEnd()){
										extend++;
										edited.setEnd(maxPeak.getEnd());
										edited.setName(edited.getName()+".adj5p");
									}
									else{
										if(maxPeak.getEnd()<gene.getEnd()){
											if(maxPeak.getEnd()>gene.getStart()){
												trim++;
												edited.setEnd(maxPeak.getEnd());
												edited.setName(edited.getName()+".adj5p");
											}
											else{
												logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
											}
										}
									}
								}
								else{
									if(maxPeak.getStart()<gene.getStart()){
										extend++;
										edited.setStart(maxPeak.getStart());
										edited.setName(edited.getName()+".adj5p");
									}
									else{
										if(maxPeak.getStart()>gene.getStart()){
											if(maxPeak.getStart()<gene.getEnd()){
												trim++; 
												edited.setStart(maxPeak.getStart());
												edited.setName(edited.getName()+".adj5p");
											}
											else{
												logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
											}
										}
									}
								}
							}
							else{
								list.add(0, Double.NaN);
							}
						}
						else{
							list.add(0, Double.NaN);
						}
						geneToMaxZScoreMap.put(gene, list);
						//logger.info(edited.toBED());
						bwBed.write(edited.toBED()+"\n");
					}
				}
				bw.write("Genes with trimmed 5p ends: "+trim+"\n");
				bw.write("Genes with extended 5p ends: "+extend+"\n");
				bwBed.close();
				bw.close();
			}
		}
		
		return geneToMaxZScoreMap;
	}
	
	/**
	 * Trims or extends transcripts at 5' and/or 3' end to the best peak.
	 * @param geneTo5pPeakMap
	 * @param geneTo3pPeakMap
	 * @param outputName
	 * @throws IOException
	 */
	private void trimAndExtendAllIsoforms(Map<Gene,List<Annotation>> geneToPeakMap,String outputName,Map<Gene,List<Double>> mapp,String whichEnd) throws IOException{
		
		BufferedWriter bwBed = new BufferedWriter(new FileWriter(outputName+".trimmed.all.bed"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".trimmed.all.summary.txt"));
		
		int count=0;
		if(whichEnd.equalsIgnoreCase("3p")){
			logger.info("Trimming all isoform with 3p ends");
			int trim =0;
			int extend = 0;
			//For each chromosome in the annotation set
			for(String chr:annotations.keySet()){
				//For each gene
				for(Gene gene:annotations.get(chr)){
					
					count++;
					boolean geneTrimmed = false;
					boolean peakTrimmed = false;
					Collection<Gene> editedGenes = new TreeSet<Gene>();
					
					//if there is a 3' peak for this gene
					if(geneToPeakMap.containsKey(gene)){
						
						geneTrimmed = false;
						peakTrimmed = false;
						
						double prevZScore = mapp.get(gene).get(0);
						List<Annotation> peaks = geneToPeakMap.get(gene);
						//Removes all peaks with z-score less than 10% of the max z-score
						peaks = cleanPeaks(peaks,prevZScore);
						//While peaks is not empty
						while(!peaks.isEmpty()){
							peakTrimmed = false;
							//Get the next highest peak
							Annotation peak3 = pop(peaks,prevZScore);
							peaks.remove(peak3);
								
							Gene edited = gene.copy();
							if(gene.isNegativeStrand()){ 
								if(peak3.getStart()<edited.getStart() && (edited.getEnd()-peak3.getStart())>200){
									extend++;
									edited.setStart(peak3.getStart());
									edited.setName(edited.getName()+".adj3p");
								}
								else{
									if(peak3.getStart()>edited.getStart()){
										if(peak3.getStart()<edited.getEnd() && (edited.getEnd()-peak3.getStart())>200){
											trim++;
											edited.setStart(peak3.getStart());
											edited.setName(edited.getName()+".adj3p");
										}
										else{
											logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
										}
									}
								}
							}
							else{
								if(peak3.getEnd()>edited.getEnd() && (peak3.getEnd()-edited.getStart())>200){
									extend++;
									edited.setEnd(peak3.getEnd());
									edited.setName(edited.getName()+".adj3p");
								}
								else{
									if(peak3.getEnd()<edited.getEnd()){
										if(peak3.getEnd()>edited.getStart() && (peak3.getEnd()-edited.getStart())>200){
											trim++;
											edited.setEnd(peak3.getEnd());
											edited.setName(edited.getName()+".adj3p");
										}
										else{
											logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
										}
									}
								}							
							}
							bwBed.write(edited.toBED()+"\n");
						}
					}
					//logger.info(edited.toBED());
					if(count%500==0){
						logger.info("Trimming all isoforms done for "+count+" genes");
					}
				}
			}
			bw.write("Genes with trimmed 3p ends: "+trim+"\n");
			bw.write("Genes with extended 3p ends: "+extend+"\n");
			bwBed.close();
			bw.close();

		}
		else if(whichEnd.equalsIgnoreCase("5p")){
			logger.info("Trimming all isoforms with 5p ends");
			int trim = 0;
			int extend=0;
			//ONLY TRIM COMPLETE TRANSCRIPTS
			//annotations = complete;
			//For each chromosome in the annotation set
			for(String chr:annotations.keySet()){
				//For each gene
				for(Gene gene:annotations.get(chr)){					
					//if there is a 5' peak for this gene
					if(geneToPeakMap.containsKey(gene)){
						
						double prevZScore = mapp.get(gene).get(0);
						List<Annotation> peaks = geneToPeakMap.get(gene);
						//Removes all peaks with z-score less than 10% of the max z-score
						peaks = cleanPeaks(peaks,prevZScore);
						//While peaks is not empty
						while(!peaks.isEmpty()){
							//Get the next highest peak
							Annotation peak5 = pop(peaks,prevZScore);
							peaks.remove(peak5);
							//A new gene copy
							Gene edited = gene.copy();
							if(gene.isNegativeStrand()){ 
								if(peak5.getEnd()>gene.getEnd()){
									extend++;
									edited.setEnd(peak5.getEnd());
									edited.setName(edited.getName()+".adj5p");
								}
								else{
									if(peak5.getEnd()<gene.getEnd()){
										if(peak5.getEnd()>gene.getStart()){
											trim++;
											edited.setEnd(peak5.getEnd());
											edited.setName(edited.getName()+".adj5p");
										}
										else{
											logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
										}
									}
								}
							}
							else{
								if(peak5.getStart()<gene.getStart()){
									extend++;
									edited.setStart(peak5.getStart());
									edited.setName(edited.getName()+".adj5p");
								}
								else{
									if(peak5.getStart()>gene.getStart()){
										if(peak5.getStart()<gene.getEnd()){
											trim++; 
											edited.setStart(peak5.getStart());
											edited.setName(edited.getName()+".adj5p");
										}
										else{
											logger.info("For "+gene.getName()+" the 5' end is downstream of 3' end.");
										}
									}
								}
							}
							bwBed.write(edited.toBED()+"\n");
						}
					}
				}
			}
			bw.write("Genes with trimmed 5p ends: "+trim+"\n");
			bw.write("Genes with extended 5p ends: "+extend+"\n");
			bwBed.close();
			bw.close();

		}				
	}
	
}

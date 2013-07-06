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
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.AlignmentModel.AlignmentCount;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readFilters.SplicedReadFilter;

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
	IsoformMap isoformMap;
	Map<String, IntervalTree<Gene>> intervalTrees;
	
	static final String usage = "Usage: AddEndRNASeqToScripture -task <task name> "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-3p <3P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\n\t\t-5p <5P Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\n\t\t-full <full length Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-annotations <Reconstruction bed file. [BED by default]> "+
			
			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\n\t\t-window <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends. Default = 5bp> "+
			"\n\t\t-extension <Specifies the size of the fixed window use to scan the genome. We recommend a size of 2-5bp for better resolution of gene ends> "+
			"\n\t\t-strand <Specifies the mate that is in the direction of trnascription> "+

			"\n";
	
	/**
	 * Instantiates the class and calls the function that finds complete transcripts and trims the reconstructions
	 * @param bamFile5p
	 * @param bamFile3p
	 * @param str
	 * @param annotationFile
	 * @param outputName
	 * @param windowS
	 * @param fullBam
	 * @param ext
	 * @throws IOException
	 */
	public AddEndRNASeqToScripture(File bamFile5p,File bamFile3p,TranscriptionRead str,String annotationFile,String outputName,int windowS,File fullBam,int ext) throws IOException{
		
		model5p=new AlignmentModel(bamFile5p.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,str,false); 
		
		TranscriptionRead oppStrand = TranscriptionRead.UNSTRANDED;
		if(str == (TranscriptionRead.FIRST_OF_PAIR))
			oppStrand = TranscriptionRead.SECOND_OF_PAIR;
		else 
			if(str == (TranscriptionRead.SECOND_OF_PAIR))
				oppStrand = TranscriptionRead.FIRST_OF_PAIR;
		if(oppStrand==TranscriptionRead.FIRST_OF_PAIR)
			logger.info("Opp is first");
		else
			if(oppStrand==TranscriptionRead.SECOND_OF_PAIR)
				logger.info("Opp is second");
		strand = str;
		model3p=new AlignmentModel(bamFile3p.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		//Read annotation file
		//annotationParser = new BEDFileParser(annotationFile);	
				
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		initiateIntervalTrees(annotations);
		
		windowSize = windowS;
		extension = ext;
		
		isoformMap = buildIsoformMap(annotations);
				
		//numberOfIsoformsPerGene(outputName,fullBam);
		findCompleteTranscripts(outputName);
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
		strand = str;
		//Read annotation file
		//annotationParser = new BEDFileParser(annotationFile);	
				
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		windowSize = windowS;
		extension = ext;
		
		isoformMap = buildIsoformMap(annotations);
				
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
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputName+".isoforms.info"));
		BufferedWriter bwCov = new BufferedWriter(new FileWriter(outputName+".isoforms.coverage.info"));
		model=new AlignmentModel(fullBam.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand);
		
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
					avg+=new ScanStatisticScore(model,isoform).getAverageCoverage(model);
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
		
		Map<Gene,List<Annotation>> geneTo5pPeakMap = new HashMap<Gene,List<Annotation>>();
		Map<Gene,List<Annotation>> geneTo3pPeakMap = new HashMap<Gene,List<Annotation>>();
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
			logger.info("Processing "+chr);
			//If 5' or 3' end RNA-seq does not have data for it, dont run
			if(model5p.getRefSequenceLambda(chr)==0.0 
					|| model3p.getRefSequenceLambda(chr)==0.0){
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
					if(bestExtension>extension){
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
					chrToGenesMap = new HashMap<String,Collection<Gene>>();
					g = new ArrayList<Gene>();
					g.add(ge);
					chrToGenesMap.put(gene.getChr(), g);
					space = new TranscriptomeSpace(chrToGenesMap);
					
					double[] nulls5p = get5pNullDistribution(ge,space);
					
					start = 0;
					end =0;
					bestExtension = getDistanceToClosestSameOrientation3pGene(gene);
					logger.info(gene.getName()+" getDistanceToClosestSameOrientation3pGene "+bestExtension);
					if(bestExtension>extension){
						bestExtension = extension;
					}
					if(gene.isNegativeStrand()){
						start = bestExtension;
						//logger.info("Start: "+start+" End: "+end);
					}
					else{
						end = bestExtension;
					}
					
					Gene gs = gene.copy();
					//EXPAND IS A STRAND-INDEPENDENT FUNCTION
					gs.expand(start, end);
					chrToGenesMap = new HashMap<String,Collection<Gene>>();
					g = new ArrayList<Gene>();
					g.add(gs);
					chrToGenesMap.put(gene.getChr(), g);
					space = new TranscriptomeSpace(chrToGenesMap);
					
					double[] nulls3p = get3pNullDistribution(gs,space);
					//logger.info("5p null : "+nulls5p[0]+" "+nulls5p[1]);
					//logger.info("3p null : "+nulls3p[0]+" "+nulls3p[1]);
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
//						Collection<Annotation> peaks = new ArrayList<Annotation>();
						//Compute the z-score and take the minimum of the z-scores
//						Iterator<? extends Window> witer = model5p.getCoordinateSpace().getWindowIterator(windowSize, chr, isoform.getStart(), isoform.getEnd(), 0);
						//Iterator<? extends Window> witer = model5p.getCoordinateSpace().getWindowIterator(windowSize, chr, gene.getStart(), gene.getEnd(), 0);
						
						/*
						 * ITERATE IN THE TRANSCRIPTOME SPACE ONLY
						 */
/*						int start = 0;
						int end =0;
						//logger.info("Gene "+gene.toBED());
						//logger.info("Start: "+start+" End: "+end);
						if(gene.isNegativeStrand()){
							end = extension;
						}
						else{
							start = extension;
						}
						//logger.info("Start: "+start+" End: "+end);
						Gene ge = gene.copy();
						//EXPAND IS A STRAND-INDEPENDENT FUNCTION
						ge.expand(start, end);*/
						//logger.info("After expansion Gene "+ge.toBED());
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
//							if(windowCount5p>0)
//								logger.info(window.toUCSC()+" "+windowCount5p);
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
/* 						start = 0;
						end =0;
						if(gene.isNegativeStrand()){
							start = extension;
							//logger.info("Start: "+start+" End: "+end);
						}
						else{
							end = extension;
						}
						//logger.info("3p: Start: "+start+" End: "+end);
						Gene gs = gene.copy();
						//EXPAND IS A STRAND-INDEPENDENT FUNCTION
						gs.expand(start, end);*/
						//logger.info("After 3p expansion "+gs.toBED());
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
							//logger.info("Count = "+windowCount3p+" zscore = "+zscore3p);
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
								edited.setName(edited.getName()+".trimmed3p");
							}
							else{
								if(maxPeak.getStart()>edited.getStart()){
									if(maxPeak.getStart()<edited.getEnd()){
										trim3p++;
										edited.setStart(maxPeak.getStart());
										edited.setName(edited.getName()+".trimmed3p");
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
								edited.setName(edited.getName()+".trimmed3p");
							}
							else{
								if(maxPeak.getEnd()<edited.getEnd()){
									if(maxPeak.getEnd()>edited.getStart()){
										trim3p++;
										edited.setEnd(maxPeak.getEnd());
										edited.setName(edited.getName()+".trimmed3p");
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
								extend5p++;
								edited.setEnd(peak5.getEnd());
								geneTrimmed=true;
								peakTrimmed=true;
								edited.setName(edited.getName()+".trimmed5p");
							}
							else{
								if(peak5.getEnd()<gene.getEnd()){
									if(peak5.getEnd()>gene.getStart()){
										trim5p++;
										edited.setEnd(peak5.getEnd());
										geneTrimmed=true;
										peakTrimmed=true;
										edited.setName(edited.getName()+".trimmed5p");
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
								edited.setName(edited.getName()+".trimmed5p");
							}
							else{
								if(peak5.getStart()>gene.getStart()){
									if(peak5.getStart()<gene.getEnd()){
										trim5p++; 
										edited.setStart(peak5.getStart());
										geneTrimmed=true;
										peakTrimmed=true;
										edited.setName(edited.getName()+".trimmed5p");
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
					
					//For each 3' peak
					for(Annotation peak3:geneTo3pPeakMap.get(gene)){
						
						//For each edited gene
						for(Gene edited:editedGenes){
							Gene editedAgain = edited.copy();
							if(gene.isNegativeStrand()){ 
								if(peak3.getStart()<edited.getStart()){
									extend3p++;
									editedAgain.setStart(peak3.getStart());
									editedAgain.setName(edited.getName()+".trimmed3p");
								}
								else{
									if(peak3.getStart()>edited.getStart()){
										if(peak3.getStart()<edited.getEnd()){
											trim3p++;
											editedAgain.setStart(peak3.getStart());
											editedAgain.setName(edited.getName()+".trimmed3p");
										}
										else{
											logger.info("For "+gene.getName()+" the 3' end is upstream of 5' end.");
										}
									}
								}
							}
							else{
								if(peak3.getEnd()>edited.getEnd()){
									extend3p++;
									editedAgain.setEnd(peak3.getEnd());
									editedAgain.setName(edited.getName()+".trimmed3p");
								}
								else{
									if(peak3.getEnd()<edited.getEnd()){
										if(peak3.getEnd()>edited.getStart()){
											trim3p++;
											editedAgain.setEnd(peak3.getEnd());
											editedAgain.setName(edited.getName()+".trimmed3p");
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
		Iterator<Alignment> readiter = model5p.getOverlappingReads(window, false);
		//for all reads in the window
		while(readiter.hasNext()){
			//AlignmentCount read = readiter.next();
			Alignment read = readiter.next();
			if(passesStartChecks(read,window,orientation)){
				windowCount += read.getWeight();
			}
		}
		return windowCount;
	}
	
	private double get3pWindowCount(Annotation window,Strand orientation){
		double windowCount = 0.0;
		//Get the reads in the window
		window.setOrientation(orientation);
		//Iterator<AlignmentCount> readiter = model3p.getOverlappingReadCountsStranded(window, false);
		Iterator<Alignment> readiter = model3p.getOverlappingReads(window, false);
		//for all reads in the window
		while(readiter.hasNext()){
			//AlignmentCount read = readiter.next();
			Alignment read = readiter.next();
			if(passesEndChecks(read,window,orientation)){
				windowCount += read.getWeight();
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
					if(passesStartChecks(read,window,annotation.getOrientation())){
						windowCount += read.getWeight();
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
					//Orientation for the 3p models will be set already on opposite mate of 5p models 
					//So sae checks as 5p models
					if(passesEndChecks(read,window,annotation.getOrientation())){
						windowCount += read.getWeight();
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
	private boolean passesStartChecks(Alignment read,Annotation window,Strand orientation){
		
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
			//Check if read is the correct read
			//if read starts in window
			if(((strand==(TranscriptionRead.FIRST_OF_PAIR) && align.getIsFirstMate()) || 
					(strand==(TranscriptionRead.SECOND_OF_PAIR) && !align.getIsFirstMate())) 
					//We used readEndFallsInWindow if the read would be in the same orientation as the gene
					//Here is it oppostite so
					//	&& (readEndFallsInWindow(read,window))
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
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"dowork");
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
		//new AddEndRNASeqToScripture(strand,argMap.getMandatory("annotations"),argMap.getOutput(),argMap.getInteger("window", DEFAULT_WINDOW_SIZE),new File(argMap.getMandatory("full")),argMap.getInteger("extension", DEFAULT_EXTENSION));

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
				if(max != null) {
					closestUpstream = max.getValue().getStart() - gene.getEnd();
				}
			}
			else{
				Node<Gene> max  = chrTree.max(gene.getStart()-1, gene.getStart() );
				if(max != null) {
					closestUpstream = gene.getStart() - max.getValue().getEnd();
				}
			}
			
		}
		return closestUpstream;
		
	}
	
	public int getDistanceToClosestSameOrientation3pGene(Gene gene) {
		IntervalTree<Gene> chrTree = intervalTrees.get(gene.getChr());
		int closestUpstream = Integer.MAX_VALUE;
		if(chrTree != null) {
			if(!gene.isNegativeStrand()){
				Node<Gene> max  = chrTree.max(gene.getEnd()+1, gene.getEnd()+2 );
				if(max != null && !max.getValue().overlaps(gene) && max.getValue().getOrientation()==gene.getOrientation()) {
					closestUpstream = max.getValue().getStart() - gene.getEnd();
					logger.info(gene.toBED());
					logger.info(max.getValue().toBED());
				}
			}
			else{
				Node<Gene> max  = chrTree.max(gene.getStart()-1, gene.getStart() );
				if(max != null && !max.getValue().overlaps(gene) && max.getValue().getOrientation()==gene.getOrientation()) {
					closestUpstream = gene.getStart() - max.getValue().getEnd();
					logger.info(gene.toBED());
					logger.info(max.getValue().toBED());
				}
			}
			
		}
		return closestUpstream;
		
	}
	
	/*static String usage=" args[0]=bam file \n\t args[1]=annotation file \n\t args[2]: Output name"
			+"\n\t args[3]= window size \n\t args[4] transcription strand";*/
	
}

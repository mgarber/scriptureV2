package nextgen.core.scripture.statistics;

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

import nextgen.core.annotation.Gene;
import nextgen.core.model.AlignmentModel;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

/**
 * @author skadri
 *
 */
public class NormalizedGeneEndCoverageForLengthClasses {

	
	static final String usage = "Usage: NormalizedGeneCoverage -task <task name> "+
			"\n\tcoverage: Computes the normalized coverage along the normalized length of the gene" +
			"\n\t\t-annotations <Specific regions to segment [BED by default]> "+
			"\n\t\t-alignment <Alignment (mapped to genome)> "+
			"\n\t\t-bins <Number of bins [100 by default]>"+
			"\n\t\t-windowSize <Step size/window size [By default: GeneLength/#bins]"+
			"\n\t\t-out <Output Filename>"+
			"\n\t\t-stranded"+
			"\n\t\t-lengthClasses <Number of length classes to make. DEFAULT:1>";
	
	static Logger logger = Logger.getLogger(NormalizedGeneEndCoverageForLengthClasses.class.getName());
	int numBins;
	double[] coverage;
	int numLenClasses;
	AlignmentModel model;
	
	public NormalizedGeneEndCoverageForLengthClasses(String[] args)throws IOException{
	
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"coverage");
		
		/*
		 * Read the annotation file
		 */
		String annotationsFile = argMap.getMandatory("annotations");
		/*
		 * Check the format of the annotation files and call the GTF or BED parser accordingly
		 */
		Map<String,Collection<Gene>> annotations =  BEDFileParser.loadDataByChr(new File(annotationsFile));
		
		numBins = argMap.isPresent("bins")? argMap.getInteger("bins") : 100;
		coverage = new double[numBins];
		numLenClasses = argMap.getInteger("lengthClasses", 1);
		/*
		 * Read the name of the alignment file
		 */
		String alignmentFile = argMap.getMandatory("alignment");
		
		boolean isStranded = argMap.isPresent("stranded");
		boolean isSecondRead = argMap.isPresent("isSecondRead");
		
		/* 
		 * ITERATE OVER THE GENES, GET LENGTHS AND THEN SPLIT INTO SAID NUMBER OF CLASSES
		 */
		Map<Gene,Double> geneLengths = new HashMap<Gene,Double>();
		for(String chr:annotations.keySet()){
			
			for(Gene gene:annotations.get(chr)){
					geneLengths.put(gene, new Double(gene.getSize()));
			}
		}
		double pct=1.0/(double)numLenClasses;
		List<Gene>[] geneClasses = new List[numLenClasses];
		for(int i=0;i<numLenClasses;i++){
			double len1 = Statistics.quantile(l2a(geneLengths.values()), pct*(i));
			double len2 = Statistics.quantile(l2a(geneLengths.values()), pct*(i+1));
			logger.info("Length class "+(i+1)+" : "+len1+"-"+len2);
			geneClasses[i] = getGenesOfLength(geneLengths,len1,len2);
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(argMap.getOutput()+".class"+(i+1)+".bed"));
			//Write a bed file of these genes
			for(Gene g:geneClasses[i]){
				bw.write(g.toBED()+"\n");
			}
			bw.close();
		}
		
		
		/*
		 * Iterate through each list of gene class
		 */
		for(int k=0;k<numLenClasses;k++){
			logger.info("Processing class: "+(k+1));
			double[] scores;
			int counter = 0;
			
			if(!isStranded){	
				/*
				 * Initialize the data model using the alignment file
				 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
				*/
				model = new AlignmentModel();
				/*
				 * For each gene
				 */
				for(Gene annotation: geneClasses[k]){
					if(libDataModel.hasDataForChromosome(annotation.getChr())){
						double max = 0.0;
						System.out.print(annotation.getName()+"\t");
						double[] c = new double[numBins];
						int annotationLength = annotation.getTranscriptLength();
						//System.out.println(Math.round((float)((annotationLength)/(float)numBins)));
						int step = argMap.isPresent("windowSize")? argMap.getInteger("windowSize") : (int)((float)(annotationLength)/(float)numBins);
						//Step = gene length/ # bins
						int threshold = ((step)*numBins);
						//System.out.println(" threshold: "+(step)*numBins);
						//System.out.print(annotation.getName()+" "+annotation.getSize()+" "+annotation.getTranscriptLength());
						if(annotationLength<threshold || threshold<numBins){
							//DONT PROCESS
							//System.out.println(" does not pass.");
						}
						else{
							//System.out.println(" passes");
							//From gene end, for each bin of size step
							if(annotation.getOrientation().equals("+")){
								for(int i=0;i<numBins;i++){
									int j=i+1;
									Gene subannotation = null;
									int relativeStart = annotation.getTranscriptLength()-(j*step);
									int relativeEnd = annotation.getTranscriptLength()-(i*step);
									if (relativeStart< 0){
										relativeStart = 0;
									}
									
									subannotation = annotation.trim(relativeStart, relativeEnd);
									//System.out.println("Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
									scores = libDataModel.scoreGene(subannotation);
									//COVERAGE
									c[numBins-i-1] = scores[3];
									if(scores[3]>max)
										max = scores[3];
									//System.out.println("Step = "+step);
									//System.out.println("Start = "+relativeStart+" End = "+relativeEnd);
								}
								
							}
							else{
								for(int i=0;i<numBins;i++){
									int j=i+1;
									Gene subannotation = null;
									int relativeStart = i*step;
									int relativeEnd = j*step;
									if(relativeEnd>annotation.getTranscriptLength()){
										relativeEnd = annotation.getTranscriptLength();
									}
									subannotation = annotation.trim(relativeStart, relativeEnd);
									//System.out.println(annotation.getName()+" Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
									//System.out.println(annotation.getName()+" "+annotation.getSize()+" Start: "+relativeStart+" End: "+relativeEnd +" index = "+(numBins-i-1));
									scores = libDataModel.scoreGene(subannotation);
									c[numBins-i-1] = scores[3];
									if(scores[3]>max)
										max = scores[3];
								}
								//System.out.println("Step: "+step+" "+annotation.getName()+": Start = "+relativeStart+" End = "+relativeEnd+" AnnotationLength = "+subannotation.getTranscriptLength());
								
								//System.out.println(coverage[i]);	
							}
							if(max>0){
								for(int i=0;i<numBins;i++){
									c[i] = c[i]/max;
									coverage[i] += c[i];
									//System.out.print(c[i]+"\t");
								}
								//System.out.println();
								counter++;
							}	
						}
					}
				}
			}
			else{
				String sizes = argMap.get("sizeFile");
				AlignmentDataModel alignmentsP=new GenericAlignmentDataModel(alignmentFile, sizes, false, 5,false,true);
				AlignmentDataModelStats alignmentDataP = new AlignmentDataModelStats(alignmentsP);
				ContinuousDataAlignmentModel libDataModelP = new ContinuousDataAlignmentModel(alignmentDataP);
				alignmentsP.setPositiveStranded();
				
				AlignmentDataModel alignmentsN=new GenericAlignmentDataModel(alignmentFile, sizes, false, 5,false,true);
				AlignmentDataModelStats alignmentDataN = new AlignmentDataModelStats(alignmentsN);
				ContinuousDataAlignmentModel libDataModelN = new ContinuousDataAlignmentModel(alignmentDataN);
				alignmentsN.setNegativeStranded();

				if(isSecondRead){
					alignmentsP.setSecondRead();
					alignmentsN.setSecondRead();
				}else{
					alignmentsP.setFirstRead();
					alignmentsN.setFirstRead();
				}
				for(Gene annotation: geneClasses[k]){
					if(libDataModelP.hasDataForChromosome(annotation.getChr())|| libDataModelN.hasDataForChromosome(annotation.getChr())){
						double max = 0.0;
						double[] c = new double[numBins];
						int annotationLength = annotation.getTranscriptLength();
						//System.out.println(Math.round((float)((annotationLength)/(float)numBins)));
						int step = argMap.isPresent("windowSize")? argMap.getInteger("windowSize") : (int)((float)(annotationLength)/(float)numBins);
						//Step = gene length/ # bins
						int threshold = ((step)*numBins);
						//System.out.println(" threshold: "+(step)*numBins);
						//System.out.print(annotation.getName()+" "+annotation.getSize()+" "+annotation.getTranscriptLength());
						if(annotationLength<threshold || threshold<numBins){
							//DONT PROCESS
							//System.out.println(" does not pass.");
						}
						else{
							//System.out.println(" passes");
							//From gene end, for each bin of size step
							if(annotation.getOrientation().equals("+")){
								for(int i=0;i<numBins;i++){
									int j=i+1;
									Gene subannotation = null;
									int relativeStart = annotation.getTranscriptLength()-(j*step);
									int relativeEnd = annotation.getTranscriptLength()-(i*step);
									if (relativeStart< 0){
										relativeStart = 0;
									}
									
									subannotation = annotation.trim(relativeStart, relativeEnd);
									//System.out.println("Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
									scores = DGE.scoreStrandedGene(libDataModelP,libDataModelN,subannotation);
									c[numBins-i-1] = scores[3];
									if(scores[3]>max)
										max = scores[3];
									//System.out.println("Step = "+step);
									//System.out.println("Start = "+relativeStart+" End = "+relativeEnd);
								}
								
							}
							else{
								for(int i=0;i<numBins;i++){
									int j=i+1;
									Gene subannotation = null;
									int relativeStart = i*step;
									int relativeEnd = j*step;
									if(relativeEnd>annotation.getTranscriptLength()){
										relativeEnd = annotation.getTranscriptLength();
									}
									subannotation = annotation.trim(relativeStart, relativeEnd);
									//System.out.println(annotation.getName()+" Start: "+relativeStart+" End: "+relativeEnd +" Length: "+subannotation.getSize()+" index = "+(numBins-i-1));
									//System.out.println(annotation.getName()+" "+annotation.getSize()+" Start: "+relativeStart+" End: "+relativeEnd +" index = "+(numBins-i-1));
									scores = DGE.scoreStrandedGene(libDataModelP,libDataModelN,subannotation);
									c[numBins-i-1] = scores[3];
									if(scores[3]>max)
										max = scores[3];
								}
								//System.out.println("Step: "+step+" "+annotation.getName()+": Start = "+relativeStart+" End = "+relativeEnd+" AnnotationLength = "+subannotation.getTranscriptLength());
								
								//System.out.println(coverage[i]);	
							}
							if(max>0){
							for(int i=0;i<numBins;i++){
								c[i] = c[i]/max;
								coverage[i] += c[i];
							}
							counter++;
							}
						}
					}
				}
			}	
			for(int i=0;i<numBins;i++){
				coverage[i] =coverage[i]/(double)counter;
			}
			/*
			 * Output Filename
			 */
			//String outputFileName = argMap.get("out");
			BufferedWriter outBw =new BufferedWriter(new FileWriter(argMap.getOutput()+".class"+(k+1)+".coverage"));
			/*
			 * Write the coverage one in one row
			 */
			//BufferedWriter outBw = new BufferedWriter(new FileWriter(outputFileName));
			for(int i=0;i<coverage.length;i++){
				outBw.write(coverage[i]+"\n");
			}
			outBw.close();
			//writeOutputToFile(outputFileName,coverage);
		}
	}
	
	/**
	 * Returns a list of genes in the specified length range
	 * @param geneLengths
	 * @param len1
	 * @param len2
	 * @return
	 */
	private List<Gene> getGenesOfLength(Map<Gene,Double> geneLengths,double len1,double len2){
		
		List<Gene> genes = new ArrayList<Gene>();
		for(Gene g:geneLengths.keySet()){
			if(g.getTranscriptLength()>=len1 && g.getTranscriptLength()<=len2){
				genes.add(g);
			}
		}
		return genes;
	}
	
	/**
	 * This function writes the output to two files
	 * @param outFileName Name of the output file provided by user
	 * @throws IOException
	 */
	private static void writeOutputToFile (String outFileName,double[] coverage) throws IOException{
		
		/*
		 * Write the coverage in a row
		 */
		BufferedWriter outBw = new BufferedWriter(new FileWriter(outFileName));
		for(int i=0;i<coverage.length;i++){
			outBw.write(coverage[i]+"\t");
		}
		outBw.newLine();
		outBw.close();
	}
	
	private double[] l2a(Collection<Double> collection){
		double[] rtrn=new double[collection.size()];
	
		int i=0;
		for(Double val: collection){rtrn[i++]=val;}
	
		return rtrn;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		new NormalizedGeneEndCoverageForLengthClasses(args);
	}

}

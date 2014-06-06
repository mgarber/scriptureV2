package nextgen.core.scripture.statistics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.general.CloseableFilterIterator;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.PairedEndFilter;
import nextgen.core.scripture.BuildScriptureCoordinateSpace;
import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.error.ParseException;
import broad.core.math.ScanStatistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

public class TestScanStatistics {
	
	static Logger logger = Logger.getLogger(TestScanStatistics.class.getName());
	Map<String,Collection<Gene>> annotations;
	private AlignmentModel model;
	private static double DEFAULT_ALPHA = 0.05;
	private CoordinateSpace space;
	int counter = 1000;
	
	static final String usage = "Usage: CalculateSignific -task doWork "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-in <Reconstruction bed file. [BED by default]> "+
			"\n\t\t-maxSize <Maximum insert size. Default=500bp> "+
			"\n";

	public TestScanStatistics(String annotationFile,File bamFile,TranscriptionRead strand) throws IOException{
		
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
//		model.addFilter(new PairedAndProperFilter());
		dowork(annotations);
	}
	
	/**
	 * This function will go through each gene and calculate the number of scan statistics for each gene in various ways
	 * @return
	 * @throws IOException 
	 */
	private void dowork(Map<String,Collection<Gene>> geneMap){
		
		for(String chr:geneMap.keySet()){
//			try{
			//MAKE A MAP OF GENE TO ISOFORMS
//			Map<Gene,Set<Gene>> isoformMap = BuildScriptureCoordinateSpace.getIsoformMap(geneMap.get(chr));			
			//For each gene
			for(Gene gene:geneMap.get(chr)){			
				logger.info(gene.getName()+"Gene length = "+gene.size());
				//1. Calculate scores using normal reads
				logger.info("Calculate scores using normal reads");
//				double[] scores = getScoresReads(gene);
//				logger.info(scores[0]+"\t"+scores[1]);
				//2. Calculate scores using fragments - without correction
				logger.info("Calculate scores using fragments - without correction");
				double[] scores = getScoresFragments(gene);
				logger.info(scores[0]+"\t"+scores[1]);
				
				double s;
				double globalPairedlambda = 5.609947485195715E7/3.095693983E9;
				//3A. Calculate scores using alpha = (w)*lambdaf (lambdaf = lambda+f)
/*				logger.info("Calculate scores using alpha = (w)*lambdaf (lambdaf = lambda+f)");
				double s = getScoresFragmentsLambda(scores[0],(scores[2]+518.0)/model.getGlobalLength(),gene.getSize());
				logger.info(s);
				//3B. Calculate scores using alpha = (w+lf)*lambdaf (lambdaf = lambda+f)change in gene size
				logger.info("Calculate scores using alpha = (w+lf)*lambdaf (lambdaf = lambda+f) change in gene size");
				s = getScoresFragmentsLambda(scores[0],(scores[2]+518.0)/model.getGlobalLength(),(gene.getSize()+518));
				logger.info(s);
*/				//4A. Calculate scores using alpha = (w)*lambdaf (lambdaf = lambda*f) NO change in gene size
				logger.info("Calculate scores using alpha = (w)*lambdaf (lambdaf = lambda*f)NO change in gene size");
				s = getScoresFragmentsLambda(scores[0],(scores[2]*270.0)/model.getGlobalLength(),gene.getSize());
				logger.info(s);
				//4B. Calculate scores using alpha = (w*lf)*lambdaf (lambdaf = lambda*f) change in gene size
				logger.info("Calculate scores using alpha = (w*lf)*lambdaf (lambdaf = lambda*f)change in gene size");
				s = getScoresFragmentsLambda(scores[0],(scores[2]*270.0)/model.getGlobalLength(),(gene.getSize()+270));
				logger.info(s);
				//4C. Calculate scores using alpha = (w*lf)*lambdaf (lambdaf = lambda*f) change in gene size
/*				logger.info("Calculate scores using alpha = (w/lf)*lambdaf (lambdaf = lambda*f) change in gene size divide by insert size");
				s = getScoresFragmentsLambda(scores[0],(scores[2]*270.0)/model.getGlobalLength(),(gene.getSize()/518));
				logger.info(s);
*/			}
//			} 
		}
		
	}
	
	/**
	 * Returns paired end counts and scan p-value for the specified gene
	 * @param gene
	 * @return
	 */
	private double[] getScoresReads(Annotation gene){
		
		logger.info("Lambda = "+model.getGlobalLambda());
		double[] scores = new double[2];
		//Get all reads overlapping the transcript
		scores[0] = model.getCount(gene, true);
		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		
		return scores;
	}
	
	
	/**
	 * Returns paired end counts and scan p-value for the specified gene
	 * @param gene
	 * @return
	 */
	private double[] getScoresFragments(Annotation gene){
		
		double globalPairedLambda = model.getGlobalPairedFragments()/model.getGlobalLength();
		
		logger.info("Lambda = "+globalPairedLambda);
		double[] scores = new double[3];
		scores[0] = 0.0;
		//Get all reads overlapping the transcript
		//get Alignments over the whole region
		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(gene,true), new PairedEndFilter());
		//For each read,
		while(iter.hasNext()){					
			Alignment read = iter.next();
			boolean countRead = true;
			for(Annotation mate:read.getReadAlignments(space)){
				if(!BuildScriptureCoordinateSpace.compatible(gene,mate)){
					//logger.debug("Read "+mate.toUCSC()+" is not compatible with isoform with "+isoform.getExons().length);
					countRead=false;
					break;
				}
			}
			//For the assembly, we need to treat each read separately	
			if(countRead){
				scores[0] += read.getWeight();
			}	
		}
		iter.close();
		//logger.info("Count = "+scores[0]+" Int version "+new Double(scores[0]).intValue()+" global paired lambda = "+globalPairedLambda+" gene size = "+model.getCoordinateSpace().getSize(gene)+ " or "+gene.size()+" global length = "+model.getGlobalLength()+" global lambda = "+model.getGlobalLambda());
		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), globalPairedLambda,gene.size(), model.getGlobalLength());
		//logger.info("Pvalue = "+scores[1]);
//		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		scores[2] = model.getGlobalPairedFragments();
		return scores;
	}
	

	/**
	 * Returns paired end counts and scan p-value for the specified gene
	 * @param gene
	 * @return
	 */
	private double getScoresFragmentsLambda(double count,double lambda,int geneSize){
		
		logger.info("Lambda = "+lambda);
		//logger.info("Count = "+scores[0]+" Int version "+new Double(scores[0]).intValue()+" global paired lambda = "+globalPairedLambda+" gene size = "+model.getCoordinateSpace().getSize(gene)+ " or "+gene.size()+" global length = "+model.getGlobalLength()+" global lambda = "+model.getGlobalLambda());
		double score = ScanStatistics.calculatePVal(new Double(count).intValue(), lambda, geneSize, model.getGlobalLength());
		//logger.info("Pvalue = "+scores[1]);
//		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		
		return score;
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
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"doWork");
		TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand").equalsIgnoreCase("first")){
			//System.out.println("First read");
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(argMap.get("strand").equalsIgnoreCase("second")){
			//System.out.println("Second read");
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}
		else
			System.out.println("no strand");
		
		double G = 3.095693983E9;
		double lambdaR = 0.026025680826696183;
		double lambdaP = 0.01812177662263368;
		System.out.println(ScanStatistics.calculatePVal(373, lambdaR,3705.0, G));
		System.out.println(ScanStatistics.calculatePVal(208, lambdaP,3705.0, G));
		System.out.println(ScanStatistics.calculatePVal(44, lambdaP,777, G));
		System.out.println(ScanStatistics.calculatePVal(549, lambdaP,4261, G));
		
		
		Gene g = new Gene("chr1",100,5000,"gene",Strand.NEGATIVE);
		System.out.println(g.toBED());
		Gene gs = g.trimGene(g.length()-1000, g.length());
		System.out.println(gs.toBED());

		//new TestScanStatistics(argMap.getInput(),new File(argMap.getMandatory("alignment")),strand);
		
		
	}
}

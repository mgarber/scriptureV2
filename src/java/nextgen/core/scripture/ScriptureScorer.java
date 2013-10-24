package nextgen.core.scripture;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readFilters.MappingQualityFilter;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

/**
 * @author skadri
 *
 */
public class ScriptureScorer {
	
	static Logger logger = Logger.getLogger(ScriptureScorer.class.getName());
	static int DEFAULT_MIN_MAPPING_QUALITY = -1;
	static final int COUNT_SCORE = 0;
	static final int RPKM_SCORE = 1;
	static final int RPK_SCORE = 2;
	static final int PVAL_SCORE = 6;
	private static final TranscriptionRead DEFAULT_TXN_READ =  TranscriptionRead.UNSTRANDED;
	private TranscriptionRead strand;
	private String outName;
	
	Map<String,Collection<Gene>> annotations;
	
//	boolean useConstituentExons;
//	boolean useConstituentIntrons;
	int minimumMappingQuality;
	boolean weighReadsFlag;
	boolean removePCRDuplicatesFlag; 
	boolean normalizedOutput;
	
	public ScriptureScorer(File bamFile,ArgumentMap argMap) throws IOException{
		
		init(argMap);
		AlignmentModel model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,true);
		
	}
	
	
	public ScriptureScorer(List<String> alignmentFiles,String annotationFile,ArgumentMap argMap) throws IOException{
				
		init(argMap);
		/*
		 * Check the format of the annotation file and call the GTF or BED parser accordingly
		 */
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		List<String> rows = new ArrayList<String>();
		
		//TODO: MAKE GENES
		
		// HashMap of gene to rowName
		HashMap<Gene, String> duplicateNameMap = new HashMap<Gene, String>();
		// HashMap of gene name to the number of duplicates
		HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
		
		Collection<Gene> annotationCollection = new ArrayList<Gene>();
		for(Collection<Gene> chrAnnotations : annotations.values()) {
			annotationCollection.addAll(chrAnnotations);
		}
		
		int duplicates=0;
		for(String chr:annotations.keySet()){
			for(Gene gene:annotations.get(chr)){
				if(!rows.contains(gene.getName())) {
					String name = gene.getName();
					if(duplicateNameMap.containsKey(gene)){
						logger.info("Entry for "+name+" already exists");
					}
					else{
						rows.add(name);
						duplicateMap.put(name, 1);
						duplicateNameMap.put(gene, name);
					}
				}
				// If the gene name has another annotation
				else {
					
					if(duplicateNameMap.containsKey(gene)){
						logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
					}
					else{
						//Row name is now the geneName appended with the duplicate number
						duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
						String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
						rows.add(name);
						duplicateNameMap.put(gene, name);
						//logger.warn("Duplicated annotation : "+ gene.toBED());
						duplicates++;
					}
				}
			}
		}
			 
		/*
		 * The columns of the matrix will be each alignment file name
		 */
		List<String> cols = new ArrayList<String>();
		for(String name: alignmentFiles){
			cols.add(name);
		}
		
		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders countsMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders rpkmMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders rpkMatrix = new MatrixWithHeaders(rows,cols);

		//TODO: Add this
//		File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;

		//FOR EACH ALIGNMENT FILE
		for(int i=0;i<alignmentFiles.size();i++){
			
			String outputname = alignmentFiles.get(i)+"."+outName+".bed";
			logger.info("Processing : "+alignmentFiles.get(i));
			//Map<String, Integer> maskFileData = ContinuousDataAlignmentModel.parseMaskFiles(maskFiles);
			
			/*
			 * Initialize the data models for all alignment files
			 * @param: <alignment flieName> <coordinate space> <read filters> <read or create paired end bam file> <load as fragments>
			 */
			boolean pairedFlag = !argMap.isPresent("singleEnd");
			logger.info("Paired flag is "+pairedFlag);
			CoordinateSpace space = new TranscriptomeSpace(annotations);
			AlignmentModel model = new AlignmentModel(alignmentFiles.get(i), space, new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true);
			
			//Add read filters
			model.addFilter(new MappingQualityFilter(minimumMappingQuality,minimumMappingQuality));
			
			Map<Gene, double[]> scores=runScore(model);
			for(String chr : annotations.keySet()){
				for(Gene g:annotations.get(chr)){
					//Might be redundant. Check if Scan statistic score class sets this
					g.setBedScore(scores.get(g)[PVAL_SCORE]);
					g.setExtraFields(scores.get(g));
					countsMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scores.get(g)[COUNT_SCORE]);
					rpkmMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scores.get(g)[RPKM_SCORE]);
					rpkMatrix.set(duplicateNameMap.get(g), alignmentFiles.get(i), scores.get(g)[RPK_SCORE]);
				}
			}
			
			//WRITE THE BED FILE
			//BuildScriptureCoordinateSpace.write(outputname, annotations);
		}
		
		
		BufferedWriter outBw = new BufferedWriter(new FileWriter(outName+".counts.txt"));
		countsMatrix.write(outBw);
		outBw.close();
		outBw = new BufferedWriter(new FileWriter(outName+".rpkm.txt"));
		rpkmMatrix.write(outBw);
		outBw.close();
		outBw = new BufferedWriter(new FileWriter(outName+".rpk.txt"));
		rpkMatrix.write(outBw);
		outBw.close();
		
		BufferedWriter bedBw = new BufferedWriter(new FileWriter(outName+".geneNameMap.bed"));
		for(Gene gene: annotationCollection){
			bedBw.write(gene.toBED());
			bedBw.append("\t"+duplicateNameMap.get(gene));
			bedBw.newLine();
		}
		bedBw.close();
		
		if(normalizedOutput){
			writeNormalizedMatrix(outName, rpkMatrix);
		}
		
	}

	/**
	 * Initializes and sets the parameters
	 * @param argMap
	 */
	public void init(ArgumentMap argMap){
		
		//Output name
		outName = argMap.getOutput();

		//Transcription strand
		strand = DEFAULT_TXN_READ;
		if(argMap.get("strand").equalsIgnoreCase("first")){
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(argMap.get("strand").equalsIgnoreCase("second")){
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}
		else
			logger.info("no strand");
		
//		useConstituentExons = argMap.containsKey("useConstituentExons");
//		useConstituentIntrons = argMap.containsKey("useConstituentIntrons");
		/*
		 * Parameters if not provided are set to defaults:
		 * maxPrematureEnd : default 0
		 * max3Pextension : default 0
		 * window: default 500
		 * minMappingQuality : default -1
		 */
		minimumMappingQuality = argMap.isPresent("minMappingQuality")? argMap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;

		/*
		 * FLAG for WEIGHING READS BY NH FLAG
		 * TRUE by default
		 * Convert string to boolean
		 */
		weighReadsFlag = argMap.isPresent("weighReadsFlag")? (Boolean.parseBoolean(argMap.get("weighReadsFlag"))): true;
		
		/*
		 * FLAG for REMOVING PCR DUPLICATES 
		 * FALSE by default
		 * Convert string to boolean
		 */
		//TODO: Add this read filter
		removePCRDuplicatesFlag = argMap.isPresent("removePCRDuplicatesFlag")? (Boolean.parseBoolean(argMap.get("removePCRDuplicatesFlag"))): false;

		/*
		 * FLAG to return a normalized matrix 
		 * FALSE by default
		 * Convert string to boolean
		 */
		normalizedOutput = argMap.isPresent("normalizedOutput")? (Boolean.parseBoolean(argMap.get("normalizedOutput"))): false;
	}
	/**
	 * Writes the normalized matrix
	 * @param fileName
	 * @param resultMatrix
	 * @throws IOException
	 */
	public static void writeNormalizedMatrix(String fileName, MatrixWithHeaders resultMatrix) throws IOException{
		
		Map<String,Double> factors = calculateNormalizationFactors(resultMatrix);
		MatrixWithHeaders normalizedMatrix = resultMatrix.multiplyColumnsWithConstants(factors);
		String normalizedFileName = fileName+".normalized";
		String normalizedFactorFile = fileName+".coverage";
		BufferedWriter outBw = new BufferedWriter(new FileWriter(normalizedFileName));
		normalizedMatrix.write(outBw);
		outBw.close();
		
		BufferedWriter fBw = new BufferedWriter(new FileWriter(normalizedFactorFile));
		for(String ss:factors.keySet()){
			fBw.write(ss+"\t"+((double)1.0/factors.get(ss))+"\n");
		}
		fBw.close();
	}

	
	public static Map<String,Double> calculateNormalizationFactors(MatrixWithHeaders mat){
		
		int M = mat.getNumberColumns();
		int N = mat.getNumberRows();
		
		Map<String,Double> factors = new HashMap<String,Double>();
		double[] means = new double[N];
		//CALCULATE GEOMETRIC MEAN FOR ALL N GENES
		for(int i=0;i<N;i++){
			means[i] = 1.0;
			for(int j=0;j<M;j++){
				means[i] = means[i] * mat.get(i, j);
			}
			means[i] = Math.pow(means[i],((double)1.0/(double)M));
		}
		//For each sample
		for(String c:mat.getColumnNames()){
			
			double[] col = new double[N];
			for(int i=0;i<N;i++){
				col[i] = mat.get(i, c)/means[i]; 
			}
			factors.put(c, (double)1.0/Statistics.median(col));
		}
		return factors;
	}
	
	/**
	 * This function converts an array of double to List of type Double
	 * @param double[]
	 * @return
	 */
	public static List<Double> array2List(double[] array){
		List<Double> lt=new ArrayList<Double>();

		for(int i=0;i<array.length;i++){
			lt.add(array[i]);
		}
	
		return lt;
	}
	
	/**
	 * Retures scores for the annotations
	 * @param annotations
	 * @param save
	 * @param maskFileData
	 * @param alignments
	 * @param scores
	 * @return
	 * @throws IOException
	 */
	private Map<Gene, double[]> runScore(AlignmentModel model) throws IOException {
		
		Map<Gene, double[]> scores=new TreeMap<Gene, double[]>();
		
		if(strand.equals(TranscriptionRead.UNSTRANDED)){
			for(String chr:annotations.keySet()){
				if(!model.containsReference(chr)){

				}
			}
		}
		for(String chr : annotations.keySet()) {
			logger.info("processing " + chr);			
			Collection<Gene> chrAnnotations = annotations.get(chr);
			scores.putAll(scoreGenes(model, chrAnnotations));
		}
		return scores;
	}
	
	public static void main(String[] args)throws IOException{
		 
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"score");
		
		if(argMap.containsKey("alignment")){
			new ScriptureScorer(new File(argMap.getMandatory("alignment")),argMap);
		}
		else{
			/*
			 * Read the names of the alignment file into an array
			 */
			BufferedReader br = new BufferedReader(new FileReader(argMap.getMandatory("alignments")));
			List<String> alignmentFiles = new ArrayList<String>();
			String s;
			while((s = br.readLine())!= null){
				alignmentFiles.add(s);
			}
			new ScriptureScorer(alignmentFiles,argMap.getMandatory("annotations"),argMap);
		}
	}
	
	public Map<Gene, double[]> scoreGenes(AlignmentModel model,Collection<Gene> genes){
		Map<Gene, double[]> rtrn=new TreeMap<Gene, double[]>();
		
		for(Gene gene: genes){
			if(!model.containsReference(gene.getChr())){
				double[] s = new double[10];
				//p-value
				s[6] = 1.0;
				rtrn.put(gene, s);
			}
			else{
				rtrn.put(gene, new ScanStatisticScore(model, gene).getScores());
			}
		}
		return rtrn;
	}
	
	static final String usage = "Usage: ScriptureScorer -task score "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Single bam file> "+
			"\n\t\tOR"+
			"\n\t\t-alignments <List of bam files to score from> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+		
			"\n\t\t-strand <VALUES: first, second, unstranded. Specifies the mate that is in the direction of transcription DEFAULT: Unstranded> "+
			"\n\t\t-annotations <Annotations bed file. [BED by default]> "+

			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\n\t\t-singleEnd Only if data is single end"+
			"\n";

}

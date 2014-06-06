package nextgen.core.scripture;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.model.AlignmentModel;

import org.apache.commons.collections15.Predicate;
import org.broad.igv.Globals;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

/**
 * This class takes an alignment file 
 * @author skadri
 *
 */
public class SignificanceCalculator {
	
	TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
	
	public SignificanceCalculator(File bamFile,ArgumentMap argMap) throws IOException{
		
		//Matrix contains the count of the sample of interest. 
		MatrixWithHeaders matrix = new MatrixWithHeaders(argMap.getInput());
		//Matrix contains the lengths of the genes. 
		MatrixWithHeaders lengths = new MatrixWithHeaders(argMap.getMandatory("lengths"));
		
		if(argMap.get("strand").equalsIgnoreCase("first")){
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(argMap.get("strand").equalsIgnoreCase("second")){
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}

		boolean pairedFlag = !argMap.isPresent("singleEnd");
		System.out.println("Paired flag is "+pairedFlag);
		
		List<String> sigNames = new ArrayList<String>();
		AlignmentModel model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),pairedFlag,strand,true);
		
		//Column number in count matrix to be used. Starts at 0
		int columnNumber = argMap.getInteger("column",0);
		System.out.println("Using column number "+columnNumber+" equivalent to "+matrix.getColoumnName(columnNumber));
		
		List<String> sigRows = new ArrayList<String>();
		Map<String,Double> pvalues = new HashMap<String,Double>();
		
		//FOr each row in lengths file
		for(int i=0;i<lengths.rowDimension();i++){
			String rowName = lengths.getRowName(i);
			int count = 0;
			if(matrix.containsRow(rowName)){
				count = new Double(matrix.get(rowName, columnNumber)).intValue();
			}
			double pval = AlignmentDataModelStats.calculatePVal(count, model.getGlobalLambda(), lengths.get(rowName, 0),new Double(model.getGlobalLength()));
			
			pvalues.put(rowName, pval);
			
			if(pval<0.05){
				sigRows.add(rowName);
				sigNames.add(rowName);
			}
		}
		
		MatrixWithHeaders outputMatrix = new MatrixWithHeaders(matrix, sigRows, matrix.getColumnNames());
		
		outputMatrix.write(argMap.getOutput());
		
		FileWriter writer=new FileWriter(argMap.getOutput()+".names");
				
		for(String sig:sigNames){
			writer.write(sig+"\n");
		}
		writer.close();
		
		writer = new FileWriter(argMap.getOutput()+".pvalues");
		for(String p:pvalues.keySet()){
			writer.write(p+"\t"+pvalues.get(p)+"\n");
		}
		writer.close();
	}
	

	
	public static void main(String[] args)throws IOException{
		 
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"score");
		
		new SignificanceCalculator(new File(argMap.getMandatory("alignment")),argMap);
	}
	
	static final String usage = "Usage: SignificanceCalculator -task score "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			
			"\n\n\t\t-alignment <Single bam file> "+
			"\n\t\t-in <Matrix of counts> "+
			"\n\t\t-lengths <Lengths matrix file. First column has lengths> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-column <Column number in the file to use counts from. Starts at 0. Default is 0> "+
			"\n\t\t-strand <VALUES: first, second, unstranded. Specifies the mate that is in the direction of transcription DEFAULT: Unstranded> "+
			"\n\t\t-singleEnd <Specify is single end alignment > "+
			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\n";

}



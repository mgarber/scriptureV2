package broad.pda.seq.rap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.rpc.ServiceException;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;

import uk.ac.ebi.webservices.axis1.EmbossMatcherClient;
import uk.ac.ebi.webservices.axis1.stubs.emboss.matcher.InputParameters;

/**
 * This class calls the EMBOSS MATCHER program to two sets of nucleotide sequences
 * @author skadri
 *
 */
public class Homology extends CommandLineProgram {

	static Logger logger = Logger.getLogger(Homology.class.getName());
	
	@Option(doc="Input1 FASTA", shortName="I")
    String INPUT1;
	
	@Option(doc="Input2 FASTA", shortName="R")
    String INPUT2;
	
	@Option(doc="Email", shortName="E")
    String EMAIL;
	
	@Option(doc="output", shortName="O")
    String OUTPUT;
	
	Map<String,Double> scores;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.exit(new Homology().instanceMain(args));

	}

	@Override
	protected int doWork() {
		try {

			//Initialize all result matrices
			scores = new HashMap<String,Double>();
			
			// prepare service client
			EmbossMatcherClient client = new EmbossMatcherClient();
			InputParameters parameters = new InputParameters();

			List<Sequence> aseq=new FastaSequenceIO(INPUT1).loadAll();
			List<Sequence> bseq=new FastaSequenceIO(INPUT2).loadAll();
			
			for(Sequence s1:aseq){
				String seq1 = s1.getSequenceBases();
				//System.out.println("Seq1:"+seq1);
				for(Sequence s2:bseq){
					String seq2 = s2.getSequenceBases();
					
					parameters.setAsequence(seq1);
					parameters.setBsequence(seq2);
					
					int result = runMatcher(client,parameters);
					if(result!=0){
						logger.error("Program is exiting");
						return 1;
					}
					else{
						//TODO
					}
				}
			}			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return 0;
	}

	private int runMatcher(EmbossMatcherClient client,InputParameters parameters) throws ServiceException, IOException, InterruptedException{
		
		String jobid = client.runApp(EMAIL,"SomeTitle", parameters);
		logger.info("JobID:"+jobid);
		String status = client.checkStatus(jobid);
		//While running
		while (status.equals("RUNNING")) {
			Thread.sleep(1000);
			System.out.print(".");
			status = client.checkStatus(jobid);
		}
		System.out.println();
		
		if (status.equals("ERROR")) {
			logger.error("exception: An error occurred when attempting to get the job status.");
			return 1;
		}

		if (status.equals("FAILURE")) {
			logger.error("exception: The job failed.");
			return 1;
		}

		if (status.equals("NOT_FOUND")) {
			logger.error("exception: The job cannot be found.");
			return 1;
		}
		
		
		//If run successfully, get result
		// get system properties
		String tmpdir = System.getProperty("java.io.tmpdir");
		logger.info("Temp directory is "+tmpdir);
		String lineseparator = System.getProperty("line.separator");

		// get results
		String[] resultfiles = client.getResults(jobid, tmpdir + "/" + jobid, null);
		for(String s:resultfiles){
			logger.info("Result file:"+s);
		}
		
		// out
		BufferedReader reader = new BufferedReader(new FileReader(resultfiles[0]));
		logger.info("Result file 0 is \n"+resultfiles[0]);
		StringBuilder result = new StringBuilder();
		String line = "";
		while ((line = reader.readLine()) != null) {
			result.append(line);
			result.append(lineseparator);
		}
		logger.info("Result is:"+result.toString());
		reader.close();
		
		//TODO: Delete the files
		
		return 0;
	}
}

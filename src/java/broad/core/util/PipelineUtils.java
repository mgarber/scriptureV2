package broad.core.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class PipelineUtils {

	static int waitTime=60000; //1 minute
	
	public static int bsubProcess(Runtime run, String command) throws IOException, InterruptedException {
		//System.err.println("BSUB: " + command);
		Process p=run.exec(command);
		p.waitFor();
		int exitVal=p.exitValue();
		InputStream is = p.getInputStream();
		writeError(is);
		is.close();
		is = p.getErrorStream();
		writeError(is);
		is.close();
		//System.err.println("END BSUB " );
		
		return exitVal;
	}
	
	public static int bsubProcess(Runtime run, String jobID, String command, String output) throws IOException, InterruptedException {
		return bsubProcess(run, jobID, command, output, "week");
	}
	public static int bsubProcess(Runtime run, String jobID, String command, String output, String queue) throws IOException, InterruptedException {
		String fullCommand="bsub -R rusage[mem=5] -q " + queue+" -J "+jobID+ " -o "+output+" "+command;
		return bsubProcess(run, fullCommand);
	}

	public static int bsubProcess(Runtime run, String jobID, String command, String output, String queue, int memory) throws IOException, InterruptedException {
		String fullCommand="bsub -M " + memory + " -q " + queue+" -J "+jobID+ " -o "+output+" "+command;
		return bsubProcess(run, fullCommand);
	}
	
	public static int bsubMultiProcess(Runtime run, String jobID, String command, String output, String queue) throws IOException, InterruptedException {
		String fullCommand="bsub -R rusage[mem=20] -n 4,8"+ " -q " + queue+" -J "+jobID+ " -o "+output+" "+command;
		System.out.println(fullCommand);
		return bsubProcess(run, fullCommand);
	}
	/**
	 * Submit one command to hour queue with no memory requirement
	 * @param run
	 * @param jobID
	 * @param command
	 * @param output
	 * @return
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static int bsubSmallProcess(Runtime run, String jobID, String command, String output) throws IOException, InterruptedException {
		String fullCommand = "bsub -q hour -J " + jobID + " -o " + output + " " + command;
		return bsubProcess(run, fullCommand);
	}

	/**
	 * Submit one command to hour queue with 4GB requirement
	 * @param run
	 * @param jobID
	 * @param command
	 * @param output
	 * @return
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static int bsubMediumProcess(Runtime run, String jobID, String command, String output) throws IOException, InterruptedException {
		String fullCommand = "bsub -q hour -R rusage[mem=4] -J " + jobID + " -o " + output + " " + command;
		return bsubProcess(run, fullCommand);
	}

	public static int bsubProcess(Runtime run, String jobID, String[] commands, String output, String queue) throws IOException, InterruptedException {
		String fullCommand="bsub -R rusage[mem=5] -q "+queue+" -J "+jobID+ " -o "+output;
		
		fullCommand+=" "+commands[0];
		for(int i=1; i<commands.length; i++){fullCommand+="; "+commands[i];}
		
		System.out.println(fullCommand);
		return bsubProcess(run, fullCommand);
	}
	
	
	public static int bsubProcess(Runtime run, String jobID, String[] commands, String output) throws IOException, InterruptedException {
		return bsubProcess(run, jobID, commands, output, "week");
	}
	
	/**
	 * Submit several commands to hour queue with no memory requirement
	 * @param run
	 * @param jobID
	 * @param commands
	 * @param output
	 * @return
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static int bsubSmallProcess(Runtime run, String jobID, String[] commands, String output) throws IOException, InterruptedException {
		String fullCommand="bsub -q hour -J "+jobID+ " -o "+output;
		
		fullCommand+=" "+commands[0];
		for(int i=1; i<commands.length; i++){fullCommand+="; "+commands[i];}
		
		System.out.println(fullCommand);
		return bsubProcess(run, fullCommand);
	}

	
	public static void writeError(InputStream errorStream) throws IOException {
		BufferedReader reader=	new BufferedReader(new InputStreamReader(errorStream));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
		}
		System.err.println();
	}
	
	public static void waitForJobs(String jobID, Runtime run) throws IOException, InterruptedException {	
		waitForJobs(jobID, run, waitTime, true);
	}

	
	public static void waitForJobs(String jobID, Runtime run, boolean verbose) throws IOException, InterruptedException {	
		waitForJobs(jobID, run, waitTime, verbose);
	}
	
	public static void waitForJobs(String jobID, Runtime run, int waitTime, boolean verbose) throws IOException, InterruptedException {	
		int running=5;
		while(running>1){
			Thread.sleep(waitTime);
			Process blatProc = run.exec("bjobs -J "+jobID);
			BufferedReader out = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));
			running=parseReply(out);
			if(verbose) System.err.println("Checking status: "+running+" jobs still running");
		}
		System.err.println("done");
		PipelineUtils.checkForLSFFailures(jobID, run);
	}
	

	
    public static void waitForAllJobs(ArrayList<String> jobIDs, Runtime run) throws InterruptedException, IOException {
    	if(jobIDs.isEmpty()) return;
    	boolean running = true;
    	boolean found = false;
    	int failedJobs = 0;
		while(running) {
			Thread.sleep(10000);
	       	found = false;  
	       	ArrayList<String> stillRunning = new ArrayList<String>();
		    for(String jobID : jobIDs) {
		    	Process blatProc = run.exec("bjobs -J " + jobID);
		    	BufferedReader out = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));
		    	int r = parseReply(out);
		    	try {
		    		PipelineUtils.checkForLSFFailures(jobID, run);
		    	} catch (IllegalArgumentException e) {
		    		failedJobs++;
		    	}
		    	if(r>1) {
		    		// At least one job is still running
		    		found = true;
		    		stillRunning.add(jobID);
		    	} 
		    }
		    jobIDs.clear();
		    jobIDs.addAll(stillRunning);
		    if(!found) running = false;
		}
	       
		if(failedJobs > 0) {
			throw new IllegalArgumentException(failedJobs + " jobs failed.");
		}
		
  		
    }

    
    public static void bkillAll(ArrayList<String> jobIDs, Runtime run) throws IOException {
    	for(String jobID : jobIDs) {
    		Process kill = run.exec("bkill " + jobID);
    	}
    }
    
    /**
     * Wait for a subset of jobs to finish
     * @param jobIDs the set of job IDs to wait for
     * @param numJobs the number of successful jobs needed
     * @param run
     * @throws InterruptedException
     * @throws IOException
     * @return the list of jobs still running
     */
    public static ArrayList<String> waitForEnoughJobs(ArrayList<String> jobIDs, int numJobs, Runtime run) throws InterruptedException, IOException {
    	boolean running = true;
    	boolean found = false;
    	int successfulJobs = 0;
    	int failedJobs = 0;
    	int failMax = jobIDs.size() - numJobs;
		while(running) {
			Thread.sleep(10000);
	       	found = false;  
	       	ArrayList<String> stillRunning = new ArrayList<String>();
		    for(String jobID : jobIDs) {
		    	Process blatProc = run.exec("bjobs -J " + jobID);
		    	BufferedReader out = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));
		    	int r = parseReply(out);
		    	try {
		    		PipelineUtils.checkForLSFFailures(jobID, run);
		    		successfulJobs++;
		    		if(successfulJobs >= numJobs) {
		    			// Enough jobs have completed successfully
		    			jobIDs.remove(jobID);
		    			return jobIDs;
		    		}
		    	} catch (IllegalArgumentException e) {
		    		failedJobs++;
		    		if(failedJobs > failMax) {
		    			throw new IllegalArgumentException(failedJobs + " failed, max allowable was " + failMax + ".");
		    		}
		    	}
		    	if(r>1) {
		    		// At least one job is still running
		    		found = true;
		    		stillRunning.add(jobID);
		    	} 
		    }
		    jobIDs.clear();
		    jobIDs.addAll(stillRunning);
		    if(!found) running = false;
		}
	       
		if(successfulJobs < numJobs) {
			throw new IllegalArgumentException("Only " + successfulJobs + " finished, needed " + numJobs + ".");
		}
				
		return jobIDs;
		
    }


    
	
	private static int parseReply(BufferedReader lsfOut) throws IOException {
		String line = null;
		String [] lineInfo = null;
		int i=0;
		while((line = lsfOut.readLine()) != null) {i++;}
		return i;
	}
	
	public static void checkForLSFFailures(String jobID, Runtime run) throws IOException {
		Process proc = run.exec("bjobs -a -J "+jobID);
		BufferedReader out = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String line = null;
		int i=0;
		while((line = out.readLine()) != null) {
			//split the line
			String[] tokens=line.split(" ");
			String completionStatus=tokens[2];
			if(completionStatus.equalsIgnoreCase("EXIT") || line.contains("job killed")){
				System.err.println("WARN: Job "+tokens[0]+" FAILED"); 
				throw new IllegalArgumentException("Job " + jobID + " failed");
			}
			//System.err.println(i);
			//System.err.println(line);
			i++;
		}
	}

	public static String getJobID() {
		String jobID="U"+System.nanoTime();
		return jobID;
	}
}

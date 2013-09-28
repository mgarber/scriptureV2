package nextgen.core.job;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;
import java.util.TreeMap;


import org.apache.log4j.Logger;

import broad.core.parser.StringParser;


/**
 * @author prussell
 *
 */
public class LSFJob implements Job {
	
	private Runtime runtime;
	private String jobId;
	private String cmmd;
	private String outFile;
	private String queue;
	private int memory;
	private static Logger logger = Logger.getLogger(LSFJob.class.getName());
	static int waitTime=60000; //1 minute
	
	/**
	 * @param jobID Job ID
	 */
	public LSFJob(String jobID) {
		this(jobID, null, null, 4);
	}

	/**
	 * @param jobID Job ID
	 * @param command Command
	 * @param outputFile LSF output file
	 */
	public LSFJob(String jobID, String command, String outputFile) {
		this(jobID, command, outputFile, 4);
	}
	
	/**
	 * @param jobID Job ID
	 * @param command Command
	 * @param outputFile LSF output file
	 * @param memoryRequest Memory request in Gb
	 */
	public LSFJob(String jobID, String command, String outputFile, int memoryRequest) {
		this(jobID, command, outputFile, "week", memoryRequest);
	}
	
	/**
	 * @param jobID Job ID
	 * @param command Command
	 * @param outputFile LSF output file
	 * @param queueName LSF queue name
	 */
	public LSFJob(String jobID, String command, String outputFile, String queueName) {
		this(jobID, command, outputFile, queueName, 4);
	}
	
	/**
	 * @param jobID Job ID
	 * @param command Command
	 * @param outputFile LSF output file
	 * @param queueName LSF queue name
	 * @param memoryRequest Memory request in Gb
	 */
	public LSFJob(String jobID, String command, String outputFile, String queueName, int memoryRequest) {
		this(Runtime.getRuntime(), jobID, command, outputFile, queueName, memoryRequest);
	}

	/**
	 * @param run Runtime
	 * @param command Command
	 */
	public LSFJob(Runtime run, String command) {
		this(run, command, "week", 4);
	}
	
	/**
	 * @param run Runtime
	 * @param command Command
	 * @param queueName LSF queue name
	 * @param memoryRequest Memory request in Gb
	 */
	public LSFJob(Runtime run, String command, String queueName, int memoryRequest) {
		runtime = run;
		jobId = generateJobID();
		cmmd = command;
		outFile = jobId + ".bsub.out";
		queue = queueName;
		memory = memoryRequest;
	}

	/**
	 * @param run Runtime
	 * @param jobID Job ID
	 * @param command Command
	 * @param outputFile LSF output file
	 * @param queueName LSF queue name
	 */
	public LSFJob(Runtime run, String jobID, String command, String outputFile, String queueName) {
		this(run, jobID, command, outputFile, queueName, 4);
	}

	
	/**
	 * @param run Runtime
	 * @param jobID Job ID
	 * @param command Command
	 * @param outputFile LSF output file
	 * @param queueName LSF queue name
	 * @param memoryRequest Memory request in Gb
	 */
	public LSFJob(Runtime run, String jobID, String command, String outputFile, String queueName, int memoryRequest) {
		runtime = run;
		jobId = jobID;
		cmmd = command;
		outFile = outputFile;
		queue = queueName;
		memory = memoryRequest;
	}
	
	@Override
	public String getID() {
		return jobId;
	}
	
	@Override
	public void submit() throws IOException, InterruptedException {
		@SuppressWarnings("unused")
		int exitVal = bsubProcess(runtime, jobId, cmmd, outFile, queue, memory);
	}

	@Override
	public void waitFor() throws IOException, InterruptedException {
		waitFor(waitTime);
	}

	@Override
	public void waitFor(int interval) throws IOException, InterruptedException {
		int running = 5;
		while(running>1){
			Thread.sleep(interval);
			Process blatProc = runtime.exec("bjobs -J " + jobId);
			BufferedReader out = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));
			running=LSFJob.parseReply(out);
		}
		LSFJob job = new LSFJob(jobId);
		if(job.failed()) {
			logger.warn("Job " + jobId + " failed.");
		}
	}
	
	@Override
	public boolean isPending() throws IOException, InterruptedException {
		return lsfStatus().equals("PEND");
	}

	@Override
	public boolean isRunning() throws IOException, InterruptedException {
		return lsfStatus().equals("RUN");
	}

	@Override
	public boolean completed() throws IOException, InterruptedException {
		return succeeded() || failed();
	}

	@Override
	public boolean succeeded() throws IOException, InterruptedException {
		return lsfStatus().equals("DONE");
	}

	@Override
	public boolean isSuspended() throws IOException, InterruptedException {
		return lsfStatus().contains("SUSP");
	}

	@Override
	public boolean failed() throws IOException {
		Process proc = runtime.exec("bjobs -a -J "+jobId);
		BufferedReader out = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String line = null;
		while((line = out.readLine()) != null) {
			//split the line
			String[] tokens=line.split(" ");
			String completionStatus=tokens[2];
			if(completionStatus.equalsIgnoreCase("EXIT") || line.contains("job killed")){
				logger.warn("Job "+ tokens[0] + " (" + jobId + ") FAILED");
				return true;
			}
		}
		return false;
	}

	@Override
	public void kill() throws IOException {
		@SuppressWarnings("unused")
		Process kill = runtime.exec("bkill " + jobId);
	}

	private String lsfStatus() throws IOException, InterruptedException {
		Map<String, String> allJobsStatus = allJobsStatus(runtime);
		if(!allJobsStatus.containsKey(jobId)) {
			throw new IllegalArgumentException("Job " + jobId + " does not have a job status.");
		}
		return allJobsStatus.get(jobId);
	}
	
	
	/**
	 * @return A unique job ID based on system time
	 */
	public static String generateJobID() {
		String jobID="U"+System.nanoTime();
		return jobID;
	}

	private static int parseReply(BufferedReader lsfOut) throws IOException {
		@SuppressWarnings("unused")
		String line = null;
		int i=0;
		while((line = lsfOut.readLine()) != null) {i++;}
		return i;
	}


	private static Map<String, String> allJobsStatus(Runtime run) throws IOException, InterruptedException {
		Process blatProc = run.exec("bjobs -aw");
		blatProc.waitFor();
		BufferedReader b = new BufferedReader(new InputStreamReader(blatProc.getInputStream()));
		StringParser s = new StringParser();
		Map<String, String> rtrn = new TreeMap<String, String>();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			String jobID = s.asString(6);
			String status = s.asString(2);
			rtrn.put(jobID, status);
		}
		return rtrn;
	}

	private static void writeError(InputStream errorStream) throws IOException {
		BufferedReader reader=	new BufferedReader(new InputStreamReader(errorStream));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
		}
		System.err.println();
	}


	private static int bsubProcess(Runtime run, String jobID, String command, String output, String queue, int memory) throws IOException, InterruptedException {
		String fullCommand="bsub -R rusage[mem=" + memory + "] -q " + queue+" -J "+jobID+ " -o "+output+" "+command;
		return bsubProcess(run, fullCommand);
	}

	private static int bsubProcess(Runtime run, String command) throws IOException, InterruptedException {
		//System.err.println("BSUB: " + command);
		Process p=run.exec(command);
		p.waitFor();
		int exitVal=p.exitValue();
		InputStream is = p.getInputStream();
		LSFJob.writeError(is);
		is.close();
		is = p.getErrorStream();
		LSFJob.writeError(is);
		is.close();
		//System.err.println("END BSUB " );
		
		return exitVal;
	}



}

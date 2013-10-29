package nextgen.core.job;

import java.io.File;
import java.io.IOException;
import java.util.List;

import nextgen.core.pipeline.OGSUtils;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.JobInfo;
import org.ggf.drmaa.JobTemplate;
import org.ggf.drmaa.Session;

/**
 * @author prussell
 *
 */
public class OGSJob implements Job {
	
	private Session session;
	private JobTemplate jobTemplate;
	private String jobID;
	File scriptFile;
	private static Logger logger = Logger.getLogger(OGSJob.class.getName());
	private boolean deleteScriptAfterSubmitting;
	private String command;
	
	/**
	 * @param drmaaSession DRMAA session
	 * @param cmmd Command to submit
	 * @throws DrmaaException
	 * @throws IOException
	 */
	public OGSJob(Session drmaaSession, String cmmd) throws DrmaaException, IOException {
		this(drmaaSession, OGSUtils.createScriptFile(cmmd), true);
		command = cmmd;
	}
	
	/**
	 * @param drmaaSession DRMAA session
	 * @param cmmd Command to submit
	 * @param jobName A simple name for the job e.g. to identify it in output files
	 * @throws DrmaaException
	 * @throws IOException
	 */
	public OGSJob(Session drmaaSession, String cmmd, String jobName) throws DrmaaException, IOException {
		this(drmaaSession, OGSUtils.createScriptFile(cmmd), true, jobName);
		command = cmmd;
	}


	/**
	 * @param drmaaSession DRMAA session
	 * @param cmmd Command to submit
	 * @param deleteScriptFileAfterSubmitting Delete the arbitrarily named script file after submitting job
	 * @throws DrmaaException
	 * @throws IOException
	 */
	public OGSJob(Session drmaaSession, String cmmd, boolean deleteScriptFileAfterSubmitting) throws DrmaaException, IOException {
		this(drmaaSession, OGSUtils.createScriptFile(cmmd), deleteScriptFileAfterSubmitting);
		command = cmmd;
	}

	/**
	 * @param drmaaSession DRMAA session
	 * @param cmmd Command to submit
	 * @param scriptFileName Script file to create with the command
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 * @throws IOException
	 */
	public OGSJob(Session drmaaSession, String cmmd, String scriptFileName, boolean deleteScriptFileAfterSubmitting) throws DrmaaException, IOException {
		this(drmaaSession, OGSUtils.createScriptFile(cmmd, scriptFileName), deleteScriptFileAfterSubmitting, null);
		command = cmmd;
	}

	/**
	 * @param drmaaSession DRMAA session
	 * @param cmmd Command to submit
	 * @param scriptFileName Script file to create with the command
	 * @param jobName A simple name for the job e.g. to identify it in output files
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 * @throws IOException
	 */
	public OGSJob(Session drmaaSession, String cmmd, boolean deleteScriptFileAfterSubmitting, String jobName) throws DrmaaException, IOException {
		this(drmaaSession, OGSUtils.createScriptFile(cmmd, null), deleteScriptFileAfterSubmitting, jobName);
		command = cmmd;
	}

	/**
	 * @param drmaaSession DRMAA session
	 * @param cmmd Command to submit
	 * @param scriptFileName Script file to create with the command
	 * @param jobName A simple name for the job e.g. to identify it in output files
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 * @throws IOException
	 */
	public OGSJob(Session drmaaSession, String cmmd, String scriptFileName, boolean deleteScriptFileAfterSubmitting, String jobName) throws DrmaaException, IOException {
		this(drmaaSession, OGSUtils.createScriptFile(cmmd, scriptFileName), deleteScriptFileAfterSubmitting, jobName);
		command = cmmd;
	}

	/**
	 * @param drmaaSession DRMAA session
	 * @param script Script file to submit
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 */
	public OGSJob(Session drmaaSession, File script, boolean deleteScriptFileAfterSubmitting) throws DrmaaException {
		this(drmaaSession, script, null, deleteScriptFileAfterSubmitting, null);
	}

	/**
	 * @param drmaaSession DRMAA session
	 * @param script Script file to submit
	 * @param jobName A simple name for the job e.g. to identify it in output files
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 */
	public OGSJob(Session drmaaSession, File script, boolean deleteScriptFileAfterSubmitting, String jobName) throws DrmaaException {
		this(drmaaSession, script, null, deleteScriptFileAfterSubmitting, jobName);
	}
	
	/**
	 * @param drmaaSession DRMAA session
	 * @param script Script file to submit
	 * @param args Arguments to script
	 * @param jobName A simple name for the job e.g. to identify it in output files
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 */
	public OGSJob(Session drmaaSession, File script, List<String> args, boolean deleteScriptFileAfterSubmitting, String jobName) throws DrmaaException {
		session = drmaaSession;
		
		// Set up the job template
		jobTemplate = session.createJobTemplate();
		String workingDirectory = System.getProperty("user.dir");
		jobTemplate.setWorkingDirectory(workingDirectory);
		String time = Long.valueOf(System.currentTimeMillis()).toString();
		String uniqueName = jobName == null ? time : jobName + "_" + time;
		jobTemplate.setJobName(uniqueName);
		jobTemplate.setErrorPath(":" + uniqueName + ".err");
		jobTemplate.setOutputPath(":" + uniqueName + ".out");
		
		scriptFile = script;
		jobTemplate.setRemoteCommand(scriptFile.getAbsolutePath());
		if(args != null) jobTemplate.setArgs(args);
		deleteScriptAfterSubmitting = deleteScriptFileAfterSubmitting;
	}
	
	/**
	 * @return The job template (will be null after job has been submitted)
	 */
	public JobTemplate getJobTemplate() {
		return jobTemplate;
	}
	
	@Override
	public String getID() {
		try {
			return jobID;
		} catch (NullPointerException e) {
			throw new IllegalStateException("OGS job doesn't have ID until submitted");
		}
	}

	@Override
	public void submit() throws IOException, InterruptedException, DrmaaException {
		jobID = session.runJob(jobTemplate);
		logger.info("Job submitted with system-assigned ID " + jobID + " and name " + jobTemplate.getJobName() + ".");
		if(command != null) {
			logger.info("Command: " + command);
		}
		session.deleteJobTemplate(jobTemplate);
		if(deleteScriptAfterSubmitting) OGSUtils.deleteScriptFile(scriptFile);
	}

	@Override
	public void waitFor() throws IOException, InterruptedException, DrmaaException {
		logger.info("Waiting for job " + jobID + "...");
		@SuppressWarnings("unused")
		JobInfo info = session.wait(jobID, Session.TIMEOUT_WAIT_FOREVER);
		logger.info("Done waiting for job.");
	}

	@Override
	public void waitFor(int interval) throws IOException, InterruptedException, DrmaaException {
		waitFor();
	}

	@Override
	public boolean isPending() throws IOException, InterruptedException, DrmaaException {
		int status = session.getJobProgramStatus(jobID);
		if(status == Session.QUEUED_ACTIVE) return true;
		if(status == Session.SYSTEM_ON_HOLD) return true;
		if(status == Session.USER_ON_HOLD) return true;
		if(status == Session.USER_SYSTEM_ON_HOLD) return true;
		return false;
	}

	@Override
	public boolean isRunning() throws IOException, InterruptedException, DrmaaException {
		int status = session.getJobProgramStatus(jobID);
		return status == Session.RUNNING;
	}

	@Override
	public boolean completed() throws IOException, InterruptedException, DrmaaException {
		return succeeded() || failed();
	}

	@Override
	public boolean succeeded() throws IOException, InterruptedException, DrmaaException {
		int status = session.getJobProgramStatus(jobID);
		return status == Session.DONE;
	}

	@Override
	public boolean failed() throws IOException, DrmaaException {
		int status = session.getJobProgramStatus(jobID);
		return status == Session.FAILED;
	}

	@Override
	public void kill() throws IOException, DrmaaException {
		session.control(jobID, Session.TERMINATE);
	}

	@Override
	public boolean isSuspended() throws IOException, InterruptedException, DrmaaException {
		int status = session.getJobProgramStatus(jobID);
		if(status == Session.SYSTEM_SUSPENDED) return true;
		if(status == Session.USER_SUSPENDED) return true;
		return false;
	}
	
}

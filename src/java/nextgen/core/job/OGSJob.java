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
	
	/**
	 * @param drmaaSession DRMAA session
	 * @param cmmd Command to submit
	 * @throws DrmaaException
	 * @throws IOException
	 */
	public OGSJob(Session drmaaSession, String cmmd) throws DrmaaException, IOException {
		this(drmaaSession, OGSUtils.createScriptFile(cmmd), true);
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
		this(drmaaSession, OGSUtils.createScriptFile(cmmd, scriptFileName), deleteScriptFileAfterSubmitting);
	}
	
	/**
	 * @param drmaaSession DRMAA session
	 * @param script Script file to submit
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 */
	public OGSJob(Session drmaaSession, File script, boolean deleteScriptFileAfterSubmitting) throws DrmaaException {
		this(drmaaSession, script, null, false);
	}
	
	/**
	 * @param drmaaSession DRMAA session
	 * @param script Script file to submit
	 * @param args Arguments to script
	 * @param deleteScriptFileAfterSubmitting Delete the script file after submitting job
	 * @throws DrmaaException
	 */
	public OGSJob(Session drmaaSession, File script, List<String> args, boolean deleteScriptFileAfterSubmitting) throws DrmaaException {
		session = drmaaSession;
		jobTemplate = session.createJobTemplate();
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
		return jobID;
	}

	@Override
	public void submit() throws IOException, InterruptedException, DrmaaException {
		jobID = session.runJob(jobTemplate);
		logger.info("Job submitted with ID " + jobID);
		session.deleteJobTemplate(jobTemplate);
		if(deleteScriptAfterSubmitting) OGSUtils.deleteScriptFile(scriptFile);
	}

	@Override
	public void waitFor() throws IOException, InterruptedException, DrmaaException {
		@SuppressWarnings("unused")
		JobInfo info = session.wait(jobID, Session.TIMEOUT_WAIT_FOREVER);
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

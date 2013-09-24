package nextgen.core.job;

import java.io.IOException;

import org.ggf.drmaa.DrmaaException;

/**
 * @author prussell
 *
 */
public interface Job {
	
	/**
	 * @return Job ID
	 */
	public String getID();
	
	/**
	 * Submit the job
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public void submit() throws IOException, InterruptedException, DrmaaException;
	
	/**
	 * Wait for the job to complete successfully or fail
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public void waitFor() throws IOException, InterruptedException, DrmaaException;
	
	/**
	 * Wait for the job, repeatedly checking at a specified time interval
	 * @param interval Time interval to check job status
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public void waitFor(int interval) throws IOException, InterruptedException, DrmaaException;
	
	/**
	 * @return Whether the job is actively pending
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public boolean isPending() throws IOException, InterruptedException, DrmaaException;
	
	/**
	 * @return Whether the job is actively pending
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public boolean isSuspended() throws IOException, InterruptedException, DrmaaException;

	/**
	 * @return Whether the job is actively running
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public boolean isRunning() throws IOException, InterruptedException, DrmaaException;
	
	/**
	 * @return Whether the job is done (includes aborted, completed successfully, failed)
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public boolean completed() throws IOException, InterruptedException, DrmaaException;
	
	/**
	 * @return Whether the job has run and completed successfully
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	public boolean succeeded() throws IOException, InterruptedException, DrmaaException;
	
	/**
	 * @return Whether the job has run and failed, or been aborted before running
	 * @throws IOException
	 * @throws DrmaaException 
	 */
	public boolean failed() throws IOException, DrmaaException;
	
	/**
	 * Abort or kill the job
	 * @throws IOException
	 * @throws DrmaaException 
	 */
	public void kill() throws IOException, DrmaaException;
	
}

package nextgen.core.pipeline;

import java.io.IOException;
import java.util.Collection;

import org.ggf.drmaa.DrmaaException;

/**
 * @author prussell
 *
 */
public class JobUtils {
	
	/**
	 * Wait for all jobs to complete (normally or fail)
	 * @param jobs Jobs to wait for
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public static void waitForAll(Collection<Job> jobs) throws IOException, InterruptedException, DrmaaException {
		for(Job job : jobs) {
			job.waitFor();
		}
	}
	
	/**
	 * Wait for all jobs to complete (normally or fail), repeatedly checking at a specified time interval
	 * @param jobs Jobs to wait for
	 * @param interval Time interval
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public static void waitForAll(Collection<Job> jobs, int interval) throws IOException, InterruptedException, DrmaaException {
		for(Job job : jobs) {
			job.waitFor(interval);
		}
	}
	
	/**
	 * @param jobs Jobs to check
	 * @return Whether all jobs have completed (normally or failed)
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public static boolean allCompleted(Collection<Job> jobs) throws IOException, InterruptedException, DrmaaException {
		for(Job job : jobs) {
			if(!job.completed()) return false;
		}
		return true;
	}
	
	/**
	 * @param jobs Jobs to check
	 * @return Whether all jobs have completed successfully
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException
	 */
	public static boolean allSucceeded(Collection<Job> jobs) throws IOException, InterruptedException, DrmaaException {
		for(Job job : jobs) {
			if(!job.succeeded()) return false;
		}
		return true;
	}
	
	/**
	 * @param jobs Jobs to kill
	 * @throws IOException
	 * @throws DrmaaException
	 */
	public static void killAll(Collection<Job> jobs) throws IOException, DrmaaException {
		for(Job job : jobs) {
			job.kill();
		}
	}
	
}

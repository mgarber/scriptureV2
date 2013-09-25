package nextgen.core.job;

import java.io.IOException;
import java.util.Collection;

public class JobUtils {
	
	public static void waitForAll(Collection<Job> jobs) throws IOException, InterruptedException {
		for(Job job : jobs) {
			job.waitFor();
		}
	}
	
	public static void waitForAll(Collection<Job> jobs, int interval) throws IOException, InterruptedException {
		for(Job job : jobs) {
			job.waitFor(interval);
		}
	}
	
	public static boolean allCompleted(Collection<Job> jobs) throws IOException, InterruptedException {
		for(Job job : jobs) {
			if(!job.completed()) return false;
		}
		return true;
	}
	
	public static boolean allSucceeded(Collection<Job> jobs) throws IOException, InterruptedException {
		for(Job job : jobs) {
			if(!job.succeeded()) return false;
		}
		return true;
	}
	
	public static void killAll(Collection<Job> jobs) throws IOException {
		for(Job job : jobs) {
			job.kill();
		}
	}
	
}

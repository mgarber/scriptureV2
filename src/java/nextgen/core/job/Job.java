package nextgen.core.job;

import java.io.IOException;
import java.util.Collection;

public interface Job {
	
	public String getID();
	public void submit() throws IOException, InterruptedException;
	public void waitFor() throws IOException, InterruptedException;
	public void waitFor(int interval) throws IOException, InterruptedException;
	public boolean isPending() throws IOException, InterruptedException;
	public boolean isRunning() throws IOException, InterruptedException;
	public boolean completed() throws IOException, InterruptedException;
	public boolean succeeded() throws IOException, InterruptedException;
	public boolean failed() throws IOException;
	public void kill() throws IOException;
	
}

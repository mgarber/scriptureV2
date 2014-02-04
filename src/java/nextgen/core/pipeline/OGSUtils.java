package nextgen.core.pipeline;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;
import org.ggf.drmaa.SessionFactory;

/**
 * @author prussell
 *
 */
public class OGSUtils {
	
	private static Logger logger = Logger.getLogger(OGSUtils.class.getName());
	
	/**
	 * Create a script file to run one command, with an arbitrary name based on the system time
	 * @param cmmd Command for the script to run
	 * @return The file
	 * @throws IOException 
	 */
	public static File createScriptFile(String cmmd) throws IOException {
		return createScriptFile(cmmd, null);
	}
	
	/**
	 * Create a script file to run one command
	 * @param cmmd Command for the script to run
	 * @param fileName The file name
	 * @return The file
	 * @throws IOException
	 */
	public static File createScriptFile(String cmmd, String fileName) throws IOException {
		String f = fileName == null ? "script_" + System.currentTimeMillis() : fileName;
		FileWriter w = new FileWriter(f);
		w.write(cmmd + "\n");
		w.close();
		File file = new File(f);
		file.setExecutable(true, false);
		return file;
	}
	
	/**
	 * Delete a file
	 * @param file File to delete
	 */
	public static void deleteScriptFile(File file) {
		boolean deleted = file.delete();
		if(!deleted) {
			logger.warn("Couldn't delete file " + file.getName() + ".");
		}
	}

	/**
	 * Get a DRMAA session
	 * Warning: there should only be one active DRMAA session at a time
	 * @return A DRMAA session
	 * @throws DrmaaException 
	 */
	public static Session getDrmaaSession() throws DrmaaException {
		Scheduler.logger.warn("Starting a DRMAA session.");
		Scheduler.logger.warn("There should only be one active DRMAA session at a time.");
		Session rtrn = SessionFactory.getFactory().getSession();
		rtrn.init(null);
		OGSUtils.attachShutDownHook(rtrn);
		return rtrn;
	}

	/**
	 * Add a shutdown hook that will close the DRMAA session upon program termination
	 * This should be added to every main method that will use a Scheduler
	 * @param drmaaSession 
	 * @throws DrmaaException
	 */
	public static void attachShutDownHook(final Session drmaaSession) throws DrmaaException {
		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				try {
					Scheduler.logger.info("Ending DRMAA session");
					drmaaSession.exit();
				} catch (DrmaaException e) {
					e.printStackTrace();
					Scheduler.logger.warn("DRMAA session might not be closed");
				}
			}
		});
		Scheduler.logger.info("Attached shutdown hook to close DRMAA session upon JVM exit.");
	}
	
}

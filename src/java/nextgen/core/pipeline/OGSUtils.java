package nextgen.core.pipeline;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author prussell
 *
 */
public class OGSUtils {
	
	/**
	 * Create a script file to run one command, with an arbitrary name based on the system time
	 * @param cmmd Command for the script to run
	 * @return The file
	 * @throws IOException 
	 */
	public static File createScriptFile(String cmmd) throws IOException {
		String fileName = "script_" + System.currentTimeMillis();
		return createScriptFile(cmmd, fileName);
	}
	
	/**
	 * Create a script file to run one command
	 * @param cmmd Command for the script to run
	 * @param fileName The file name
	 * @return The file
	 * @throws IOException
	 */
	public static File createScriptFile(String cmmd, String fileName) throws IOException {
		FileWriter w = new FileWriter(fileName);
		w.write(cmmd + "\n");
		w.close();
		File file = new File(fileName);
		file.setExecutable(true, false);
		return file;
	}
	
	/**
	 * Delete a file
	 * @param file File to delete
	 */
	public static void deleteScriptFile(File file) {
		@SuppressWarnings("unused")
		boolean deleted = file.delete();
	}
	
}

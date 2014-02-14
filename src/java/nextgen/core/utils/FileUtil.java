/**
 * 
 */
package nextgen.core.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author prussell
 *
 */
public final class FileUtil {

	/**
	 * Get the lines of a file as a list
	 * @param fileName File name
	 * @return The lines of the file in a list of strings
	 * @throws IOException
	 */
	public static List<String> fileLinesAsList(String fileName) throws IOException {
		List<String> rtrn = new ArrayList<String>();
		FileReader r = new FileReader(fileName);
		BufferedReader b = new BufferedReader(r);
		while(b.ready()) {
			String line = b.readLine();
			if(line.length() > 0) {
				rtrn.add(line);
			}
		}
		r.close();
		b.close();
		return rtrn;
	}

}

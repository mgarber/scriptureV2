package nextgen.core.readers;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

import org.apache.log4j.Logger;

import broad.core.parser.StringParser;

import nextgen.core.annotation.PSLRecord;

/**
 * Iterate through a PSL file
 * @author prussell
 *
 */
public class PSLReader implements Iterator<PSLRecord> {
	
	private BufferedReader reader;
	private static Logger logger = Logger.getLogger(PSLReader.class.getName());
	
	public PSLReader(String pslFile) throws FileNotFoundException {
		FileReader r = new FileReader(pslFile);
		reader = new BufferedReader(r);
	}
	
	public void close() {
		try {
			reader.close();
		} catch (IOException e) {
			logger.warn("Couldn't close buffered reader");
			e.printStackTrace();
		}
	}
	
	@Override
	public boolean hasNext() {
		try {
			boolean rtrn = reader.ready();
			if(!rtrn) {
				close();
			}
			return rtrn;
		} catch(IOException e) {
			logger.warn("Couldn't check if buffered reader is ready. Returning false for hasNext().");
			e.printStackTrace();
			close();
			return false;
		}
	}

	@Override
	public PSLRecord next() {
		if(!hasNext()) {
			return null;
		}
		try {
			while(true) {
				String line = reader.readLine();
				StringParser s = new StringParser();
				s.parse(line);
				String[] tokens = s.getStringArray();
				if(tokens == null) {
					continue;
				}
				if(tokens.length != 21) {
					continue;
				}
				return new PSLRecord.Factory().create(tokens);
			}
		} catch (IOException e) {
			logger.warn("Couldn't read from file");
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not implemented");
	}

}

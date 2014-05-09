package nextgen.core.berkeleydb;

import java.io.File;

import org.apache.log4j.Logger;

import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;

/*
 * http://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/persist_first.html
 */
/**
 * A database environment
 * @author prussell
 *
 */
public class DatabaseEnvironment {
	
	private Environment environment;
	private static Logger logger = Logger.getLogger(DatabaseEnvironment.class.getName());
	
	public DatabaseEnvironment() {}
	
	/**
	 * Set up the environment
	 * @param home Environment home directory
	 * @param readOnly Whether the environment should be read only
	 * @param transactional Whether the environment should be transactional (see berkeley db documentation)
	 */
	public void setup(File home, boolean readOnly, boolean transactional) {
		
		EnvironmentConfig config = new EnvironmentConfig();
		config.setReadOnly(readOnly);
		config.setAllowCreate(!readOnly);
		config.setTransactional(transactional);
		
		environment = new Environment(home, config);
		
		printCurrentProperties();
	}
	
	public Environment getEnvironment() {
		return environment;
	}
	
	/**
	 * Close the environment
	 */
	public void close() {
		if(environment != null) {
			environment.close();
		}
	}
	
	public void printCurrentProperties() {
		EnvironmentConfig config = environment.getConfig();
		String r = config.getReadOnly() ? "is" : "is not";
		String t = config.getTransactional() ? "is" : "is not";
		logger.info("");
		logger.info("Cache size is " + config.getCacheSize() + ". Database " + r + " read only. Database " + t + " transactional.");
	}
	
}

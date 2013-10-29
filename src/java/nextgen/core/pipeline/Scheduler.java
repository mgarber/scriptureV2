package nextgen.core.pipeline;

import org.apache.log4j.Logger;

/**
 * @author prussell
 *
 */
public enum Scheduler {
	
	/**
	 * Load Sharing Facility
	 */
	LSF,
	/**
	 * Open Grid Scheduler (formerly Sun Grid Engine)
	 */
	OGS;
	
	/**
	 * DRMAA session
	 * Only one session should be open at a time
	 * The session must be closed upon exiting a program
	 * Add a call to attachShutDownHook() to all main methods using a Scheduler - this will close the session upon exit
	 */
	
	static Logger logger = Logger.getLogger(Scheduler.class.getName());
	
	/**
	 * Get a comma separated list of all the scheduler names
	 * @return List of all scheduler names
	 */
	public static String getCommaSeparatedList() {
		if(Scheduler.values().length == 0) {
			throw new IllegalStateException("No values");
		}
		String rtrn = Scheduler.values()[0].toString();
		for(int i = 1; i < Scheduler.values().length; i++) {
			rtrn += ", " + Scheduler.values()[i].toString();
		}
		return rtrn;
	}
	
	@Override
	public String toString() {
		switch(this) {
		case LSF: return "LSF";
		case OGS: return "OGS";
		default: throw new IllegalArgumentException("Not implemented.");
		}
	}
	
	/**
	 * @param name Name
	 * @return Scheduler
	 */
	public static Scheduler fromString(String name) {
		if(name.equals("LSF")) {
			return LSF;
		}
		if(name.equals("OGS")) {
			return OGS;
		}
		throw new IllegalArgumentException("Scheduler name " + name + " not recognized.");
	}
	
}

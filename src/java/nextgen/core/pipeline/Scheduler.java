package nextgen.core.pipeline;

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

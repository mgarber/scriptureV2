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
	 * Sun Grid Engine
	 */
	SGE;
	
	@Override
	public String toString() {
		switch(this) {
		case LSF: return "LSF";
		case SGE: return "SGE";
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
		if(name.equals("SGE")) {
			return SGE;
		}
		throw new IllegalArgumentException("Scheduler name " + name + " not recognized.");
	}
	
}

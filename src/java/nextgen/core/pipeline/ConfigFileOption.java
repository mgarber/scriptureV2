package nextgen.core.pipeline;

import java.util.TreeSet;

/**
 * @author prussell
 * An option in a config file
 */
public class ConfigFileOption {
	
	private String flag;
	private boolean required;
	private TreeSet<Integer> numVals;
	private boolean repeatable;
	private String defaultVal;
	
	/**
	 * @param optionFlag Flag
	 * @param numValues Number of values including flag
	 * @param fewerValuesOK Can have any number of values between 1 and the specified number
	 * @param isRepeatable Whether the option can be specified multiple times in the file
	 * @param isRequired Whether the option is required for the section it's in
	 */
	public ConfigFileOption(String optionFlag, int numValues, boolean fewerValuesOK, boolean isRepeatable, boolean isRequired) {
		this(optionFlag, numValues, fewerValuesOK, isRepeatable, isRequired, null);
	}
	
	/**
	 * @param optionFlag Flag
	 * @param numValues Number of values including flag
	 * @param fewerValuesOK Can have any number of values between 1 and the specified number
	 * @param isRepeatable Whether the option can be specified multiple times in the file
	 * @param isRequired Whether the option is required for the section it's in
	 * @param defaultValue Default value not including the flag itself
	 */
	public ConfigFileOption(String optionFlag, int numValues, boolean fewerValuesOK, boolean isRepeatable, boolean isRequired, String defaultValue) {
		defaultVal = defaultValue;
		flag = optionFlag;
		required = isRequired;
		repeatable = isRepeatable;
		if(isRequired && hasDefault()) {
			throw new IllegalArgumentException("Required option " + flag + " cannot specify a default value.");
		}
		numVals = new TreeSet<Integer>();
		if(!fewerValuesOK) {
			numVals.add(Integer.valueOf(numValues));
		} else {
			for(int i = 1; i <= numValues; i++) {
				numVals.add(Integer.valueOf(i));
			}
		}
	}
	
	/**
	 * Whether the option has a default value
	 * @return True iff the option has a default value specified
	 */
	protected boolean hasDefault() {
		return defaultVal != null;
	}
	
	protected ConfigFileOptionValue getDefaultValue() {
		if(defaultVal == null) {
			return null;
		}
		return new ConfigFileOptionValue(this, flag + " " + defaultVal);
	}
	
	/**
	 * Get option flag
	 * @return Flag
	 */
	public String getName() {
		return flag;
	}
	
	protected boolean isRequired() {
		return required;
	}
	
	protected TreeSet<Integer> getAllowableNumbersOfValues() {
		return numVals;
	}
	
	protected boolean isRepeatable() {
		return repeatable;
	}
	
	@Override
	public String toString() {
		String rtrn = flag + "_";
		rtrn += numVals + "_values_";
		if(required) {
			rtrn += "required_";
		} else {
			rtrn += "not_required_";
		}
		if(repeatable) {
			rtrn += "repeatable";
		} else {
			rtrn += "not_repeatable";
		}
		return rtrn;
	}
	
	@Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) {
			return false;
		}
		ConfigFileOption c = (ConfigFileOption)o;
		return c.toString().equals(toString());
	}
	
	@Override
	public int hashCode() {
		return toString().hashCode();
	}
	
}
package nextgen.core.pipeline;

import broad.core.parser.StringParser;

/**
 * @author prussell
 *
 */
public class ConfigFileOptionValue {
	
	private ConfigFileOption option;
	private StringParser stringParser;
	private String line;
	
	protected ConfigFileOptionValue(ConfigFileOption op, String fileLine) {
		this(op, fileLine, true);
	}
	
	protected ConfigFileOptionValue(ConfigFileOption op, String fileLine, boolean validate) {
		line = fileLine;
		option = op;
		stringParser = new StringParser();
		stringParser.parse(fileLine);
		if(validate && !option.getAllowableNumbersOfValues().contains(Integer.valueOf(stringParser.getFieldCount()))) {
			throw new IllegalArgumentException("Line must have " + option.getAllowableNumbersOfValues() + " fields: " + fileLine);
		}
	}
	
	/**
	 * Get number of fields in line including the flag itself
	 * @return The number of fields in option line
	 */
	public int getActualNumValues() {
		return stringParser.getFieldCount();
	}
	
	/**
	 * Get string value for field
	 * @param fieldNumber Field number where the flag itself is field 0
	 * @return The value as a string
	 */
	public String asString(int fieldNumber) {
		return stringParser.asString(fieldNumber);
	}
	
	/**
	 * Get int value for field
	 * @param fieldNumber Field number where the flag itself is field 0
	 * @return The value as an int
	 */
	public int asInt(int fieldNumber) {
		return stringParser.asInt(fieldNumber);
	}
	
	/**
	 * Get double value for field
	 * @param fieldNumber Field number where the flag itself is field 0
	 * @return The value as a double
	 */
	public double asDouble(int fieldNumber) {
		return stringParser.asDouble(fieldNumber);
	}
	
	/**
	 * Get the actual line from the config file without the flag field
	 * @return The line minus the flag itself
	 */
	public String getLineMinusFlag() {
		return getLastFields(1);
	}
	
	/**
	 * Get the last fields of actual line from config file
	 * @param firstFieldToGet First field to include
	 * @return A space separated string starting with the first field and continuing to the end of config file line, or empty string if fewer fields
	 */
	public String getLastFields(int firstFieldToGet) {
		String rtrn = "";
		for(int i = firstFieldToGet; i < stringParser.getFieldCount(); i++) {
			rtrn += stringParser.asString(i);
			if(i < stringParser.getFieldCount() - 1) {
				rtrn += " ";
			}
		}
		return rtrn;		
	}
	
	/**
	 * Get the flag (first field on line)
	 * @return The flag
	 */
	public String getFlag() {
		return stringParser.asString(0);
	}
	
	/**
	 * Get the actual line from the config file
	 * @return Config file line specifying the option and value
	 */
	public String getFullOptionLine() {
		return line;
	}
			
	@Override
	public String toString() {
		return option.toString() + "_" + line;
	}
	
	@Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) {
			return false;
		}
		ConfigFileOptionValue c = (ConfigFileOptionValue)o;
		return c.toString().equals(toString());
	}
	
	@Override
	public int hashCode() {
		return toString().hashCode();
	}

	
}
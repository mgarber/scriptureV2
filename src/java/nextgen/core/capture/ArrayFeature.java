package nextgen.core.capture;

import nextgen.core.pipeline.ConfigFileOptionValue;

/**
 * @author prussell
 * An element of the array that can be specified in the config file
 */
public interface ArrayFeature {
	
	/**
	 * @return Name of feature
	 */
	public String name();
	
	/**
	 * @return Description of a valid config file line to specify parameters of this feature
	 */
	public String configFileLineDescription();
	
	/**
	 * @param value Config file value
	 * @return Whether the value is a valid specification of this feature
	 */
	public boolean validConfigFileValue(ConfigFileOptionValue value);
	
	/**
	 * @param value Config file value
	 */
	public void setParametersFromConfigFile(ConfigFileOptionValue value);
	

}

package nextgen.core.pipeline;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

/**
 * @author prussell
 * A section of a config file indicated by a header string
 */
public class ConfigFileSection {
	
	private String name;
	private Collection<ConfigFileOption> allowableOptions;
	private Map<String, ConfigFileOption> allowableOptionsByName;
	private boolean required;
	
	/**
	 * @param sectionName Section name used in header
	 * @param isRequired Whether the section is required
	 */
	public ConfigFileSection(String sectionName, boolean isRequired) {
		required = isRequired;
		name = sectionName;
		allowableOptions = new ArrayList<ConfigFileOption>();
		allowableOptionsByName = new TreeMap<String, ConfigFileOption>();
	}
	
	/**
	 * Whether the section is required
	 * @return True iff the section is required
	 */
	public boolean isRequired() {
		return required;
	}

	/**
	 * @param sectionName Section name
	 * @param optionCollection Collection of options to instantiate with
	 * @param isRequired Whether the section is required
	 */
	public ConfigFileSection(String sectionName, Collection<ConfigFileOption> optionCollection, boolean isRequired) {
		this(sectionName, isRequired);
		addAllowableOptions(optionCollection);
	}
	
	/**
	 * Get section name
	 * @return Section name
	 */
	public String getName() {
		return name;
	}
	
	protected Collection<ConfigFileOption> getAllowableOptions() {
		return allowableOptions;
	}
	
	protected ConfigFileOption getAllowableOption(String optionName) {
		return getAllowableOption(optionName, true);
	}
	
	protected ConfigFileOption getAllowableOption(String optionName, boolean validate) {
		if(validate && !allowableOptionsByName.containsKey(optionName)) {
			throw new IllegalArgumentException("File section " + name + " does not have allowable option " + optionName + ".");
		}
		return allowableOptionsByName.get(optionName);
	}
	
	/**
	 * Add an option
	 * @param option Option to add
	 */
	public void addAllowableOption(ConfigFileOption option) {
		allowableOptions.add(option);
		allowableOptionsByName.put(option.getName(), option);
	}
	
	/**
	 * Add a collection of options
	 * @param optionCollection Collection of options
	 */
	public void addAllowableOptions(Collection<ConfigFileOption> optionCollection) {
		for(ConfigFileOption option : optionCollection) {
			addAllowableOption(option);
		}
	}
	
	@Override
	public String toString() {
		String rtrn = name + "_";
		if(required) rtrn += "required_";
		else rtrn += "not_required_";
		for(ConfigFileOption option : allowableOptions) {
			rtrn += option.toString() + "_";
		}
		return rtrn;
	}
	
	@Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) {
			return false;
		}
		ConfigFileSection c = (ConfigFileSection)o;
		return c.toString().equals(toString());
	}
	
	@Override
	public int hashCode() {
		return toString().hashCode();
	}

	
}
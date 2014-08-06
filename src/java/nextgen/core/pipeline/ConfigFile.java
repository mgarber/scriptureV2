/**
 * 
 */
package nextgen.core.pipeline;

import broad.core.parser.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;


/**
 * @author prussell
 *
 */
public class ConfigFile {
	
	
	/**
	 * Comment character
	 */
	public static char COMMENT_CHAR = '#';
	
	/**
	 * Character that begins each section header
	 */
	public static char SECTION_HEADER_CHAR = ':';
	
	private ArrayList<ConfigFileSection> fileSections;
	private Map<String, ConfigFileSection> allowableSectionsByName;
	private Map<ConfigFileSection, Map<ConfigFileOption, Collection<ConfigFileOptionValue>>> optionsFromFile;
	private Map<ConfigFileSection, List<OptionValuePair>> orderedOptionsBySection;
	private static Logger logger = Logger.getLogger(ConfigFile.class.getName());
	
	/**
	 * Instantiate config file and parse file
	 * @param sections Sections
	 * @param fileName File name
	 * @throws IOException 
	 */
	public ConfigFile(Collection<ConfigFileSection> sections, String fileName) throws IOException {
		fileSections = new ArrayList<ConfigFileSection>();
		fileSections.addAll(sections);
		allowableSectionsByName = new TreeMap<String, ConfigFileSection>();
		for(ConfigFileSection section : fileSections) {
			allowableSectionsByName.put(section.getName(), section);
		}
		parseFile(fileName);
		validateFileAndGetDefaults();
	}
		
	/**
	 * Whether the file specifies the option in the indicated section
	 * @param section The section
	 * @param option The option
	 * @return Whether the option is specified in the section
	 */
	public boolean hasOption(ConfigFileSection section, ConfigFileOption option) {
		return optionsFromFile.get(section).containsKey(option);
	}
	
	/**
	 * Whether the file includes the section
	 * @param section The section
	 * @return Whether the section is specified in the file
	 */
	public boolean hasSection(ConfigFileSection section) {
		return optionsFromFile.containsKey(section);
	}
	
	/**
	 * Get string values for all the config file option values
	 * @param values Config file option value objects
	 * @param fieldNumber Field number to use from each value object
	 * @return This field from each option as a string
	 */
	public static Collection<String> valuesAsStrings(Collection<ConfigFileOptionValue> values, int fieldNumber) {
		Collection<String> rtrn = new ArrayList<String>();
		for(ConfigFileOptionValue value : values) {
			rtrn.add(value.asString(fieldNumber));
		}
		return rtrn;
	}
	
	/**
	 * Get all values for the option in the file
	 * @param section File section
	 * @param option The option
	 * @return The values for the option in the section or null if option is not specified in the section
	 */
	public Collection<ConfigFileOptionValue> getOptionValues(ConfigFileSection section, ConfigFileOption option) {
		if(!optionsFromFile.get(section).containsKey(option)) {
			return null;
		}
		return optionsFromFile.get(section).get(option);
	}
	
	/**
	 * Get the value for an option that is specified at most once in the file
	 * @param section File section
	 * @param option The option
	 * @return The value for the option in the section or null if option is not specified in the section
	 */
	public ConfigFileOptionValue getSingleValue(ConfigFileSection section, ConfigFileOption option) {
		if(!optionsFromFile.get(section).containsKey(option)) {
			return null;
		}
		if(optionsFromFile.get(section).get(option).size() > 1) {
			throw new IllegalArgumentException("Can't get single value for option " + option.getName() + ": option is specified more than once in config file.");
		}
		return optionsFromFile.get(section).get(option).iterator().next();
	}
	
	/**
	 * Get the only value field as a string for an option that is specified at most once in the file and has one value field not counting the flag itself
	 * @param section File section
	 * @param option The option
	 * @return The single value for the singleton option or null if file does not have the value
	 */
	public String getSingleValueString(ConfigFileSection section, ConfigFileOption option) {
		ConfigFileOptionValue value = getSingleValue(section, option);
		if(value == null) {
			return null;
		}
		if(value.getActualNumValues() != 2) {
			throw new IllegalArgumentException("Can't get single value for option " + option.getName() + ": option must specify 2 fields including flag.");
		}
		return value.asString(1);
	}
	
	/**
	 * Get the only value field as a double for an option that is specified at most once in the file and has one value field not counting the flag itself
	 * @param section File section
	 * @param option The option
	 * @return The single value for the singleton option
	 */
	public double getSingleValueDouble(ConfigFileSection section, ConfigFileOption option) {
		ConfigFileOptionValue value = getSingleValue(section, option);
		if(value == null) {
			throw new IllegalArgumentException("File does not have option " + option.getName());
		}
		if(value.getActualNumValues() != 2) {
			throw new IllegalArgumentException("Can't get single value for option " + option.getName() + ": option must specify 2 fields including flag.");
		}
		return value.asDouble(1);
	}
	
	/**
	 * Get the only value field as an int for an option that is specified at most once in the file and has one value field not counting the flag itself
	 * @param section File section
	 * @param option The option
	 * @return The single value for the singleton option
	 */
	public int getSingleValueInt(ConfigFileSection section, ConfigFileOption option) {
		ConfigFileOptionValue value = getSingleValue(section, option);
		if(value == null) {
			throw new IllegalArgumentException("File does not have option " + option.getName());
		}
		if(value.getActualNumValues() != 2) {
			throw new IllegalArgumentException("Can't get single value for option " + option.getName() + ": option must specify 2 fields including flag.");
		}
		return value.asInt(1);
	}
	
	/**
	 * Get config file options and values in the order listed in a section
	 * @param section The section
	 * @return List of options and values
	 */
	public List<OptionValuePair> getOrderedOptionsAndValues(ConfigFileSection section) {
		if(!orderedOptionsBySection.containsKey(section)) {
			throw new IllegalArgumentException("Section not found: " + section.getName());
		}
		return orderedOptionsBySection.get(section);
	}
	
	
	/**
	 * Get config file values in the order listed in a section
	 * @param section The section
	 * @return List of values
	 */
	public List<ConfigFileOptionValue> getOrderedValues(ConfigFileSection section) {
		List<OptionValuePair> pairs = getOrderedOptionsAndValues(section);
		List<ConfigFileOptionValue> rtrn = new ArrayList<ConfigFileOptionValue>();
		for(OptionValuePair pair : pairs) {
			rtrn.add(pair.value());
		}
		return rtrn;
	}
	
	private static boolean isSectionHeaderLine(String fileLine) {
		if(fileLine.length() < 1) return false;
		return fileLine.charAt(0) == SECTION_HEADER_CHAR;
	}
	
	private static boolean isComment(String fileLine) {
		if(fileLine.length() < 1) return false;
		return fileLine.charAt(0) == COMMENT_CHAR;
	}
	
	@SuppressWarnings("unused")
	private ConfigFileSection getSection(String sectionName) {
		return getSection(sectionName, true);
	}
	
	private ConfigFileSection getSection(String sectionName, boolean validate) {
		if(validate && !allowableSectionsByName.containsKey(sectionName)) {
			exitWithMessageAndHelpMenu("Section " + sectionName + " not recognized.");
		}
		return allowableSectionsByName.get(sectionName);
	}
	
	private void parseFile(String file) throws IOException {
		logger.info("Reading config file " + file + "...");
		try {
			optionsFromFile = new HashMap<ConfigFileSection, Map<ConfigFileOption, Collection<ConfigFileOptionValue>>>();
			orderedOptionsBySection = new HashMap<ConfigFileSection, List<OptionValuePair>>();
	
			FileReader r = new FileReader(file);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			boolean started = false;
			Map<ConfigFileOption, Collection<ConfigFileOptionValue>> currentSectionOptions = new HashMap<ConfigFileOption, Collection<ConfigFileOptionValue>>();
			List<OptionValuePair> currentSectionOptionValuePairs = new ArrayList<OptionValuePair>();
			ArrayList<String> currentSectionName = new ArrayList<String>();
			
			while(b.ready()) {
				String line = b.readLine();
				s.parse(line);
				// Skip blank lines
				if(s.getFieldCount() == 0 || line.isEmpty()) {
					s.clear();
					continue;
				}
				// Skip comments
				if(isComment(line)) {
					s.clear();
					continue;
				}
				if(isSectionHeaderLine(line)) {
					// Encounter new section header
					String newSectionName = line.substring(1);
					logger.info("Reading section " + newSectionName);
					if(started) {
						// We have already been reading a previous section
						Map<ConfigFileOption, Collection<ConfigFileOptionValue>> tmpOptions = new HashMap<ConfigFileOption, Collection<ConfigFileOptionValue>>();
						List<OptionValuePair> tmpOptionValuePairs = new ArrayList<OptionValuePair>();
						// Store the options from the previous section
						tmpOptions.putAll(currentSectionOptions);
						tmpOptionValuePairs.addAll(currentSectionOptionValuePairs);
						currentSectionOptions.clear();
						currentSectionOptionValuePairs.clear();
						ConfigFileSection section = allowableSectionsByName.get(currentSectionName.get(0));
						if(optionsFromFile.containsKey(section)) {
							exitWithMessageAndHelpMenu("Section is specified more than once: " + section.getName());
						}
						currentSectionName.clear();
						currentSectionName.add(newSectionName);
						optionsFromFile.put(section, tmpOptions);
						orderedOptionsBySection.put(section, tmpOptionValuePairs);
					} else {
						// This is the first section header
						currentSectionName.clear();
						currentSectionName.add(newSectionName);
						started = true;
					}
					continue;
				}
				if(currentSectionName.isEmpty()) {
					// This is an option line but there was no section header
					exitWithMessageAndHelpMenu("First line of file must be section header.");
				}
				// Make option and value from line
				String optionName = s.asString(0);
				ConfigFileSection section = getSection(currentSectionName.get(0), true);
				ConfigFileOption option = section.getAllowableOption(optionName);
				if(!currentSectionOptions.containsKey(option)) {
					Collection<ConfigFileOptionValue> tmpCol = new ArrayList<ConfigFileOptionValue>();
					currentSectionOptions.put(option, tmpCol);
				} else {
					if(!option.isRepeatable()) {
						exitWithMessageAndHelpMenu("Option " + optionName + " in section " + section.getName() + " is not repeatable.");
					}
				}
				// Store option and value
				ConfigFileOptionValue val = new ConfigFileOptionValue(option, line, false);
				currentSectionOptions.get(option).add(val);
				currentSectionOptionValuePairs.add(new OptionValuePair(option, val));
			}
			
			// Store the last section
			if(started) {
				Map<ConfigFileOption, Collection<ConfigFileOptionValue>> tmpOptions = new HashMap<ConfigFileOption, Collection<ConfigFileOptionValue>>();
				List<OptionValuePair> tmpOptionValuePairs = new ArrayList<OptionValuePair>();
				tmpOptions.putAll(currentSectionOptions);
				tmpOptionValuePairs.addAll(currentSectionOptionValuePairs);
				currentSectionOptions.clear();
				currentSectionOptionValuePairs.clear();
				ConfigFileSection section = allowableSectionsByName.get(currentSectionName.get(0));
				if(optionsFromFile.containsKey(section)) {
					exitWithMessageAndHelpMenu("Section is specified more than once: " + section.getName());
				}
				optionsFromFile.put(section, tmpOptions);
				orderedOptionsBySection.put(section, tmpOptionValuePairs);
			} 
	
			r.close();
			b.close();
			logger.info("Successfully read config file.");
		} catch(IllegalArgumentException e) {
			e.printStackTrace();
			exitWithHelpMenu();
		}
	}
	
	private void validateFileAndGetDefaults() {
		
		// Make sure sections and options are recognized
		for(ConfigFileSection section : optionsFromFile.keySet()) {
			String sectionName = section.getName();
			if(!allowableSectionsByName.containsKey(sectionName)) {
				exitWithMessageAndHelpMenu("Section name " + sectionName + " not recognized.");
			}
			for(ConfigFileOption option : optionsFromFile.get(section).keySet()) {
				if(!allowableSectionsByName.get(sectionName).getAllowableOptions().contains(option)) {
					exitWithMessageAndHelpMenu("Option name " + option.getName() + " in section " + sectionName + " not recognized.");
				}
				if(optionsFromFile.get(section).get(option).size() > 1 && !option.isRepeatable()) {
					exitWithMessageAndHelpMenu("Option " + option.getName() + " in section " + sectionName + " is not repeatable.");
				}
			}
		}
		
		// Look for required sections and options
		for(String sectionName : allowableSectionsByName.keySet()) {
			ConfigFileSection section = allowableSectionsByName.get(sectionName);
			if(!optionsFromFile.containsKey(section)) {
				if(section.isRequired()) {
					exitWithMessageAndHelpMenu("Required section " + sectionName + " is missing.");
				} else {
					continue;
					//optionsFromFile.put(section, new HashMap<ConfigFileOption, Collection<ConfigFileOptionValue>>());
				}
			}
			for(ConfigFileOption option : section.getAllowableOptions()) {
				if(!optionsFromFile.get(section).containsKey(option)) {
					if(option.isRequired()) {
						exitWithMessageAndHelpMenu("Required option " + option.getName() + " in section " + sectionName + " is missing.");
					} else {
						Collection<ConfigFileOptionValue> tmpCol = new ArrayList<ConfigFileOptionValue>();
						if(option.getDefaultValue() != null) {
							tmpCol.add(option.getDefaultValue());
							optionsFromFile.get(section).put(option, tmpCol);
						}
					}
				}
				if(optionsFromFile.get(section).containsKey(option)) {
					for(ConfigFileOptionValue value : optionsFromFile.get(section).get(option)) {
						if(value == null) {
							continue;
						}
						if(!option.getAllowableNumbersOfValues().contains(Integer.valueOf(value.getActualNumValues()))) {
							exitWithMessageAndHelpMenu("Option " + option.getName() + " in section " + sectionName + " must have " + option.getAllowableNumbersOfValues() + " fields including the flag itself. (Line = " + value.getFullOptionLine() + ")");
						}
					}
				}
			}
		}
		
	}
	
	private String getHelpMenu() {

		String helpMenu = "\n************************************************************\n";
		helpMenu += "************         CONFIG FILE GUIDE         *************\n";
		helpMenu += "************************************************************\n";
		
		helpMenu += "\n" + COMMENT_CHAR + "<Comment>\n";
		helpMenu += COMMENT_CHAR + "<Comment>\n";
		helpMenu += SECTION_HEADER_CHAR + "<Section_name>\n";
		helpMenu += SECTION_HEADER_CHAR + "<Section_name>\n";
		
		for(String sectionName : allowableSectionsByName.keySet()) {
			helpMenu += "\n";
			ConfigFileSection section = allowableSectionsByName.get(sectionName);
			helpMenu += SECTION_HEADER_CHAR + sectionName;
			if(section.isRequired()) {
				helpMenu += "\t(Required)";
			}
			helpMenu += "\n";
			for(ConfigFileOption option : section.getAllowableOptions()) {
				for(Integer numVals : option.getAllowableNumbersOfValues()) {
					String helpMenuLine = option.getName() + "\t";
					for(int i = 2; i <= numVals.intValue(); i++) {
						helpMenuLine += "<value" + Integer.valueOf(i-1).toString() + ">\t";
					}
					if(option.isRequired()) {
						helpMenuLine += "(Required)\t";
					} 
					if(option.isRepeatable()) {
						helpMenuLine += "(Repeatable)\t";
					}
					if(option.hasDefault()) {
						helpMenuLine += "(Default = " + option.getDefaultValue().getLineMinusFlag() + ")\t";
					}
					helpMenuLine += "\n";
					helpMenu += helpMenuLine;
				}
			}
		}
		
		helpMenu += "\n************************************************************\n";
		helpMenu += "************************************************************\n";
		
		return helpMenu;
		
	}
	
	/**
	 * Print help menu and exit
	 */
	public void exitWithHelpMenu() {
		
		String message = getHelpMenu();
		message += "\n\n";
		
		System.err.print(message);
		System.exit(0);
		
	}

	
	/**
	 * @param errorMessage Error message to print
	 */
	private void exitWithMessageAndHelpMenu(String errorMessage) {
		
		String message = getHelpMenu();
		message += "\nInvalid config file:\n";
		message += errorMessage + "\n\n";
		
		System.err.print(message);
		System.exit(0);
		
	}
	
	
	public class OptionValuePair {
		
		private ConfigFileOption option;
		private ConfigFileOptionValue value;
		
		public OptionValuePair(ConfigFileOption op, ConfigFileOptionValue val) {
			option = op;
			value = val;
		}
		
		public ConfigFileOption option() {return option;}
		public ConfigFileOptionValue value() {return value;}
		
	}
	
	
	
}

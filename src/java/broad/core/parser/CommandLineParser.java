/**
 * 
 */
package broad.core.parser;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Stack;
import java.util.TreeSet;


/**
 * @author prussell
 * Full-featured command line parser
 */
public final class CommandLineParser {

	
	private boolean isParsed;
	private ArrayList<String> programDescription;

	private HashMap<String,String> stringArgDescriptions;
	private HashMap<String,String> intArgDescriptions;
	private HashMap<String,String> floatArgDescriptions;
	private HashMap<String,String> doubleArgDescriptions;
	private HashMap<String,String> boolArgDescriptions;	
	
	private HashMap<String,String> stringArgDefaults;
	private HashMap<String,Integer> intArgDefaults;
	private HashMap<String,Float> floatArgDefaults;
	private HashMap<String,Double> doubleArgDefaults;
	private HashMap<String,Boolean> boolArgDefaults;	
	
	private HashSet<String> requiredArgs;
	private HashMap<String,String> commandLineValues;

	
	/**
	 * 
	 */
	public CommandLineParser() {
		isParsed = false;
		
		stringArgDescriptions = new HashMap<String,String>();
		intArgDescriptions = new HashMap<String,String>();
		floatArgDescriptions = new HashMap<String,String>();
		doubleArgDescriptions = new HashMap<String,String>();
		boolArgDescriptions = new HashMap<String,String>();	
		
		stringArgDefaults = new HashMap<String,String>();
		intArgDefaults = new HashMap<String,Integer>();
		floatArgDefaults = new HashMap<String,Float>();
		doubleArgDefaults = new HashMap<String,Double>();
		boolArgDefaults = new HashMap<String,Boolean>();	
		
		programDescription = new ArrayList<String>();
		requiredArgs = new HashSet<String>();
		commandLineValues = new HashMap<String,String>();

	}
	
	/**
	 * Sets program description to be printed as part of help menu
	 * @param description The program description
	 */
	public void setProgramDescription(String description) {
		programDescription.add(description);
	}
	
	/**
	 * Adds new string argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addStringArg(String flag, String description, boolean required) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		stringArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
	}
	
	/**
	 * Adds new string argument to set of arguments and stores default
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 * @param def default value
	 */
	public void addStringArg(String flag, String description, boolean required, String def) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		stringArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
		stringArgDefaults.put(flag, def);
	}

	/**
	 * Adds new int argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */	
	public void addIntegerArg(String flag, String description, boolean required) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		intArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
	}

	/**
	 * Adds new int argument to set of arguments and stores default
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 * @param def default value
	 */	
	public void addIntegerArg(String flag, String description, boolean required, Integer def) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		intArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
		intArgDefaults.put(flag, def);
	}

	/**
	 * Adds new float argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addFloatArg(String flag, String description, boolean required) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		floatArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);		
		if(required) requiredArgs.add(flag);
	}
	
	/**
	 * Adds new float argument to set of arguments and stores default
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 * @param def default value
	 */
	public void addFloatArg(String flag, String description, boolean required, Float def) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		floatArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);		
		if(required) requiredArgs.add(flag);
		floatArgDefaults.put(flag, def);
	}
	
	/**
	 * Adds new double argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addDoubleArg(String flag, String description, boolean required) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		doubleArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
	}

	/**
	 * Adds new double argument to set of arguments and stores default
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 * @param def default value
	 */
	public void addDoubleArg(String flag, String description, boolean required, Double def) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		doubleArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
		doubleArgDefaults.put(flag, def);
	}

	/**
	 * Adds new boolean argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addBooleanArg(String flag, String description, boolean required) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		boolArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
	}
	
	/**
	 * Adds new boolean argument to set of arguments and stores default
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 * @param def default value
	 */
	public void addBooleanArg(String flag, String description, boolean required, Boolean def) {
		enforceUniqueFlag(flag);
		enforceUniqueDescription(description);
		boolArgDescriptions.put(flag, description);
		commandLineValues.put(flag, null);
		if(required) requiredArgs.add(flag);
		boolArgDefaults.put(flag, def);
	}
	
	/**
	 * Parse command arguments
	 * If command line is not in proper form, prints help menu and exits
	 * If a required argument is missing, prints help menu and exits
	 * @param args the command line arguments passed to a main program
	 */
	public void parse(String[] args) {
		
		isParsed = false;
				
		commandLineValues.clear();
		int i=0;
		while(i < args.length) {
			
			// Stop when output redirection is encountered
			if(args[i].contentEquals(">") || args[i].contentEquals(">&") || args[i].contentEquals(">!") || args[i].contentEquals(">&!") || args[i].contentEquals("|") || args[i].contentEquals(">>") || args[i].contentEquals(">>&")) break;
			
			// A flag shouldn't be the last item
			if(args.length == i+1) {
				printHelpMessage();
				System.exit(-1);
			}
			
			// Make sure flag exists
			// Can't see same flag twice
			// Next item should not be a flag
			if(!hasFlag(args[i]) || commandLineValues.containsKey(args[i]) || hasFlag(args[i+1])) {
				printHelpMessage();
				System.exit(-1);
			}
			
			// Add entries to map
			commandLineValues.put(args[i], args[i+1]);
			
			// Skip to next flag
			i += 2;
						
		}
		
		// Make sure all required arguments have been provided
		for(String req : requiredArgs) {
			if(!commandLineValues.containsKey(req)) {
				System.err.println("\n------------------------------------------------------");
				System.err.println("Invalid command line: argument " + req + " is required");
				System.err.println("------------------------------------------------------\n");
				printHelpMessage();
				System.exit(-1);
			}
		}
		
		isParsed = true;
		
	}
	
	/**
	 * Get value of String parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return String specified on command line or null if parameter was not specified
	 */
	public String getStringArg(String flag) {
		
		// Make sure command line has been parsed
		if(!isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()");
		}
		
		// Make sure parameter type is correct
		if(!stringArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get String value for non-String parameter " + flag);
		}
		
		if(commandLineValues.get(flag) == null) {
			if(stringArgDefaults.containsKey(flag)) return stringArgDefaults.get(flag);
		}
		
		return commandLineValues.get(flag);
	}

	/**
	 * Get value of int parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Integer specified on command line or null if parameter was not specified
	 */
	public int getIntArg(String flag) {
		
		// Make sure command line has been parsed
		if(!isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!intArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Integer value for non-Integer parameter " + flag); 
		}
		
		if(commandLineValues.get(flag) == null) {
			if(intArgDefaults.containsKey(flag)) return intArgDefaults.get(flag);
		}
		
		return Integer.valueOf(commandLineValues.get(flag),10).intValue();
	}

	/**
	 * Get value of float parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Float specified on command line or null if parameter was not specified
	 */
	public float getFloatArg(String flag) {
		
		// Make sure command line has been parsed
		if(!isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!floatArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Float value for non-Float parameter " + flag); 
		}
		
		if(commandLineValues.get(flag) == null) {
			if(floatArgDefaults.containsKey(flag)) return floatArgDefaults.get(flag).floatValue();
		}
		
		return Float.valueOf(commandLineValues.get(flag)).floatValue();
	}

	/**
	 * Get value of double parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Double specified on command line or null if parameter was not specified
	 */
	public double getDoubleArg(String flag) {
		
		// Make sure command line has been parsed
		if(!isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!doubleArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Double value for non-Double parameter " + flag); 
		}
		
		if(commandLineValues.get(flag) == null) {
			if(doubleArgDefaults.containsKey(flag)) return doubleArgDefaults.get(flag);
		}
		
		return Double.valueOf(commandLineValues.get(flag)).doubleValue();
	}

	/**
	 * Get value of boolean parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Boolean specified on command line or null if parameter was not specified
	 */
	public boolean getBooleanArg(String flag) {
		
		// Make sure command line has been parsed
		if(!isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!boolArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Boolean value for non-Boolean parameter " + flag); 
		}
		
		if(commandLineValues.get(flag) == null) {
			if(boolArgDefaults.containsKey(flag)) return boolArgDefaults.get(flag);
		}
		
		return Boolean.valueOf(commandLineValues.get(flag)).booleanValue();
	}
	
	
	
	/**
	 * Prints program description plus argument flags and descriptions
	 */
	private void printHelpMessage() {
		System.err.println();
		if(!programDescription.isEmpty()) {
			for(String s : programDescription) System.err.println(s + "\n");
			System.err.println();
		}
		
		TreeSet<String> args = new TreeSet<String>();
		
		for(String key : stringArgDescriptions.keySet()) {
			String msg = key + " <String>\t" + stringArgDescriptions.get(key);
			if(requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + stringArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : intArgDescriptions.keySet()) {
			String msg = key + " <int>\t" + intArgDescriptions.get(key);
			if(requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + intArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : floatArgDescriptions.keySet()) {
			String msg = key + " <float>\t" + floatArgDescriptions.get(key);
			if(requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + floatArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : doubleArgDescriptions.keySet()) {
			String msg = key + " <double>\t" + doubleArgDescriptions.get(key);
			if(requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + doubleArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : boolArgDescriptions.keySet()) {
			String msg = key + " <boolean>\t" + boolArgDescriptions.get(key);
			if(requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + boolArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}

		for(String s : args) {
			System.err.println(s);
		}
		System.err.println();
	
	}
	
	/**
	 * Checks if flag has already been added
	 * @param flag
	 * @return true if and only if flag has already been used
	 */
	private boolean hasFlag(String flag) {
		return (stringArgDescriptions.containsKey(flag) || intArgDescriptions.containsKey(flag) || floatArgDescriptions.containsKey(flag) || doubleArgDescriptions.containsKey(flag) || boolArgDescriptions.containsKey(flag));
	}
	
	/**
	 * Checks if description has already been added
	 * @param description
	 * @return true if and only if description has already been used
	 */
	private boolean hasDescription(String description) {
		return (stringArgDescriptions.containsValue(description) || intArgDescriptions.containsValue(description) || floatArgDescriptions.containsValue(description) || doubleArgDescriptions.containsValue(description) || boolArgDescriptions.containsValue(description));
	}
	
	/**
	 * Causes client to crash with error message if argument flag has already been used
	 * @param flag
	 */
	private void enforceUniqueFlag(String flag) {
		if(hasFlag(flag)) {
			System.err.println("Flag " + flag + " has already been used."); 
			System.exit(-1);
		}
	}
	
	/**
	 * Causes client to crash with error message if argument description has already been used
	 * @param description
	 */
	private void enforceUniqueDescription(String description) {
		if(hasDescription(description)) {
			System.err.println("Description " + description + " has already been used."); 
			System.exit(-1);
		}
	}
	
	
	
}

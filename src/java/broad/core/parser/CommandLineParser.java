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
		this.isParsed = false;
		
		this.stringArgDescriptions = new HashMap<String,String>();
		this.intArgDescriptions = new HashMap<String,String>();
		this.floatArgDescriptions = new HashMap<String,String>();
		this.doubleArgDescriptions = new HashMap<String,String>();
		this.boolArgDescriptions = new HashMap<String,String>();	
		
		this.stringArgDefaults = new HashMap<String,String>();
		this.intArgDefaults = new HashMap<String,Integer>();
		this.floatArgDefaults = new HashMap<String,Float>();
		this.doubleArgDefaults = new HashMap<String,Double>();
		this.boolArgDefaults = new HashMap<String,Boolean>();	
		
		this.programDescription = new ArrayList<String>();
		this.requiredArgs = new HashSet<String>();
		this.commandLineValues = new HashMap<String,String>();

	}
	
	/**
	 * Sets program description to be printed as part of help menu
	 * @param description The program description
	 */
	public void setProgramDescription(String description) {
		this.programDescription.add(description);
	}
	
	/**
	 * Adds new string argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addStringArg(String flag, String description, boolean required) {
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.stringArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
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
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.stringArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
		this.stringArgDefaults.put(flag, def);
	}

	/**
	 * Adds new int argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */	
	public void addIntegerArg(String flag, String description, boolean required) {
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.intArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
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
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.intArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
		this.intArgDefaults.put(flag, def);
	}

	/**
	 * Adds new float argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addFloatArg(String flag, String description, boolean required) {
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.floatArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);		
		if(required) this.requiredArgs.add(flag);
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
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.floatArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);		
		if(required) this.requiredArgs.add(flag);
		this.floatArgDefaults.put(flag, def);
	}
	
	/**
	 * Adds new double argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addDoubleArg(String flag, String description, boolean required) {
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.doubleArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
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
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.doubleArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
		this.doubleArgDefaults.put(flag, def);
	}

	/**
	 * Adds new boolean argument to set of arguments
	 * Client will crash if the same argument flag or description has already been added
	 * @param flag the command line flag for the argument
	 * @param description the description of the argument
	 * @param required whether parameter is required
	 */
	public void addBooleanArg(String flag, String description, boolean required) {
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.boolArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
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
		this.enforceUniqueFlag(flag);
		this.enforceUniqueDescription(description);
		this.boolArgDescriptions.put(flag, description);
		this.commandLineValues.put(flag, null);
		if(required) this.requiredArgs.add(flag);
		this.boolArgDefaults.put(flag, def);
	}
	
	/**
	 * Parse command arguments
	 * If command line is not in proper form, prints help menu and exits
	 * If a required argument is missing, prints help menu and exits
	 * @param args the command line arguments passed to a main program
	 */
	public void parse(String[] args) {
		
		this.isParsed = false;
				
		this.commandLineValues.clear();
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
			if(!this.hasFlag(args[i]) || this.commandLineValues.containsKey(args[i]) || this.hasFlag(args[i+1])) {
				printHelpMessage();
				System.exit(-1);
			}
			
			// Add entries to map
			this.commandLineValues.put(args[i], args[i+1]);
			
			// Skip to next flag
			i += 2;
						
		}
		
		// Make sure all required arguments have been provided
		for(String req : this.requiredArgs) {
			if(!this.commandLineValues.containsKey(req)) {
				System.err.println("\n------------------------------------------------------");
				System.err.println("Invalid command line: argument " + req + " is required");
				System.err.println("------------------------------------------------------\n");
				printHelpMessage();
				System.exit(-1);
			}
		}
		
		this.isParsed = true;
		
	}
	
	/**
	 * Get value of String parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return String specified on command line or null if parameter was not specified
	 */
	public String getStringArg(String flag) {
		
		// Make sure command line has been parsed
		if(!this.isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()");
		}
		
		// Make sure parameter type is correct
		if(!this.stringArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get String value for non-String parameter " + flag);
		}
		
		if(this.commandLineValues.get(flag) == null) {
			if(this.stringArgDefaults.containsKey(flag)) return this.stringArgDefaults.get(flag);
		}
		
		return this.commandLineValues.get(flag);
	}

	/**
	 * Get value of Integer parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Integer specified on command line or null if parameter was not specified
	 */
	public Integer getIntegerArg(String flag) {
		
		// Make sure command line has been parsed
		if(!this.isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!this.intArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Integer value for non-Integer parameter " + flag); 
		}
		
		if(this.commandLineValues.get(flag) == null) {
			if(this.intArgDefaults.containsKey(flag)) return this.intArgDefaults.get(flag);
		}
		
		return Integer.valueOf(this.commandLineValues.get(flag),10);
	}

	/**
	 * Get value of Float parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Float specified on command line or null if parameter was not specified
	 */
	public Float getFloatArg(String flag) {
		
		// Make sure command line has been parsed
		if(!this.isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!this.floatArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Float value for non-Float parameter " + flag); 
		}
		
		if(this.commandLineValues.get(flag) == null) {
			if(this.floatArgDefaults.containsKey(flag)) return this.floatArgDefaults.get(flag);
		}
		
		return Float.valueOf(this.commandLineValues.get(flag));
	}

	/**
	 * Get value of Double parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Double specified on command line or null if parameter was not specified
	 */
	public Double getDoubleArg(String flag) {
		
		// Make sure command line has been parsed
		if(!this.isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!this.doubleArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Double value for non-Double parameter " + flag); 
		}
		
		if(this.commandLineValues.get(flag) == null) {
			if(this.doubleArgDefaults.containsKey(flag)) return this.doubleArgDefaults.get(flag);
		}
		
		return Double.valueOf(this.commandLineValues.get(flag));
	}

	/**
	 * Get value of Boolean parameter specified by flag
	 * @param flag The command line flag for the argument
	 * @return Boolean specified on command line or null if parameter was not specified
	 */
	public Boolean getBooleanArg(String flag) {
		
		// Make sure command line has been parsed
		if(!this.isParsed) {
			throw new IllegalStateException("Cannot get parameter value without first calling method parse()"); 
		}
		
		// Make sure parameter type is correct
		if(!this.boolArgDescriptions.containsKey(flag)) {
			throw new IllegalArgumentException("Trying to get Boolean value for non-Boolean parameter " + flag); 
		}
		
		if(this.commandLineValues.get(flag) == null) {
			if(this.boolArgDefaults.containsKey(flag)) return this.boolArgDefaults.get(flag);
		}
		
		return Boolean.valueOf(this.commandLineValues.get(flag));
	}
	
	
	
	/**
	 * Prints program description plus argument flags and descriptions
	 */
	private void printHelpMessage() {
		System.err.println();
		if(!this.programDescription.isEmpty()) {
			for(String s : this.programDescription) System.err.println(s + "\n");
			System.err.println();
		}
		
		TreeSet<String> args = new TreeSet<String>();
		
		for(String key : this.stringArgDescriptions.keySet()) {
			String msg = key + " <String>\t" + this.stringArgDescriptions.get(key);
			if(this.requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + this.stringArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : this.intArgDescriptions.keySet()) {
			String msg = key + " <int>\t" + this.intArgDescriptions.get(key);
			if(this.requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + this.intArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : this.floatArgDescriptions.keySet()) {
			String msg = key + " <float>\t" + this.floatArgDescriptions.get(key);
			if(this.requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + this.floatArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : this.doubleArgDescriptions.keySet()) {
			String msg = key + " <double>\t" + this.doubleArgDescriptions.get(key);
			if(this.requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + this.doubleArgDefaults.get(key) + ")\n";
			args.add(msg); 
		}
		for(String key : this.boolArgDescriptions.keySet()) {
			String msg = key + " <boolean>\t" + this.boolArgDescriptions.get(key);
			if(this.requiredArgs.contains(key)) msg += " (required)\n";
			else msg += " (default=" + this.boolArgDefaults.get(key) + ")\n";
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
		return (this.stringArgDescriptions.containsKey(flag) || this.intArgDescriptions.containsKey(flag) || this.floatArgDescriptions.containsKey(flag) || this.doubleArgDescriptions.containsKey(flag) || this.boolArgDescriptions.containsKey(flag));
	}
	
	/**
	 * Checks if description has already been added
	 * @param description
	 * @return true if and only if description has already been used
	 */
	private boolean hasDescription(String description) {
		return (this.stringArgDescriptions.containsValue(description) || this.intArgDescriptions.containsValue(description) || this.floatArgDescriptions.containsValue(description) || this.doubleArgDescriptions.containsValue(description) || this.boolArgDescriptions.containsValue(description));
	}
	
	/**
	 * Causes client to crash with error message if argument flag has already been used
	 * @param flag
	 */
	private void enforceUniqueFlag(String flag) {
		if(this.hasFlag(flag)) {
			System.err.println("Flag " + flag + " has already been used."); 
			System.exit(-1);
		}
	}
	
	/**
	 * Causes client to crash with error message if argument description has already been used
	 * @param description
	 */
	private void enforceUniqueDescription(String description) {
		if(this.hasDescription(description)) {
			System.err.println("Description " + description + " has already been used."); 
			System.exit(-1);
		}
	}
	
	
	
}

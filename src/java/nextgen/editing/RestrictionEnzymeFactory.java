/**
 * 
 */
package nextgen.editing;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

import broad.core.parser.StringParser;


/**
 * Factory class for restriction enzymes
 * All new restriction enzymes should be added here
 * @author prussell
 *
 */
public class RestrictionEnzymeFactory {
	
	public static Logger logger = Logger.getLogger(RestrictionEnzymeFactory.class.getName());

	/**
	 * Get restriction enzyme by name
	 * @param enzymeName Enzyme name
	 * @return Restriction enzyme object representing the activity of this enzyme
	 */
	public static SingleCleavageTypeIIRestrictionEnzyme getRestrictionEnzyme(String enzymeName) {
		RestrictionEnzymeName r = RestrictionEnzymeName.fromString(enzymeName);
		switch(r) {
		case BSAI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BSAI, "GGTCTC", 7, -5);
		case SAPI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.SAPI, "GCTCTTC", 8, -4);
		case BSMBI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BSMBI, "CGTCTC", 7, -5);
		case BBSI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BBSI, "GAAGAC", 8, -6);
		case BTSI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BTSI, "GCAGTG", 8, 0);
		case BMRI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BMRI, "ACTGGG", 11, -4);
		case BGLII:
			return new SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeName.BGLII, "AGATCT", 1, 1);
		case SALI:
			return new SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeName.SALI, "GTCGAC", 1, 1);
		case SCAI:
			return new SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeName.SCAI, "AGTACT", 3, 3);
		case STUI:
			return new SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeName.STUI, "AGGCCT", 3, 3);
		case SRFI:
			return new SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeName.SRFI, "GCCCGGGC", 4, 4);
		case PSII:
			return new SingleCleavageTypeIIRestrictionEnzyme(RestrictionEnzymeName.PSII, "TTATAA", 3, 3);
		default:
			throw new IllegalStateException("Enzyme " + r.toString() + " not supported.");
		}
		
	}
	
	/**
	 * Get restriction enzyme by name
	 * @param enzymeName Enzyme name
	 * @return Restriction enzyme object representing the activity of this enzyme
	 */
	public static TypeIISRestrictionEnzyme getTypeIISRestrictionEnzyme(String enzymeName) {
		RestrictionEnzymeName r = RestrictionEnzymeName.fromString(enzymeName);
		switch(r) {
		case BSAI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BSAI, "GGTCTC", 7, -5);
		case SAPI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.SAPI, "GCTCTTC", 8, -4);
		case BSMBI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BSMBI, "CGTCTC", 7, -5);
		case BBSI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BBSI, "GAAGAC", 8, -6);
		case BTSI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BTSI, "GCAGTG", 8, 0);
		case BMRI:
			return new TypeIISRestrictionEnzyme(RestrictionEnzymeName.BMRI, "ACTGGG", 11, -4);
		default:
			throw new IllegalStateException("Enzyme " + r.toString() + " not supported as type IIS.");
		}
		
	}
	

	
	/**
	 * Get collection of restriction enzymes specified in file, one name per line
	 * @param file File name
	 * @return Collection of restriction enzymes named in the file
	 * @throws IOException
	 */
	public static Collection<SingleCleavageTypeIIRestrictionEnzyme> readFromFile(String file) throws IOException {
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		Collection<SingleCleavageTypeIIRestrictionEnzyme> rtrn = new ArrayList<SingleCleavageTypeIIRestrictionEnzyme>();
		while(b.ready()) {
			String line = b.readLine();
			rtrn.add(getRestrictionEnzyme(line));
		}
		r.close();
		b.close();
		return rtrn;
	}
	
	/**
	 * Get collection of restriction enzymes specified in file, one name per line
	 * @param file File name
	 * @return Collection of restriction enzymes named in the file
	 * @throws IOException
	 */
	public static Collection<TypeIISRestrictionEnzyme> readFromFileAsTypeIIS(String file) throws IOException {
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		Collection<TypeIISRestrictionEnzyme> rtrn = new ArrayList<TypeIISRestrictionEnzyme>();
		while(b.ready()) {
			String line = b.readLine();
			rtrn.add(getTypeIISRestrictionEnzyme(line));
		}
		r.close();
		b.close();
		return rtrn;
	}
	

	/**
	 * Get a set of enzyme pairs from a file
	 * @param listFileEnzymePairs The file where each line is <left_enzyme_name> <right_enzyme_name>
	 * @return Collection of enzyme pairs specified in the file
	 * @throws IOException
	 */
	public static Collection<RestrictionEnzymePair> readPairsFromFile(String listFileEnzymePairs) throws IOException {
		FileReader r = new FileReader(listFileEnzymePairs);
		BufferedReader b = new BufferedReader(r);
		Collection<RestrictionEnzymePair> rtrn = new ArrayList<RestrictionEnzymePair>();
		StringParser s = new StringParser();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			String e1 = s.asString(0);
			String e2 = s.asString(1);
			RestrictionEnzymePair pair = new RestrictionEnzymePair(getRestrictionEnzyme(e1), getRestrictionEnzyme(e2));
			rtrn.add(pair);
		}
		r.close();
		b.close();
		return rtrn;
	}



	/**
	 * Restriction enzyme name
	 * @author prussell
	 *
	 */
	public enum RestrictionEnzymeName {
		
		/**
		 * BglII
		 */
		BGLII,
		
		/**
		 * SalI
		 */
		SALI,
		
		/**
		 * ScaI
		 */
		SCAI,
		
		/**
		 * StuI
		 */
		STUI,
		
		/**
		 * SrfI
		 */
		SRFI,
		
		/**
		 * PsiI
		 */
		PSII,
		
		/**
		 * BsaI
		 */
		BSAI,
		
		/**
		 * SapI
		 */
		SAPI,
		
		/**
		 * BsmBI
		 */
		BSMBI,
		
		/**
		 * BbsI
		 */
		BBSI,
		
		/**
		 * BtsI
		 */
		BTSI,
				
		/**
		 * BmRI
		 */
		BMRI;
		
		/**
		 * Get a comma separated list of all available names
		 * @return List of available names
		 */
		public static String commaSeparatedList() {
			String rtrn = "";
			for(int i = 0; i < values().length; i++) {
				if(i > 0) {
					rtrn += ", ";
				}
				rtrn += values()[i];
			}
			return rtrn;
		}
		
		/**
		 * Create from string
		 * @param name Name
		 * @return Object corresponding to this name
		 */
		public static RestrictionEnzymeName fromString(String name) {
			for(RestrictionEnzymeName r : values()) {
				if(name.equals(r.toString())) {
					return r;
				}
			}
			throw new IllegalArgumentException("Enzyme " + name + " not recognized. Options: " + commaSeparatedList());
		}
		
		@Override
		public String toString() {
			switch(this) {
			case BSAI:
				return "BsaI";
			case SAPI:
				return "SapI";
			case BSMBI:
				return "BsmBI";
			case BBSI:
				return "BbsI";
			case BTSI:
				return "BtsI";
			case BMRI:
				return "BmRI";
			case BGLII:
				return "BglII";
			case PSII:
				return "PsiI";
			case SALI:
				return "SalI";
			case SCAI:
				return "ScaI";
			case SRFI:
				return "SrfI";
			case STUI:
				return "StuI";
			default:
				throw new IllegalStateException("Element not supported.");
			}
		}
		
	}


}

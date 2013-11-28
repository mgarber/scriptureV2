/**
 * 
 */
package nextgen.synbio;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

/**
 * @author prussell
 *
 */
public class RestrictionEnzymeFactory {
	
	/**
	 * Get restriction enzyme by name
	 * @param enzymeName Enzyme name
	 * @return Restriction enzyme object representing the activity of this enzyme
	 */
	public static TypeIISRestrictionEnzyme getRestrictionEnzyme(String enzymeName) {
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
			throw new IllegalStateException("Enzyme " + r.toString() + " not supported.");
		}
		
	}
	
	
	/**
	 * Get collection of restriction enzymes specified in file, one name per line
	 * @param file File name
	 * @return Collection of restriction enzymes named in the file
	 * @throws IOException
	 */
	public static Collection<TypeIISRestrictionEnzyme> readFromFile(String file) throws IOException {
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		Collection<TypeIISRestrictionEnzyme> rtrn = new ArrayList<TypeIISRestrictionEnzyme>();
		while(b.ready()) {
			String line = b.readLine();
			rtrn.add(getRestrictionEnzyme(line));
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
			default:
				throw new IllegalStateException("Element not supported.");
			}
		}
		
	}

}

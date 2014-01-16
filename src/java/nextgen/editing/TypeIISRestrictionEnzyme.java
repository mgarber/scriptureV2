/**
 * 
 */
package nextgen.editing;


import org.apache.log4j.Logger;


/**
 * Single cleavage type II restriction enzyme that cleaves 5' of the recognition site on the bottom strand
 * @author prussell
 */
public class TypeIISRestrictionEnzyme  extends SingleCleavageTypeIIRestrictionEnzyme {
	
	public static Logger logger = Logger.getLogger(TypeIISRestrictionEnzyme.class.getName());

	
	/**
	 * @param enzymeName Enzyme name
	 * @param topStrandRecognitionSequence Recognition sequence on top strand
	 * @param topStrandCleavagePosition Cleavage position on top strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 * @param bottomStrandCleavagePosition Cleavage position on bottom strand relative to cleavage just 5' of the first position of the recognition sequence (e.g. 5' ---G   CWGC--- 3' would be 1.)
	 */
	public TypeIISRestrictionEnzyme(RestrictionEnzymeFactory.RestrictionEnzymeName enzymeName, String topStrandRecognitionSequence, int topStrandCleavagePosition, int bottomStrandCleavagePosition) {
		super(enzymeName, topStrandRecognitionSequence, topStrandCleavagePosition, bottomStrandCleavagePosition);
		if(bottomStrandCleavagePosition > 0) {
			throw new IllegalArgumentException("Enzymes that cleave within or 3' of the recognition site on the bottom strand not supported (you provided bottom strand cleavage postion = " + bottomStrandCleavagePosition + ")");
		}
	}


}

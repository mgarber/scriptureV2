package broad.core.primer3;

import java.lang.reflect.Field;

/*
* $Id: Primer3Configuration.java 35971 2007-02-12 16:01:50Z mgarber $
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2004 by the
* Broad Institute/Massachusetts Institute of Technology.
* All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support
* whatsoever. Neither the Broad Institute nor MIT can be responsible for its
* use, misuse, or functionality.
*
* $Header$
*/

/**
 * A dumb data object that represents configuration information
 * used to drive primer3.
 * <br><br>
 * @author Andrew R. Zimmer
 */
public class Primer3Configuration implements Cloneable{

    public boolean canViolateConstraints;
    
    public double DNAConcentration;
    
    public double optimalMeltingTemp;

    public double minMeltingTemp;

    public double maxMeltingTemp;

    public boolean interpretBasesLiberally;

    public int minQualityScore;

    public double minGCContent;

    public double maxGCContent;

    public boolean useGCClamp;

    public double saltConcentration;

    public int numNBasesAccepted;

    public double selfAnyAlignScore;

    public double selfEndAlignScore;

    public int optimalPrimerSize;

    public int minPrimerSize;

    public int maxPrimerSize;

    public int maxNumPrimersToReturn;

    private Primer3PrimerPickingMode mPrimerPickingMode;

    public int minProductSize;

    public int maxProductSize;
    
    public int optimalProductSize;

    public double maxMeltingTempDiff;

    public double overMeltingTempPenaltyWeight;

    public double underMeltingTempPenaltyWeight;

    public double overGCContentPenaltyWeight;

    public double underGCContentPenaltyWeight;

    public double overLengthPenaltyWeight;

    public double underLengthPenaltyWeight;
    
    public double underProductSizePenaltyWeight;
    
    public double overProductSizePenaltyWeight;

    private boolean mPickLeftPrimers;

    private boolean mPickRightPrimers;

    private boolean mPickPrimerPair;

    public int primerWindowSize;
    
    public String missprimingLibraryFile;

	public float maxMisspriming;
	
	public int maxPolyX;
	
	public int primerNumReturn;

	public int primerMin3PrimeOverlapOfJunction=4;

	public int primerMin5PrimerOverlapOfJunction=7;

	
	public void setPrimerNumReturn(int num){this.primerNumReturn=num;}
	public int getPrimerNumReturn(){return this.primerNumReturn;}
	
    /**
     * Sets the primer3 primer picking mode
     * @param mode the mode
     */
    public void setPrimerPickingMode(Primer3PrimerPickingMode mode) {
        // 3 booleans here are required by native C code.

        if (Primer3PrimerPickingMode.FORWARD_ONLY.equals(mode)) {
            mPickLeftPrimers = true;
            mPickRightPrimers = false;
            mPickPrimerPair = false;
        }
        else if (Primer3PrimerPickingMode.REVERSE_ONLY.equals(mode)) {
            mPickLeftPrimers = false;
            mPickRightPrimers = true;
            mPickPrimerPair = false;
        }
        else if (Primer3PrimerPickingMode.PCR.equals(mode)) {
            mPickLeftPrimers = false;
            mPickRightPrimers = false;
            mPickPrimerPair = true;
        }
        else {
            throw new RuntimeException("primer picking mode must be one of FORWARD_ONLY, REVERSE_ONLY, or PCR");
        }

        mPrimerPickingMode = mode;
    }
    
    public Primer3Configuration copy() throws IllegalArgumentException, IllegalAccessException {
    	Primer3Configuration p3c = new Primer3Configuration();
    	Field [] fields = Primer3Configuration.class.getDeclaredFields();
    	
    	for(int i = 0; i < fields.length; i++) {
    		fields[i].set(p3c, fields[i].get(this));
    	}
    	/*
    	p3c.canViolateConstraints = canViolateConstraints;
    	p3c.DNAConcentration      = DNAConcentration;
    	p3c.interpretBasesLiberally = interpretBasesLiberally;
    	p3c.maxGCContent            = maxGCContent;
    	p3c.maxMeltingTemp          = maxMeltingTemp;
    	p3c.maxMeltingTempDiff      = maxMeltingTempDiff;
    	p3c.maxNumPrimersToReturn   = maxNumPrimersToReturn;
    	p3c.maxPrimerSize	        = maxPrimerSize;
    	p3c.maxPrimersToRetain      = maxPrimersToRetain;
    	p3c.maxProductSize          = maxProductSize;
    	p3c.minGCContent            = minGCContent;
    	p3c.minMeltingTemp          = minMeltingTemp;
    	p3c.minPrimerSize           = minPrimerSize;
    	p3c.minProductSize          = minProductSize;
    	p3c.minQualityScore         = minQualityScore;
    	p3c.mPickLeftPrimers        = mPickLeftPrimers;
    	p3c.mPickPrimerPair         = mPickPrimerPair;
    	p3c.mPickRightPrimers       = mPickRightPrimers;
    	p3c.mPrimerPickingMode      = mPrimerPickingMode;
    	*/
    	return p3c;
    }
    
    public Primer3PrimerPickingMode getPrimerPickingMode() {
        return mPrimerPickingMode;
    }
	public void setPrimerMin3PrimeOverlapOfJunction(int min3) {
		this.primerMin3PrimeOverlapOfJunction=min3;
	}
	
	public void setPrimerMin5PrimeOverlapOfJunction(int min5) {
		this.primerMin5PrimerOverlapOfJunction=min5;
	}

}
package broad.core.primer3;

/*
* $Id: Primer3ConfigurationFactory.java 35971 2007-02-12 16:01:50Z mgarber $
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
 * A factory that produces configurations for primer
 * picking.
 * <br><br>
 * @author Andrew R. Zimmer
 */
public class Primer3ConfigurationFactory {
	public static final String HUMAN_LONG_REP="/seq/mgarber/tools/RepeatMasker/Libraries/20060314/homo_sapiens/longlib";
    /**
     * Returns the primer3 configuration used for
     * finishing primer walks
     * @return the primer3 configuration used for
     * finishing primer walks
     */
    public static Primer3Configuration getFinishingPrimerWalkConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 55.0;
        config.maxMeltingTemp = 60.0;
        config.minMeltingTemp = 53.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 20;
        config.minGCContent = 20.0;
        config.maxGCContent = 80.0;
        config.useGCClamp = false;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 1200;
        config.selfEndAlignScore = 800;
        config.optimalPrimerSize = 25;
        config.minPrimerSize = 22;
        config.maxPrimerSize = 28;
        config.maxNumPrimersToReturn = 5000;
        config.canViolateConstraints = true;
        config.primerWindowSize = 600;

        return config;
    }
    /**
     * Returns the primer3 configuration used for
     * short range PCR
     * @return the primer3 configuration used for
     * finishing primer walks
     */
    public static Primer3Configuration getLongRangePCRConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 62.0;
        config.minMeltingTemp = 58.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 15.0;
        config.maxGCContent = 85.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 23;
        config.minPrimerSize = 19;
        config.maxPrimerSize = 30;
        config.maxProductSize = 500;
        config.minProductSize = 300;
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;

        return config;
    }
    
    /**
     * Returns the primer3 configuration used for
     * short range PCR
     * @return the primer3 configuration used for
     * finishing primer walks
     * use human mispriming library
-primer size: 18-30
-primer GC: 40-60
-max self complem.: 5
-max 3' self complem.: 2
-max polyX: 3 
     */
    public static Primer3Configuration getOptimalPCRConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 62.0;
        config.minMeltingTemp = 58.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 40.0;
        config.maxGCContent = 60.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 5.0;
        config.selfEndAlignScore = 2.0;
        config.optimalPrimerSize = 23;
        config.minPrimerSize = 20;
        config.maxPrimerSize = 30;
        config.maxProductSize = 11000;
        config.minProductSize = 5500;
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.maxPolyX = 3;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;

        return config;
    }    
    
    /**
     * Returns the primer3 configuration used for
     * short range PCR
     * @return the primer3 configuration used for
     * finishing primer walks
     */
    public static Primer3Configuration getShortRangePCRConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 63.0;
        config.minMeltingTemp = 58.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 10.0;
        config.maxGCContent = 90.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 18;
        config.minPrimerSize = 17;
        config.maxPrimerSize = 21;
        config.maxProductSize = 234;
        config.minProductSize = 100;
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;

        return config;
    }
    
    public static Primer3Configuration getShortRangeShortReadPCRConfiguration() {
    	Primer3Configuration config = getShortRangePCRConfiguration();
    	config.maxProductSize = 130;
    	config.minProductSize = 60;
    	//config.useGCClamp = false;
        config.minMeltingTemp = 57.0;
        //config.maxMisspriming = 5.0f;
    	return config;
    }
    
    /**
     * Returns a less stringent version of the short range PCR configuration
     * 
     */
    public static Primer3Configuration getShortRangeShortReadPCRConfiguration2() {
    	Primer3Configuration config = getShortRangeShortReadPCRConfiguration();
    	config.minProductSize = 50;
    	config.minPrimerSize = 16;
    	config.minMeltingTemp = 56.0;

    	return config;
    }

    /**
     * Returns a less stringent version of the short range PCR configuration 2
     * Basically it reduces the minumum product size to 50
     * 
     */
    public static Primer3Configuration getShortRangeShortReadPCRConfiguration3() {
    	Primer3Configuration config = getShortRangeShortReadPCRConfiguration2();
    	config.useGCClamp = false;
    	config.minPrimerSize = 15;
    	config.maxMisspriming = 5.0f;
    	config.minMeltingTemp = 55.0;
    	return config;
    }
    
    /**
     * Returns a less stringent version of the short range PCR configuration
     * 
     */
    public static Primer3Configuration getShortRangePCRConfiguration2() {
    	Primer3Configuration config = getShortRangePCRConfiguration();
    	config.minProductSize = 50;
    	config.minPrimerSize = 16;
    	config.minMeltingTemp = 56.0;

    	return config;
    }

    /**
     * Returns a less stringent version of the short range PCR configuration 2
     * Basically it reduces the minumum product size to 50
     * 
     */
    public static Primer3Configuration getShortRangePCRConfiguration3() {
    	Primer3Configuration config = getShortRangePCRConfiguration2();
    	config.useGCClamp = false;
    	config.minPrimerSize = 15;
    	config.maxMisspriming = 5.0f;
    	return config;
    }
    
 
    /**
     * Configuration used for running unit tests.  This
     * configuration <b>should never change</b>, as the
     * unit tests rely on it.  Even if the configuration becomes
     * "unreasonable", or "unrealistic" over time, its purpose
     * is to stay fixed so that it is possible to run unit tests
     * reliably.
     * @return the Primer3Configuration used in unit tests
     */
    public static Primer3Configuration getFinPrimerWalkUnitTestConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 55.0;
        config.maxMeltingTemp = 60.0;
        config.minMeltingTemp = 53.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 20;
        config.minGCContent = 20.0;
        config.maxGCContent = 80.0;
        config.useGCClamp = false;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 1200;
        config.selfEndAlignScore = 800;
        config.optimalPrimerSize = 25;
        config.minPrimerSize = 22;
        config.maxPrimerSize = 28;
        config.maxNumPrimersToReturn = 50;
        config.primerWindowSize = 600;

        return config;
    }

    
    /*Returns the primer3 configuration used for qPCR*/
    public static Primer3Configuration getJenRTPCRConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 61.0;
        config.minMeltingTemp = 59.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 10.0;
        config.maxGCContent = 90.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 17;
        config.maxPrimerSize = 24;
        config.maxProductSize = 120;
        config.minProductSize = 65;
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;
        config.maxPolyX=3;
        config.primerNumReturn=20;
        config.primerMin3PrimeOverlapOfJunction=4;
        config.primerMin5PrimerOverlapOfJunction=7;
        return config;
    } 
    
    /*Returns the primer3 configuration used for qPCR*/
    public static Primer3Configuration getqRTPCRConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 61.0;
        config.minMeltingTemp = 59.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 10.0;
        config.maxGCContent = 90.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 18;
        config.maxPrimerSize = 21;
        
        config.maxProductSize = 120;
        config.minProductSize = 65;
                
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;
        
        return config;
    }
    
    /*Returns the primer3 configuration used for ChIP-qPCR*/
    public static Primer3Configuration getChIPqPCRConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 61.0;
        config.minMeltingTemp = 59.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 10.0;
        config.maxGCContent = 90.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 18;
        config.maxPrimerSize = 21;
        
        config.maxProductSize = 140;
        config.minProductSize = 100;
                
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;
        
        return config;
    }
    
    
    /*Returns the primer3 configuration used for qPCR*/
    public static Primer3Configuration getPATConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 61.0;
        config.minMeltingTemp = 59.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 10.0;
        config.maxGCContent = 90.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 18;
        config.maxPrimerSize = 21;
        
        config.maxProductSize = 120;
        config.minProductSize = 22;
                
        
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;
        
        return config;
    }

    
    /* Returns the primer3 configuration used for RNA qPCR design for
     * assaying RAP efficiency.
     * Jesse Engreitz
     * March 6, 2012
     */
    public static Primer3Configuration getRAPqPCRConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 63.0;
        config.minMeltingTemp = 57.0;
        config.maxMeltingTempDiff = 3.0;
        
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 20.0;
        config.maxGCContent = 80.0;
        config.useGCClamp = false;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 8;
        config.selfEndAlignScore = 3;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 18;
        config.maxPrimerSize = 27;
        
        config.maxProductSize = 100;
        config.minProductSize = 60;
                
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;
        
        return config;
    }

    /**
     * The Primer3 configuration settings used for RNA qPCR
     * Added by Pam Russell May 11 2012
     * @return The Primer3 configuration settings used for RNA qPCR
     */
    public static Primer3Configuration getQpcrConfiguration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 63.0;
        config.minMeltingTemp = 57.0;
        config.maxMeltingTempDiff = 3.0;
        
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 20.0;
        config.maxGCContent = 80.0;
        config.useGCClamp = false;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 8;
        config.selfEndAlignScore = 3;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 18;
        config.maxPrimerSize = 27;
        
        config.maxProductSize = 100;
        config.minProductSize = 60;
                
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        
        return config;
    }

    
    public static Primer3Configuration getTest1Configuration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 55.0;
        config.maxMeltingTemp = 63.0;
        config.minMeltingTemp = 53.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 20;
        config.minGCContent = 20.0;
        config.maxGCContent = 80.0;
        config.useGCClamp = false;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 1200;
        config.selfEndAlignScore = 800;
        config.optimalPrimerSize = 25;
        config.minPrimerSize = 22;
        config.maxPrimerSize = 28;
        config.maxNumPrimersToReturn = 5000;
        config.canViolateConstraints = true;
        config.primerWindowSize = 600;

        return config;
    }

    public static Primer3Configuration getTest2Configuration() {
        Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 55.0;
        config.maxMeltingTemp = 60.0;
        config.minMeltingTemp = 53.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 20;
        config.minGCContent = 40.0;
        config.maxGCContent = 60.0;
        config.useGCClamp = false;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 1200;
        config.selfEndAlignScore = 800;
        config.optimalPrimerSize = 25;
        config.minPrimerSize = 22;
        config.maxPrimerSize = 28;
        config.maxNumPrimersToReturn = 5000;
        config.canViolateConstraints = true;
        config.primerWindowSize = 600;

        return config;
    }
	public static Primer3Configuration getSyntheticConfiguration() {
		Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 61.0;
        config.minMeltingTemp = 59.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 40.0;
        config.maxGCContent = 60.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 20;
        config.maxPrimerSize = 20;
        
        config.maxProductSize = 500;
        config.minProductSize = 21;
                
        
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = false;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0.01;
        config.underLengthPenaltyWeight=0.01;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;
        
        return config;
	}
	
	public static Primer3Configuration getSyntheticConfiguration2() {
		Primer3Configuration config = new Primer3Configuration();
        config.optimalMeltingTemp = 60.0;
        config.maxMeltingTemp = 61.0;
        config.minMeltingTemp = 59.0;
        config.interpretBasesLiberally = true;
        config.minQualityScore = 0;
        config.minGCContent = 40.0;
        config.maxGCContent = 60.0;
        config.useGCClamp = true;
        config.saltConcentration = 50.0;
        config.DNAConcentration = 50.0;
        config.numNBasesAccepted = 0;
        config.selfAnyAlignScore = 0;
        config.selfEndAlignScore = 0;
        config.optimalPrimerSize = 20;
        config.minPrimerSize = 20;
        config.maxPrimerSize = 34;
        
        config.maxProductSize = 900;
        config.minProductSize = 35;
                
        
        config.maxNumPrimersToReturn = 1;
        config.canViolateConstraints = true;
        config.primerWindowSize = 0;
        config.overLengthPenaltyWeight=0;
        config.underLengthPenaltyWeight=0;
        config.underMeltingTempPenaltyWeight=0.05;
        config.overMeltingTempPenaltyWeight=0.05;
        config.missprimingLibraryFile=HUMAN_LONG_REP;
        
        return config;
	}
	
}

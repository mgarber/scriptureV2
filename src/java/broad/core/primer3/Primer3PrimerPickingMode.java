package broad.core.primer3;

/*
* $Id: Primer3PrimerPickingMode.java 23619 2006-02-01 15:50:00Z mgarber $
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
 * A type-safe enum of the different "modes"
 * of operation of primer3.
 * <br><br>
 * Valid modes are: <br>
 * FORWARD_ONLY: instructs the primer picker to pick
 * only primers in the forward direction.<br>
 * REVERSE_ONLY: instructs the primer picker to
 * pick only primers in the reverse direction.<br>
 * PCR: instructs the primer picker to pick
 * paired primers for use in PCR experiments.
 */
public class Primer3PrimerPickingMode {

    private String mModeName;

    /**
     * Only picks forward primers
     */
    public static final Primer3PrimerPickingMode FORWARD_ONLY = new Primer3PrimerPickingMode("FORWARD_ONLY");

    /**
     * Only picks reverse primers
     */
    public static final Primer3PrimerPickingMode REVERSE_ONLY = new Primer3PrimerPickingMode("REVERSE_ONLY");

    /**
     * Picks paired primers for PCR.
     */
    public static final Primer3PrimerPickingMode PCR = new Primer3PrimerPickingMode("PCR");

    /**
     * Creates a new Primer3PrimerPickingMode
     * with the mode name
     * @param modeName the name of the mode
     */
    private Primer3PrimerPickingMode(String modeName) {
        if (modeName == null) {
            throw new IllegalArgumentException("modeName must be non-null in Primer3PrimerPickingMode.Primer3PrimerPickingMode");
        }

        mModeName = modeName;
    }

    public boolean equals(Object obj) {
        if (obj instanceof Primer3PrimerPickingMode) {
            Primer3PrimerPickingMode otherMode = (Primer3PrimerPickingMode)obj;
            return otherMode.mModeName.equals(mModeName);
        }
        return false;
    }

    public String toString() {
        return mModeName;
    }
}
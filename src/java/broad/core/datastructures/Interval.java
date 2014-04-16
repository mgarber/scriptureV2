/*
 * $Id: Interval.java 74570 2008-10-30 20:23:05Z tsharpe $
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or functionality.
 */
package broad.core.datastructures;

import com.sleepycat.persist.model.Persistent;

/**
 * Semi-open interval on the integer number line.
 * Turf covered runs from the start value inclusive, up to, but not including, the end value.
 *
 * @author tsharpe
 * @version $Revision: 74570 $
 */
public interface Interval
{
    // bit-wise definitions from which the other constants are composed
    static final int HAS_LESSER_PART = 1;
    static final int HAS_OVERLAPPING_PART = 2;
    static final int HAS_GREATER_PART = 4;

    static final int IS_ADJACENT_AND_EMPTY = 0;
    static final int IS_STRICTLY_LESS = HAS_LESSER_PART; // 1
    static final int IS_SUBSET = HAS_OVERLAPPING_PART; // 2
    static final int IS_LEFT_OVERHANGING_OVERLAPPER = HAS_LESSER_PART | HAS_OVERLAPPING_PART; // 3
    static final int IS_STRICTLY_GREATER = HAS_GREATER_PART; // 4
    // there is no value that equals 5, since that would imply overhanging on left and right without overlapping
    static final int IS_RIGHT_OVERHANGING_OVERLAPPER = HAS_GREATER_PART | HAS_OVERLAPPING_PART; // 6
    static final int IS_SUPERSET = HAS_LESSER_PART | HAS_OVERLAPPING_PART | HAS_GREATER_PART; // 7

    /**
     * Returns the starting point of the interval.
     * @return The start.
     */
    int getStart();

    /**
     * Returns the ending point of the interval.
     * The interval is not regarded as including this point.
     * @return The end.
     */
    int getEnd();

    /**
     * End - start.
     */
    int getLength();

    /**
     * Returns a constant that describes the relationship of this interval
     * to a specified interval with regard to position on the number line.
     * @param interval The interval to compare this one to.
     * @return One of the IS_* constants defined above.
     */
    int getRelationship( Interval interval );

    /**
     * Returns true if this interval ends where the specified interval starts,
     * or vice versa.
     * @param interval The interval to compare this one to.
     * @return True, if adjacent.
     */
    boolean isAdjacent( Interval interval );

    /**
     * A perfectly trivial implementation of the Interval interface.
     */
    @Persistent
    public static class Impl
        implements Interval,java.io.Serializable
    {
    	
    	/**
    	 * For Berkeley DB only
    	 * Do not use this constructor
    	 */
    	public Impl() {
    		mStart = -1;
    		mEnd = -1;
    	}
    	
        public Impl( int start, int end )
        {
            if ( end < start )
            {
                throw new IllegalArgumentException("End must be at least as large as start. End=" + end + " Start=" + start);
            }
            mStart = start;
            mEnd = end;
        }

        public int getStart()
        {
            return mStart;
        }

        public int getEnd()
        {
            return mEnd;
        }

        public int getLength()
        {
            return mEnd - mStart;
        }

        public int getRelationship( Interval interval )
        {
            if ( interval == null )
            {
                throw new IllegalArgumentException("interval cannot be null");
            }
            int result = 0;
            if ( mStart < interval.getStart() )
                result = HAS_LESSER_PART;
            if ( mEnd > interval.getEnd() )
                result |= HAS_GREATER_PART;
            if ( mStart < interval.getEnd() && interval.getStart() < mEnd )
                result |= HAS_OVERLAPPING_PART;
            return result;
        }

        public boolean isAdjacent( Interval interval )
        {
            if ( interval == null )
            {
                throw new IllegalArgumentException("interval cannot be null");
            }
            return mStart == interval.getEnd() || mEnd == interval.getStart();
        }

        private final int mStart;
        private final int mEnd;
    }
}

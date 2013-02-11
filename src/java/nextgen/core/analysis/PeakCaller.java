package nextgen.core.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.CoordinateSpace;

/**
 * @author prussell
 *
 */
public interface PeakCaller {

	/**
	 * Score windows
	 * @param windows Windows to score
	 */
	public void scoreWindows(Collection<Annotation> windows);
	
	/**
	 * Filter windows that do not contain peaks
	 * @param windows Windows to filter
	 * @return Filtered collection of windows
	 * @throws IOException
	 */
	public Collection<Annotation> filterWindows(Collection<Annotation> windows) throws IOException;
	
	/**
	 * Trim coarse peak to optimal coordinates
	 * @param peak The peak
	 * @return The trimmed peak
	 */
	public Annotation trimPeak(Annotation peak);
	
	/**
	 * Merge a collection of peaks
	 * @param peaks Peaks to merge
	 * @return Collection of merged peaks
	 */
	public Collection<Annotation> mergePeaks(Collection<Annotation> peaks);

	/**
	 * Write a set of peaks to a file
	 * @param windows Regions to write
	 * @param out Output file name
	 * @throws IOException
	 */
	public void writeResults(Collection<Annotation> windows, String out) throws IOException;
	
	/**
	 * Write one peak to a file stream
	 * @param windows Region to write
	 * @param writer File writer
	 * @throws IOException
	 */
	public void writeResult(Collection<Annotation> windows, FileWriter writer) throws IOException;

	/**
	 * Set the coordinate space
	 * @param space The coordinate space to use
	 */
	public void setCoordinateSpace(CoordinateSpace space);
	
	/**
	 * Get the coordinate space
	 * @return The coordinate space
	 */
	public CoordinateSpace getCoordinateSpace();

}


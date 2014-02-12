package nextgen.sequentialbarcode.readlayout;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import nextgen.core.alignment.SmithWatermanAlignment;

import org.apache.log4j.Logger;

/**
 * A sequence of read elements that should be present in a read
 * Can have other stuff between them
 * @author prussell
 *
 */
public class ReadLayout {
	
	private ArrayList<ReadSequenceElement> elements;
	private int readLen;
	public static Logger logger = Logger.getLogger(ReadLayout.class.getName());
	
	/**
	 * @param elementSequence Sequence of elements expected in each read
	 * @param readLength Read length
	 */
	public ReadLayout(ArrayList<ReadSequenceElement> elementSequence, int readLength) {
		// Check that all lengths add up to at most read length
		int totalLen = 0;
		for(ReadSequenceElement elt : elementSequence) {
			totalLen += elt.getLength();
		}
		if(totalLen > readLength) {
			throw new IllegalArgumentException("Total length of read sequence elements (" + totalLen + ") must be at most declared read length (" + readLength + ").");
		}
		elements = elementSequence;
		readLen = readLength;
	}
	
	/**
	 * Get names of read elements expected in each read
	 * @return Ordered list of read element names
	 */
	public ArrayList<String> getElementNames() {
		ArrayList<String> rtrn = new ArrayList<String>();
		for(ReadSequenceElement elt : elements) {
			rtrn.add(elt.elementName());
		}
		return rtrn;
	}
	
	/**
	 * Get read elements expected in each read
	 * @return Ordered list of read elements
	 */
	public ArrayList<ReadSequenceElement> getElements() {
		return elements;
	}
	
	/**
	 * Get a list of matched elements in the read sequence
	 * Only returns an object if the whole layout matches the read sequence
	 * Each item on list corresponds to one read sequence element in the layout
	 * List item is an ordered list of matched elements for that read element
	 * If whole layout does not match read sequence, returns null
	 * @param readSequence Read sequence to search for matches to this layout
	 * @return List of lists of matched elements from this layout that appear in the read
	 */
	public List<List<ReadSequenceElement>> getMatchedElements(String readSequence) {
		
		List<List<ReadSequenceElement>> rtrn = new ArrayList<List<ReadSequenceElement>>();
		for(int i = 0; i < elements.size(); i++) {
			rtrn.add(new ArrayList<ReadSequenceElement>());
		}
		
		// Check that the read sequence has the read length required by this layout
		if(readSequence.length() != readLen) {
			logger.debug("WRONG_LENGTH\tRead layout " + toString() + " does not match read " + readSequence + " because lengths are different: " + readLen + ", " + readSequence.length());
			return null;
		}
		
		// Look for all the elements in order; can have other stuff between them
		Iterator<ReadSequenceElement> elementIter = elements.iterator();
		
		// For repeatable elements, save the first occurrences of their "next" element so can keep looking up until next element
		Map<ReadSequenceElement, Integer> stopSignalPos = new HashMap<ReadSequenceElement, Integer>();
		for(ReadSequenceElement elt : elements) {
			if(elt.isRepeatable()) {
				int posNext = SmithWatermanAlignment.ungappedMatchStartOnFirstSequence(readSequence, elt.getStopSignalForRepeatable(), SmithWatermanAlignment.DEFAULT_MATCH_SCORE, SmithWatermanAlignment.DEFAULT_MISMATCH_SCORE, (float)0.9);
				if(posNext != -1) {
					stopSignalPos.put(elt, Integer.valueOf(posNext));
					logger.debug("STOP_SIGNAL\t for element " + elt.getId() + " is at position " + posNext);
				}
			}
		}
		
		// Get the current element and look ahead to the next element
		// If current element is repeatable, will use next element to know when to stop looking for current element
		ReadSequenceElement currElt = elementIter.next();
		ReadSequenceElement nextElt = null;
		if(elementIter.hasNext()) {
			nextElt = elementIter.next();
		}
		int currStart = 0;
		
		// Make sure all elements have been found at least once in the specified order
		boolean[] found = new boolean[elements.size()];
		
		while(currStart < readLen) {
			logger.debug("");
			logger.debug("CURRENT_START\t" + currStart);
			logger.debug("CURRENT_ELEMENT\t" + currElt.elementName());
			logger.debug("NEXT_ELEMENT\t" + (nextElt == null ? null : nextElt.elementName()));
			// If too far along in the read and have not found everything required, return null
			if(currStart + currElt.getLength() > readLen) {
				logger.debug("NO_MATCH_FOR_LAYOUT\tNo match for element " + currElt.elementName() + " in read " + readSequence);
				return null;
			}
			// If current element is repeatable, look for next element at this position
			if(currElt.isRepeatable() && !(currStart + nextElt.getLength() > readLen)) {
				boolean lookNext = false;
				if(!stopSignalPos.containsKey(currElt)) {
					lookNext = true;
				} else if(stopSignalPos.get(currElt).intValue() == currStart) {
					lookNext = true;
				}
				if(lookNext) {
					logger.debug("LOOKING_FOR_NEXT_ELT\tLooking for " + nextElt.elementName() + " at position " + currStart);
					if(nextElt != null && nextElt.matchesSubstringOf(readSequence, currStart)) {
						if(!found[elements.indexOf(currElt)]) {
							// Next element was found before any instance of current element
							logger.debug("FOUND_NEXT_BEFORE_CURRENT\tFound match for next element " + nextElt.elementName() + " before any instance of " + currElt.elementName());
							return null;
						}
						logger.debug("FOUND_NEXT_OK\tFound match for next element " + nextElt.elementName() + " at start position " + currStart + " of read " + readSequence);
						found[elements.indexOf(nextElt)] = true;
						rtrn.get(elements.indexOf(nextElt)).add(nextElt.matchedElement(readSequence.substring(currStart, currStart + nextElt.getLength())));
						logger.debug("NUM_MATCHES\tThere are " + rtrn.get(elements.indexOf(nextElt)).size() + " matches for this element");
						currStart += nextElt.getLength();
						if(!elementIter.hasNext()) {
							// We have found a match for the last element; return
							logger.debug("MATCHED_LAYOUT\tFound match for entire read layout");
							return rtrn;
						} 
						// Now look for the element after "nextElt"
						currElt = elementIter.next();
						logger.debug("NEW_CURR_ELT\t" + currElt.elementName());
						if(elementIter.hasNext()) {
							nextElt = elementIter.next();
							logger.debug("NEW_NEXT_ELT\t" + nextElt.elementName());
						} else {
							logger.debug("NEW_NEXT_ELT\tnull");
							nextElt = null;
						}
						continue;
					}
				}
				logger.debug("NOT_LOOKING_FOR_NEXT_ELT\tNot looking for " + nextElt.elementName() + " at position " + currStart);
			}
			// Look for current element
			if(currElt.matchesSubstringOf(readSequence, currStart)) {
				// Found an instance of current element
				logger.debug("MATCHED_CURRENT_ELEMENT\tFound match for element " + currElt.elementName() + " at start position " + currStart + " of read " + readSequence);
				found[elements.indexOf(currElt)] = true;
				// Add to return data structure
				rtrn.get(elements.indexOf(currElt)).add(currElt.matchedElement(readSequence.substring(currStart, currStart + currElt.getLength())));
				logger.debug("NUM_MATCHES\tThere are " + rtrn.get(elements.indexOf(currElt)).size() + " matches for this element");
				// Change current position to end of element
				currStart += currElt.getLength();
				// Now will look for the next element unless the current element is repeatable
				if(!currElt.isRepeatable()) {
					currElt = nextElt;
					nextElt = elementIter.hasNext() ? elementIter.next() : null;
					logger.debug("NEW_CURR_ELT\t" + currElt.elementName());
					logger.debug("NEW_NEXT_ELT\t" + nextElt.elementName());
				}
				continue;
			}
			logger.debug("No match for element " + currElt.elementName() + " at start position " + currStart + " of read " + readSequence);
			currStart++;
		}
		return null;
	}
	
	public String toString() {
		Iterator<ReadSequenceElement> iter = elements.iterator();
		String rtrn = iter.next().elementName();
		while(iter.hasNext()) {
			rtrn += "_" + iter.next().elementName();
		}
		return rtrn;
	}
	
}

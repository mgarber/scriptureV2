package nextgen.editing.crispr;

import java.io.IOException;

import org.apache.commons.collections15.Predicate;

public interface GuideRNAPredicate extends Predicate<GuideRNA> {
	
	/**
	 * Get the name of this predicate
	 * @return Predicate name
	 */
	public String getPredicateName();
	
	/**
	 * Get a short explanation (no spaces) of why the predicate evaluates to false, e.g. for inclusion in a bed name
	 * @return Short string explanation of false value
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	public String getShortFailureMessage(GuideRNA g) throws IOException, InterruptedException;

}

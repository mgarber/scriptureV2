package nextgen.core.berkeleydb;

import com.sleepycat.persist.evolve.Deleter;
import com.sleepycat.persist.evolve.Mutations;

/**
 * Factory for Mutations objects: collections of mutations for configuring class evolution
 * @author prussell
 *
 */
public class MutationsFactory {
	
	/**
	 * Get deleter for a field
	 * @param declaringClass Class name
	 * @param declaringClassVersion Class version
	 * @param fieldName Name of field to delete
	 * @return A Mutations object containing this deleter
	 */
	public static Mutations getDeleter(String declaringClass, int declaringClassVersion, String fieldName) {
		Mutations rtrn = new Mutations();
		addDeleter(rtrn, declaringClass, declaringClassVersion, fieldName);
		return rtrn;
	}
	
	/**
	 * Add a deleter to an existing Mutations object
	 * @param mutations Mutations object to add to
	 * @param declaringClass Class name
	 * @param declaringClassVersion Class version
	 * @param fieldName Name of field to delete
	 */
	public static void addDeleter(Mutations mutations, String declaringClass, int declaringClassVersion, String fieldName) {
		mutations.addDeleter(new Deleter(declaringClass, declaringClassVersion, fieldName));
	}
	
}

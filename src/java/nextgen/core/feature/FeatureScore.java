package nextgen.core.feature;

import nextgen.core.annotation.Annotation;

/**
 * Represents a generic scoring object
 * @author mguttman
 *
 * @param <T>
 */
public interface FeatureScore<T> {

	/**
	 * Returns the underlying feature
	 */
	T getFeature();
	
}

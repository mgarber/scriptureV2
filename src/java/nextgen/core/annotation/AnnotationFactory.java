package nextgen.core.annotation;

import broad.core.error.ParseException;

public interface AnnotationFactory<T extends Annotation> {
	T create(String[] rawFields) throws ParseException;
}

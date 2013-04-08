package broad.core.annotation;

//Testing GitHub2

import broad.core.error.ParseException;

public interface AnnotationFactory<T extends LightweightGenomicAnnotation> extends nextgen.core.annotation.AnnotationFactory<T> {
	T create(String [] rawFields) throws ParseException;
	T create(GenomicAnnotation baseAnnotation);
	T create(String name);
}

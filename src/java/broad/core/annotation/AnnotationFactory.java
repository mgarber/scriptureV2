package broad.core.annotation;

//Testing GitHub8

import broad.core.error.ParseException;

public interface AnnotationFactory<T extends LightweightGenomicAnnotation> extends nextgen.core.general.TabbedReader.Factory<T> {
	T create(String [] rawFields) throws ParseException;
	T create(GenomicAnnotation baseAnnotation);
	T create(String name);
}

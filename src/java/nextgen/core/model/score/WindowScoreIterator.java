package nextgen.core.model.score;

import java.util.Iterator;

import nextgen.core.annotation.Annotation;
import nextgen.core.feature.Window;
import net.sf.samtools.util.CloseableIterator;

public class WindowScoreIterator<T extends WindowScore> implements CloseableIterator<T> {

	Iterator<? extends Annotation> itr;
	WindowProcessor<T> processor;
	T previousScore = null;

	public WindowScoreIterator(Iterator<? extends Annotation> windowIterator, WindowProcessor<T> processor, Annotation region){
		this.itr = windowIterator;
		this.processor = processor;
		processor.initRegion(region);
	}

	@Override
	public boolean hasNext() {
		return this.itr.hasNext();
	}

	@Override
	public T next() {
		Annotation w = this.itr.next();
		T score= processor.processWindow(w, previousScore);
		//T score= processor.processWindow(w);
		previousScore=score;
		return score;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void close() {
		processor.finishedRegion();
	}
}

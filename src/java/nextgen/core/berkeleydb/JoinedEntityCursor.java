package nextgen.core.berkeleydb;

import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;

import com.sleepycat.persist.EntityCursor;

public class JoinedEntityCursor<E> implements Iterator<E> {

	private List<EntityCursor<E>> cursors;
	private EntityCursor<E> currentCursor;
	private Iterator<E> currentIterator;
	private int numCursors;
	private int currentCursorNum;
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(JoinedEntityCursor.class.getName());
	
	public JoinedEntityCursor (List<EntityCursor<E>> entityCursors) {
		cursors = entityCursors;
		numCursors = cursors.size();
		currentCursorNum = 0;
		currentCursor = cursors.get(currentCursorNum);
		currentIterator = currentCursor.iterator();
	}
	
	public void close() {
		for(EntityCursor<E> cursor : cursors) {
			cursor.close();
		}
	}
	
	@Override
	public boolean hasNext() {
		boolean currHasNext = currentIterator.hasNext();
		if(currHasNext) {
			return true;
		}
		currentCursor.close();
		currentCursorNum++;
		if(currentCursorNum >= numCursors) {
			return false;
		}
		currentCursor = cursors.get(currentCursorNum);
		currentIterator = currentCursor.iterator();
		return currentIterator.hasNext();
	}

	@Override
	public E next() {
		return currentIterator.next();
	}

	@Override
	public void remove() {
		currentIterator.remove();
	}

}

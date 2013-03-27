package xp.test.Converter;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.SingleEndAlignment;

/**
 *  Created on 2013-3-11  
 */
public class SAMRecordToAlignmentIterator implements CloseableIterator<Alignment> {
	private CloseableIterator<SAMRecord> iter;
	
    
	public SAMRecordToAlignmentIterator(CloseableIterator<SAMRecord> iter) {
		super();
		this.iter = iter;
	}

	@Override
	
	public boolean hasNext() {
		// TODO Auto-generated method stub
		return iter.hasNext();
	}

	@Override
	public Alignment next() {
		// TODO Auto-generated method stub
		SAMRecord a=iter.next();
		Alignment b = new SingleEndAlignment(a);
		return b;
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void close() {
		// TODO Auto-generated method stub
		iter.close();
	}

}

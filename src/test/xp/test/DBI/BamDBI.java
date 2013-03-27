package xp.test.DBI;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import xp.test.Converter.BamToSingleEndIterator;
import xp.test.Converter.SAMRecordToAlignmentIterator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import nextgen.core.alignment.Alignment;
import nextgen.core.coordinatesystem.GenomicSpace;

/**
 *  Created on 2013-3-11  
 */
public class BamDBI implements AlignmentDBI<Alignment> {

	private SAMFileReader data;
	private Integer globalCount;
    ArrayList<String> chrs;
    HashMap<String,Long> chrToLength;
	
    
    
    
    public BamDBI(SAMFileReader data, ArrayList<String> chrs,
			HashMap<String, Long> chrToLength) {
		super();
		this.data = data;
		this.chrs = chrs;
		this.chrToLength = chrToLength;
	}
    
    public BamDBI(SAMFileReader data, GenomicSpace gs)
    {
    	
    	chrs=(ArrayList<String>) gs.getReferenceNames();
    	chrToLength = new HashMap<String,Long>();
    	for(String chr: chrs)
    	{
    		chrToLength.put(chr, gs.getLength(chr));
    	}
    	
    }

	@Override
	public Iterator<Alignment> query(String chr, int start, int end) {
		// TODO Auto-generated method stub
		return new SAMRecordToAlignmentIterator(data.query(chr, start, end, false));
	}

	@Override
	public Iterator<Alignment> iterate() {
		// TODO Auto-generated method stub
		return new SAMRecordToAlignmentIterator(data.iterator());
	}

	@Override
	public int getCount(String chr, int start, int stop) {
		// TODO Auto-generated method stub
		
		return 0;
	}

	@Override
	public double getGlobalCount() {
		// TODO Auto-generated method stub
		if (globalCount!=null) return globalCount;
		countGlobal();
		return globalCount;
	}
	private void countGlobal()
	{
		int i=0;
        Iterator<SAMRecord> iter=data.iterator();
        while(iter.hasNext())
        {
        	i+=1;iter.next();
        }
        globalCount=i; //what if we count twice? does the iterator reset?
	}

	@Override
	public double getLocalCount(String chr, int start, int stop) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public List<String> getChrs() {
		// TODO Auto-generated method stub
		return chrs;
	}

	@Override
	public long getChrLength(String chr) {
		// TODO Auto-generated method stub
		return chrToLength.get(chr);
	}
	
	

}

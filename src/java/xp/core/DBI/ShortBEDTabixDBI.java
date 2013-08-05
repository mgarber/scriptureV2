package xp.test.DBI;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;


import broad.core.annotation.ShortBED;

import net.sf.samtools.TabixReader;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;

/**
 *  Created on 2013-3-4  
 */
public class ShortBEDTabixDBI implements AlignmentDBI<ShortBED> {
	static Logger logger = Logger.getLogger(ShortBEDTabixDBI.class.getName());
	private TabixReader db;
	private HashMap<String,Long> chrToLen = new HashMap<String,Long>();
	private boolean hasGlobalCount = false;
	private double globalCount;
     
	public ShortBEDTabixDBI(String fileDb) throws IOException
	{
		Path path=Paths.get(fileDb);
		this.db=new TabixReader(path);
	}
	
	
	public ShortBEDTabixDBI(String fileDb, String chromSizeFileName) throws IOException
	{
		this(fileDb);
		readChrLength(chromSizeFileName);	
	}
	
	public ShortBEDTabixDBI(String fileDb, GenomicSpace chromSize) throws IOException
	{
		this(fileDb);
		setChrLength(chromSize);	
	}
	@Override
    public IteratorShortBED query(String chr,int beg, int end)
    {
    	Iterator<String> iter=queryString(chr,beg,end);
    	return new IteratorShortBED(iter);
    	
    }
	public IteratorShortBED query(String chr)
    {
		int beg=0;
		int end=(int)getChrLength(chr);
    	Iterator<String> iter=queryString(chr,beg,end);
    	return new IteratorShortBED(iter);
    	
    }
    private class IteratorShortBED implements Iterator<ShortBED>
    {
    	private Iterator<String>  iter;
    	public IteratorShortBED(Iterator<String> _iter)
    	{
    		this.iter=_iter;
    	}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return iter.hasNext();
		}

		@Override
		public ShortBED next() {
			// TODO Auto-generated method stub
			String s = iter.next();
			String[] x = s.split("\t");
			return  new ShortBED("noname",x[0],Integer.valueOf(x[1]), Integer.valueOf(x[2]));
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub
//		  throw new NotImplementedException();	
		}
    	
    }
	public Iterator<String> queryString(String chr, int beg, int end) {
		return db.query(chr, beg, end);
	}
	
	public Iterator<ShortBED> query(Annotation annotation)
	{
		IteratorShortBED _iter=this.query(annotation.getChr(),annotation.getStart(),annotation.getEnd());
		return _iter;
	}

	@Override
	public int getCount(String chr, int beg, int end) {
		Iterator<String> iter=this.queryString(chr, beg, end);
		int i=0;
		while(iter.hasNext()) {iter.next();i++;}
		return i;
	}
	
	public double getGlobalCount()
	{
		if (hasGlobalCount)
			return globalCount;
		Iterator<String> iter=db.iterator();
		int i=0;
		logger.info("Counting Global Stat...");
		while(iter.hasNext()) 
			
		{
			iter.next();
			if (i%1000000==0)
				logger.info("counting "+i+" entries");
			i++;
		}
		hasGlobalCount=true;
		globalCount=(double)i;
		return globalCount;
	}
	public String readHeader()
	{
	 //TO DO
		return null;
	}
	@Override
	public double getLocalCount(String chr, int start, int end) {
		// TODO Auto-generated method stub
		return (double)getCount(chr,start,end);
	}
	@Override
	public List<String> getChrs()
	{
		Iterator<String> iter=chrToLen.keySet().iterator();
		List<String> retv=new ArrayList<String>();
		while(iter.hasNext())
		{
			retv.add(iter.next());
		}
		return retv;
	}
	@Override
	public long getChrLength(String chr) {
		// TODO Auto-generated method stub
		if (chrToLen.containsKey(chr)) return chrToLen.get(chr);
		else return 0;
	}
	public void readChrLength(String chrSizeFileName) throws IOException
	{
		File chrSizeFile= new File(chrSizeFileName);
		readChrLength(chrSizeFile);
	}
	public void readChrLength(File chrSizeFile) throws IOException
	{
		if (!chrSizeFile.exists())
		{
			throw new IOException("Can't find the "+chrSizeFile.toString());
			
		}
		FileInputStream fis = new FileInputStream(chrSizeFile);
		DataInputStream in = new DataInputStream(fis);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		  //Read File Line By Line
		  while ((strLine = br.readLine()) != null)   {
		  // Print the content on the console
		  String[] x=strLine.split("\t");
		  chrToLen.put(x[0].trim(), Long.valueOf(x[1]));
		  }
		br.close();
	}
	public void setChrLength(GenomicSpace g)
	{
		Collection<String> chrs=g.getReferenceNames();
		for (String i: chrs)
		{
			chrToLen.put(i, g.getLength(i));
			
		}
	}

	@Override
	public Iterator<ShortBED> iterate()
	{
		Iterator<String> iter=db.iterator();
		IteratorShortBED iterShortBED=new IteratorShortBED(iter);
		return iterShortBED;
	}
}

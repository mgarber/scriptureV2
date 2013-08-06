package xp.core.Utils;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;

import net.sf.samtools.util.SortingCollection;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;

import org.apache.log4j.Logger;

import broad.core.annotation.ShortBED;

/**
 *  Created on 2013-3-7  
 */
public class AlignmentUtils {
	static Logger logger = Logger.getLogger(AlignmentUtils.class.getName());
	
	private static int MAX_NUM=500000;
    public static int getMAX_NUM() {
		return MAX_NUM;
	}


	public static void setMAX_NUM(int mAX_NUM) {
		MAX_NUM = mAX_NUM;
	}


	private static File tmpDir = new File(System.getProperty("java.io.tmpdir"));
 

	public static SortingCollection<ShortBED>  toShortBEDSortingCollection(Iterator<? extends Annotation> iter)
 {
	 
	 SortingCollection<ShortBED> sortingArray = SortingCollection.newInstance(ShortBED.class, new ShortBEDCodec(), new ShortBEDComparator(), MAX_NUM , tmpDir);
	 
	 while(iter.hasNext())
	 {
		 Annotation a=iter.next();
		 //ShortBED b=new ShortBED("n",a.getChr(),a.getStart(),a.getEnd());
		 //sortingArray.add(b);
		 sortingArray.add((ShortBED)a);
	 }
	 return sortingArray;
	 
	 
 }
	public static SortingCollection<ShortBED>  toShortBEDSortingCollection(Collection<? extends Alignment> alignmentCollection)
	{
		SortingCollection<ShortBED> sortingArray = SortingCollection.newInstance(ShortBED.class, new ShortBEDCodec(), new ShortBEDComparator(), MAX_NUM , tmpDir);
		for(Alignment a: alignmentCollection)
		{
			sortingArray.add((ShortBED)a);
			
		}
		return sortingArray;
		
	}
	
	
	
	
	
}

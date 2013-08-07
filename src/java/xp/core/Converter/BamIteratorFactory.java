package xp.core.Converter;


import java.io.File;
import java.util.Iterator;

import nextgen.core.alignment.Alignment;

import org.apache.log4j.Logger;

/*
 *  Interface for iterate bam as single end or paired end  
 *  
 *  for example:
 *      do_somthing_with_alignment_iterator(BamIteratorFactory.makeIterator(bamfile,"s"))
 *      do_somthing_with_alignment_iterator(BamIteratorFactory.makeIterator(bamfile,"p"))
 *      
 *      
 *      it is an friendly interface to repeat your same procedure
 *      between treat bam as single end and treat bam as paired end
 *      and compare them.
 *      all you need to do is change the factory parameters.
 * 
 * @author  zhuxp
 * @version Experimental
 *      
 */


public class BamIteratorFactory


{
	private static final Logger logger = Logger
			.getLogger(BamIteratorFactory.class);
	public static Iterator<Alignment> makeIterator(File bamFile)
	{
		String t=bamFile.getName();
		String[] a=t.split("\\.");
		//logger.info(bamFile.getName());
		//logger.info(Arrays.toString(a));
		if (a[a.length-1].equalsIgnoreCase("bam"))
		{
	        return makeIterator(bamFile,"s");
		}
		else  
		if (a[a.length-1].equalsIgnoreCase("bamlist") || a[a.length-1].equalsIgnoreCase("bams"))
		{
			logger.info("database is bamlist format");
			//return new BamlistToAlignmentIterator(bamFile.getAbsolutePath());
		}
		return null;
	}
	
	public static Iterator<Alignment> makeIterator(File bamFile, String mode)
	{
		
		return makeIterator(bamFile,mode,5000); //default maxAllowableInsert is 5000;
	}
	public static Iterator<Alignment> makeIterator(File bamFile, String mode, int maxAllowableInsert)
	{
		
		
		if ( mode.equals("p") || mode.equals("P") || mode.equals("Paired") || mode.equals("paired") )
		{
			return new BamToPairedEndIterator(bamFile,maxAllowableInsert);
		}
		
		else if (mode.equals("s") || mode.equals("S") || mode.equals("Single") || mode.equals("single"))
		{
			return new BamToSingleEndIterator(bamFile);
		}
		else
		{
			return new BamToSingleEndIterator(bamFile);
		}
		
		//return new BamToSingleEndIterator(bamFile);
		
	}
	
	
	

	
	
}

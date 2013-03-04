package nextgen.core.readers;



import java.io.File;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.SingleEndAlignment;

/**
 * 
 * 
 * Reading a BAM file and iterate as SingleEndAlignment.
 * Interface for program to use Alignment Iterator
 * 
 * @EXAMPLE
 * {@code
 *
 * BamToSingleEndIterator iter = new BamToSingleEndIterator(File bamfile);
 * while (iter.hasNext()) {
 *		    Alignment  record= iter.next();
 *		    System.out.println(record);
 *		   
 *		}	
 *		iter.close();
 *}
 *
 * @author zhuxp
 * @version Experimental
 * 
 * 
 * 
 */

public class BamToSingleEndIterator implements AlignmentIterator
{
	private SAMFileReader reader;
	private SAMFileHeader header;
	private SAMRecordIterator iter;
	private Alignment nextAlignment;
	public BamToSingleEndIterator(File bamFile)
	{
		reader = new SAMFileReader(bamFile);
		header = reader.getFileHeader();
		iter = reader.iterator();
		nextAlignment=_next();
	}
	
	public Alignment next()
	{  
		Alignment lastAlignment=(Alignment) nextAlignment;
		nextAlignment=_next();
		return lastAlignment;
	}
	
	
	public boolean hasNext()
	{
		
		if (nextAlignment != null)
		{
			return true;
		}
		else
		{
			return false;
		}
	
	}
	
	private Alignment _next()
	{
	
		SAMRecord record;
		while (iter.hasNext())
		{
		try {
			record=iter.next();
			}
			catch (SAMFormatException e)
			{
				System.err.println("catch SAMFormatException");   	
				continue;
			};
		if(record.getReadUnmappedFlag()) continue; 
		SingleEndAlignment retv=new SingleEndAlignment(record);
		return retv;
	
		}
		return null;
		
	
	}
	public void close()
	{
		reader.close();
	}
	
}
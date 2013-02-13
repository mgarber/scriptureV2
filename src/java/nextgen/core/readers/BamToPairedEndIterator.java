package nextgen.core.readers;

import java.io.File;
import java.util.Collection;
import java.util.TreeMap;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.PairedEndAlignment;
import nextgen.core.alignment.PairedEndAlignment.TranscriptionRead;
import nextgen.core.alignment.SingleEndAlignment;
/**
 * Revised from PairedEndWriter
 * Goal : Iterate the paired end reads in bam file.
 * Reading a BAM file and iterate the PairEndAlignment.
 * 
 * @EXAMPLE
 * {@code
 *
 * BamToPairedEndIterator iter = new BamToPairedEndIterator(File bamfile);
 * while (iter.hasNext()) {
 *		    Alignment  record= iter.next();
 *		    System.out.println(record);
 *		  if ( record instanceof PairedEndAlignment)
 *		  {
 *			    System.out.println("First Mate:\t"+((PairedEndAlignment) record).getFirstMate());
 *			    System.out.println("Second Mate:\t"+((PairedEndAlignment) record).getSecondMate());
 *		  }
 *			    System.out.println();
 *		   
 *		}	
 *		iter.close();
 *}
 *
 * @author zhuxp
 * @version Experimental
 * @see BamToPairedEndIteratorTest
 * 
 * 
 */
public class BamToPairedEndIterator {


	
	
	/**
	 * Read Bam File and Iterate the Paired End Reads in alignment format
	 *
	 */
	
	//static Logger logger = Logger.getLogger(BamToPairedEndIterator.class.getName());
		

	private SAMFileReader reader;
	private SAMFileHeader header;
	private SAMRecordIterator iter;
	private TreeMap<String,Alignment> bufferCollection; 
	private String flag;
	
	
	
	private Alignment nextAlignment;
	private TranscriptionRead txnRead;
	
	
	
	//static public final String mateLineFlag="mateLine";

	public BamToPairedEndIterator(File bamFile)
	{
		
		this(bamFile,TranscriptionRead.FIRST_OF_PAIR,"a");
	}
	
	 public BamToPairedEndIterator(File bamFile, TranscriptionRead txnRead)
	 {
		 this(bamFile,txnRead,"a");
	 }	
	/**
	 * @param bamFile Input SAM or BAM file to extract header and/or reads from.
	 * @param 
	 */
	public BamToPairedEndIterator(File bamFile, TranscriptionRead txnRead, String flag) {
		/**
		 * @param bamFile Input SAM or BAM file to extract header and/or reads from.
		 * @param  txnRead  if first read or second read is in direction of trasncription.
		 *         default: first read. TranscriptionRead.FIRST_OF_PAIR
		 * @param flag  "a" to report paired end reads and unpaired end reads
		 *        "p" to only report the paired end reads.
		 * @generator an PairedEndAlignment iterator.
		 *      
		 */
		
		this.flag=flag;
		reader = new SAMFileReader(bamFile);
		header = reader.getFileHeader();
		
		//set header to pairedEnd
		//header.setAttribute(mateLineFlag, "mergedPairedEndFormat");
		iter = reader.iterator();
		bufferCollection = new TreeMap<String,Alignment>();
		this.txnRead=txnRead;
		nextAlignment=_next();
		
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
	
	public Alignment next()
	/*
	 *  iterator the paired end first, then iterator the rest.
	 *  
	 *   report the single read or read with unmapped mate depends on the flag of iterator 
	 *   for now flag=="a" , report all 
	 *   
	 *   @return Alignment object  ( could be SingleEndAlignment or PairedEndAlignment )
	 *   
	 *        
	 */
	{  
	 Alignment lastAlignment=(Alignment) nextAlignment;
	
		nextAlignment=_next();
		return lastAlignment;
		
	}
	private Alignment _next()
		{
		   /*
		    *  private class used to iterate the bam file
		    *    
		    *      
		    */
			while (iter.hasNext())
			{
			SAMRecord record;
			try {	
			record=iter.next();
			}
			catch (SAMFormatException e)
			{
				System.err.println("catch SAMFormatException");
			continue;	
				
			};
			//ignore the unmapped 
			if(record.getReadUnmappedFlag()) continue;
			
			//this is the second read and reverse it 
			if (txnRead.equals(TranscriptionRead.FIRST_OF_PAIR) &&  record.getSecondOfPairFlag() )
			
			{
			
				record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
		       
			}
			
			// this is the first read and reverse it
			if (txnRead.equals(TranscriptionRead.SECOND_OF_PAIR) && record.getFirstOfPairFlag())
			
			{
				record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
			}
			
			
			
			//If the read is not paired or the mate is unmapped, write it as it is
			if(!record.getReadPairedFlag() || record.getMateUnmappedFlag()){
				//mate unmapped so just write it
				record.setMateUnmappedFlag(true);
				SingleEndAlignment s1= new SingleEndAlignment(record);
				// if flag is a, report the mate unmmaped read as single end. 
				if (flag.equals("a"))
				{
				return s1;
				}
				else
				{
					continue;
				}
		    }
			
		
			SingleEndAlignment s1= new SingleEndAlignment(record);
			String name=s1.getName();
			if (bufferCollection.containsKey(name))
			  {
				SingleEndAlignment s2 = (SingleEndAlignment) bufferCollection.get(name);
				
				if (! s2.getChromosome().equals(s1.getChromosome()))
						{
					     if (flag.equals("a"))
					     {
				           return s1; //treat as single end . if mate chromsome is different.	
					     }
					     else
					     {
					    	 bufferCollection.put(name, s1); //put it into buffer and continue...
					    	 continue;
					     }
					     }
				bufferCollection.remove(name);
				PairedEndAlignment p = new PairedEndAlignment(s2,s1);
				return p;
			  }
			  else
			  {
				bufferCollection.put(name, s1);
			  }
			
			}
			/* 
			 *  iterate the single end data ( equals to write the remains) 
			 *  only if the flag is "a"
			 *  
			 */
           if (flag.equals("a"))
        	   
           {
        	   
			while (!bufferCollection.isEmpty() )
			{
				String name=bufferCollection.firstKey();
				Alignment s1=bufferCollection.get(name);
				bufferCollection.remove(name);
				return s1;
				
			}
           }
			return null;
			
		}
	
	
	public void close() {
		reader.close();
		
	}
	public int getBufferSize()
	{
		return bufferCollection.size();
	}
	
}

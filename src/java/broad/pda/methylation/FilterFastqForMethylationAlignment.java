package broad.pda.methylation;

import java.io.File;
import java.util.Iterator;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.fastq.*;


import org.apache.commons.lang3.StringUtils;

public class FilterFastqForMethylationAlignment extends CommandLineProgram {
    private static final Log log = Log.getInstance(FilterFastqForMethylationAlignment.class);
	
    @Usage
    public String USAGE = "Filters a FASTQ file for methylation alignment.";
    
    @Option(doc = "FASTQ file 1", shortName = "1")
	public File INPUT_1;
    
    @Option(doc = "FASTQ file 2", shortName = "2")
	public File INPUT_2;
    
    @Option(doc = "Output FASTQ base", shortName = "O")
    public String OUTPUT;
    
    @Option(doc = "Filtered reads FASTQ base", shortName = "FR")
    public File FILTERED_READS;  
    
    @Override
    protected int doWork() {
    	try {
			IoUtil.assertFileIsReadable(INPUT_1);
			IoUtil.assertFileIsReadable(INPUT_2);
			
			FastqReader reader1 = new FastqReader(INPUT_1);
			FastqReader reader2 = new FastqReader(INPUT_2);
			
			FastqWriter writer1 = new FastqWriterFactory().newWriter(new File(OUTPUT + "_1.fq"));
			FastqWriter writer2 = new FastqWriterFactory().newWriter(new File(OUTPUT + "_2.fq"));
			
			FastqWriter filterWriter1 = new FastqWriterFactory().newWriter(new File(FILTERED_READS + "_1.fq"));
			FastqWriter filterWriter2 = new FastqWriterFactory().newWriter(new File(FILTERED_READS + "_2.fq"));			

			int countFilteredN = 0;
			int countFilteredGC = 0;
			int counter = 0;
			
			while (reader1.hasNext()) {
				FastqRecord rec1 = reader1.next();
				FastqRecord rec2 = reader2.next();
				

	            final String frec1Name = getReadName(rec1.getReadHeader());
	            final String frec2Name = getReadName(rec2.getReadHeader());
	            final String baseName = getBaseName(frec1Name, frec2Name, reader1, reader2);  // errors out if the names are incorrect
				
				boolean filter = false;
				
				// Count number of N bases and filter those with 2+
				if (StringUtils.countMatches(rec1.getReadString(), "N") > 1 || StringUtils.countMatches(rec2.getReadString(), "N") > 1) {
					filter = true;
					countFilteredN++;
				}
				
				// Check to see if the read might be unconverted.
				if (MethylationUtils.guessUnconverted(rec1.getReadString()) || MethylationUtils.guessUnconverted(rec2.getReadString())) {
					filter = true;
					countFilteredGC++;
				}

				if (!filter) {
					writer1.write(rec1);
					writer2.write(rec2);
				} else {
					filterWriter1.write(rec1);
					filterWriter2.write(rec2);
				}
				
				counter++;
				if (counter % 1000000 == 0) log.info(counter + " read pairs completed.");
			}

			log.info(countFilteredN + " reads were filtered due to too many N bases.");
			log.info(countFilteredGC + " reads were filtered due to appearing to be unconverted.");
			log.info((counter - countFilteredN - countFilteredGC) + " reads were retained.");
			
			writer1.close();
			writer2.close();
			filterWriter1.close();
			filterWriter2.close();
    		
    	} catch (Exception e) {
    		log.error(e);
    		return 1;
    	}
    	return 0;
    }
    
    
   
    
    
    // COPIED FROM PICARD FastqToSam.jar
    
    // Read names cannot contain blanks
    private String getReadName(final String fastaqHeader) {
        final int idx = fastaqHeader.indexOf(" ");
        return (idx == -1) ? fastaqHeader : fastaqHeader.substring(0,idx); 
    }
    
    
    /** Returns read baseName and asserts correct pair read name format:
     * <ul>
     * <li> Paired reads must either have the exact same read names or they must contain at least one "/"
     * <li> and the First pair read name must end with "/1" and second pair read name ends with "/2"
     * <li> The baseName (read name part before the /) must be the same for both read names
     * <li> If the read names are exactly the same but end in "/2" or "/1" then an exception will be thrown 
     * </ul>
     */
    String getBaseName(final String readName1, final String readName2, final FastqReader freader1, final FastqReader freader2) {
        String [] toks = getReadNameTokens(readName1, 1, freader1);
        final String baseName1 = toks[0] ;
        final String num1 = toks[1] ;

        toks = getReadNameTokens(readName2, 2, freader2);
        final String baseName2 = toks[0] ;
        final String num2 = toks[1];

        if (!baseName1.equals(baseName2)) {
            throw new PicardException(String.format("In paired mode, read name 1 (%s) does not match read name 2 (%s)", baseName1,baseName2));
        }

        final boolean num1Blank = StringUtil.isBlank(num1);
        final boolean num2Blank = StringUtil.isBlank(num2);
        if (num1Blank || num2Blank) {
            if(!num1Blank) throw new PicardException("Pair 1 number is missing (" +readName1+ "). Both pair numbers must be present or neither.");       //num1 != blank and num2   == blank
            else if(!num2Blank) throw new PicardException("Pair 2 number is missing (" +readName2+ "). Both pair numbers must be present or neither."); //num1 == blank and num =2 != blank 
        } else {
            if (!num1.equals("1")) throw new PicardException("Pair 1 number must be 1 ("+readName1+")");
            if (!num2.equals("2")) throw new PicardException("Pair 2 number must be 2 ("+readName2+")");
        }

        return baseName1 ;
    }


    /** Breaks up read name into baseName and number separated by the last / */
    private String [] getReadNameTokens(final String readName, final int pairNum, final FastqReader freader) {
        if(readName.equals("")) throw new PicardException("Pair read name "+pairNum+" cannot be empty: "+readName);

        final int idx = readName.lastIndexOf("/");
        final String result[] = new String[2];

        if (idx == -1) {
            result[0] = readName;
            result[1] = null;
        } else {
            result[1] = readName.substring(idx+1, readName.length()); // should be a 1 or 2
            
            if(!result[1].equals("1") && !result[1].equals("2")) {    //if not a 1 or 2 then names must be identical
                result[0] = readName;
                result[1] = null;
            }
            else {
                result[0] = readName.substring(0,idx); // baseName
            }
        }

        return result ;
    }
  
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new FilterFastqForMethylationAlignment().instanceMain(args));
	}
}

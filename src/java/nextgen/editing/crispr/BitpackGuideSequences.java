package nextgen.editing.crispr;

import java.io.*;
import java.util.*;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationFileReader;
import nextgen.core.annotation.BasicAnnotation;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.seq.rap.CRISPRDesigner;

/**
 * @author engreitz
 * Use DumpAllGuideSequences to make a BED file first, then bitpack with this class.
 * Order is very important - do not change the BED files or bitpack after running this class.
 * Does not deal with N's ... these get encoded as A's
 */
public class BitpackGuideSequences extends CommandLineProgram {
    private static final Log log = Log.getInstance(BitpackGuideSequences.class);
    
    final static byte A = 0x0;
    final static byte C = 0x1;
    final static byte G = 0x2;
    final static byte T = 0x3;
    
	@Option(doc="Bed file containing guide locations and sequences")
	public File BED;
	
	@Option(doc="Output file with 20mer guide sequences, one per line")
	public File OUTPUT;
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new BitpackGuideSequences().instanceMain(args));
	}
	
    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(BED);
        IoUtil.assertFileIsWritable(OUTPUT);
        
        try {
        	
        	CloseableIterator<Annotation> itr = AnnotationFileReader.read(BED, Annotation.class, new BasicAnnotation.Factory());
           	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(OUTPUT));
			
        	while (itr.hasNext()) {
        		Annotation next = itr.next();
        		byte[] bytes = sequenceToBits(next.getName().substring(0,20));
        		//String reconstructed = bytesToSequence(bytes, 20);
        		//System.out.println("OK");
        		//log.info(next.getName() + " --> bytes --> " + reconstructed);
        		out.write(bytes);
        	}
 
			out.close();
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }
    
    
    public static byte[] sequenceToBits(String seq) {
    	seq = seq.toUpperCase();
    	byte[] bytes = new byte[(seq.length()*2 + 7) / 8];   
    	for (int i = 0; i < seq.length()*2; i += 2) {
    		if (seq.charAt(i/2) == 'A') {
    			bytes[i/8] |= A << (i%8);
    		} else if (seq.charAt(i/2) == 'C') {
    			bytes[i/8] |= C << (i%8);
    		} else if (seq.charAt(i/2) == 'G') {
    			bytes[i/8] |= G << (i%8);
    		} else if (seq.charAt(i/2) == 'T') {
    			bytes[i/8] |= T << (i%8);
    		}
    	}
    	return bytes;
    }
    
    
    public static String bytesToSequence(byte[] bytes, int n) {
    	char[] chars = new char[n];
    	for (int i = 0; i < n*2; i += 2) {
    		if (((byte) (bytes[i/8] >> (i%8)) & 0b00000011) == A) chars[i/2] = 'A';
    		else if (((byte) (bytes[i/8] >> (i%8)) & 0b00000011) == C) chars[i/2] = 'C';
    		else if (((byte) (bytes[i/8] >> (i%8)) & 0b00000011) == G) chars[i/2] = 'G';
    		else if (((byte) (bytes[i/8] >> (i%8)) & 0b00000011) == T) chars[i/2] = 'T';
    	}
    	return String.valueOf(chars);
    }
    
    public static void printByte(byte b) {
		System.out.println(String.format("%8s", Integer.toBinaryString(b & 0xFF)).replace(' ', '0'));
    }
    
	
}

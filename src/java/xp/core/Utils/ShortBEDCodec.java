package xp.core.Utils;
import java.io.InputStream;
import java.io.OutputStream;

import org.apache.log4j.Logger;
import broad.core.annotation.ShortBED;

import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.RuntimeEOFException;
import net.sf.samtools.util.SortingCollection;
import nextgen.core.alignment.Alignment;

public class ShortBEDCodec implements SortingCollection.Codec<ShortBED>
{
	
	static Logger logger = Logger.getLogger(ShortBEDCodec.class.getName());
	private final BinaryCodec binaryCodec = new BinaryCodec();
	public ShortBED decode() {
		String chr;
		try {
		 chr=this.binaryCodec.readNullTerminatedString();
        }
        catch (RuntimeEOFException e) {   //VERY IMPORTANT , TO TELL THAT THE FILE IS ENDING
            return null;
        }
		int start=this.binaryCodec.readInt();
		int stop=this.binaryCodec.readInt();
		ShortBED a= new ShortBED("noname",chr,start,stop);
		return  a;
	}

	@Override
	public void encode(ShortBED a) {
		//String name=a.getName();
		String chr=a.getReferenceName();
		int start=a.getStart();
		int end=a.getEnd();
		//this.binaryCodec.writeString(name, false, true); 
		this.binaryCodec.writeString(chr, false, true); // write null at end of string 
		this.binaryCodec.writeInt(start);
		this.binaryCodec.writeInt(end);
		
		
	}

	@Override
	public void setInputStream(InputStream is) {
		this.binaryCodec.setInputStream(is);
	}

	@Override
	public void setOutputStream(OutputStream os) {
		this.binaryCodec.setOutputStream(os);
	}
	@Override
	public ShortBEDCodec clone()
	{  
		return new ShortBEDCodec();
	}
	
}
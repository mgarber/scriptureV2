package xp.core.Utils;
/**
 *  Created on 2013-3-7  
 */

import java.io.InputStream;
import java.io.OutputStream;

import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.RuntimeEOFException;
import net.sf.samtools.util.SortingCollection;

import org.apache.log4j.Logger;
public class JieCodec implements SortingCollection.Codec<JieCode>{
	static Logger logger = Logger.getLogger(JieCodec.class.getName());
	private final BinaryCodec binaryCodec = new BinaryCodec();
	@Override
	public JieCode decode() {
		// TODO Auto-generated method stub
		int tid;
		try {
		 tid=this.binaryCodec.readInt();
        }
        catch (RuntimeEOFException e) {   //VERY IMPORTANT , TO TELL THAT THE FILE IS ENDING
            return null;
        }
		int pos=this.binaryCodec.readInt();
		boolean isStart=this.binaryCodec.readBoolean();
		int classIndex=this.binaryCodec.readInt();
		int USB=this.binaryCodec.readInt();
		JieCode a= new JieCode(tid,pos,isStart,classIndex,USB);
		return  a;
	}

	@Override
	public void encode(JieCode a) {
		// TODO Auto-generated method stub
		this.binaryCodec.writeInt(a.getTid());
		this.binaryCodec.writeInt(a.getPos());
		this.binaryCodec.writeBoolean(a.isStart());
		this.binaryCodec.writeInt(a.getClassIndex());
		this.binaryCodec.writeInt(a.getUSB());
	}

	@Override
	public void setInputStream(InputStream is) {
		// TODO Auto-generated method stub
		this.binaryCodec.setInputStream(is);
		
	}

	@Override
	public void setOutputStream(OutputStream os) {
		// TODO Auto-generated method stub
		this.binaryCodec.setOutputStream(os);
		
	}
	public JieCodec clone()
	{
		return new JieCodec();
	}

}


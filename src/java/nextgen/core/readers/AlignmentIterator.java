package nextgen.core.readers;




//import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;

public interface AlignmentIterator 
{
	public Alignment next();
	
	public boolean hasNext();
	
}
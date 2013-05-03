package nextgen.core.alignment;

import java.util.ArrayList;

import org.apache.commons.collections15.Predicate;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.model.AlignmentModel;
import junit.framework.TestCase;

public class PairedEndAlignmentTest extends TestCase{

	public void testInstantiateAlignmentModel() {
		
		//NO transcription read NO fragment boolean
		System.out.println("****************************************");
		System.out.println("NO TXN READ NO FRAG BOOLEAN");
		System.out.println("****************************************");
		AlignmentModel model = new AlignmentModel("test.bam", null, new ArrayList<Predicate<Alignment>>(),true);
		
		CloseableIterator<Alignment> iter = model.getReadIterator();
		while(iter.hasNext()){
			Alignment read = iter.next();
			checkInstance(read);
		}
		iter.close();
		//NO transcription read NO fragment boolean
		System.out.println("****************************************");
		System.out.println("TXN READ FIRST NO FRAG BOOLEAN");
		System.out.println("****************************************");
		model = new AlignmentModel("test.bam", null, new ArrayList<Predicate<Alignment>>(),true,TranscriptionRead.FIRST_OF_PAIR);
		
		iter = model.getReadIterator();
		while(iter.hasNext()){
			Alignment read = iter.next();
			checkInstance(read);
		}
		iter.close();
		//NO transcription read NO fragment boolean
		System.out.println("****************************************");
		System.out.println("TXN READ SECOND NO FRAG BOOLEAN");
		System.out.println("****************************************");
		model = new AlignmentModel("test.bam", null, new ArrayList<Predicate<Alignment>>(),true,TranscriptionRead.SECOND_OF_PAIR);
		
		iter = model.getReadIterator();
		while(iter.hasNext()){
			Alignment read = iter.next();
			checkInstance(read);
		}
		iter.close();
		//NO transcription read NO fragment boolean
		System.out.println("****************************************");
		System.out.println("TXN READ SECOND FRAG BOOLEAN");
		System.out.println("****************************************");
		model = new AlignmentModel("test.bam", null, new ArrayList<Predicate<Alignment>>(),true,TranscriptionRead.SECOND_OF_PAIR,false);
		
		iter = model.getReadIterator();
		while(iter.hasNext()){
			Alignment read = iter.next();
			checkInstance(read);
		}

	}
	
	private void checkInstance(Alignment read){
		if(SingleEndAlignment.class.isInstance(read)){
			System.out.println(read.getName()+" is a single end alignment with orientation "+ read.getOrientation());
		}
		else if(PairedReadAlignment.class.isInstance(read)){
			System.out.println(read.getName()+" is a paired read alignment with orientation "+ read.getOrientation());
		}
		else if(FragmentAlignment.class.isInstance(read)){
			System.out.println(read.getName()+" is a fragment alignment with orientation "+ read.getOrientation());
		}
	}

}

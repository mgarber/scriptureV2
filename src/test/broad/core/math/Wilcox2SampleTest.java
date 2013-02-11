package broad.core.math;

import broad.core.math.Wilcox2Sample;
import junit.framework.TestCase;


public class Wilcox2SampleTest extends TestCase {
	
	public void testSimple() {
		double[] x={1,200,3,400,5,600,7,800,9,1000};
		double[] y={100, 12, 300,14, 500, 6, 700, 8, 900, 20};
		Wilcox2Sample W2S=new Wilcox2Sample(x,y);
		double p = 0.7988593499960496;
		double z = -0.2548235957188128;
		System.err.println(W2S.getP()+" "+W2S.getZ());
		assertEquals("Bad pval",p,W2S.getP());
		assertEquals("Bad z",z, W2S.getZ());
	}
	
	public void testSimpleWithDupplicatedVals() {
		double[] x={1,200,200,3,400,5,600,7,800,9,1000,1000};
		double[] y={100, 200,12, 300,14, 500, 6, 700, 8, 900, 20,200};
		Wilcox2Sample W2S=new Wilcox2Sample(x,y);
		double p = 0.5829195472688133;
		double z = -0.5491251783869153;
		System.err.println(W2S.getP()+" "+W2S.getZ());
		assertEquals("Bad pval",p,W2S.getP());
		assertEquals("Bad z",z, W2S.getZ());
	}

}

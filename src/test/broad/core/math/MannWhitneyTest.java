package broad.core.math;

import broad.core.math.MannWhitney;
import junit.framework.TestCase;


public class MannWhitneyTest  extends TestCase{
	
	public void testWithTies() {
		double[] x={1,2,2,3,4,5,6,7,8,9,10};
		double[] y={10, 10, 11, 12, 13,14, 15, 16, 17, 18, 19, 20};
		double pval = 6.321587792512684E-5;
		double z    = -4.000473456828313;
		MannWhitney MW=new MannWhitney(x,y);
		//System.err.println("P: " + MW.getP()  + " Z: " + MW.getZ() );
		assertEquals("Pvalue test failed",pval, MW.getP());
		assertEquals("Z score failed", z, MW.getZ());		
		
	}
	
	public  void testSimple(){
		double[] x={1,2,3,4,5,6,7,8,9,10};
		double[] y={11, 12, 13,14, 15, 16, 17, 18, 19, 20};
		double pval = 1.5705228423079642E-4;
		double z    = -3.779644730092272;
		MannWhitney MW=new MannWhitney(x,y);
		assertEquals("Pvalue test failed",pval, MW.getP());
		assertEquals("Z score failed", z, MW.getZ());
	}
	

}

package broad.core.math;

import broad.core.math.MathUtil;
import junit.framework.TestCase;

public class MathUtilTest extends TestCase {

	public void testRootFinding() {
		double [] coefs = { -27,-72,-6, 1};
		double [] expected = {12.122893784632387, -5.734509942225077, -0.3883838424073199};
				double [][] roots = MathUtil.roots(coefs);
		for (int i = 0; i < roots.length; i++) {
			System.out.println("Root " + i + " " + roots[i][0] + " + i" + roots[i][1]);
			assertEquals("Root "+i +" does not match", expected[i] , roots[i][0]);
		}
	}
}

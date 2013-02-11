package broad.core.siphy;

import java.text.DecimalFormat;
import java.util.List;

import broad.core.siphy.TreeScalerIO;
import broad.core.siphy.TreeScalerIO.ScaledWindow;

import junit.framework.TestCase;

public class TestTreeScalerIO extends TestCase {

	public void testStitchingOmegaWindows() throws Exception {
		TreeScalerIO tsio = new TreeScalerIO();
		
		String [] rawData1 = {"0", "0.5","1"};
		ScaledWindow sw1 = new ScaledWindow(rawData1, "1", 1, 0);		
		tsio.addWindow(sw1);
		
		String [] rawData2 = {"1", "1","1"};
		ScaledWindow sw2 = new ScaledWindow(rawData2, "1", 1, 0);	
		tsio.addWindow(sw2);
		
		tsio.stitch(2,0);
		
		List<ScaledWindow> scalings = tsio.getScaledWindows();
		
		assertEquals("Wrong number of stitched windows",1,scalings.size());
		ScaledWindow testW = scalings.get(0);
		assertEquals("Wrong start", 0, testW.getStart());
		assertEquals("Wrong end", 2, testW.getEnd());
		assertEquals("Wrong averaged omega",0.75,testW.getOmega());
		
		tsio.clear();
		
		String [] rawData3 = {"2","0.75","1.75"};
		String [] rawData4 = {"4","0.75","1.75"};
		
		tsio.addWindow(new ScaledWindow(rawData1, "1", 1, 0));
		tsio.addWindow(new ScaledWindow(rawData2, "1", 1, 0));
		tsio.addWindow(new ScaledWindow(rawData3, "1", 1, 0));
		tsio.addWindow(new ScaledWindow(rawData4, "1", 1, 0));
		
		tsio.stitch(2,0);
		scalings = tsio.getScaledWindows();
		assertEquals("Wrong number of stitched windows",2,scalings.size());
		
		testW = scalings.get(0);
		assertEquals("Wrong start", 0, testW.getStart());
		assertEquals("Wrong end", 2, testW.getEnd());
		assertEquals("Wrong averaged omega",0.75,testW.getOmega());
		
		testW = scalings.get(1);
		assertEquals("Wrong start", 1, testW.getStart());
		assertEquals("Wrong end", 3, testW.getEnd());
		assertEquals("Wrong averaged omega",0.875,testW.getOmega());
			
		tsio.clear();
		
		tsio.addWindow(new ScaledWindow(rawData1, "1", 1, 0));
		tsio.addWindow(new ScaledWindow(rawData2, "1", 1, 0));
		tsio.addWindow(new ScaledWindow(rawData3, "1", 1, 0));
		tsio.addWindow(new ScaledWindow(rawData4, "1", 1, 0));
		
		tsio.stitch(3,0);
		scalings = tsio.getScaledWindows();
		assertEquals("Wrong number of stitched windows",1,scalings.size());
		
		testW = scalings.get(0);
		System.out.println(testW.toFullString(new DecimalFormat("##0.####")));
		assertEquals("Wrong start", 0, testW.getStart());
		assertEquals("Wrong end", 3, testW.getEnd());
		assertEquals("Wrong averaged omega",0.75,testW.getOmega());
		System.out.println(scalings);

	}

}

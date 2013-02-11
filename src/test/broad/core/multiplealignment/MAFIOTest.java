package broad.core.multiplealignment;

import java.net.URL;

import broad.core.multiplealignment.MAFIO;

import junit.framework.TestCase;

public class MAFIOTest extends TestCase {
	public void testCreateIndex() throws java.io.IOException {
		URL mafTestURL = getClass().getResource("test.maf");
		MAFIO mafIO = new MAFIO();
		System.out.println("test.maf URL: " + mafTestURL.toExternalForm() + " test.maf file " + mafTestURL.getFile());
	}
}

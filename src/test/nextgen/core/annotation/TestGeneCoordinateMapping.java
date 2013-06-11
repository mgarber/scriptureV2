package nextgen.core.annotation;


import nextgen.core.annotation.Gene;

public class TestGeneCoordinateMapping  extends junit.framework.TestCase{
	public void testPlusMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		Gene g = new Gene("chr1", 1000000, 1005200, "test_gene", 0, "+", exonStarts, exonEnds);
		
		int testPositionShort = 900000;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000001;
		assertEquals(1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000050;
		assertEquals(50, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001250;
		assertEquals(100, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001249;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001350;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001351;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002000;
		assertEquals(100+100, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002100;
		assertEquals(100+100+100, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005000;
		assertEquals(400, g.genomicToTranscriptPosition(testPositionShort));		

		testPositionShort = 1005199;
		assertEquals(400+199, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000000;
		assertEquals(0, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
	}
	
	public void testMinusMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		Gene g = new Gene("chr1", 1000000, 1005200, "test_gene", 0, "-", exonStarts, exonEnds);
		int transcriptLength = 0;
		for (int i = 0 ; i < exonStarts.length; i++) {
			transcriptLength += exonEnds[i] - exonStarts[i];
		}
		int testPositionShort = 900000;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000001;
		assertEquals(transcriptLength - 2, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000050;
		assertEquals(transcriptLength - 1 - 50, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001250;
		assertEquals(200+200+100-1, g.genomicToTranscriptPosition(testPositionShort));
		
		
		testPositionShort = 1001249;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001350;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1001351;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002000;
		assertEquals(200+200-1, g.genomicToTranscriptPosition(testPositionShort));

		testPositionShort = 1002001;
		assertEquals(200+199-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002100;
		assertEquals(200+100-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1002200;
		assertEquals(-1, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1005001;
		assertEquals(199-1, g.genomicToTranscriptPosition(testPositionShort));		

		testPositionShort = 1005199;
		assertEquals(0, g.genomicToTranscriptPosition(testPositionShort));
		
		testPositionShort = 1000000;
		assertEquals(200+200+100 + 100-1, g.genomicToTranscriptPosition(testPositionShort));
		
	}
	public void testPositivePositionMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		Gene g = new Gene("chr1", 1000000, 1005200, "test_gene", 0, "+", exonStarts, exonEnds);
		
		assertEquals(1000000, g.transcriptToGenomicPosition(0));
		
		int testPositionShort = 1000001;
		assertEquals(1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1000050;
		assertEquals(50, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));


		testPositionShort = 1002001;
		assertEquals(100+100+1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1002100;
		assertEquals(100+100+100, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005001;
		assertEquals(100+100+200+1, g.genomicToTranscriptPosition(testPositionShort));		
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005199;
		assertEquals(400+199, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		
	}	

	
	public void testNegativePositionMapping() throws Exception {
		int [] exonStarts = {1000000, 1001250,1002000, 1005000};
		int [] exonEnds = {1000100, 1001350,1002200, 1005200};
		Gene g = new Gene("chr1", 1000000, 1005200, "test_gene", 0, "-", exonStarts, exonEnds);
		
		assertEquals(1005200-1, g.transcriptToGenomicPosition(0));
		assertEquals(0, g.genomicToTranscriptPosition(1005200-1));
		
		assertEquals(200, g.genomicToTranscriptPosition(1002199));
		assertEquals(1002199, g.transcriptToGenomicPosition(200));

		int testPositionShort = 1000001;
		assertEquals(200+200+100+99-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1000050;
		assertEquals(200+200+100+50-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));


		testPositionShort = 1002001;
		assertEquals(200+199-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1002100;
		assertEquals(200+100-1, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005001;
		assertEquals(199-1, g.genomicToTranscriptPosition(testPositionShort));		
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		testPositionShort = 1005199;
		assertEquals(0, g.genomicToTranscriptPosition(testPositionShort));
		assertEquals(testPositionShort, g.transcriptToGenomicPosition(g.genomicToTranscriptPosition(testPositionShort)));
		
		
	}
	
	
}

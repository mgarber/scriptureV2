package broad.core.multiplealignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.RandomAccessFile;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;


import Jama.Matrix;

import junit.framework.TestCase;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;
import broad.core.multiplealignment.MAFAlignment;
import broad.core.multiplealignment.MAFIO;
import broad.core.multiplealignment.MultipleAlignment;
import broad.core.multiplealignment.MultipleAlignmentIOFactory;
import broad.core.multiplealignment.MAFAlignment.MAFMultipleAlignmentBlock;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;

public class MAFMultipleAlignmentTest extends TestCase {
	MAFAlignment test;
	
	protected void setUp() throws Exception{
		URL mafTestURL = getClass().getResource("test.maf");
		String idxFileName = mafTestURL.getFile() + ".index";
		System.out.println("test.maf index file: " + idxFileName);
		File idxFile = new File(idxFileName);
		if(!idxFile.exists()) {
			test = new MAFAlignment();
			test.createIndex(mafTestURL.getFile());
			test.writeIndex(idxFileName);
		}
		test = new MAFAlignment(idxFileName);
	}
	
	public void testAppend() throws Exception {
		
	}
	
	
	public void testConcatenation2() throws Exception {
		URL mafTestURL = getClass().getResource("test2.maf");	
		MAFIO mafio = new MAFIO();
		MAFAlignment test2 = mafio.load(mafTestURL.getFile());
		
		MAFAlignment problemAln = test2.getSubAlignment(128765340, 128765506, false);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
		problemAln.write(bw);
		bw.flush();
		
	}
	
	public void testLoadOnlyRequestedSeqs() throws Exception {
		URL mafTestURL = getClass().getResource("test.maf");	
		MAFIO mafio = new MAFIO();
		ArrayList<String> seqsToLoad = new ArrayList<String>(2);
		seqsToLoad.add("hg18");
		seqsToLoad.add("canFam2");
		test = mafio.load(mafTestURL.getFile(), seqsToLoad);
		assertEquals("Expecting hg18 and canFam2 but got: " + test.getAlignedSequenceIds(), 2, test.getAlignedSequenceIds().size());
		assertEquals(2, test.getAlignedSequences().size());
		
		seqsToLoad.add("panTro2");
		test = mafio.load(mafTestURL.getFile(), seqsToLoad);
		assertEquals(3, test.getAlignedSequenceIds().size());
		assertEquals(3, test.getAlignedSequences().size());
	}
	
	public void testCountAligningSequences() throws Exception{
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		
		MultipleAlignment combined = new MultipleAlignment();
		combined.setReferenceId(test.getReferenceId());
		
		MAFAlignment subAln = test.getSubAlignment(111323, 111333, false);
		combined.append(subAln);
		AlignedSequence seq = subAln.getAlignedSequence("hg18");
		assertEquals("TGAACAAGAT", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("dasNov1");
		assertEquals("TGAACAAGAA", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("xenTro2");
		assertEquals("-GAACAAGAT", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("oryLat1");
		assertEquals("----------", seq.getSequenceBases());
		assertEquals(8, subAln.getNumberOfAligningSequences());
		assertEquals(test.getAlignedSequenceIds().size(), subAln.getAlignedSequenceIds().size());
	}
	
	public void testCompressAppendedAlignment() throws Exception{
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		
		MultipleAlignment combined = new MultipleAlignment();
		combined.setReferenceId(test.getReferenceId());
		
		MAFAlignment subAln = test.getSubAlignment(111323, 111333, false);
		combined.append(subAln);
		AlignedSequence seq = subAln.getAlignedSequence("hg18");
		assertEquals("TGAACAAGAT", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("dasNov1");
		assertEquals("TGAACAAGAA", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("xenTro2");
		assertEquals("-GAACAAGAT", seq.getSequenceBases());
		
		//MultipleAlignment ma = test.getSubAlignment(111323, 111333, false);
		subAln = test.getSubAlignment(111486, 111494, false);
		combined.append(subAln);
		seq = subAln.getAlignedSequence("hg18");
		assertEquals("CGGAGGAA", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("dasNov1");
		assertEquals("CAGAGGAA", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("ornAna1");
		assertEquals("CAGAAGAA", seq.getSequenceBases());	
		
		subAln = test.getSubAlignment( 119880 + 519, 120475 + 5, false);
		combined.append(subAln);
		seq = subAln.getAlignedSequence("hg18");
		System.out.println("Precompression seq " + seq.getSequenceBases());
		assertEquals("c-ataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("rheMac2");
		assertEquals("ccataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		
		seq = combined.getReference();
		System.out.println("combined sequences: " + combined.getAlignedSequenceIds());
		assertNotNull("Combined reference sequence is null", seq);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
		combined.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));
		combined.write(bw);
		bw.flush();
		assertEquals("TGAACAAGAT" + "CGGAGGAA" + "c-ataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		
		combined.compress();
		seq = combined.getReference();
		assertEquals("TGAACAAGAT" + "CGGAGGAA" + "cataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
	}
	
	public void testGetSubAlignment() throws Exception {
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		
		MAFAlignment subAln = test.getSubAlignment(111323, 111333, false);
		AlignedSequence seq = subAln.getAlignedSequence("hg18");
		assertEquals("TGAACAAGAT", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("dasNov1");
		assertEquals("TGAACAAGAA", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("xenTro2");
		assertEquals("-GAACAAGAT", seq.getSequenceBases());
		
		//MultipleAlignment ma = test.getSubAlignment(111323, 111333, false);
		subAln = test.getSubAlignment(111486, 111494, false);
		seq = subAln.getAlignedSequence("hg18");
		assertEquals("CGGAGGAA", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("dasNov1");
		assertEquals("CAGAGGAA", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("ornAna1");
		assertEquals("CAGAAGAA", seq.getSequenceBases());	
		
		subAln = test.getSubAlignment( 119880 + 519, 120475 + 5, false);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
		test.write(bw);
		bw.flush();
		seq = subAln.getAlignedSequence("hg18");
		assertEquals("c-ataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("rheMac2");
		assertEquals("ccataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		System.out.println(seq.getSequenceBases());	
		
	}
	
	public void testCompressSubalignment() throws Exception {
		URL mafTestURL = getClass().getResource("test.maf");
		
		test.load(mafTestURL.getFile(), 119880 + 519, 120475 + 5);
		test.compress();
		AlignedSequence seq = test.getAlignedSequence("hg18");
		System.out.println(seq.getSequenceBases());	
		assertEquals("cataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		seq = test.getAlignedSequence("rheMac2");
		assertEquals("cataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		System.out.println(seq.getSequenceBases());	
		
		test.load(mafTestURL.getFile());
		MAFAlignment subAln = test.getSubAlignment( 119880 + 519, 120475 + 5, false);
		subAln.compress();
		seq = subAln.getAlignedSequence("hg18");
		System.out.println(seq.getSequenceBases());	
		assertEquals("cataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		seq = subAln.getAlignedSequence("rheMac2");
		assertEquals("cataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
	}

	public void testBlockParsing() throws Exception {
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		MAFMultipleAlignmentBlock first = test.getBlockIterator().next();
		assertEquals(64264, first.getReferenceStart());
		assertEquals(64264 + 553, first.getReferenceEnd());
		assertEquals(553, first.getReference().getUngappedLength());
	}
	
	public void testCompressOneBlock() throws Exception {
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		MAFMultipleAlignmentBlock first = test.getBlockIterator().next();
		int firstStart = first.getReferenceStart();
		int firstEnd   = first.getReferenceEnd();
		String firstAlignedSequence = first.getReference().getSequenceBases();
		System.out.println("first uncompressed: " + firstAlignedSequence);
		
		assertEquals(firstEnd - firstStart, first.getReferenceEnd() - first.getReferenceStart());
		assertTrue("gapped aligned sequence length should be larger than alignment length but alignment length is " + (first.getReferenceEnd() - first.getReferenceStart()) + " and first sequence length is " + firstAlignedSequence.length() , firstAlignedSequence.length() > first.getReferenceEnd() - first.getReferenceStart());
		assertTrue("First block sequence must contain gaps", firstAlignedSequence.contains("-"));
		
		first.compress();
		int firstNewStart = first.getReferenceStart();
		int firstNewEnd   = first.getReferenceEnd();
		String firstCompressedBases = first.getReference().getSequenceBases();
		System.out.println("first compressed: " + firstCompressedBases);
		
		assertEquals(firstStart, firstNewStart);
		assertEquals(firstEnd, firstNewEnd);
		assertTrue("First block compressed sequence must NOT contain gaps", !firstCompressedBases.contains("-"));
		assertEquals(first.length(), firstCompressedBases.length());
		
		
	}
	public void testCompress() throws Exception{
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		test.compress();
		
		AlignedSequence ref = test.getReference();
		String refBases = ref.getSequenceBases();
		assertTrue("Comressed reference should not contain a gap symbol", !refBases.contains("-"));
		assertEquals(test.length(), ref.getEnd() - ref.getStart());
	}

	public void testClosestOffset() throws Exception {
		assertEquals(117271,test.getClosestOffset(111294));
		assertEquals(117757,test.getClosestOffset(111295));
		assertEquals(117757,test.getClosestOffset(111295));  
		assertEquals(169333,test.getClosestOffset(125817));
	}
	
	public void testLoadSubalignment() throws IOException, ParseException {
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile(),111323,111333);
		AlignedSequence seq = test.getAlignedSequence("hg18");
		assertEquals("TGAACAAGAT", seq.getSequenceBases());
		seq = test.getAlignedSequence("dasNov1");
		assertEquals("TGAACAAGAA", seq.getSequenceBases());
		seq = test.getAlignedSequence("xenTro2");
		assertEquals("-GAACAAGAT", seq.getSequenceBases());
		
		//MultipleAlignment ma = test.getSubAlignment(111323, 111333, false);
		test.load(mafTestURL.getFile(), 111486, 111494);
		seq = test.getAlignedSequence("hg18");
		assertEquals("CGGAGGAA", seq.getSequenceBases());
		seq = test.getAlignedSequence("dasNov1");
		assertEquals("CAGAGGAA", seq.getSequenceBases());
		seq = test.getAlignedSequence("ornAna1");
		assertEquals("CAGAAGAA", seq.getSequenceBases());	
		
		test.load(mafTestURL.getFile(), 119880 + 520, 119880 + 530);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
		test.write(bw);
		bw.flush();
		seq = test.getAlignedSequence("hg18");
		System.out.println(seq.getSequenceBases());
		
		test.load(mafTestURL.getFile(), 119880 + 519, 120475 + 5);
		test.write(bw);
		bw.flush();
		seq = test.getAlignedSequence("hg18");
		assertEquals("c-ataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		seq = test.getAlignedSequence("rheMac2");
		assertEquals("ccataatccNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNtcgag", seq.getSequenceBases());
		System.out.println(seq.getSequenceBases());		
		
		test.load(mafTestURL.getFile(), 71257 + 330, 71257 + 331);
		seq = test.getAlignedSequence("hg18");
		assertEquals("A", seq.getSequenceBases());
	}
	
	public void testGetColumn() throws Exception{
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		test.encode();
		AlignedSequence seq = test.getAlignedSequence("hg18");
		//System.out.println(seq.getSequenceBases());
		Map<String, Short> col0 = test.getColumn(0);
		assertEquals(4, (short) col0.get("hg18"));
		
		Map<String, Short> col64264 = test.getColumn(64264);
		assertEquals(1, (short) col64264.get("hg18"));
		assertEquals(1, (short) col64264.get("rheMac2"));
		assertEquals(4, (short) col64264.get("panTro2"));
		
		Map<String, Short> col64265 = test.getColumn(64265);
		assertEquals(1, (short) col64265.get("hg18"));
		assertEquals(1, (short) col64265.get("rheMac2"));
		assertEquals(4, (short) col64265.get("panTro2"));
		
		Map<String, Short> col64266 = test.getColumn(64266);
		assertEquals(3, (short) col64266.get("hg18"));
		assertEquals(1, (short) col64266.get("rheMac2"));
		assertEquals(4, (short) col64266.get("panTro2"));
		
		
		Map<String, Short> col111493 = test.getColumn(111493);
		assertEquals(0, (short) col111493.get("hg18"));
		assertEquals(0, (short) col111493.get("rheMac2"));
		assertEquals(0, (short) col111493.get("panTro2"));
		assertEquals(0, (short) col111493.get("dasNov1"));
		assertEquals(0, (short) col111493.get("ornAna1"));
		assertEquals(1, (short) col111493.get("galGal3"));
		assertEquals(0, (short) col111493.get("anoCar1"));
		assertEquals(0, (short) col111493.get("oryLat1"));
		assertEquals(4, (short) col111493.get("xenTro2"));
		
		Map<String, Short> col111521 = test.getColumn(111521);
		assertEquals(0, (short) col111521.get("hg18"));
		assertEquals(2, (short) col111521.get("rheMac2"));
		assertEquals(0, (short) col111521.get("panTro2"));
		assertEquals(0, (short) col111521.get("dasNov1"));
		assertEquals(0, (short) col111521.get("ornAna1"));
		assertEquals(0, (short) col111521.get("galGal3"));
		assertEquals(0, (short) col111521.get("anoCar1"));
		assertEquals(2, (short) col111521.get("oryLat1"));
		assertEquals(4, (short) col111521.get("xenTro2"));
	}
	
	public void testGetColumnsAsVector() throws Exception{
		URL mafTestURL = getClass().getResource("test.maf");		
		test.load(mafTestURL.getFile());
		test.encode();
		AlignedSequence seq = test.getAlignedSequence("hg18");
		Map<String, Matrix> col = test.getColumnsAsVector(0, 10);
		assertEquals(21, col.keySet().size());
		//System.out.println("seqId list: " + col.keySet());
		Iterator<String> seqIt = col.keySet().iterator();
		while(seqIt.hasNext()) {
			String seqId = seqIt.next();
			Matrix seqData = col.get(seqId);
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 10; j++) {
					assertEquals(0d,seqData.get(i, j));
				}
			}
		}

		
		col = test.getColumnsAsVector(64264, 4);
		assertEquals(21, col.keySet().size());

		assertEquals(1d, col.get("hg18").get(1, 0));
		assertEquals(0d, col.get("hg18").get(2, 0));
		assertEquals(1d, col.get("rheMac2").get(1, 0));
		assertEquals(0d, col.get("panTro2").get(1,0));
		
		assertEquals(1d, col.get("hg18").get(1, 1));
		assertEquals(0d, col.get("hg18").get(2, 1));
		assertEquals(1d, col.get("rheMac2").get(1, 1));
		assertEquals(0d, col.get("panTro2").get(1,1));
		
		assertEquals(1d, col.get("hg18").get(3, 2));
		assertEquals(0d, col.get("hg18").get(2, 2));
		assertEquals(1d, col.get("rheMac2").get(1, 2));
		assertEquals(0d, col.get("panTro2").get(1,2));
		
		
		col = test.getColumnsAsVector(111493, 29);
		assertEquals(21, col.keySet().size());
		
		col.get("hg18").print(3, 0);
		assertEquals(1d, col.get("hg18").get(0, 0));
		assertEquals(1d, col.get("rheMac2").get(0, 0));
		assertEquals(1d, col.get("panTro2").get(0, 0));
		assertEquals(1d, col.get("dasNov1").get(0, 0));
		assertEquals(1d, col.get("ornAna1").get(0, 0));
		assertEquals(1d, col.get("galGal3").get(1,0));
		assertEquals(0d, col.get("galGal3").get(0,0));
		assertEquals(1d, col.get("anoCar1").get(0,0));
		assertEquals(1d, col.get("oryLat1").get(0, 0));
		assertEquals(0d, col.get("xenTro2").get(0, 0));
		
		assertEquals(1d, col.get("hg18").get(0, 28));
		assertEquals(1d, col.get("rheMac2").get(2, 28));
		assertEquals(0d, col.get("rheMac2").get(0, 28));
		assertEquals(1d, col.get("panTro2").get(0, 28));
		assertEquals(0d, col.get("panTro2").get(1, 28));
		assertEquals(1d, col.get("dasNov1").get(0, 28));
		assertEquals(1d, col.get("ornAna1").get(0, 28));
		assertEquals(1d, col.get("galGal3").get(0,28));
		assertEquals(0d, col.get("galGal3").get(1,28));
		assertEquals(1d, col.get("anoCar1").get(0,28));
		assertEquals(1d, col.get("oryLat1").get(2, 28));
		assertEquals(0d, col.get("oryLat1").get(0, 28));
		assertEquals(0d, col.get("xenTro2").get(0, 28));
	}
	
	public void testCreateIndex() throws java.io.IOException, ParseException{
		URL mafTestURL = getClass().getResource("test.maf");
		System.out.print("opening " + mafTestURL.getFile());
		MAFAlignment testAln = new MAFAlignment();
		System.out.print(" ... creating index ");
		testAln.createIndex(mafTestURL.getFile());
		System.out.println(" .Done");
		IntervalTree<Long> index = testAln.getIndex();
		RandomAccessFile raf = new RandomAccessFile(mafTestURL.getFile(), "r");
		try {
			Iterator<Node<Long>> it = index.iterator();
			while(it.hasNext()) {
				Node<Long> alnIndex = it.next();
				long offset = alnIndex.getValue();
				raf.seek(offset);
				String alignmentLine = raf.readLine();
				assertEquals("Alignment line must start with a ", "a", alignmentLine.substring(0,1));
			}
		} finally {
			raf.close();
		}
	}
	



}

package broad.core.motif;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.stat.clustering.Cluster;



import broad.core.error.ParseException;
import broad.core.motif.PositionWeightMatrix;
import broad.core.motif.PositionWeightMatrixIO;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;

import junit.framework.TestCase;

public class PositionWeightMatrixTest extends TestCase {
	 PositionWeightMatrix pwm1;
	 PositionWeightMatrix pwm2;
	 PositionWeightMatrix pwm3;
	 PositionWeightMatrix pwm4;
	 PositionWeightMatrix pwm5;
	 PositionWeightMatrix p53;
	 List<PositionWeightMatrix> toClusterList;
	 PositionWeightMatrix incompatibleWithClusters;
	 Sequence mdm2;
	 
	 void createPWMs() throws IOException, ParseException {
		 InputStream is = getClass().getResourceAsStream("test_5.pwm");
		 PositionWeightMatrixIO io = new PositionWeightMatrixIO();
		 io.load(is);
		 is.close();
		 assertEquals("Wrong number of matrices in test_5.pwm or error in pwmIO",5,io.getMatrices().size());
		 pwm1 = io.getMatrices().get(0);
		 pwm2 = io.getMatrices().get(1);
		 pwm3 = io.getMatrices().get(2);
		 pwm4 = io.getMatrices().get(3);
		 pwm5 = io.getMatrices().get(4);	
		 assertEquals("YCAGARGGAKTADTGTTCTGYAGATCTATTGATCTTCT",pwm1.getConsensus());
		 assertEquals("ARGCMAGTKYAGAWMYWTGGCAV",pwm5.getConsensus());
		 

	 }

	 void createClusterList() throws IOException, ParseException {
		 InputStream is = getClass().getResourceAsStream("test.clustering.pwm");
		 PositionWeightMatrixIO io = new PositionWeightMatrixIO();
		 io.load(is);
		 is.close();
		 int numToCluster = io.getMatrices().size();
		 toClusterList = new ArrayList<PositionWeightMatrix>(numToCluster);
		 assertEquals("Wrong number of matrices in test.clustering.pwm or error in pwmIO",6, numToCluster);
		 for (int i = 0; i < numToCluster - 1; i++) {
			 toClusterList.add(io.getMatrices().get(i));
		 }
		 incompatibleWithClusters = io.getMatrices().get(numToCluster - 1);
	 }
	 void createP53PWM() throws IOException, ParseException {
		 InputStream is = getClass().getResourceAsStream("p53.pwm");
		 PositionWeightMatrixIO io = new PositionWeightMatrixIO();
		 io.load(is);
		 is.close();	
		 p53 = io.getMatrices().get(0);
		 //assertEquals("BAD CONSENSUS","VBGRACATGYCCGGGCATGT", p53.getConsensus());
		     assertEquals("BAD CONSENSUS","RRACATGYCCGGGCATGTYB", p53.getConsensus());
		 for(int i = 0; i < p53.size(); i++) {
			 System.out.println(p53.get(i).getInformationContent()+"  ");
		 }
	 }
	 
	 void loadMDM2Sec() throws Exception {
		 FastaSequenceIO fsio = new FastaSequenceIO();
		 InputStream is = getClass().getResourceAsStream("mdm2.fa");
		 mdm2 = fsio.loadAll(is).get(0);
		 is.close();
	 }
	 
	 public void testInformationContent() throws Exception {
		 createPWMs();
		 assertEquals("Column one IC is wrong", 1.0073126603161544, pwm1.get(0).getInformationContent());
	 }
	 
	 public void testReverseComplement() throws Exception {
		 createPWMs();
		 assertEquals("AGAAGATCAATAGATCTRCAGAACAHTAMTCCYTCTGR", pwm1.reverseComplement().getConsensus());
		 assertEquals("GCTSRRAKAGSCTYRGMKGVC", pwm2.reverseComplement().getConsensus());
	 }
	 
	 public void testSetLeftAndRightStartOfInformation() throws Exception {
		 createPWMs();
		 pwm1.setAlignPos(7);
		 for(int i = 0; i < pwm1.size(); i++) {
			 System.out.println(i + "\t" + pwm1.get(i).getInformationContent() + "\t" + (i < pwm1.getLeftHighInfoStart() || i > pwm1.getRightHighInfoStart()  ? "NO" : "YES"));
		 }
		 //assertEquals("pwm1 start is wrong",0,pwm1.getLeftHighInfoStart());
	 }
	 
	 public void testTrimmByIC() throws Exception {
		 createPWMs();
		 pwm1.setAlignPos(6);
		 
		 PositionWeightMatrix trimmed  = pwm1.trimByInformationContent(2);
		 assertEquals("Trimmed size and left/right info start should be the same",trimmed.size(), pwm1.getRightHighInfoStart() - pwm1.getLeftHighInfoStart() + 1);
		 for(int i = 0; i < trimmed.size(); i++) {
			 assertEquals("IC at original position " + (i +  pwm1.getLeftHighInfoStart()) + " does not match corresp trimmed position", trimmed.get(i).getInformationContent(), pwm1.get(i + pwm1.getLeftHighInfoStart()).getInformationContent());
		 }
		 
		 pwm2.setAlignPos(7);
		 
		 trimmed  = pwm2.trimByInformationContent(2);
		 assertEquals("Trimmed size and left/right info start should be the same",trimmed.size(), pwm2.getRightHighInfoStart() - pwm2.getLeftHighInfoStart() + 1);
		 for(int i = 0; i < trimmed.size(); i++) {
			 assertEquals("IC at original position " + (i +  pwm2.getLeftHighInfoStart()) + " does not match corresp trimmed position", trimmed.get(i).getInformationContent(), pwm2.get(i + pwm2.getLeftHighInfoStart()).getInformationContent());
		 }
		 
		 trimmed = pwm1.trimByInformationContent(1.285);
		 assertEquals("Trimmed size by IC=0.6 was not correct",36,trimmed.size());
	 }
	 
	 public void testMatch() throws Exception {
		 int numOfMatches = 5;
		 double [] bg = {0.3, 0.2, 0.2, 0.3};
		 createP53PWM();
		 PositionWeightMatrix background = new PositionWeightMatrix("BG");
		 for(int i = 0; i < p53.size(); i++) {
			 background.addColumn(bg);
		 }
		 loadMDM2Sec();
		 
		 WindowSlider slider = mdm2.getSlider(p53.size(), p53.size() - 1);
		 double[] matches = new double[numOfMatches];
		 double[] matchPositions = new double[numOfMatches];
		 while(slider.hasNext()) {
			 SequenceRegion window = slider.next();
			 char[] windowChrs = window.getSequenceBases().toCharArray();
			 window.reverse();
			 char[] reverseWindowChrs = window.getSequenceBases().toCharArray();
			 
			 double score = Math.max(p53.getLogLikelihood(windowChrs) - background.getLogLikelihood(windowChrs), 
					 p53.getLogLikelihood(reverseWindowChrs) - background.getLogLikelihood(reverseWindowChrs))  ;
			 System.out.println(window.getStart() + " of " + mdm2.getLength());
			 if(score > 0) {
				 for(int i = 0; i < numOfMatches; i++) {
					 if (matches[i] < score ) {
						 matches[i] = score;
						 matchPositions[i] = window.getStart();
						 break;
					 }
				 }
			 }
		 }
		 
		 
		 for(int i = 0; i < numOfMatches; i++) {
			 System.out.println("pos " + matchPositions[i] + " " + matches[i]);
		 }
	 }
	 
	 public void testNextKmer() throws Exception {
		 int [] kmer = {0,0,0,0};
		 
		 
		 int [] next = PositionWeightMatrix.getNextKmer(kmer);
		 assertEquals("First entry is wrong",0,next[0]);
		 assertEquals("Second entry is wrong",0,next[1]);
		 assertEquals("Third entry is wrong",0,next[2]);
		 assertEquals("Fourth entry is wrong",1,next[3]);
		 
		 next = PositionWeightMatrix.getNextKmer(next);
		 assertEquals("First entry is wrong",0,next[0]);
		 assertEquals("Second entry is wrong",0,next[1]);
		 assertEquals("Third entry is wrong",0,next[2]);
		 assertEquals("Fourth entry is wrong",2,next[3]);
		 
		 int [] kmer2 = {0,0,0,3};
		 next = PositionWeightMatrix.getNextKmer(kmer2);
		 assertEquals("First entry is wrong",0,next[0]);
		 assertEquals("Second entry is wrong",0,next[1]);
		 assertEquals("Third entry is wrong",1,next[2]);
		 assertEquals("Fourth entry is wrong",0,next[3]);
		 
		 next = PositionWeightMatrix.getNextKmer(next);
		 assertEquals("First entry is wrong",0,next[0]);
		 assertEquals("Second entry is wrong",0,next[1]);
		 assertEquals("Third entry is wrong",1,next[2]);
		 assertEquals("Fourth entry is wrong",1,next[3]);
		 
		 int [] kmer3 = {0,3,3,2};
		 next = PositionWeightMatrix.getNextKmer(kmer3);
		 assertEquals("First entry is wrong",0,next[0]);
		 assertEquals("Second entry is wrong",3,next[1]);
		 assertEquals("Third entry is wrong",3,next[2]);
		 assertEquals("Fourth entry is wrong",3,next[3]);
		 
		 next = PositionWeightMatrix.getNextKmer(next);
		 assertEquals("First entry is wrong",1,next[0]);
		 assertEquals("Second entry is wrong",0,next[1]);
		 assertEquals("Third entry is wrong",0,next[2]);
		 assertEquals("Fourth entry is wrong",0,next[3]);
		 
		 next = PositionWeightMatrix.getNextKmer(next);
		 assertEquals("First entry is wrong",1,next[0]);
		 assertEquals("Second entry is wrong",0,next[1]);
		 assertEquals("Third entry is wrong",0,next[2]);
		 assertEquals("Fourth entry is wrong",1,next[3]);
		 
		 int [] kmer4 = {3,3,3,3};
		 next = PositionWeightMatrix.getNextKmer(kmer4);
		 assertEquals("First entry is wrong",0,next[0]);
		 assertEquals("Second entry is wrong",0,next[1]);
		 assertEquals("Third entry is wrong",0,next[2]);
		 assertEquals("Fourth entry is wrong",0,next[3]);
		 
		 
	 }
	 
	 public void testKLDistance() throws IOException, ParseException {
		 createClusterList();
		 KMeansPlusPlusClusterer<PositionWeightMatrix> clusterer = new KMeansPlusPlusClusterer<PositionWeightMatrix>(new Random());
		 List<Cluster<PositionWeightMatrix>> cluster =  clusterer.cluster(toClusterList, 2, 100);
		 BufferedWriter out = new BufferedWriter(new PrintWriter(System.err));
		 cluster.get(0).getCenter().write(out, NumberFormat.getNumberInstance());
		 cluster.get(1).getCenter().write(out, NumberFormat.getNumberInstance());
		 out.flush();
		 //System.err.print(cluster);
		 System.err.println(toClusterList.get(0).kullbackLeiber(toClusterList.get(1)));
	 }

}

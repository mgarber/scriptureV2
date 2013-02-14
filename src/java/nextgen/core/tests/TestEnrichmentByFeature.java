/**
 * 
 */
package nextgen.core.tests;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import java.util.Iterator;


import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.ProperPairFilter;
import nextgen.core.alignment.Alignment;

/**
 * @author prussell
 * Compare an RNA seq sample to background expression
 * Compare read count ratio in individual features to overall ratio
 */
public class TestEnrichmentByFeature {
	
	private ScanStatisticDataAlignmentModel backgroundData;
	private ScanStatisticDataAlignmentModel signalData;
	private ScanStatisticDataAlignmentModel backgroundGenomicData;
	private ScanStatisticDataAlignmentModel signalGenomicData;
	private double globalRatio;
	private Collection<Gene> genes;
	static Logger logger = Logger.getLogger(TestEnrichmentByFeature.class.getName());
	private TranscriptomeSpace transcriptomeSpace;
	private GenomicSpace genomeSpace;
	private static int DEFAULT_MAX_FRAGMENT_LENGTH = 2000;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	private static boolean DEFAULT_FULLY_CONTAINED = false;
	private boolean fullyContained;
	
	private TestEnrichmentByFeature(String backgroundBamFile, String signalBamFile, String bedFile, String chrSizeFile) throws IOException {
		this(backgroundBamFile, signalBamFile, bedFile, chrSizeFile, DEFAULT_MAX_FRAGMENT_LENGTH, DEFAULT_MAX_GENOMIC_SPAN, DEFAULT_FULLY_CONTAINED);
	}
	
	private TestEnrichmentByFeature(String backgroundBamFile, String signalBamFile, String bedFile, String chrSizeFile, int maxFragmentLength, int maxGenomicSpan, boolean fullyContained) throws IOException {
		this.fullyContained = fullyContained;
		this.genomeSpace = new GenomicSpace(chrSizeFile);
		this.transcriptomeSpace = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File (bedFile)));
		this.backgroundData = new ScanStatisticDataAlignmentModel(backgroundBamFile, this.transcriptomeSpace);
		this.backgroundGenomicData = new ScanStatisticDataAlignmentModel(backgroundBamFile, genomeSpace);
		this.signalGenomicData = new ScanStatisticDataAlignmentModel(signalBamFile, genomeSpace);
		this.signalData = new ScanStatisticDataAlignmentModel(signalBamFile, this.transcriptomeSpace);		
		this.genes = BEDFileParser.loadData(new File(bedFile));
		double signalNumReads = this.signalData.getGlobalNumReads();
		double backgroundNumReads = this.backgroundData.getGlobalNumReads();
		this.globalRatio = signalNumReads/backgroundNumReads;
		logger.info("Signal has " + signalNumReads + " total reads. Background has " + backgroundNumReads + " total reads.");
		logger.info("The global read count ratio is " + this.globalRatio);
		// Add read filters
		this.backgroundData.addFilter(new ProperPairFilter());
		this.backgroundData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
		this.backgroundData.addFilter(new FragmentLengthFilter(this.transcriptomeSpace, maxFragmentLength));
		this.signalData.addFilter(new ProperPairFilter());
		this.signalData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
		this.signalData.addFilter(new FragmentLengthFilter(this.transcriptomeSpace, maxFragmentLength));
		
		//this.signalGenomicData.addFilter(new ProperPairFilter());
		//this.signalGenomicData.addFilter(new GenomicSpanFilter(maxGenomicSpan));

		//this.backgroundGenomicData.addFilter(new ProperPairFilter());
		//this.backgroundGenomicData.addFilter(new GenomicSpanFilter(maxGenomicSpan));
	}
	
	/**
	 * Get the ratio of total counts among all the genes in collection
	 * @param geneSet The genes
	 * @return The overall ratio of counts
	 * @throws IOException
	 */
	private double getRatio(Collection<? extends Annotation> geneSet) throws IOException	{
		double backgroundCount = 0;
		double signalCount = 0;
		for(Annotation gene : geneSet) {
			backgroundCount += this.backgroundData.getCount(gene, this.fullyContained);
			signalCount += this.signalData.getCount(gene, this.fullyContained);
		}
		return signalCount/backgroundCount;
	}
		
	/**
	 * Get the ratio of counts for the transcribed regions of the gene
	 * @param gene The gene
	 * @return The ratio of counts
	 * @throws IOException
	 */
	private double getRatio(Annotation gene) throws IOException {
		return this.signalData.getCount(gene, this.fullyContained)/this.backgroundData.getCount(gene, this.fullyContained);
	}
	
	/**
	 * Get the ratio of counts normalized by the global ratio
	 * @param gene The gene
	 * @return The normalized ratio
	 * @throws IOException
	 */
	public double getEnrichment(Annotation gene) throws IOException {
		return getRatio(gene)/this.globalRatio;
	}
	
	/**
	 * Calculate enrichment off of counts
	 * @param signalCount Signal count
	 * @param backgroundCount Background count
	 * @return The ratio divided by global count ratio
	 */
	public double calculateEnrichment(double signalCount, double backgroundCount) {
		double ratio = signalCount / backgroundCount;
		return ratio / this.globalRatio;
	}

	/**
	 * Get the ratio of total counts among all the genes in collection normalized by the global ratio
	 * @param geneSet The genes
	 * @return The normalized overall ratio
	 * @throws IOException
	 */
	private double getEnrichment(Collection<? extends Annotation> geneSet) throws IOException {
		return getRatio(geneSet)/this.globalRatio;
	}
	
	
	/**
	 * Write global normalized ratios for each feature type
	 * @param outPerGeneFile Output file for enrichments per feature per gene
	 * @param outGlobalFile Output file for global enrichments
	 * @throws IOException
	 */
	private void computeAndWriteEnrichments(String outPerGeneFile, String outGlobalFile) throws IOException {
		FileWriter globalWriter = new FileWriter(outGlobalFile);
		FileWriter perGeneWriter = new FileWriter(outPerGeneFile);
		double allTranscriptsSignalCount = 0;
		double allTranscriptsBackgroundCount = 0;
		double allCDSsSignalCount = 0;
		double allCDSsBackgroundCount = 0;
		double all5UTRsSignalCount = 0;
		double all5UTRsBackgroundCount = 0;
		double all3UTRsSignalCount = 0;
		double all3UTRsBackgroundCount = 0;
		double allIntronsSignalCount = 0;
		double allIntronsBackgroundCount = 0;
		int numExpressed = 0;
		
		String header = "gene\t";
		header += "signal_count_transcript\tbackground_count_transcript\tenrichment_transcript\t";
		header += "signal_count_cds\tbackground_count_cds\tenrichment_cds\t";
		header += "signal_count_utr5t\tbackground_count_utr5\tenrichment_utr5\t";
		header += "signal_count_utr3\tbackground_count_utr3\tenrichment_utr3\t";
		header += "signal_count_introns\tbackground_count_introns\tenrichment_introns\n";
		perGeneWriter.write(header);
		
		logger.info("Calculating enrichments for " + this.genes.size() + " genes...");
		int numDone = 0;

		for(Gene gene : this.genes) {

			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Finished " + numDone + " genes of which " + numExpressed + " are significantly expressed.");
			}
			
			String geneName = gene.getName();
			if(!this.backgroundGenomicData.isExpressed(gene, 0.05)) {
				continue;
			}
			
			numExpressed++;

			Gene cds = gene.getCDS();
			Gene utr5 = gene.get5UTRGene();
			Gene utr3 = gene.get3UTRGene();
			Gene introns = gene.getIntrons();
			
			String transcriptEnrichment = "NA";
			String cdsEnrichment = "NA";
			String utr5Enrichment = "NA";
			String utr3Enrichment = "NA";
			String intronEnrichment = "NA";
			
			String transcriptSignalCountString = "NA";
			String cdsSignalCountString = "NA";
			String utr5SignalCountString = "NA";
			String utr3SignalCountString = "NA";
			String intronSignalCountString = "NA";
			
			String transcriptBackgroundCountString = "NA";
			String cdsBackgroundCountString = "NA";
			String utr5BackgroundCountString = "NA";
			String utr3BackgroundCountString = "NA";
			String intronBackgroundCountString = "NA";
						
			
			try {
				double transcriptSignalCount = this.signalData.getCount(gene, this.fullyContained);
				double transcriptBackgroundCount = this.backgroundData.getCount(gene, this.fullyContained);
				transcriptSignalCountString = Double.valueOf(transcriptSignalCount).toString();
				transcriptBackgroundCountString = Double.valueOf(transcriptBackgroundCount).toString();
				allTranscriptsSignalCount += transcriptSignalCount;
				allTranscriptsBackgroundCount += transcriptBackgroundCount;
				transcriptEnrichment = Double.valueOf(calculateEnrichment(transcriptSignalCount,transcriptBackgroundCount)).toString();
			} catch(Exception e) {}
			
			try {
				double cdsSignalCount = this.signalData.getCount(cds, this.fullyContained);
				double cdsBackgroundCount = this.backgroundData.getCount(cds, this.fullyContained);
				cdsSignalCountString = Double.valueOf(cdsSignalCount).toString();
				cdsBackgroundCountString = Double.valueOf(cdsBackgroundCount).toString();
				allCDSsSignalCount += cdsSignalCount;
				allCDSsBackgroundCount += cdsBackgroundCount;
				cdsEnrichment = Double.valueOf(calculateEnrichment(cdsSignalCount,cdsBackgroundCount)).toString();
			} catch(Exception e) {}
			
			if(cds != null) {
				if(cds.getSize() > 2) {
					try {
						double utr5SignalCount = this.signalData.getCount(utr5, this.fullyContained);
						double utr5BackgroundCount = this.backgroundData.getCount(utr5, this.fullyContained);
						utr5SignalCountString = Double.valueOf(utr5SignalCount).toString();
						utr5BackgroundCountString = Double.valueOf(utr5BackgroundCount).toString();
						all5UTRsSignalCount += utr5SignalCount;
						all5UTRsBackgroundCount += utr5BackgroundCount;
						utr5Enrichment = Double.valueOf(calculateEnrichment(utr5SignalCount,utr5BackgroundCount)).toString();
					} catch(Exception e) {}
					try {
						double utr3SignalCount = this.signalData.getCount(utr3, this.fullyContained);
						double utr3BackgroundCount = this.backgroundData.getCount(utr3, this.fullyContained);
						utr3SignalCountString = Double.valueOf(utr3SignalCount).toString();
						utr3BackgroundCountString = Double.valueOf(utr3BackgroundCount).toString();
						all3UTRsSignalCount += utr3SignalCount;
						all3UTRsBackgroundCount += utr3BackgroundCount;
						utr3Enrichment = Double.valueOf(calculateEnrichment(utr3SignalCount,utr3BackgroundCount)).toString();
					} catch(Exception e) {}
				}
			}
			
			if(introns != null) {
				double intronSignalCount = 0;
				double intronBackgroundCount = 0;
				for(Annotation intron : introns.getExonSet()) {
					//logger.info("Intron " + intron + " for gene " + gene);
					try {
						intronSignalCount += this.signalGenomicData.getCount(intron, true);
						double count = this.backgroundGenomicData.getCount(intron, true);
						intronBackgroundCount += count;
						//logger.info("Found " + count);
					
						/*
						if (count > 0.0) {
							Iterator<Alignment> itr = this.backgroundGenomicData.getOverlappingReads(intron, true);
							while  (itr.hasNext()) {
								Alignment next = itr.next();
								System.out.println(next.getReadName());
								System.out.println(next.getFragment(genomeSpace).iterator().next().toUCSC());
							}
						}*/
						
					} catch (Exception e) {
						logger.error(e.toString());					
					}
				}
				intronSignalCountString = Double.valueOf(intronSignalCount).toString();
				intronBackgroundCountString = Double.valueOf(intronBackgroundCount).toString();
				allIntronsSignalCount += intronSignalCount;
				allIntronsBackgroundCount += intronBackgroundCount;
				intronEnrichment = Double.valueOf(calculateEnrichment(intronSignalCount,intronBackgroundCount)).toString();
			}
			
			String line = geneName + "\t";
			line += transcriptSignalCountString + "\t" + transcriptBackgroundCountString + "\t" + transcriptEnrichment + "\t";
			line += cdsSignalCountString + "\t" + cdsBackgroundCountString + "\t" + cdsEnrichment + "\t";
			line += utr5SignalCountString + "\t" + utr5BackgroundCountString + "\t" + utr5Enrichment + "\t";
			line += utr3SignalCountString + "\t" + utr3BackgroundCountString + "\t" + utr3Enrichment + "\t";
			line += intronSignalCountString + "\t" + intronBackgroundCountString + "\t" + intronEnrichment + "\t";
			perGeneWriter.write(line + "\n");

		}
		
		globalWriter.write("Features\tSignalCount\tBackgroundCount\tEnrichment\n");
		globalWriter.write("Transcripts\t" + allTranscriptsSignalCount + "\t" + allTranscriptsBackgroundCount + "\t" + calculateEnrichment(allTranscriptsSignalCount, allTranscriptsBackgroundCount) + "\n");
		globalWriter.write("CDSs\t" + allCDSsSignalCount + "\t" + allCDSsBackgroundCount + "\t" + calculateEnrichment(allCDSsSignalCount, allCDSsBackgroundCount) + "\n");
		globalWriter.write("5'UTRs\t" + all5UTRsSignalCount + "\t" + all5UTRsBackgroundCount + "\t" + calculateEnrichment(all5UTRsSignalCount, all5UTRsBackgroundCount) + "\n");
		globalWriter.write("3'UTRs\t" + all3UTRsSignalCount + "\t" + all3UTRsBackgroundCount + "\t" + calculateEnrichment(all3UTRsSignalCount, all3UTRsBackgroundCount) + "\n");
		globalWriter.write("Introns\t" + allIntronsSignalCount + "\t" + allIntronsBackgroundCount + "\t" + calculateEnrichment(allIntronsSignalCount, allIntronsBackgroundCount) + "\n");
		globalWriter.close();
		perGeneWriter.close();
		logger.info("Wrote enrichments per gene and feature to file " + outPerGeneFile + ".");
		logger.info("Wrote global enrichments to file " + outGlobalFile + ".");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Background bam file", true);
		p.addStringArg("-s", "Signal bam file", true);
		p.addStringArg("-g", "Genes bed file", true);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addStringArg("-ot", "Output file for table of all genes", true);
		p.addStringArg("-og", "Output file for global enrichments by feature type", true);
		p.addIntArg("-mf", "Max fragment length", false, Integer.valueOf(DEFAULT_MAX_FRAGMENT_LENGTH));
		p.addIntArg("-mg", "Max genomic span", false, Integer.valueOf(DEFAULT_MAX_GENOMIC_SPAN));
		p.addBooleanArg("-fc", "Fully contained reads in features", false, Boolean.valueOf(DEFAULT_FULLY_CONTAINED));
		p.parse(args);
		String backgroundFile = p.getStringArg("-b");
		String signalFile = p.getStringArg("-s");
		String geneBedFile = p.getStringArg("-g");
		String chrSizeFile = p.getStringArg("-c");
		String outTable = p.getStringArg("-ot");
		String outGlobal = p.getStringArg("-og");
		int maxFragmentLength = p.getIntArg("-mf");
		int maxGenomicSpan = p.getIntArg("-mg");
		boolean fullyContained = p.getBooleanArg("-fc");
		
		TestEnrichmentByFeature ebf = new TestEnrichmentByFeature(backgroundFile, signalFile, geneBedFile, chrSizeFile, maxFragmentLength, maxGenomicSpan, fullyContained);
		
		ebf.computeAndWriteEnrichments(outTable, outGlobal);
		logger.info("All done.");

	}

}

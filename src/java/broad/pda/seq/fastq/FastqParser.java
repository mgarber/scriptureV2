package broad.pda.seq.fastq;

import java.io.*;
import java.util.*;

import org.apache.log4j.Logger;


public class FastqParser implements Iterator<FastqSequence>{
	static Logger logger = Logger.getLogger(FastqParser.class.getName());
	Collection<FastqSequence> sequences;
	File fastqFile;
	BufferedReader reader;
	private boolean done;
	private int numberOfSeq;
	
	/**
	 * @deprecated This constructure is highly discouraged as it opens a reader. Use the empty constructor instead and set the 
	 * file to read, then call start.
	 * @param fastqFile
	 * @throws FileNotFoundException
	 */
			
	public FastqParser(File fastqFile) throws FileNotFoundException{
		this.fastqFile=fastqFile;
		reader=new BufferedReader(new InputStreamReader(new FileInputStream(fastqFile)));
		done=false;
	}
	
	public FastqParser() {
		super();
	}
	
	public void start(File fastqParser) throws FileNotFoundException {
		this.fastqFile = fastqParser;
		reader=new BufferedReader(new InputStreamReader(new FileInputStream(fastqFile)));
		done=false;
	}
	
	public void start (BufferedReader br) {
		reader=br;
		done=false;
	}
	
	public void convertToFasta(String save)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(this.fastqFile)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
        		String sequence=seq.getSequence();
        		//int num=10;
        		//String lastNBps=getLastBps(sequence, num);
    			//String polyN=polyN("T", num);
    			//if(lastNBps.equalsIgnoreCase(polyN)){System.err.println(sequence);}
    			writer.write(seq.toFasta());
        	}
        }
        writer.close();
	}
	
	public void convertToNumberedFasta(String save)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(this.fastqFile)));
		String nextLine;
    	int i=1;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
        		String sequence=seq.getSequence();
    			writer.write(seq.toFasta(i));
    			i++;
        	}
        }
        writer.close();
	}
	
	public Collection<FastqSequence> parse(File file)throws IOException{
		Collection rtrn=new ArrayList();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine,thirdLine, fourthLine);
        		String sequence=seq.getSequence();
        		//int num=10;
        		//String lastNBps=getLastBps(sequence, num);
    			//String polyN=polyN("T", num);
    			//if(lastNBps.equalsIgnoreCase(polyN)){System.err.println(sequence);}
    			
        		rtrn.add(seq);
        	}
        	
        	
        }
        return rtrn;
	}
	
	
	private String getLastBps(String sequence, int num){
		return sequence.substring(sequence.toCharArray().length-num);
	}
	
	private String polyN(String letter, int num){
		String rtrn="";
		for(int i=0; i<num; i++){rtrn=rtrn+letter;}
		return rtrn;
	}
	
	public Collection<FastqSequence> getSequences() throws IOException{
		if(this.sequences!=null){
		return this.sequences;
		}
		else{this.sequences=this.parse(fastqFile); return this.sequences;}
	}

	public File[] writeChunks(String save, int chunkSize) throws IOException {
		Collection<File> rtrn=new TreeSet<File>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(this.fastqFile)));
		String nextLine;
		FileWriter writer=new FileWriter(save+".0.fq");
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		if(i%chunkSize ==0){
        			writer.close();
        			rtrn.add(new File(save+"."+(i/chunkSize)+".fq")); 
        			writer=new FileWriter(save+"."+(i/chunkSize)+".fq"); 
        			System.err.println(i);
        		}
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
        		//String sequence=seq.getSequence();
        		writer.write(seq.toFastq()+"\n");
        		i++;
        		if(i % 100000 == 0){System.err.println("Iterating.. "+ i);}
        	}
        	
        	
        }
        this.numberOfSeq=i;
        File[] files=new File[rtrn.size()];
        int counter=0;
        for(File file: rtrn){files[counter++]=file;}
		return files;
	}
	
	public int getNumberOfSequences(){
		if(this.numberOfSeq>0){return this.numberOfSeq;}
		else{
			try{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(this.fastqFile)));
			String nextLine;
			int i=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	        	if(nextLine.startsWith("@")){i++;}
	        }
	        this.numberOfSeq=i;
	        reader.close();
			}catch(IOException ex){}
			return this.numberOfSeq;
		}
	}

	public boolean hasNext() {
		return !done;
	}

	public FastqSequence next() {
		FastqSequence seq = null;
		try{
			String nextLine;
			
	        if ((nextLine = reader.readLine()) != null) {
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
				String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
	        }
	        done = nextLine == null ;
		}catch(Exception ex){
			this.done=true; 
			logger.error("Exception thrown while reading fastq file",ex);
		}

		return seq;
	}

	public void close() throws IOException{
		reader.close();
	}

	public void remove() {}
	
}

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;

public class predictBatchBLAST {

	
	public static HashMap<String, String> getFasta(String path) throws IOException{
		HashMap<String, String> fasta_dict = new HashMap<>();
		String prot_id="";
		String prot_seq="";
		
		BufferedReader br = new BufferedReader(new FileReader(path));	
		String line = br.readLine();
		 while (true) {
			 if(line==null)//for terminating
				 break;
			 if(line.startsWith(">")==true){
				 prot_seq="";
				 StringTokenizer st1 = new StringTokenizer(line, "|");//split
			    	prot_id=st1.nextToken();
			    	prot_id=st1.nextToken();
			    	fasta_dict.put(prot_id, line);
					 line = br.readLine();
			    	 while(line.startsWith(">")==false){//until > add sequence 
				    		prot_seq+=line;
					    	line = br.readLine();
					    	if(line==null)
								 break;
				    	}
			 }			
		    	

			    }	 	
		 br.close();
		return fasta_dict;
		
	}
	
	public static void main(String[] args, Vector<String> ecnums, long time,String ROOTPATH, List<String> test_ids, String fastaFile) throws IOException, InterruptedException {
		int k = 5;
		int evalue = 20;
		
		   
	        String method = "blast";
	        HashMap<String, String> fasta_dict = new HashMap<>();
			

		
		for (int i =0; i< ecnums.size(); i++) {  
    		HashMap<String, List<List<String>>> simHashHash = new HashMap<>();
    		
			String path = ROOTPATH+ File.separator+ecnums.get(i)+File.separator + method;
			String testpath = ROOTPATH + File.separator+ "testResult"+ File.separator+ time + File.separator+ ecnums.get(i)+File.separator+method;

		File workdir = new File(testpath);
		workdir.mkdirs();	
		List<String> pos = Files.readAllLines(Paths.get(ROOTPATH + File.separator +ecnums.get(i) + "/positive.ids"));
		List<String> neg =  Files.readAllLines(Paths.get(ROOTPATH + File.separator +ecnums.get(i) + "/negative.ids"));
		String predFile = workdir + File.separator+ecnums.get(i)+".preds";
		String confFile = workdir + File.separator+ecnums.get(i)+".confs";

		String posPredFile = path + File.separator+"ppreds.txt";
		String negPredFile = path + File.separator+"npreds.txt";

		 String[] cmdArray = new String[11];

		   // first argument is the program we want to open
		 	
		   cmdArray[0] = "lib/ncbi-blast-2.7.1+/bin/blastp";

		   cmdArray[1] = "-query";

		   cmdArray[2] = fastaFile;
		   
		   cmdArray[3] = "-db";
		   
		   cmdArray[4] = ROOTPATH+ "/" + ecnums.get(i)+ File.separator+ "blast" + File.separator+  ecnums.get(i)+".blastdb";
		   
		   cmdArray[5] = "-outfmt"; 
		   
		   cmdArray[6] = "6"; 
 
		   cmdArray[7] = "-out"; 
		   
		   cmdArray[8] = workdir + File.separator+ "blast.out" ; 
		
		   cmdArray[9] = "-evalue"; 
		   
		   cmdArray[10] = String.valueOf(evalue); 
		   
		   String cmd = "lib/ncbi-blast-2.7.1+/bin/blastp -query " +fastaFile +" -db "
			 		+ ROOTPATH+ "/" + ecnums.get(i)+ File.separator+ "blast" + File.separator+  ecnums.get(i)+".blastdb"
	 				+ " -outfmt 6 -out " +workdir + File.separator+"blast.out -evalue " +String.valueOf(evalue);
	         ProcessBuilder pb = new ProcessBuilder(cmd.split(" ")); 
	        
	         Process process = pb.start();
	         
	        
	         try {
	             int exitValue = process.waitFor();
	         } catch (InterruptedException e) {
	             // TODO Auto-generated catch block
	        	 System.out.print("blastp is not working!");
	             e.printStackTrace();
	         }
		   
	         Blast blast = new Blast();
        	 Vector<Double> preds = new Vector<>();

	         List<String> blastLines = Files.readAllLines(Paths.get(cmdArray[8]));
	         if(blastLines.size()==0){
	        	 preds.add(0.0);
	         }
	         else{
	        	 simHashHash = blast.parseTabBlast(cmdArray[8]);
	 	 		preds = new Vector<>();
	 			double pred;

	 	 		for(int m = 0; m < test_ids.size(); m++){
	 				pred = blast.blastknn(simHashHash.get(test_ids.get(m)), pos, neg, k);
	 				preds.add(pred);

	 			}
	         }
	 				
			PrintWriter final_file = new PrintWriter(predFile, "UTF-8");
			for(int a=0; a<preds.size(); a++)
			final_file.println(preds.get(a));
			final_file.close();
		
		utils u = new utils();	
		u.calculateConfidence(posPredFile, negPredFile, predFile, confFile);		

		}

	}

}
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
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

public class predictBatchSPMAP {

		public static void main(String[] args, Vector<String> ecnums, long time, String ROOTPATH, List<String> test_ids, String fastaFile) throws IOException, InterruptedException {
		int sigTh  = -15;
		int subseqlen = 5;
		String predictProg = "lib/svmlight/svm_classify";
	    String method = "spmap";	
		for (int i =0; i< ecnums.size(); i++) {  
			File workdir = new File(ROOTPATH + File.separator+ "testResult"+ File.separator+ time + File.separator+ ecnums.get(i));
			workdir.mkdirs();
			String path = ROOTPATH+ File.separator+ecnums.get(i)+File.separator + method;
			String testpath = workdir+File.separator+method;
	        String modelfile = path+File.separator+ "model.svm";

	
			workdir = new File(testpath);
			workdir.mkdirs();	
			
			String predFile = testpath + File.separator+ecnums.get(i)+".preds";
			String confFile = testpath + File.separator+ecnums.get(i)+".confs";
	        String batchVect = testpath + File.separator+"test.vec";

			String posPredFile = path + File.separator+"ppreds.txt";
			String negPredFile = path + File.separator+"npreds.txt";

			seq2vectPSSMtest seqtest = new seq2vectPSSMtest();

			seqtest.calculateVectors(sigTh,subseqlen,ecnums.get(i),test_ids,fastaFile,time, ROOTPATH);
			
			
			String cmd = predictProg + " " + batchVect + " " +  modelfile + " " + predFile;
	         ProcessBuilder pb = new ProcessBuilder(cmd.split(" ")); 
	         Process process = pb.start();
	         
	        
	         try {
	             int exitValue = process.waitFor();
	         } catch (InterruptedException e) {
	             // TODO Auto-generated catch block
	        	 System.out.print("Svm_classify is not working!");
	             e.printStackTrace();
	         }
	        
	        utils u = new utils();	
			u.calculateConfidence(posPredFile, negPredFile, predFile, confFile);
		}
	}
}
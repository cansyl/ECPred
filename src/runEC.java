import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;

public class runEC {

	public HashMap<String, Vector<Vector<String>>> predictions (String[] args, String ROOTPATH , Vector<String> ecnums, long time, HashMap<String,Vector<Vector<String>>> predictions, List<String> idlist, String fastaFile) throws IOException, InterruptedException{
		predictBatchSPMAP obj = new predictBatchSPMAP();
		obj.main(args, ecnums, time,ROOTPATH, idlist, fastaFile);
		
		predictBatchBLAST obj2 = new predictBatchBLAST();
		obj2.main(args, ecnums, time,ROOTPATH,idlist,fastaFile);
		
		predictBatchPEPSTATS obj3 = new predictBatchPEPSTATS();

		obj3.main(args, ecnums, time,ROOTPATH,fastaFile);
		
		
		BufferedReader br = new BufferedReader(new FileReader(ROOTPATH.substring(0, ROOTPATH.length()-3)+"/subclasses/thresholds.txt"));	
		String line;
		HashMap<String, Double> thresholds = new HashMap<>();
		while((line = br.readLine())!=null){
			StringTokenizer st1 = new StringTokenizer(line, "\t");//split
	    	String ECclass=st1.nextToken();
	    	String thres = st1.nextToken();
	    	Double threshold = Double.parseDouble(thres);	    	
	    	thresholds.put(ECclass, threshold);
		}
		
		for (int i =0; i< ecnums.size(); i++) {  
			String workdir = ROOTPATH+File.separator+ecnums.get(i);
			String testDir = ROOTPATH+File.separator+"testResult"+File.separator+ time+File.separator+ecnums.get(i) +File.separator+"preds";
			final File folder1 = new File(testDir);
			folder1.mkdirs();
			for (final File fileEntry : folder1.listFiles()) 
				fileEntry.delete();	
			List<String> spreds = Files.readAllLines(Paths.get(ROOTPATH+File.separator+"testResult" + File.separator + time +File.separator+ecnums.get(i) + File.separator + "spmap" + File.separator+ecnums.get(i)+".confs"));
			List<String> bpreds = Files.readAllLines(Paths.get(ROOTPATH+File.separator+"testResult" + File.separator + time +File.separator+ecnums.get(i) + File.separator + "blast" + File.separator+ecnums.get(i)+".confs"));
			List<String> ppreds = Files.readAllLines(Paths.get(ROOTPATH+File.separator+"testResult" + File.separator + time+File.separator+ecnums.get(i) + File.separator + "pepstats" + File.separator+ecnums.get(i)+".confs"));
			List<String> weights = Files.readAllLines(Paths.get(ROOTPATH + File.separator +ecnums.get(i) + File.separator+ "weights.txt"));

			Vector<String> combined = new Vector<>();
			for(int j = 0 ; j < spreds.size(); j++){
				Double comb = Double.parseDouble(spreds.get(j)) * Double.parseDouble(weights.get(0)) + Double.parseDouble(bpreds.get(j)) * Double.parseDouble(weights.get(1)) 
						+ Double.parseDouble(ppreds.get(j)) * Double.parseDouble(weights.get(2)) ;
				combined.add(String.valueOf(comb));
			}
	
			DecimalFormat df = new DecimalFormat();
			df.setMaximumFractionDigits(2);
			BufferedWriter final_file = new BufferedWriter(new FileWriter(ROOTPATH+File.separator+"testResult"+File.separator + time+File.separator+ecnums.get(i) +File.separator + ecnums.get(i)+"_preds.txt",false));
			for(int j = 0 ; j < idlist.size(); j++){
				final_file.write(df.format(Double.parseDouble(combined.get(j))) + "\n");
			}
			final_file.close();
			for(int j = 0 ; j < idlist.size(); j++){
				BufferedWriter final_file1 = new BufferedWriter(new FileWriter(testDir+File.separator+  idlist.get(j) + ".preds",true));
				final_file1.write(ecnums.get(i) + "\t" + spreds.get(j) + "\t" + bpreds.get(j) + "\t" + ppreds.get(j) + "\t" + combined.get(j) + "\n");
				final_file1.close();
			}		
		}
		if(ecnums.contains("1.-.-.-")){
			Vector<Vector<String>> allPreds = new Vector<>();
			for (int i =0; i< ecnums.size(); i++) {  
				String workdir = ROOTPATH+File.separator+ecnums.get(i);
				Vector<String> pred = new Vector<>();
				  br = new BufferedReader(new FileReader(ROOTPATH+File.separator+"testResult"+File.separator + time +File.separator+ecnums.get(i)+File.separator + ecnums.get(i)+"_preds.txt"));	
					while((line = br.readLine())!=null){
						pred.add(line);
					}
				allPreds.add( pred);

			}
			
			for(int i=0 ; i < idlist.size(); i++){
				Vector<Vector<String>> predswithScore = new Vector<>();
				Vector<String> preds = new Vector<>();
				double maxPred =0.0;
				String mainClass = null ;
				for(int j = 0 ; j < allPreds.size(); j++){
					if(Double.parseDouble(allPreds.get(j).get(i))> maxPred){
						maxPred = Double.parseDouble(allPreds.get(j).get(i));
						mainClass = String.valueOf((j+1)) + ".-.-.-";
					}	
				}
				if(maxPred<0.4){
					preds.add("non");
					preds.add("");
					predswithScore.add(preds);
					predictions.put(idlist.get(i) , predswithScore);
				}
				else if (maxPred >= thresholds.get(mainClass)){
					preds.add(mainClass);
					preds.add(String.valueOf(maxPred));
					predswithScore.add(preds);
					predictions.put(idlist.get(i) , predswithScore);
				}
				else{
					preds.add("nop");
					preds.add("");
					predswithScore.add(preds);
					predictions.put(idlist.get(i) , predswithScore);
				}
					
			}
		}
		
		else{
			Vector<Vector<String>> predswithScore = new Vector<>();
			double maxPred = 0.0;
			String predClass ="";
			Vector<String> preds = new Vector<>();
			for (int i =0; i< ecnums.size(); i++) {  
				List<String> pred = Files.readAllLines(Paths.get(ROOTPATH+File.separator+"testResult"+File.separator + time+File.separator+ecnums.get(i) +File.separator + ecnums.get(i)+"_preds.txt"));
				if(Double.parseDouble(pred.get(0)) >= maxPred){
					maxPred = Double.parseDouble(pred.get(0));
					predClass = ecnums.get(i);
				}
			}
			if (maxPred >= thresholds.get(predClass)){
				preds.add(predClass);
				preds.add(String.valueOf(maxPred));
				predictions.get(idlist.get(0)).add(preds);
			}
			else{
				preds.add("nop");
				preds.add("");
				predictions.get(idlist.get(0)).add(preds);
			}			
		}
		return predictions;
		
	}

}
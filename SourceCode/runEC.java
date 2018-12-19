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

	public HashMap<String, HashMap<String, Double>> predictions(String[] args, String ROOTPATH, Vector<String> ecnums, long time, HashMap<String, HashMap<String, Double>> predictions, List<String> idlist, String fastaFile, String tempDir, String method) throws IOException, InterruptedException{	
	if (method.equals("spmap"))
	    {
	      predictBatchSPMAP obj = new predictBatchSPMAP();
	      predictBatchSPMAP.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
	    }
	    else if (method.equals("blast"))
	    {
	      predictBatchBLAST obj2 = new predictBatchBLAST();
	      predictBatchBLAST.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
	    }
	    else if (method.equals("pepstats"))
	    {
	      predictBatchPEPSTATS obj3 = new predictBatchPEPSTATS();
	      predictBatchPEPSTATS.main(args, ecnums, time, ROOTPATH, fastaFile, tempDir);
	    }
	    else if (method.equals("weighted"))
	    {
	      predictBatchSPMAP obj = new predictBatchSPMAP();
	      predictBatchSPMAP.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
	      
	      predictBatchBLAST obj2 = new predictBatchBLAST();
	      predictBatchBLAST.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
	      
	      predictBatchPEPSTATS obj3 = new predictBatchPEPSTATS();
	      predictBatchPEPSTATS.main(args, ecnums, time, ROOTPATH, fastaFile, tempDir);
	    }
		HashMap<String, Double> NHpredictions = new HashMap<>();
	    BufferedReader br = new BufferedReader(new FileReader(ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/subclasses/thresholds.txt"));
	    
	    HashMap<String, Double> thresholds = new HashMap();
	    String line;
	    while ((line = br.readLine()) != null)
	    {
	      StringTokenizer st1 = new StringTokenizer(line, "\t");
	      String ECclass = st1.nextToken();
	      String thres = st1.nextToken();
	      Double threshold = Double.valueOf(Double.parseDouble(thres));
	      thresholds.put(ECclass, threshold);
	    }
	    for (int i = 0; i < ecnums.size(); i++)
	    {
	      String workdir = ROOTPATH + File.separator + (String)ecnums.get(i);
	      String testDir = tempDir + File.separator + "testResult" + File.separator + time + File.separator + (String)ecnums.get(i) + File.separator + "preds";
	      File folder1 = new File(testDir);
	      folder1.mkdirs();

	      for (final File fileEntry : folder1.listFiles()) 
				fileEntry.delete();
	      Vector<String> combined = new Vector();
	      if (method.equals("spmap"))
	      {
	    	  List<String> spreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + "spmap" + File.separator + (String)ecnums.get(i) + ".confs", new String[0]));
	        for (int j = 0; j < spreds.size(); j++)
	        {
	          Double comb = Double.valueOf(Double.parseDouble(spreds.get(j)));
	          combined.add(String.valueOf(comb));
	        }
	      }
	      else if (method.equals("blast"))
	      {
	    	  List<String> bpreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + "blast" + File.separator + (String)ecnums.get(i) + ".confs", new String[0]));
	        for (int j = 0; j < bpreds.size(); j++)
	        {
	          Double comb = Double.valueOf(Double.parseDouble(bpreds.get(j)));
	          combined.add(String.valueOf(comb));
	        }
	      }
	      else if (method.equals("pepstats"))
	      {
	    	  List<String> ppreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + "pepstats" + File.separator + (String)ecnums.get(i) + ".confs", new String[0]));
	        for (int j = 0; j < ppreds.size(); j++)
	        {
	          Double comb = Double.valueOf(Double.parseDouble(ppreds.get(j)));
	          combined.add(String.valueOf(comb));
	        }
	      }
	      else if (method.equals("weighted"))
	      {
	    	  List<String> spreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + "spmap" + File.separator + ecnums.get(i) + ".confs", new String[0]));
	    	  List<String> bpreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + "blast" + File.separator + ecnums.get(i) + ".confs", new String[0]));
	    	  List<String> ppreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + "pepstats" + File.separator + ecnums.get(i) + ".confs", new String[0]));
	        List<String> weights = Files.readAllLines(Paths.get(ROOTPATH + File.separator + ecnums.get(i) + File.separator + "weights.txt", new String[0]));
	      
	        for (int j = 0; j < spreds.size(); j++)
	        {
	          Double comb = Double.valueOf(Double.parseDouble(spreds.get(j)) * Double.parseDouble(weights.get(0)) + Double.parseDouble(bpreds.get(j)) * Double.parseDouble(weights.get(1)) + 
	            Double.parseDouble(ppreds.get(j)) * Double.parseDouble(weights.get(2)));
	          combined.add(String.valueOf(comb));
	        }
	        for (int j = 0; j < idlist.size(); j++)
	        {
	          BufferedWriter final_file1 = new BufferedWriter(new FileWriter(testDir + File.separator + (String)idlist.get(j) + ".preds", true));
	          final_file1.write(ecnums.get(i) + "\t" + spreds.get(j) + "\t" + bpreds.get(j) + "\t" + ppreds.get(j) + "\t" + combined.get(j) + "\n");
	          final_file1.close();
	        }
	      }
	      DecimalFormat df = new DecimalFormat();
	      df.setMaximumFractionDigits(2);
	      BufferedWriter final_file = new BufferedWriter(new FileWriter(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + (String)ecnums.get(i) + "_preds.txt", false));
	      for (int j = 0; j < idlist.size(); j++) {
	        final_file.write(df.format(Double.parseDouble(combined.get(j))) + "\n");
	      }
	      final_file.close();
	    }
		if(ecnums.contains("1.-.-.-")){
			Vector<Vector<String>> allPreds = new Vector<>();
			for (int i =0; i< ecnums.size(); i++) {  
				String workdir = ROOTPATH+File.separator+ecnums.get(i);
				Vector<String> pred = new Vector<>();
				br = new BufferedReader(new FileReader(tempDir + File.separator + "testResult" + File.separator + time + File.separator + (String)ecnums.get(i) + File.separator + (String)ecnums.get(i) + "_preds.txt"));					while((line = br.readLine())!=null){
						pred.add(line);
					}
				allPreds.add( pred);

			}
			
			for(int i=0 ; i < idlist.size(); i++){
				NHpredictions = new HashMap<>();
				Vector<Vector<String>> predswithScore = new Vector<>();
				Vector<String> preds = new Vector<>();
				double maxPred =0.0;
				String mainClass = null ;
				for(int j = 0 ; j < allPreds.size(); j++){
					if(Double.parseDouble(allPreds.get(j).get(i)) >= thresholds.get(String.valueOf((j+1)) + ".-.-.-")) {
						NHpredictions.put(String.valueOf((j+1)) + ".-.-.-", Double.valueOf(allPreds.get(j).get(i)));
					}						
					predictions.put(idlist.get(i) , NHpredictions);
					if(Double.parseDouble(allPreds.get(j).get(i))> maxPred){
						maxPred = Double.parseDouble(allPreds.get(j).get(i));
						mainClass = String.valueOf((j+1)) + ".-.-.-";
					}	
				}
				if(maxPred<0.4){
					predswithScore.add(preds);
					NHpredictions.put("non", maxPred);
					predictions.put(idlist.get(i) , NHpredictions);
				}
				else if (maxPred >= thresholds.get(mainClass)){
					preds.add(mainClass);
					preds.add(String.valueOf(maxPred));
					predswithScore.add(preds);
					predictions.put(idlist.get(i) , NHpredictions);
				}
				else{
					NHpredictions.put("nop", Double.valueOf("0"));
					predictions.put(idlist.get(i) , NHpredictions);
				}
					
			}
		}
		
		else{
			double maxPred = 0.0;
			String predClass ="";
			Vector<String> preds = new Vector<>();
			for (int i = 0; i< ecnums.size(); i++) {  
				List<String> pred = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + (String)ecnums.get(i) + File.separator + (String)ecnums.get(i) + "_preds.txt", new String[0]));
				predClass = ecnums.get(i);
				if(Double.parseDouble(pred.get(0)) >=  thresholds.get(predClass)){
					predictions.get(idlist.get(0)).put(predClass, Double.valueOf(pred.get(0)));
				}
			}
		}
		return predictions;
		
	}

}
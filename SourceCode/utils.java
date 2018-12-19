import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.Vector;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import org.w3c.dom.views.DocumentView;

public class utils {

	public static HashMap<String, Vector<Double>> calculateMeasures(double threshold, int alpha, Vector<Double> ppreds, Vector<Double> npreds){
		double tp = 0.0;
		double fp = 0.0; 
		double tn = 0.0;
		double fn = 0.0;
		HashMap<String, Vector<Double>> result = new HashMap<>();
		Collections.sort(ppreds);
		Collections.reverse(ppreds);	
		Collections.sort(npreds);
		Collections.reverse(npreds);
		for(int i=0; i< ppreds.size(); i++){
			if(ppreds.get(i)>=threshold)
				tp+=1.0;
			else
				fn+=1.0;			
		}
		for(int i=0; i< npreds.size(); i++){
			if(npreds.get(i)<threshold)
				tn+=1.0;
			else
				fp+=1.0;
			
		}
		result.put("ppreds", ppreds);
		result.put("npreds", npreds);
		Vector<Double> v = new Vector<>();
		v.add(tp);
		result.put("tp", v);
		v= new Vector<>();
		v.add(fp);
		result.put("fp", v);
		v= new Vector<>();
		v.add(tn);
		result.put("tn", v);
		v= new Vector<>();
		v.add(fn);
		result.put("fn", v);
		v= new Vector<>();
		double recall = tp*1.0/(tp+fn);
		v.add(recall);
		result.put("recall", v);
		v= new Vector<>();
		double precision = tp*1.0/(tp+fp);
		v.add(precision);
		result.put("precision", v);
		v= new Vector<>();
		v.add(tn*1.0/(tn+fp));
		result.put("specificity", v);
		v= new Vector<>();
		v.add(tp*1.0/(tp+fn));
		result.put("sensitivity", v);
		v= new Vector<>();
		v.add((1+alpha)*(precision * recall)/(alpha*precision+recall));
		result.put("fone", v);
		Double auc = roc_score(ppreds, npreds);
		v= new Vector<>();
		v.add(auc);
		result.put("auc", v);
	      return result;
	}

	public static double roc_score_aux(HashMap<Integer, String>  vals){
		double tp = 0.0;
		double fp = 0.0;
		double roc = 0.0;
		
		Set set = vals.entrySet();
	      Iterator it = set.iterator();
	      while(it.hasNext()) {
	         Map.Entry me = (Map.Entry)it.next();
	         if(me.getValue().toString().substring(0, 1).equals("P"))
	        	 tp+=1;
	         else{
	        	 fp+=1;
	        	 roc+=tp;
	         }

	      }
		
		 if (tp == 0.0)
				 return 0;
		 if (fp == 0.0)
				 return 1 ;
		return (roc * 1.0) / (1.0 * tp * fp );
		
	}

	public static double roc_score(List<Double> ppreds, List<Double> npreds){
		HashMap<Integer, String>  vals = new HashMap<>();
		int count =0;
		for(int i=0; i<ppreds.size(); i++){
			
			String x = "P";
			x += ppreds.get(i);
			vals.put(count, x);
			count++;
		}
		for(int i=0; i<npreds.size(); i++){

			String x = "N";
			x += npreds.get(i);
			vals.put(count, x);
			count++;
		}
		System.out.println(vals);
		vals = sortByComparator(vals, true);
		return roc_score_aux(vals);
		
	}

	
	private static HashMap<Integer, String> sortByComparator(HashMap<Integer, String> unsortMap, final boolean order)
	  {
	      List<Entry<Integer, String>> list = new LinkedList<Entry<Integer, String>>(unsortMap.entrySet());

	      Collections.sort(list, new Comparator<Entry<Integer, String>>()
	      {
	          public int compare(Entry<Integer, String> o1,
	                  Entry<Integer, String> o2)
	          {
	              if (order)
	              {
	            	  final double time1 = Double.parseDouble(o1.getValue().substring(1, o1.getValue().length()));
	  	            final double time2 = Double.parseDouble(o2.getValue().substring(1,o2.getValue().length()));
	  	            return Double.compare( time2,time1);
	              }
	              else
	              {
	            	  final double time1 = Double.parseDouble(o1.getValue().substring(1, o1.getValue().length()));
	    	            final double time2 = Double.parseDouble(o2.getValue().substring(1,o2.getValue().length()));
	    	            return Double.compare( time1,time2);
	              }
	          }
	      });

	      HashMap<Integer, String> sortedMap = new LinkedHashMap<Integer, String>();
	      for (Entry<Integer, String> entry : list)
	      {
	          sortedMap.put(entry.getKey(), entry.getValue());
	      }

	      return sortedMap;
	  }

	public static void calculateConfidence(String ppf, String npf, String pf, String cf) throws IOException{
		List<String> pos = Files.readAllLines(Paths.get(ppf));
		List<String> neg = Files.readAllLines(Paths.get(npf));
		List<String> preds = Files.readAllLines(Paths.get(pf));

		 int closer2pos;
		Vector<Double> confs = new Vector<>();
		for(int i = 0; i < preds.size(); i++){
			if(preds.get(i)== "0")
				confs.add(0.5);
			closer2pos = 0;
			for(int j = 0; j < pos.size(); j++){
				if (Double.parseDouble(pos.get(j)) <= Double.parseDouble(preds.get(i)))
					break;
				closer2pos += 1;
			}
			double posconf =  ((double)(pos.size() - closer2pos) / pos.size());
			
			closer2pos = 0;
			for(int j = 0; j < neg.size(); j++){
				if (Double.parseDouble(neg.get(j)) < Double.parseDouble(preds.get(i)))
					break;
				closer2pos += 1;
			}
			double negconf =  ((double)closer2pos/ neg.size());
			
			if (posconf + negconf == 0)
					confs.add(0.5);
			else
					confs.add(posconf / (posconf + negconf));
		}
	
		PrintWriter final_file = new PrintWriter(cf, "UTF-8");
		for(int i=0; i<confs.size(); i++){
			final_file.println(confs.get(i));
		}
		final_file.close();
		
	}
	
	public static void main(String[] args) throws IOException {
	
		
	}
	
	

}
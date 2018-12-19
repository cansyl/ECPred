import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringTokenizer;
import java.util.TreeSet;
import java.util.Vector;

import javax.swing.plaf.SliderUI;

public class Blast {

	public static HashMap<String, List<List<String>>> parseTabBlast(String blastFile) throws IOException{
		
		int id1f = 0;
		int id2f = 1;
		int bitsf = 11;
		HashMap<String, List<List<String>>> result = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(blastFile));	
		String line = br.readLine();
		StringTokenizer st1;
		String queryID ;
		String[] parts = line.split("\t");
		if(parts[id1f].contains("|")){
			 st1 = new StringTokenizer(parts[id1f], "|");//split

			queryID= st1.nextToken();
			queryID = st1.nextToken();
		}
		else{
			queryID = parts[id1f];
		}
		
		StringTokenizer st2 = new StringTokenizer(parts[id2f], "|");//split
		String hitID = st2.nextToken();
		hitID = st2.nextToken();
		String prevID = queryID;
		List<List<String>> allscore = new ArrayList<List<String>>();
		ArrayList<String> score = new ArrayList<>();
		score.add(hitID);
		score.add(parts[bitsf]);
		allscore.add(score);
		result.put(queryID, allscore);
		while ((line = br.readLine())!= null){
			 parts = line.split("\t");		
			queryID = parts[id1f];
			st2 = new StringTokenizer(parts[id2f], "|");//split
			hitID = st2.nextToken();
			hitID = st2.nextToken();
			
			if (queryID.equals(prevID)==false){			
				prevID = queryID;
				score = new ArrayList<>();
				score.add(hitID);
				score.add(parts[bitsf]);
				allscore = new ArrayList<List<String>>();
				allscore.add(score);
				result.put(queryID, allscore);
			}
			else{
				
				score = new ArrayList<>();
				score.add(hitID);
				score.add(parts[bitsf]);
				result.get(queryID).add(score);					
			}
		}
	
		for (String key : result.keySet()) {		
		Collections.sort(result.get(key), new Comparator<List<String>>()
		{
		  public int compare(List<String> o1, List<String> o2)
		  {
			  final double score1 = Double.parseDouble(o1.get(1));
	            final double score2 = Double.parseDouble(o2.get(1));
	            return Double.compare( score2,score1);
		  }
		});
		}
		return result;
	}

	public static double blastknn(List<List<String>> list,  List<String> posIDset, List<String> negIDset, int knn)throws IOException{
			
		if (list == null) 
			return  0;
	
		double posSum = 0.0;
		double negSum = 0.0;
		double total = 0.0;
		int count = 0;
		for(int i=0; i<list.size(); i++){
			if (posIDset.contains(list.get(i).get(0))){
				posSum += Double.parseDouble(list.get(count).get(1));
				total += Double.parseDouble(list.get(count).get(1));
				count++;
			}
				
			else if (negIDset.contains(list.get(i).get(0))){

				negSum += Double.parseDouble(list.get(count).get(1));
				total += Double.parseDouble(list.get(count).get(1));
				count++;
			}
			if (count == knn)
				break ;	
		}
		
		if (total == 0)
		return 0;

		else
			return (posSum - negSum) / total;		
	}
		
	public static void main(String[] args) throws IOException {
		
		Vector<String> lst_positive_ids = new Vector<>();		
	}

}

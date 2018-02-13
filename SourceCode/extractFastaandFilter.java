import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

public class extractFastaandFilter {
	
	public static void createFasta( List<String> lst_training_ids,String fasta, String outf, String file) throws IOException{
		PrintWriter final_file = new PrintWriter(file+"/"+outf, "UTF-8");
		String prot_id = "";
		BufferedReader br = new BufferedReader(new FileReader(fasta));	
		HashMap<String, String> fastaArray = new HashMap<>();
		Vector<String> fastaWrite = new Vector<>();
		StringTokenizer st1;
		String line = br.readLine(), fastaString = "";
		int i=0;
		 while(line != null ){
			 fastaString = "";
			   line = line.replaceAll("/","");
			   line = line.trim().replaceAll(" +", "");
					   if(line.length()>80)
							prot_id = line.substring(1, 81);
					   else
						   prot_id  = line.substring(1, line.length());
					   
		 if(lst_training_ids.contains(prot_id)){ 
			 fastaString+=">"+prot_id+"\n";
			line = br.readLine();
		    while(line.startsWith(">")==false){//until > add sequence 
		    	 fastaString+=line;
				line = br.readLine();
				if(line==null)
					 break;
			   }
			fastaArray.put(prot_id,fastaString);
		 	}
		 else{
				line = br.readLine();
			 while(line.startsWith(">")==false){//until > add sequence 
					line = br.readLine();
					if(line==null)
						 break;
				   }
		 }
		}	
		 
		for( i=0; i<lst_training_ids.size();i++){
			fastaWrite.add(fastaArray.get(lst_training_ids.get(i)));
		}
		for( i=0; i<fastaWrite.size();i++){
			if(fastaWrite.get(i)==null)
				continue;
			final_file.println(fastaWrite.get(i));
		}
		 br.close();
		 final_file.close();			 
	}

}
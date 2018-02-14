/*ECPred: a tool for the prediction of enzymatic properties of protein sequences based on the EC Nomenclature
    Copyright (C) 2018 CanSyL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.ObjectInputStream.GetField;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.print.attribute.standard.Media;

public class ECPpred {

	
		
	public static void main(String[] args) throws IOException, InterruptedException {
		Date d1 = new Date();
		long time = d1.getTime();
		
		Vector<String> ecnums = new Vector<>() ;
		ecnums.add("1.-.-.-");ecnums.add("2.-.-.-");ecnums.add("3.-.-.-");
		ecnums.add("4.-.-.-");ecnums.add("5.-.-.-");ecnums.add("6.-.-.-");
		String ROOTPATH = "lib/EC";
		runEC runECobj = new runEC();
		if(args.length<1){
			System.out.println("Please give a fasta file as an input!");
			System.exit(0);
		}
		else if(args.length>1){
			System.out.println("Please give only a fasta file as an input!");
			System.exit(0);
		}	
		StringTokenizer st1;
		String fastaFile = args[0];
		String file_namepart1 = "",file_name;
		if(fastaFile.contains("/")){
			 st1 = new StringTokenizer(fastaFile, "/");//split
			 while(st1.hasMoreTokens())
				 file_namepart1 = st1.nextToken();
			 st1 = new StringTokenizer(file_namepart1, ".");//split
				file_name = st1.nextToken();
		}
		else{
			st1 = new StringTokenizer(fastaFile, ".");//split
			file_name = st1.nextToken();
		}
		
		String dateandtime = new SimpleDateFormat("yyyyMMdd-HHmmss").format(new Date());

		
		HashMap<String, String> protID = checkFasta(fastaFile);
		List<String> idlist = createFasta(fastaFile);

		if(idlist.size()>20){
			System.out.println("Input protein size should be less than or equal to 20!");
			System.exit(0);
		}
		HashMap<String, Vector<Vector<String>>> predictions = new HashMap<>();
		System.out.println("Main classes of input proteins are being predicted ...");
		extractFastaandFilter ext = new extractFastaandFilter();
		createFasta(idlist, fastaFile, "test.fasta", ROOTPATH+File.separator+ "testResult/" + time);
		String newfasta =  ROOTPATH+File.separator+ "testResult/" + time + File.separator + "test.fasta"; 
		predictions = runECobj.predictions(args, ROOTPATH, ecnums, time,predictions, idlist, newfasta);
	     for (Map.Entry<String, Vector<Vector<String>>> entry : predictions.entrySet()) {
	    	  if(entry.getValue().get(0).get(0).equals("non") || entry.getValue().get(0).get(0).equals("nop"))
	        		continue;
	    	  if(protID.get(entry.getKey()).length()>81)
	    		  System.out.println("Subclasses of "+protID.get(entry.getKey()).substring(1,81) + " are being predicted ...");
    		  else
    	    	  System.out.println("Subclasses of "+protID.get(entry.getKey()).substring(1,protID.get(entry.getKey()).length()) + " are being predicted ...");

	    	  for(int i = 1 ; i<4; i++){
		        		List<String> ecList = Files.readAllLines(Paths.get("lib/subclasses/"+ entry.getValue().get(i-1).get(0) + ".txt"));
		        		if(ecList.size()==0){
		        			Vector<String> preds = new Vector<>();
		        			preds.add("nop");
		    				preds.add("");
		        			predictions.get(idlist.get(0)).add(preds);
		        			break;
		        		}

		        		ecnums = new Vector<>() ;
		        		ecnums.addAll(ecList);
		        		idlist = new Vector<String>();
		        		idlist.add(entry.getKey());
		        		ext.createFasta(idlist, fastaFile, "test.fasta", ROOTPATH+File.separator+ "testResult/" + time);
		        		String indivfasta =  ROOTPATH+File.separator+ "testResult/" + time + File.separator + "test.fasta"; 
		        		predictions = runECobj.predictions(args,ROOTPATH, ecnums, time,predictions, idlist,indivfasta);
		        	
		        	 if(entry.getValue().get(i).get(0).equals("nop"))
		        		break;
		        	
		        }
	    	 }
			PrintWriter predFile = new PrintWriter("predictionResults_"+file_name+"_"+ dateandtime+".tsv", "UTF-8");
			predFile.println("Protein ID\tMain Class\tConfidence Score(max 1.0)\tSubfamily Class\tConfidence Score(max 1.0)\tSub-subfamily Class\tConfidence Score(max 1.0)\t"
					+ "Substrate Class\tConfidence Score(max 1.0)");
			idlist = createFasta(fastaFile);
		      for (int a = 0 ; a < idlist.size(); a++ ) {
		    	  for (Map.Entry<String, Vector<Vector<String>>> entry : predictions.entrySet()) {
		    		  if(idlist.get(a).equals(entry.getKey())==false)
		    			 continue;
		    		  if(protID.get(entry.getKey()).length()>81)
		    			  predFile.print(protID.get(entry.getKey()).substring(1,81));
		    		  else
		    			  predFile.print(protID.get(entry.getKey()).substring(1,protID.get(entry.getKey()).length()));
			    	  if(entry.getValue().get(0).get(0).equals("non")){
				    	  predFile.println("\tnon Enzyme\t\t\t\t\t\t\t");
				    	  continue;
			    	  }
			    	  if(entry.getValue().get(0).get(0).equals("nop")){
				    	  predFile.println("\tno Prediction\t\t\t\t\t\t\t");
				    	  continue;
			    	  }
			    	  for(int i=0; i<4; i++){
			    		  if(entry.getValue().get(i).get(0).equals("nop")){
			    			  for(int j=i; j<4; j++){
				    			  predFile.print("\t\t");
			    			  }
					    	  break;
			    		  }
				    	  predFile.print("\t"+entry.getValue().get(i).get(0)+"\t"+entry.getValue().get(i).get(1));
				    	  
			    	  }
			    	  predFile.println();
						
						}
		      }	    
	      predFile.close();
	     File deletefile = new File(ROOTPATH+File.separator+"testResult");
		 deleteDirectory(deletefile);
			   Date d2 = new Date();
				long diff = d2.getTime() - d1.getTime();
				long diffSeconds = diff / 1000 % 60;
				long diffMinutes = diff / (60 * 1000) % 60;
				long diffHours = diff / (60 * 60 * 1000) % 24;
			    long diffDays = diff / (24 * 60 * 60 * 1000);

				if(diffMinutes==0 && diffHours ==0 && diffDays==0)
					System.out.println("--- Proteins are predicted in "+diffSeconds + " seconds ---");
				else if(diffHours ==0 && diffDays==0)
					System.out.println("--- Proteins are predicted in "+diffMinutes + " minutes "+diffSeconds + " seconds ---");
				else if(diffDays==0)
					System.out.println("--- Proteins are predicted in "+diffHours + " hours "+diffMinutes + " minutes "+diffSeconds + " seconds ---");
				else
					System.out.print("--- Proteins are predicted in "+diffDays + " days "+diffHours + " hours "+diffMinutes + " minutes "+diffSeconds + " seconds ---");				

	}
	
	public static boolean deleteDirectory(File directory) {
	    if(directory.exists()){
	        File[] files = directory.listFiles();
	        if(null!=files){
	            for(int i=0; i<files.length; i++) {
	                if(files[i].isDirectory()) {
	                    deleteDirectory(files[i]);
	                }
	                else {
	                    files[i].delete();
	                }
	            }
	        }
	    }
	    return(directory.delete());
	}
	
	public static HashMap<String, String> checkFasta( String fasta) throws IOException{
		HashMap<String, String> protID_Full = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(fasta));	
		List<String> idlist = new ArrayList<String>();
		String orjinal_line,line = br.readLine(), fastaString = "";
		String prot_id;
		 while(line != null ){
			 fastaString = "";
			 if(line.startsWith(">")==false){
					System.out.println("Wrong input! Sequences should start with \">\" character.");
					System.exit(0);
				}
			 orjinal_line = line;
			 line = line.replaceAll("/","");
			   line = line.trim().replaceAll(" +", "");
			   if(line.length()>80)
					prot_id = line.substring(1, 81);
			   else
				   prot_id  = line.substring(1, line.length());
			 protID_Full.put(prot_id, orjinal_line);
			line = br.readLine();
		    while(line.startsWith(">")==false){//until > add sequence 
		    	 fastaString+=line;
				line = br.readLine();
				if(line==null)
					 break;
			   }
		    fastaString = fastaString.replaceAll("\\s+","");

		    Pattern p = Pattern.compile("[^A-Z ]", Pattern.CASE_INSENSITIVE);
		      Matcher m = p.matcher(fastaString);
		      boolean flag =false;
		      while (m.find()) {
		         System.out.println("Fasta sequence contains special character at position "  + m.start() + ": " + fastaString.charAt(m.start()) + " Your fasta: "+fastaString);
		         flag = true;
		      }
		      if(flag== true){
			      System.exit(0);
		      }
		}	
		 br.close();
		return protID_Full;
		
	}

	
	public static List<String> createFasta( String fasta) throws IOException{
		String prot_id = "";
		BufferedReader br = new BufferedReader(new FileReader(fasta));	
		List<String> idlist = new ArrayList<String>();
		String line = br.readLine(), fastaString = "";
		
		 while(line != null ){
				 fastaString = "";
			
				   line = line.replaceAll("/","");
				   line = line.trim().replaceAll(" +", "");
				   if(line.length()>80)
						prot_id = line.substring(1, 81);
				   else
					   prot_id  = line.substring(1, line.length());
				line = br.readLine();
			    while(line.startsWith(">")==false){//until > add sequence 
			    	 fastaString+=line;
					line = br.readLine();
					if(line==null)
						 break;
				   }
			    if(fastaString.length()>40)
					idlist.add(prot_id);
			 	
			}	
			 br.close();	
		 return idlist;
	}

	public static void createFasta( List<String> lst_training_ids,String fasta, String outf, String file) throws IOException{
		File workdir = new File(file);
		workdir.mkdirs();
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
			
					   fastaString+=">"+prot_id+"\n";
			 
			 
		 if(lst_training_ids.contains(prot_id)){ 
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
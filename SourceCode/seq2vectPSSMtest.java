import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.Vector;

public class seq2vectPSSMtest
{
  public static HashMap<String, Integer> blosum_dict1 = new HashMap();
  public static Vector<String> lst_positive_ids = new Vector();
  
  public static HashMap<String, String> getFastaOrg(String file)
    throws IOException
  {
    HashMap<String, String> fasta_dict = new HashMap();
    String prot_id = "";
    String prot_seq = "";
    
    BufferedReader br = new BufferedReader(new FileReader(file));
    String line = br.readLine();
    while (line != null)
    {
      if (line.startsWith(">"))
      {
        prot_seq = "";
        StringTokenizer st1 = new StringTokenizer(line, "|");
        line = line.replaceAll("/", "");
        line = line.trim().replaceAll(" +", "");
        if (line.length() > 80) {
          prot_id = line.substring(1, 81);
        } else {
          prot_id = line.substring(1, line.length());
        }
        line = br.readLine();
        while (!line.startsWith(">"))
        {
          prot_seq = prot_seq + line;
          line = br.readLine();
          if (line == null) {
            break;
          }
        }
      }
      fasta_dict.put(prot_id, prot_seq);
    }
    br.close();
    return fasta_dict;
  }
  
  public static Vector<String> extractSubsequences(String subseq, int subseq_len)
  {
    String[] subseqArray = subseq.split("");
    String lst_subeqs = "";
    Vector<String> v = new Vector();
    int count = 0;
    int total_subseqs = subseqArray.length - subseq_len + 1;
    while (count < total_subseqs)
    {
      lst_subeqs = subseqArray[count] + subseqArray[(count + 1)] + subseqArray[(count + 2)] + subseqArray[(count + 3)] + subseqArray[(count + 4)];
      v.add(lst_subeqs);
      count++;
    }
    return v;
  }
  
  public static Vector<Double> calculateSignature(String string, Vector<Double> vector, HashMap<Integer, HashMap<Integer, ArrayList<String>>> pssm, int subseqlen)
  {
    double sum = 0.0D;
    String aa_letters = "A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V";
    String[] lst_aa_letters = aa_letters.split(",");
    vector.add(Double.valueOf(0.0D));
    if ((string.contains("X")) || (string.contains("B")) || (string.contains("Z")) || (string.contains("U")) || (string.contains("O")) || (string.contains("J")))
    {
      for (int p = 1; p <= pssm.size(); p++) {
        vector.add(Double.valueOf(sum));
      }
      return vector;
    }
    for (int p = 1; p <= pssm.size(); p++)
    {
      sum = 0.0D;
      for (int i = 1; i < subseqlen + 1; i++)
      {
        double x = Double.parseDouble(((pssm.get(Integer.valueOf(p))).get(Integer.valueOf(i))).get(Arrays.asList(lst_aa_letters).indexOf(string.substring(i - 1, i)) * 2 + 1));
        sum += x;
      }
      vector.add(Double.valueOf(sum));
    }
    return vector;
  }
  
  public static float calculateVectors(int signifThreshold, int subseqlen, String ECNumber, List<String> targetList, String filename, long time, String ROOTPATH, String tempDir)
    throws IOException
  {
    int s = 0;
    Vector<Vector<Double>> result = new Vector();
    Vector<Double> vector = new Vector();
    Vector<Double> tmpvect = new Vector();
    result.add(vector);
    int numConverted = 0;
    HashMap<String, String> fasta_dict = new HashMap();
    fasta_dict = getFastaOrg(filename);
    
    List<String> lines = Files.readAllLines(Paths.get(ROOTPATH + File.separator + ECNumber + "/spmap/profile.txt", new String[0]));
    int number_of_cluster = lines.size() / subseqlen;
    HashMap<Integer, HashMap<Integer, ArrayList<String>>> pssm = new HashMap();
    ArrayList<String> aa_list = new ArrayList();
    HashMap<Integer, ArrayList<String>> cluster_ps_aa_dict2 = new HashMap();
    int count = 0;
    for (int ind2 = 1; ind2 <= number_of_cluster; ind2++)
    {
      cluster_ps_aa_dict2 = new HashMap();
      for (int k = 1; k < subseqlen + 1; k++)
      {
        aa_list = new ArrayList();
        StringTokenizer st1 = new StringTokenizer(lines.get(count++), "\t");
        String list = null;
        while (st1.hasMoreElements()) {
          list = st1.nextToken();
        }
        aa_list = new ArrayList(Arrays.asList(list.split("\\s*,\\s*")));
        cluster_ps_aa_dict2.put(Integer.valueOf(k), aa_list);
      }
      pssm.put(Integer.valueOf(ind2), cluster_ps_aa_dict2);
    }
    for (s = 0; s < targetList.size(); s++)
    {
      Vector<String> lst_subseqs = extractSubsequences(fasta_dict.get(targetList.get(s)), subseqlen);
      
      vector = new Vector();
      if (lst_subseqs.size() >= 6)
      {
        vector = calculateSignature(lst_subseqs.get(0), vector, pssm, subseqlen);
        for (int p = 1; p < lst_subseqs.size(); p++)
        {
          tmpvect = new Vector();
          tmpvect = calculateSignature(lst_subseqs.get(p), tmpvect, pssm, subseqlen);
          for (int i = 1; i <= pssm.size(); i++) {
            if ((vector.get(i)).doubleValue() < (tmpvect.get(i)).doubleValue()) {
              vector.set(i, tmpvect.get(i));
            }
          }
        }
        for (int i = 1; i <= pssm.size(); i++) {
          if ((vector.get(i)).doubleValue() < signifThreshold) {
            vector.set(i, Double.valueOf(0.0D));
          } else {
            vector.set(i, Double.valueOf(Math.exp((vector.get(i)).doubleValue() / subseqlen)));
          }
        }
        result.add(vector);
        numConverted++;
      }
    }
    writeVectors(ECNumber, result, numConverted, time, ROOTPATH, tempDir);
    return 0.0F;
  }
  
  public static void writeVectors(String ECNumber, Vector<Vector<Double>> result, int numvect, long time, String ROOTPATH, String tempDir)
    throws IOException
  {
    PrintWriter final_file = new PrintWriter(tempDir + "/testResult/" + time + File.separator + ECNumber + "/spmap/test.vec", "UTF-8");
    for (int i = 1; i <= numvect; i++)
    {
      String f_l = "1 ";
      for (int key = 1; key < (result.get(1)).size(); key++) {
        f_l = f_l + key + ":" + (result.get(i)).get(key) + " ";
      }
      final_file.println(f_l);
    }
    final_file.close();
  }
  
  public static HashMap<String, Integer> readBLOSUM62Matrix()
    throws IOException
  {
    String aa_letters = "A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V";
    String[] lst_aa_letters = aa_letters.split(",");
    HashMap<String, Integer> blosum_dict = new HashMap();
    BufferedReader br = new BufferedReader(new FileReader("blo62.csv"));
    String[] blo62_mat = null;
    
    int i = 0;
    String line;
    while ((line = br.readLine()) != null)
    {

      blo62_mat = line.split(",");
      for (int j = 0; j < blo62_mat.length; j++) {
        blosum_dict.put(lst_aa_letters[i] + "," + lst_aa_letters[j], Integer.valueOf(Integer.parseInt(blo62_mat[j])));
      }
      i++;
    }
    return blosum_dict;
  }
}
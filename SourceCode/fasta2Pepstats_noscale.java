import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

public class fasta2Pepstats_noscale
{
  public static Vector<String> parse_single(Vector<String> single)
    throws IOException
  {
    String irrel = "B,J,O,U,X,Z";
    List<String> lst_irrel = new ArrayList(Arrays.asList(irrel.split(",")));
    Vector<String> lines = new Vector();
    for (int i = 0; i < single.size(); i++) {
      if ((i != 0) && (i != 1) && (i != 8) && (i != 9) && (i != 36) && (i != 37) && (i != 47))
      {
        for (i = 2; i < 8; i++)
        {
          StringTokenizer st1 = new StringTokenizer(single.get(i), "/\t+/");
          StringTokenizer st2 = new StringTokenizer(single.get(i), "/\t+/");
          if (i == 5)
          {
            String s1 = (String)st1.nextElement();
            String[] a1 = s1.split("\\s+");
            lines.add(a1[(a1.length - 3)]);
          }
          else if (i == 6)
          {
        	String s1 = (String)st1.nextElement();
            s1 = (String)st1.nextElement();
            String[] a1 = s1.split("\\s+");
            lines.add(a1[(a1.length - 3)]);
          }
          else
          {
            while (st1.hasMoreElements())
            {
              String s1;
              String[] a1;
              String s = (String)st1.nextElement();
              String[] a = s.split("\\s+");
              lines.add(a[(a.length - 1)]);
            }
          }
        }
        for (i = 10; i < 36; i++)
        {
          String[] a = (single.get(i)).split("\\s+");
          if (!lst_irrel.contains(a[0])) {
            lines.add(a[(a.length - 1)]);
          }
        }
        for (i = 38; i < 47; i++)
        {
          String[] a = (single.get(i)).split("\\s+");
          lines.add(a[(a.length - 1)]);
        }
      }
    }
    return lines;
  }
  
  public static void parse_pepstats(String infile, String ECNumber, String type, long time, String ROOTPATH, String tempDir)
    throws IOException, InterruptedException
  {
    BufferedReader br = new BufferedReader(new FileReader(infile));
    
    int count = 0;int vecCount = 0;
    Vector<String> single = new Vector();
    Vector<Vector<String>> vects = new Vector();
    String line;
    while ((line = br.readLine()) != null)
    {

      count++;
      single.add(line);
      if (line.contains("None"))
      {
        while (count % 49 != 0)
        {
          line = br.readLine();
          count++;
        }
        count = 0;
        single = new Vector();
      }
      else if (count % 48 == 0)
      {
        vects.add(vecCount, parse_single(single));
        vecCount++;
        single = new Vector();
        count = 0;
      }
    }
    print_vector(ECNumber, vects, type, time, ROOTPATH, tempDir);
  }
  
  public static Vector<String> allFasta(String ECNumber, String type)
    throws IOException
  {
    HashMap<String, String> fasta_dict_test = new HashMap();
    Vector<String> all_ids = new Vector();
    for (String key : fasta_dict_test.keySet()) {
      all_ids.add(key);
    }
    return all_ids;
  }
  
  public static void print_vector(String ECNumber, Vector<Vector<String>> vects, String type, long time, String ROOTPATH, String tempDir)
    throws IOException, InterruptedException
  {
    int count = 0;
    PrintWriter final_file = new PrintWriter(tempDir + "/testResult/" + time + File.separator + ECNumber + "/pepstats/test.svm", "UTF-8");
    for (int i = 0; i < vects.size(); i++)
    {
      String f_l = "";
      f_l = "1 ";
      for (int j = 0; j < (vects.get(i)).size(); j++) {
        f_l = f_l + (j + 1) + ":" + (String)(vects.get(i)).get(j) + " ";
      }
      final_file.write(f_l + "\n");
    }
    final_file.close();
  }
}

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Date;
import java.util.List;
import java.util.Vector;

public class predictBatchPEPSTATS
{
  public static void main(String[] args, Vector<String> ecnums, long time, String ROOTPATH, String fastaFile, String tempDir)
    throws IOException, InterruptedException
  {
    String predictProg = ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/svmlight/svm_classify";
    String method = "pepstats";
    for (int i = 0; i < ecnums.size(); i++)
    {
      Date d1 = new Date();
      
      String path = ROOTPATH + File.separator + (String)ecnums.get(i) + File.separator + method;
      String modelfile = path + File.separator + "model.svm";
      String rangefile = path + File.separator + "rangefile";
      String testpath = tempDir + File.separator + "testResult" + File.separator + time + File.separator + (String)ecnums.get(i) + File.separator + method;
      
      String batchSVM = testpath + "/test.svm";
      String batchVect = testpath + "/test.vec";
      File workdir = new File(testpath);
      workdir.mkdirs();
      
      fasta2Pepstats_noscale fas = new fasta2Pepstats_noscale();
      
      String predFile = testpath + File.separator + (String)ecnums.get(i) + ".preds";
      String confFile = testpath + File.separator + (String)ecnums.get(i) + ".confs";
      
      String posPredFile = path + File.separator + "ppreds.txt";
      String negPredFile = path + File.separator + "npreds.txt";
      String[] cmdArray = new String[7];
      
      cmdArray[0] = (ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/EMBOSS-6.5.7/emboss/pepstats");
      
      cmdArray[1] = "-sequence";
      
      cmdArray[2] = fastaFile;
      
      cmdArray[3] = "-outfile";
      
      cmdArray[4] = (workdir + "/temp.out");
      
      cmdArray[5] = "-warning";
      
      cmdArray[6] = "FALSE";
      
      String cmd = ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/EMBOSS-6.5.7/emboss/pepstats -sequence " + fastaFile + " -outfile " + workdir + "/temp.out -warning FALSE";
      ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
      
      Process process = pb.start();
      try
      {
        process.waitFor();
        int exitVal = process.exitValue();
      }
      catch (InterruptedException e)
      {
        System.out.print("pepstats is not working!");
        e.printStackTrace();
      }
      fasta2Pepstats_noscale.parse_pepstats(cmdArray[4], ecnums.get(i), "1", time, ROOTPATH, tempDir);
      File deletefile = new File(cmdArray[4]);
      deletefile.delete();
      cmdArray = new String[1];
      PrintWriter scalefile = new PrintWriter(tempDir + File.separator + "testResult" + File.separator + time + File.separator + time + ".sh", "UTF-8");
      scalefile.println("#!/bin/bash");
      scalefile.println(ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/libsvm-3.16/svm-scale -r " + ROOTPATH + "/" + ecnums.get(i) + "/pepstats/rangefile " + batchSVM + " > " + batchVect + " 2> /dev/null");
      scalefile.close();
      cmd = "chmod +x " + tempDir + File.separator + "testResult" + File.separator + time + File.separator + time + ".sh";
      pb = new ProcessBuilder(cmd.split(" "));
      process = pb.start();
      try
      {
        process.waitFor();
        int j = process.exitValue();
      }
      catch (InterruptedException e)
      {
        System.out.print("svm-scale is not working!");
        e.printStackTrace();
      }
      cmd = tempDir + File.separator + "testResult" + File.separator + time + File.separator + time + ".sh";
      pb = new ProcessBuilder(new String[] { cmd });
      
      process = pb.start();
      try
      {
        process.waitFor();
        int k = process.exitValue();
      }
      catch (InterruptedException e)
      {
        System.out.print("Svm_classify is not working!");
        e.printStackTrace();
      }
      cmd = predictProg + " " + batchVect + " " + modelfile + " " + predFile;
      pb = new ProcessBuilder(cmd.split(" "));
      
      process = pb.start();
      try
      {
        process.waitFor();
        int m = process.exitValue();
      }
      catch (InterruptedException e)
      {
        System.out.print("Svm_classify is not working!");
        e.printStackTrace();
      }
      utils u = new utils();
      utils.calculateConfidence(posPredFile, negPredFile, predFile, confFile);
    }
  }
}
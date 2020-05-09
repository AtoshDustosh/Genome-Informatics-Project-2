package application;

import pwmcalculator.PWMmanager;

public class TFBSsearcher {
  public static final String PROMOTER = "CCCGGGGTCCAGTAGGGGGCACTCAC";
  public static final String TEST_PROMOTER = "GGGTCAGCATGGCCA";

  public static void main(String[] args) {
    PWMmanager pwmManager = new PWMmanager(PWMmanager.PFM_FILEPATH);
    String promoter = TFBSsearcher.PROMOTER;
    String tfbsString = "";

    double maxRelativeScore = 0;
    int startingPos = 0;

    double maximumScore = pwmManager.calcMaximumScore();
    double minimumScore = pwmManager.calcMinimumScore();

    System.out.println("PFM matrix\n" + pwmManager.getPFMmatrixInfo());
    System.out.println("PWM matrix\n" + pwmManager.getPWMmatrixInfo());

    System.out.println("maximum score: " + maximumScore);
    System.out.println("minimum score: " + minimumScore);

    for (int i = 0; i + pwmManager.getColumnCount() <= promoter.length(); i++) {
      String sequence = promoter.substring(i, i + pwmManager.getColumnCount());
      double absoluteScore = pwmManager.calcAbsoluteScore(sequence);
      double relativeScore = (absoluteScore - minimumScore)
          / (maximumScore - minimumScore);
      System.out.println("sequence: " + sequence);
      System.out.println("\trelative score: " + relativeScore);
      if (relativeScore > maxRelativeScore) {
        maxRelativeScore = relativeScore;
        startingPos = i;
      }
    }

    tfbsString = promoter.substring(startingPos,
        startingPos + pwmManager.getColumnCount());
    System.out.println("TFBS: " + tfbsString);

  }
}

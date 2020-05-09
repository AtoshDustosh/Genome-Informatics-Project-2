package pwmcalculator;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PWMmanager {
  public static final String PFM_FILEPATH = "src/testdata/MA0139.1.jaspar";
  public static final String PFM_TEST_FILEPATH = "src/testdata/test.jaspar";

  public static final Double BACKGROUND_FREQUENCY_A = 0.3;
  public static final Double BACKGROUND_FREQUENCY_C = 0.2;
  public static final Double BACKGROUND_FREQUENCY_G = 0.2;
  public static final Double BACKGROUND_FREQUENCY_T = 0.3;

  @SuppressWarnings("unused")
  private String pfmFileHeader = "";
  private Map<String, List<Integer>> pfmMatrix = new HashMap<>();
  private Map<String, List<Double>> pwmMatrix = new HashMap<>();

  private int columnCount = 0;
  private Map<String, Double> backgroundPMap = new HashMap<>();

  public PWMmanager(String filePath) {
    this.initialize();
    this.loadPFM(filePath);
    this.calculatePWMmatrix();
  }

  public static void main(String[] args) {
    PWMmanager pwmManager = new PWMmanager(PWMmanager.PFM_TEST_FILEPATH);
    System.out.println(pwmManager.getPFMmatrixInfo());
    System.out.println(pwmManager.getPWMmatrixInfo());
  }

  /**
   * Calculate the absolute score of a nucleotide sequence. This method
   * requires the sequence must match the column count of this PWM.
   * 
   * @param sequence nucleotide sequence
   * @return absolute score according to the PWM matrix
   */
  public double calcAbsoluteScore(String sequence) {
    double score = 0;

    // check the sequence length
    if (sequence.length() != this.columnCount) {
      System.out.println("the sequence is not fit for this PWM");
      System.exit(-1);
    }
    // check the sequence content
    String pattern = "[ACGT]+";
    Pattern regex = Pattern.compile(pattern);
    Matcher matcher = regex.matcher(sequence);
    if (matcher.matches() == false) {
      System.out.println("the sequence has invalid characters. ");
      System.exit(-1);
    }
    // calculation
    for (int i = 0; i < sequence.length(); i++) {
      String key = "" + sequence.charAt(i);
      score = score + this.valueOfPWM(key, i);
    }

    return score;
  }

  public double calcMaximumScore() {
    double score = 0;
    for (int i = 0; i < this.columnCount; i++) {
      score = score + this.maxValueOfPWMColumn(i);
    }
    return score;
  }

  public double calcMinimumScore() {
    double score = 0;
    for (int i = 0; i < this.columnCount; i++) {
      score = score + this.minValueOfPWMColumn(i);
    }
    return score;
  }

  public int getColumnCount() {
    return this.columnCount;
  }

  public String getPFMmatrixInfo() {
    String str = "";
    for (String key : this.pfmMatrix.keySet()) {
      List<Integer> row = this.pfmMatrix.get(key);
      str = str + key + " - " + row.toString() + "\n";
    }
    return str;
  }

  public String getPWMmatrixInfo() {
    String str = "";
    for (String key : this.pwmMatrix.keySet()) {
      List<Double> row = this.pwmMatrix.get(key);
      str = str + key + " - " + row.toString() + "\n";
    }
    return str;
  }

  private void initialize() {
    this.backgroundPMap.put("A", PWMmanager.BACKGROUND_FREQUENCY_A);
    this.backgroundPMap.put("C", PWMmanager.BACKGROUND_FREQUENCY_C);
    this.backgroundPMap.put("G", PWMmanager.BACKGROUND_FREQUENCY_G);
    this.backgroundPMap.put("T", PWMmanager.BACKGROUND_FREQUENCY_T);
  }

  private void loadPFM(String filePath) {
    try {
      Scanner scanner = new Scanner(new FileInputStream(filePath));

      this.pfmFileHeader = scanner.nextLine();

      String[] row1 = scanner.nextLine().split("(\\s)+");
      String[] row2 = scanner.nextLine().split("(\\s)+");
      String[] row3 = scanner.nextLine().split("(\\s)+");
      String[] row4 = scanner.nextLine().split("(\\s)+");

      if (!(row1.length == row2.length && row2.length == row3.length
          && row3.length == row4.length)) {
        System.out.println("error - input PFM matrix invalid. ");
        System.exit(-1);
      }

      List<String[]> rowsList = new ArrayList<>(
          Arrays.asList(row1, row2, row3, row4));
      for (int i = 0; i < rowsList.size(); i++) {
        String key = rowsList.get(i)[0];
        List<Integer> frequencyRow = new ArrayList<>();
        for (int j = 2; j < rowsList.get(i).length - 1; j++) {
          // the frequency occurs at [2] of the line
          frequencyRow.add(Integer.parseInt(rowsList.get(i)[j]));
        }
        this.pfmMatrix.put(key, frequencyRow);
        this.columnCount = frequencyRow.size();
      }

      scanner.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
  }

  private void calculatePWMmatrix() {
    final int columnCount = this.columnCount;
    for (String key : this.pfmMatrix.keySet()) { // for each row - A, C, G, T
      List<Double> frequencyRow = new ArrayList<>();
      for (int i = 0; i < columnCount; i++) { // for each column
        int N = 0;
        for (String tempKey : this.pfmMatrix.keySet()) { // calculate the total count of bp
          N = N + this.pfmMatrix.get(tempKey).get(i);
        }
        double rawCount = this.pfmMatrix.get(key).get(i);
        double result = rawCount + Math.sqrt(N) / 4;
        result = result / (N + Math.sqrt(N));
        result = result / this.backgroundPMap.get(key);
        result = Math.log(result) / Math.log(2);
        frequencyRow.add(result);
      }
      this.pwmMatrix.put(key, frequencyRow);
    }
  }

  private double maxValueOfPWMColumn(int i) {
    double maxValue = 0;
    for (String key : this.pwmMatrix.keySet()) {
      double value = this.pwmMatrix.get(key).get(i);
      if (value > maxValue) {
        maxValue = value;
      }
    }
    return maxValue;
  }

  private double minValueOfPWMColumn(int i) {
    double minValue = 0;
    for (String key : this.pwmMatrix.keySet()) {
      double value = this.pwmMatrix.get(key).get(i);
      if (value < minValue) {
        minValue = value;
      }
    }
    return minValue;
  }

  private double valueOfPWM(String key, int index) {
    return this.pwmMatrix.get(key).get(index);
  }

}

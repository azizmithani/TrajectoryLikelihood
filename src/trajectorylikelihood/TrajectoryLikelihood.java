/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package trajectorylikelihood;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;

/**
 *
 * @author mithani
 */
public class TrajectoryLikelihood {

    static BigDecimal[] fact;
    static BigDecimal epsilon;
    static BigDecimal natural_e;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
//        double[] jumpRates = {0.5000, 0.3000, 0.3000, 0.3000, 0.5000};
//        double[] exitRates = {2.4000, 2.2000, 2.4000, 2.6000, 2.8000, 2.6000};
//        double[] jumpRates = {0.0500, 0.0300, 0.0300, 0.0500, 0.0500, 0.0300, 0.0500, 0.0300, 0.0300, 0.0500, 0.0500, 0.0300, 0.0500, 0.0300, 0.0300, 0.0300, 0.0500, 0.0500};
//        double[] exitRates = {0.4200, 0.4000, 0.4200, 0.4400, 0.4200, 0.4000, 0.4200, 0.4000, 0.4200, 0.4400, 0.4200, 0.4000, 0.4200, 0.4000, 0.4200, 0.4400, 0.4600, 0.4400, 0.4200};


        int prec = 301;

        int n = (2 * prec) / 3;
        fact = new BigDecimal[n];
        fact[0] = new BigDecimal(1, MathContext.DECIMAL128);
        for (int i = 1; i < n; i++) {
            fact[i] = fact[i - 1].multiply(new BigDecimal(i), MathContext.DECIMAL128);
        }

        epsilon = new BigDecimal("1");
        int bits = 1000; // precision in bits = about 3.32 precision in digits
        for (int i = 0; i < bits; i++) {
            epsilon = epsilon.multiply(new BigDecimal("0.5"));
        }
        epsilon = epsilon.setScale(prec, BigDecimal.ROUND_DOWN);

        natural_e = naturalE(prec);


//        BigDecimal bd = new BigDecimal(2395376998.85156);
        System.out.println("Start: " + new Date());
        for (int i = 0; i < 10; i++) {
            //exp(bd);
        }
        System.out.println("Finish: " + new Date());


//        BigDecimal[][] rm = new BigDecimal[][]{{new BigDecimal(-3), new BigDecimal(1), new BigDecimal(2)},
//            {new BigDecimal(3), new BigDecimal(-7), new BigDecimal(4)},
//            {new BigDecimal(5), new BigDecimal(6), new BigDecimal(-11)}};
        double[][] rm = new double[][]{{-3, 1, 2},
            {3, -7, 4},
            {5, 6, -11}};

        int[] path = new int[]{1, 2, 0, 2};
        ArrayList lstPaths = new ArrayList();
        lstPaths.add(path);
        calculateProbability(lstPaths,rm);

        lstPaths = new ArrayList();
        lstPaths.add(new int[]{1,0,1,2});
        lstPaths.add(new int[]{1,2,0,2});
        lstPaths.add(new int[]{1,2,1,2});
//        System.out.print("Probability: ");
//        System.out.println(calculateProbability(lstPaths,rm));
//        double[] jumpRates = {4.0, 5.0, 2.0};
//        double[] exitRates = {7.0, 11.0, 3.0, 11.0};
        //0.5,0.3
//        double[] jumpRates = {0.5000, 0.3000, 0.3000, 0.5000, 0.5000, 0.3000, 0.3000, 0.5000, 0.5000, 0.3000, 0.5000, 0.3000, 0.3000, 0.5000, 0.5000, 0.5000, 0.3000, 0.3000, 0.3000, 0.3000, 0.5000, 0.5000, 0.5000, 0.3000, 0.3000, 0.5000};
//        double[] exitRates = {3.8000, 3.6000, 3.8000, 4.0000, 3.8000, 3.6000, 3.8000, 4.0000, 3.8000, 3.6000, 3.8000, 3.6000, 3.8000, 4.0000, 3.8000, 3.6000, 3.4000, 3.6000, 3.8000, 4.0000, 4.2000, 4.0000, 3.8000, 3.6000, 3.8000, 4.0000, 3.8000};

        int[] states = new int[]{0,1,2};
        path = new int[4];
        lstPaths = new ArrayList();
        allPaths(path, states, 0, lstPaths);
//        System.out.println(calculateProbability(lstPaths,rm));


        double evolTime = 1.0;
//        double p = trajectoryLikelihood_Miklos(jumpRates, exitRates);
//        System.out.println("Miklos's Code: " + Double.toString(p));
//        double p1 = trajectoryLikelihood(jumpRates, exitRates, evolTime);
//        System.out.println("Using Doubles: " + Double.toString(p1));
//        BigDecimal p2 = trajectoryLikelihood1(jumpRates, exitRates, evolTime);
//        System.out.println("Using BigDecimals: " + p2.toString());
//        double p4 = trajectoryLikelihood4(jumpRates, exitRates, evolTime);
//        System.out.println("Using Doubles 1: " + Double.toString(p4));
//        BigDecimal p1 = trajectoryLikelihood1(jumpRates, exitRates, evolTime);
//        System.out.println(p1.toString());


    /*
    System.out.println("BigDecimal:");
    System.out.println("Start: " + new Date());
    BigDecimal b = null;
    for (int i = 0; i <0 ; i++) {
    //            BigDecimal tl = trajectoryLikelihood4(jumpRates, exitRates, evolTime);
    BigDecimal a = new BigDecimal(1);
    b = exp(a);

    }
    System.out.println(b);
    System.out.println("Finish: " + new Date());
     */
    /*
    System.out.println("BigDecimal 1:");
    System.out.println("Start: " + new Date());
    for (int i = 0; i < 10000; i++) {
    //            BigDecimal tl = trajectoryLikelihood3(jumpRates, exitRates, evolTime);
    BigDecimal a = new BigDecimal(1.0, MathContext.DECIMAL128);
    BigDecimal b = new BigDecimal(1.7, MathContext.DECIMAL128);
    BigDecimal factor = new BigDecimal(1.0, MathContext.DECIMAL128).add(b, MathContext.DECIMAL128);

    }
    System.out.println("Finish: " + new Date());

     */


    /*        System.out.println("Start: " + new Date());
    for (int i = 0; i < 0; i++) {
    trajectoryLikelihood(jumpRates, exitRates, evolTime);
    //            double prob = trajectoryLikelihood(jumpRates, exitRates);
    }
    //        System.out.println(prob);
    System.out.println("Finish: " + new Date());
     */
    }

    static double trajectoryLikelihood_Miklos(double[] exitRate, double[] totalExitRates) {

        //DEBUGGING REPORT


        System.out.print("\nExit rates:\t");
        for (int i = 0; i < exitRate.length; i++) {
            System.out.print(exitRate[i] + ", ");
        }
        System.out.print("\nTotal exit rates:\t");
        for (int i = 0; i < totalExitRates.length; i++) {
            System.out.print(totalExitRates[i] + ", ");
        }
        System.out.println();

        double value = 1.0;
        for (int i = 0; i < exitRate.length; i++) {
            value *= exitRate[i];
        }
        // sorting the exitRates
        // lazy, bouble sorting....
        boolean sorted = false;
        while (!sorted) {
            sorted = true;
            for (int j = 1; j < totalExitRates.length; j++) {
                if (totalExitRates[j - 1] > totalExitRates[j]) {
                    double temp = totalExitRates[j - 1];
                    totalExitRates[j - 1] = totalExitRates[j];
                    totalExitRates[j] = temp;
                    sorted = false;
                }
            }
        }

        System.out.print("\nTotal exit rates after sorting:\t");
        for (int i = 0; i < totalExitRates.length; i++) {
            System.out.print(totalExitRates[i] + ", ");
        }
        System.out.println();

        double[][] coefficients = new double[totalExitRates.length][totalExitRates.length];
        double[][] tempco = new double[totalExitRates.length][totalExitRates.length];
        double[] rate = new double[totalExitRates.length];
        int[] d = new int[totalExitRates.length];
        int m = 1;
        coefficients[0][0] = 1.0;
        d[0] = 1;
        rate[0] = totalExitRates[0];

        for (int i = 1; i < totalExitRates.length; i++) {
            // updating other values
/*            if (Math.abs(exitRates[i] - rate[m - 1]) <
            Math.pow(0.0001, 1.0 / exitRates.length) ||
            Math.abs(Math.exp(-exitRates[i]) - Math.exp(-rate[m - 1])) / Math.exp(-rate[m - 1]) < 0.0001) {
            exitRates[i] = rate[m - 1];
            }
             */ for (int j = 0; j < m - 1; j++) {
                for (int k = 0; k < d[j]; k++) {
                    tempco[j][k] = 0.0;
                    double temp = 1.0 / (rate[j] - totalExitRates[i]);
                    for (int k1 = k; k1 < d[j]; k1++) {
                        tempco[j][k] += -coefficients[j][k1] * temp;
                        temp *= (k1 + 1.0) / (rate[j] - totalExitRates[i]);
                    }

                }
            }
            // if the exit rate is not new
            if (Math.abs(totalExitRates[i] - rate[m - 1]) < Math.pow(0.0001, 1.0 / totalExitRates.length) || Math.abs(Math.exp(-totalExitRates[i]) - Math.exp(-rate[m - 1])) / Math.exp(-rate[m - 1]) < 0.0001) {
                d[m - 1]++;
                for (int k = 1; k < d[m - 1]; k++) {
                    tempco[m - 1][k] = coefficients[m - 1][k - 1] / k;
                }
                tempco[m - 1][0] = 0;
                for (int j1 = 0; j1 < m - 1; j1++) {
                    double temp = 1.0 / (rate[j1] - totalExitRates[i]);
                    for (int k1 = 0; k1 < d[j1]; k1++) {
                        tempco[m - 1][0] += coefficients[j1][k1] * temp;
                        temp *= (k1 + 1.0) / (rate[j1] - totalExitRates[i]);
                    }
                }
            } else { //the exit rate is new
                // updating the old last row
                for (int k = 0; k < d[m - 1]; k++) {
                    tempco[m - 1][k] = 0.0;
                    double temp = 1.0 / (rate[m - 1] - totalExitRates[i]);
                    for (int k1 = k; k1 < d[m - 1]; k1++) {
                        tempco[m - 1][k] += -coefficients[m - 1][k1] * temp;
                        temp *= (k1 + 1.0) / (rate[m - 1] - totalExitRates[i]);
                    }
                }
                rate[m] = totalExitRates[i];
                d[m] = 1;
                m++;
                tempco[m - 1][0] = 0;
                for (int j1 = 0; j1 < m - 1; j1++) {
                    double temp = 1.0 / (rate[j1] - totalExitRates[i]);
                    for (int k1 = 0; k1 < d[j1]; k1++) {
                        tempco[m - 1][0] += coefficients[j1][k1] * temp;
                        temp *= (k1 + 1.0) / (rate[j1] - totalExitRates[i]);
                    }
                }

            }

            //updateing the coefficients
             System.out.println("Coefficients:");
            for (int j = 0; j < m; j++) {
                for (int k = 0; k < d[j]; k++) {
                    coefficients[j][k] = tempco[j][k];
                    //report for debugging reasons:
                    System.out.print(coefficients[j][k] + ",\t");
                }
                System.out.println();
            }

            System.out.println("\n\t\t\t****\t\t\t****\t\t\t****");

        }

        for (int i = 0; i < coefficients.length; i++) {
            for (int j = 0; j < coefficients[i].length; j++) {
                System.out.print(coefficients[i][j]);
                System.out.print(", ");
            }
            System.out.println("/");
        }

        double temp = 0.0;
        for (int i = 0; i < m; i++) {
            double temp1 = 0.0;
            for (int k = 0; k < d[i]; k++) {
                temp1 += coefficients[i][k];
            }
            temp += Math.exp(-rate[i]) * temp1;
        }

        System.out.print("\nTotal exit rates used actually:\t");
        for (int i = 0; i < totalExitRates.length; i++) {
            System.out.print(totalExitRates[i] + ", ");
        }
        System.out.println();
        System.out.println("Trajectory likelihood: " + (value * temp));

        return value * temp;
    }

    static public double trajectoryLikelihood(double[] rates, double[] totExitRates, double evolTime) {

        double[] zeta = totExitRates.clone();

        // Allocate space
        double[][] coefficients = new double[zeta.length][zeta.length];
        double[][] tempCoefficients = new double[zeta.length][zeta.length];
        double[] ksi = new double[zeta.length];
        int[] degeneracy = new int[zeta.length];

        int M = 0;
        ksi[0] = zeta[0];
        degeneracy[0] = 0;
        coefficients[0][0] = 1.0;

        for (int i = 1; i < zeta.length; i++) {

            int pos = find(zeta[i], ksi);
            if (pos == -1) {
                ksi[++M] = zeta[i];
                degeneracy[M] = 0;
            } else {
                degeneracy[pos]++;
            }

            for (int n = 0; n <= M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    if (ksi[n] != zeta[i]) {
                        tempCoefficients[n][k] = 0.0;
                        double factor = 1.0 / (ksi[n] - zeta[i]);
                        for (int j = k; j <= degeneracy[n]; j++) {
                            tempCoefficients[n][k] += coefficients[n][j] * factor;
                            factor *= (j + 1.0) / (ksi[n] - zeta[i]);
                        }
                        tempCoefficients[n][k] = -tempCoefficients[n][k];
                    } else { // if (ksi[n] != zeta[i])
                        if (k == 0) {
                            tempCoefficients[n][k] = 0.0;
                            for (int m = 0; m <= M; m++) {
                                if (m != n) {
                                    double factor = 1.0 / (ksi[m] - zeta[i]);
                                    for (int j = 0; j <= degeneracy[m]; j++) {
                                        tempCoefficients[n][k] += coefficients[m][j] * factor;
                                        factor *= (j + 1.0) / (ksi[m] - zeta[i]);
                                    }
                                } // end if (m != n)
                            } // for (int m = 0; m < M; m++ )
                        } else { // if ( k == 0 )
                            tempCoefficients[n][k] = coefficients[n][k - 1] / k;
                        } // end if ( k == 0 )
                    } // end if (ksi[n] != zeta[i])
                } // for (int k = 0; k < degeneracy[n]; k++)
            } //for (int n = 0; n < M; n++)

            for (int n = 0; n <= M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    coefficients[n][k] = tempCoefficients[n][k];
                }
            }

        } //for (int i = 1; i < zeta.length; i++)

        for (int i = 0; i < coefficients.length; i++) {
            for (int j = 0; j < coefficients[i].length; j++) {
                System.out.print(coefficients[i][j]);
                System.out.print(", ");
            }
            System.out.println("/");
        }

        double jump = 1.0;
        for (int i = 0; i < rates.length; i++) {
            jump *= rates[i];
        }

        double temp = 0.0;
        for (int n = 0; n <= M; n++) {
            double cT = 0.0;
            for (int k = 0; k <= degeneracy[n]; k++) {
                cT += coefficients[n][k] * Math.pow(evolTime, k);
            }
            temp += Math.exp(-ksi[n] * evolTime) * cT;
        }

        return jump * temp;

    }

    static public BigDecimal trajectoryLikelihood1(double[] rates, double[] totExitRates, double evolTime) {

        BigDecimal[] zeta = new BigDecimal[totExitRates.length];
        for (int i = 0; i < totExitRates.length; i++) {
            zeta[i] = BigDecimal.valueOf(totExitRates[i]);
        }

        // Allocate space
        BigDecimal[][] coefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[][] tempCoefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[] ksi = new BigDecimal[zeta.length];
        int[] degeneracy = new int[zeta.length];

        int M = 0;
        ksi[0] = zeta[0];
        degeneracy[0] = 0;
        coefficients[0][0] = new BigDecimal(1.0);

        for (int i = 1; i < zeta.length; i++) {

            int pos = find(zeta[i], ksi);
            if (pos == -1) {
                ksi[++M] = zeta[i];
                degeneracy[M] = 0;
            } else {
                degeneracy[pos]++;
            }

            for (int n = 0; n <= M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    if (ksi[n].compareTo(zeta[i]) != 0) {
                        tempCoefficients[n][k] = new BigDecimal(0.0, MathContext.DECIMAL128);
                        BigDecimal factor = new BigDecimal(1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128);
                        for (int j = k; j <= degeneracy[n]; j++) {
                            tempCoefficients[n][k] = tempCoefficients[n][k].add(coefficients[n][j].multiply(factor, MathContext.DECIMAL128), MathContext.DECIMAL128);
                            factor = factor.multiply(new BigDecimal(j + 1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
                        }
                        tempCoefficients[n][k] = new BigDecimal(-1.0, MathContext.DECIMAL128).multiply(tempCoefficients[n][k], MathContext.DECIMAL128);
                    } else // if (ksi[n] != zeta[i])
                    if (k == 0) {
                        tempCoefficients[n][k] = new BigDecimal(0.0, MathContext.DECIMAL128);
                        for (int m = 0; m <= M; m++) {
                            if (m != n) {
                                BigDecimal factor = new BigDecimal(1.0, MathContext.DECIMAL128).divide(ksi[m].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128);
                                for (int j = 0; j <= degeneracy[m]; j++) {
                                    tempCoefficients[n][k] = tempCoefficients[n][k].add(coefficients[m][j].multiply(factor, MathContext.DECIMAL128), MathContext.DECIMAL128);
                                    factor = factor.multiply(new BigDecimal(j + 1.0, MathContext.DECIMAL128).divide(ksi[m].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
                                }
                            } // for (int m = 0; m < M; m++ )
                        }
                    } else // if ( k == 0 )
                    {
                        tempCoefficients[n][k] = coefficients[n][k - 1].divide(new BigDecimal(k, MathContext.DECIMAL128), MathContext.DECIMAL128); //for (int n = 0; n < M; n++)
                    }
                }
            }
            for (int n = 0; n <= M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    coefficients[n][k] = tempCoefficients[n][k];
                }
            }

        } //for (int i = 1; i < zeta.length; i++)

        /*        for (int i = 0; i < coefficients.length; i++) {
        for (int j = 0; j < coefficients[i].length; j++) {
        System.out.print(coefficients[i][j] == null ? 0 : coefficients[i][j].toString());
        System.out.print(", ");
        }
        System.out.println("");
        }
         */
        BigDecimal jump = new BigDecimal(1.0, MathContext.DECIMAL128);
        for (int i = 0; i < rates.length; i++) {
            jump = jump.multiply(BigDecimal.valueOf(rates[i]), MathContext.DECIMAL128);
        }

        BigDecimal temp = new BigDecimal(0.0, MathContext.DECIMAL128);
        for (int n = 0; n <= M; n++) {
            BigDecimal cT = new BigDecimal(0.0, MathContext.DECIMAL128);
            for (int k = 0; k <= degeneracy[n]; k++) {
                cT = cT.add(coefficients[n][k].multiply(new BigDecimal(evolTime, MathContext.DECIMAL128).pow(k, MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
            }
            BigDecimal factor = exp(ksi[n].multiply(new BigDecimal(evolTime, MathContext.DECIMAL128), MathContext.DECIMAL128).negate(MathContext.DECIMAL128));
            temp = temp.add(factor.multiply(cT, MathContext.DECIMAL128));
        }

        /*        BigDecimal temp = new BigDecimal(0.0, MathContext.DECIMAL128);
        for (int n = 0; n <= M; n++) {
        BigDecimal cT = new BigDecimal(0.0, MathContext.DECIMAL128);
        for (int k = 0; k <= degeneracy[n]; k++)
        cT = cT.add(coefficients[n][k].multiply(new BigDecimal(evolTime, MathContext.DECIMAL128).pow(k, MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
        BigDecimal factor = exp(ksi[n].multiply(new BigDecimal(evolTime, MathContext.DECIMAL128), MathContext.DECIMAL128).negate(MathContext.DECIMAL128));
        temp = temp.add(factor.multiply(cT, MathContext.DECIMAL128));
        }
         */
        return jump.multiply(temp, MathContext.DECIMAL128);

    }

    static public BigDecimal trajectoryLikelihood3(double[] rates, double[] totExitRates, double evolTime) {

        // sorting the exitRates
        // lazy, bouble sorting....
        boolean sorted = false;
        while (!sorted) {
            sorted = true;
            for (int j = 1; j < totExitRates.length; j++) {
                if (totExitRates[j - 1] > totExitRates[j]) {
                    double temp = totExitRates[j - 1];
                    totExitRates[j - 1] = totExitRates[j];
                    totExitRates[j] = temp;
                    sorted = false;
                }
            }
        }

        BigDecimal[] zeta = new BigDecimal[totExitRates.length];
        for (int i = 0; i < totExitRates.length; i++) {
            zeta[i] = BigDecimal.valueOf(totExitRates[i]);
        }

        // Allocate space
        BigDecimal[][] coefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[][] tempCoefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[] ksi = new BigDecimal[zeta.length];
        int[] degeneracy = new int[zeta.length];

        int M = 0;
        ksi[0] = zeta[0];
        degeneracy[0] = 0;
        coefficients[0][0] = new BigDecimal(1.0);

        for (int i = 1; i < zeta.length; i++) {

            for (int n = 0; n < M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    tempCoefficients[n][k] = new BigDecimal(0.0, MathContext.DECIMAL128);
                    BigDecimal factor = new BigDecimal(1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128);
                    for (int j = k; j <= degeneracy[n]; j++) {
                        tempCoefficients[n][k] = tempCoefficients[n][k].add(coefficients[n][j].negate(MathContext.DECIMAL128).multiply(factor, MathContext.DECIMAL128), MathContext.DECIMAL128);
                        factor = factor.multiply(new BigDecimal(j + 1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
                    }
                }
            }

            int m = find(zeta[i], ksi);
            if (Math.abs(totExitRates[i] - ksi[M].doubleValue()) >= 0.01) {

                for (int k = 0; k <= degeneracy[M]; k++) {
                    tempCoefficients[M][k] = new BigDecimal(0.0, MathContext.DECIMAL128);
                    BigDecimal factor = new BigDecimal(1.0, MathContext.DECIMAL128).divide(ksi[M].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128);
                    for (int j = k; j <= degeneracy[M]; j++) {
                        tempCoefficients[M][k] = tempCoefficients[M][k].add(coefficients[M][j].negate(MathContext.DECIMAL128).multiply(factor, MathContext.DECIMAL128), MathContext.DECIMAL128);
                        factor = factor.multiply(new BigDecimal(j + 1.0, MathContext.DECIMAL128).divide(ksi[M].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
                    }
                }

                ksi[++M] = zeta[i];
                degeneracy[M] = 0;
                m = M;

                tempCoefficients[m][0] = new BigDecimal(0.0, MathContext.DECIMAL128);
                for (int n = 0; n < M; n++) {
                    if (ksi[n].compareTo(zeta[i]) != 0) {
                        BigDecimal factor = new BigDecimal(1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128);
                        for (int k = 0; k <= degeneracy[n]; k++) {
                            tempCoefficients[m][0] = tempCoefficients[m][0].add(coefficients[n][k].multiply(factor, MathContext.DECIMAL128), MathContext.DECIMAL128);
                            factor = factor.multiply(new BigDecimal(k + 1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
                        }
                    }
                }

            } else { // not new


                degeneracy[M]++;
                for (int k = 1; k <= degeneracy[M]; k++) {
                    tempCoefficients[M][k] = coefficients[M][k - 1].divide(new BigDecimal(k, MathContext.DECIMAL128), MathContext.DECIMAL128);
                }

                tempCoefficients[M][0] = new BigDecimal(0.0, MathContext.DECIMAL128);
                for (int n = 0; n < M; n++) {
                    BigDecimal factor = new BigDecimal(1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128);
                    for (int k = 0; k <= degeneracy[n]; k++) {
                        tempCoefficients[M][0] = tempCoefficients[M][0].add(coefficients[n][k].multiply(factor, MathContext.DECIMAL128), MathContext.DECIMAL128);
                        factor = factor.multiply(new BigDecimal(k + 1.0, MathContext.DECIMAL128).divide(ksi[n].subtract(zeta[i], MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
                    }
                }
            }

            for (int n = 0; n <= M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    coefficients[n][k] = tempCoefficients[n][k];
                //report for debugging reasons:
//                    System.out.print(coefficients[n][k].toString() + ",\t");
                }
//                System.out.println();
            }

//            System.out.println("\n\t\t\t****\t\t\t****\t\t\t****");

        } //for (int i = 1; i < zeta.length; i++)

        BigDecimal jump = new BigDecimal(1.0, MathContext.DECIMAL128);
        for (int i = 0; i < rates.length; i++) {
            jump = jump.multiply(BigDecimal.valueOf(rates[i]), MathContext.DECIMAL128);
        }

        BigDecimal temp = new BigDecimal(0.0, MathContext.DECIMAL128);
        for (int n = 0; n <= M; n++) {
            BigDecimal cT = new BigDecimal(0.0, MathContext.DECIMAL128);
            for (int k = 0; k <= degeneracy[n]; k++) {
                cT = cT.add(coefficients[n][k].multiply(new BigDecimal(evolTime, MathContext.DECIMAL128).pow(k, MathContext.DECIMAL128), MathContext.DECIMAL128), MathContext.DECIMAL128);
            }
            BigDecimal factor = exp(ksi[n].multiply(new BigDecimal(evolTime, MathContext.DECIMAL128), MathContext.DECIMAL128).negate(MathContext.DECIMAL128));
            temp = temp.add(factor.multiply(cT, MathContext.DECIMAL128));

        }
        return jump.multiply(temp, MathContext.DECIMAL128);

    }

    static public BigDecimal trajectoryLikelihood4(double[] rates, double[] totExitRates, double evolTime) {

        // sorting the exitRates
        // lazy, bouble sorting....
        boolean sorted = false;
        while (!sorted) {
            sorted = true;
            for (int j = 1; j < totExitRates.length; j++) {
                if (totExitRates[j - 1] > totExitRates[j]) {
                    double temp = totExitRates[j - 1];
                    totExitRates[j - 1] = totExitRates[j];
                    totExitRates[j] = temp;
                    sorted = false;
                }
            }
        }

        MathContext mc = MathContext.DECIMAL128;

        BigDecimal[] zeta = new BigDecimal[totExitRates.length];
        for (int i = 0; i < totExitRates.length; i++) {
            zeta[i] = BigDecimal.valueOf(totExitRates[i]);
        }

        // Allocate space
        BigDecimal[][] coefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[][] tempCoefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[] ksi = new BigDecimal[zeta.length];
        int[] degeneracy = new int[zeta.length];

        int M = 0;
        ksi[0] = zeta[0];
        degeneracy[0] = 0;
        coefficients[0][0] = new BigDecimal(1.0, mc);

        for (int i = 1; i < zeta.length; i++) {

            for (int n = 0; n < M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    tempCoefficients[n][k] = new BigDecimal(0.0, mc);
                    BigDecimal factor = new BigDecimal(1.0, mc).divide(ksi[n].subtract(zeta[i]), mc);
                    for (int j = k; j <= degeneracy[n]; j++) {
                        tempCoefficients[n][k] = tempCoefficients[n][k].add(coefficients[n][j].multiply(factor, mc));
                        factor = factor.multiply(new BigDecimal(j + 1.0, mc).divide(ksi[n].subtract(zeta[i]), mc), mc);
                    }
                    tempCoefficients[n][k] = tempCoefficients[n][k].negate();
                }
            }

//            int m = find(zeta[i], ksi);
//            if (Math.abs(totExitRates[i] - ksi[M].doubleValue()) >= 0.00001 ) {
            if (zeta[i].compareTo(ksi[M]) != 0) {

                for (int k = 0; k <= degeneracy[M]; k++) {
                    tempCoefficients[M][k] = new BigDecimal(0.0, mc);
                    BigDecimal factor = new BigDecimal(1.0, mc).divide(ksi[M].subtract(zeta[i]), mc);
                    for (int j = k; j <= degeneracy[M]; j++) {
                        tempCoefficients[M][k] = tempCoefficients[M][k].add(coefficients[M][j].multiply(factor, mc));
                        factor = factor.multiply(new BigDecimal(j + 1.0, mc).divide(ksi[M].subtract(zeta[i]), mc), mc);
                    }
                    tempCoefficients[M][k] = tempCoefficients[M][k].negate();
                }

                ksi[++M] = zeta[i];
                degeneracy[M] = 0;
//                m = M;

                tempCoefficients[M][0] = new BigDecimal(0.0, mc);
                for (int n = 0; n < M; n++) {
                    if (ksi[n].compareTo(zeta[i]) != 0) {
                        BigDecimal factor = new BigDecimal(1.0, mc).divide(ksi[n].subtract(zeta[i]), mc);
                        for (int k = 0; k <= degeneracy[n]; k++) {
                            tempCoefficients[M][0] = tempCoefficients[M][0].add(coefficients[n][k].multiply(factor, mc));
                            factor = factor.multiply(new BigDecimal(k + 1.0, mc).divide(ksi[n].subtract(zeta[i]), mc), mc);
                        }
                    }
                }
            } else { // not new
                degeneracy[M]++;
                for (int k = 1; k <= degeneracy[M]; k++) {
                    tempCoefficients[M][k] = coefficients[M][k - 1].divide(new BigDecimal(k, mc), mc);
                }

                tempCoefficients[M][0] = new BigDecimal(0.0, mc);
                for (int n = 0; n < M; n++) {
                    BigDecimal factor = new BigDecimal(1.0, mc).divide(ksi[n].subtract(zeta[i]), mc);
                    for (int k = 0; k <= degeneracy[n]; k++) {
                        tempCoefficients[M][0] = tempCoefficients[M][0].add(coefficients[n][k].multiply(factor, mc));
                        factor = factor.multiply(new BigDecimal(k + 1.0, mc).divide(ksi[n].subtract(zeta[i]), mc), mc);
                    }
                }
            }

            for (int n = 0; n <= M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    coefficients[n][k] = tempCoefficients[n][k];
                //report for debugging reasons:
//                    System.out.print(coefficients[n][k].toString() + ",\t");
                }
//                System.out.println();
            }

//            System.out.println("\n\t\t\t****\t\t\t****\t\t\t****");

        } //for (int i = 1; i < zeta.length; i++)

        BigDecimal jump = new BigDecimal(1.0, mc);
        for (int i = 0; i < rates.length; i++) {
            jump = jump.multiply(BigDecimal.valueOf(rates[i]), mc);
        }

        BigDecimal temp = new BigDecimal(0.0, mc);
        for (int n = 0; n <= M; n++) {
            BigDecimal cT = new BigDecimal(0.0, mc);
            for (int k = 0; k <= degeneracy[n]; k++) {
                cT = cT.add(coefficients[n][k].multiply(new BigDecimal(evolTime, mc).pow(k, mc), mc));
            }
            BigDecimal factor = exp(ksi[n].multiply(new BigDecimal(evolTime, mc), mc).negate());
            temp = temp.add(factor.multiply(cT, mc));

        }
        return jump.multiply(temp, mc);

    }

    static public double logsum(double a, double b) {

        if (a == 0) {
            return b;
        }
        if (b == 0) {
            return a;
        }

        return a + Math.log(1 + Math.exp(b - a));
    }

    static private int find(double value, double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            if (value == vector[i]) {
                return i;
            }
        }
        return -1;
    }

    static private int find(BigDecimal value, BigDecimal[] vector) {
        for (int i = 0; i < vector.length; i++) {
//            if (vector[i] != null && value.compareTo(vector[i]) == 0)
            if (vector[i] != null && value.subtract(vector[i]).abs().compareTo(new BigDecimal(0.5)) < 0) {
                return i;
            }
        }
        return -1;
    }

    static private BigDecimal naturalE(int prec_dig) {
        BigDecimal sum = new BigDecimal("1");
        BigDecimal fact = new BigDecimal("1");
        BigDecimal del = new BigDecimal("1");
        BigDecimal one = new BigDecimal("1");
        BigDecimal ten = new BigDecimal("10");
        int prec_bits = (prec_dig * 332) / 100;

        one = one.setScale(prec_dig, BigDecimal.ROUND_DOWN);
        for (int i = 0; i < prec_dig; i++) {
            del = del.multiply(ten);
        }
        for (int i = 1; i < prec_bits; i++) {
            fact = fact.multiply(new BigDecimal(i));
            fact = fact.setScale(prec_dig, BigDecimal.ROUND_DOWN);
            sum = sum.add(one.divide(fact, BigDecimal.ROUND_DOWN));
            if (del.compareTo(fact) < 0) {
                break;
            }
        }
        return sum.setScale(prec_dig, BigDecimal.ROUND_DOWN);
    }

    static private BigDecimal exp_series(BigDecimal x) // abs(x)<=0.5, prec digits 
    {                                   // prec digits returned 
        BigDecimal fact = new BigDecimal("1"); // factorial
        BigDecimal xp = new BigDecimal("1"); // power of x
        BigDecimal y = new BigDecimal("1"); // sum of series on x
        int n;
        int prec = 301;

        n = (2 * prec) / 3;
        for (int i = 1; i < n; i++) {
            fact = fact.multiply(new BigDecimal(i));
            fact = fact.setScale(prec, BigDecimal.ROUND_DOWN);
            xp = xp.multiply(x);
            xp = xp.setScale(prec, BigDecimal.ROUND_DOWN);
            y = y.add(xp.divide(fact, BigDecimal.ROUND_DOWN));
        }
        y = y.setScale(prec, BigDecimal.ROUND_DOWN);
//        System.out.println("exp_series x=" + x);
//        System.out.println("           y=" + y);
        return y;
    }

    static private BigDecimal exp_series(BigDecimal x, int prec) // abs(x)<=0.5, prec digits 
    {                                   // prec digits returned 

        BigDecimal xp = new BigDecimal("1"); // power of x
        BigDecimal y = new BigDecimal("1"); // sum of series on x
        int n;

        n = (2 * prec) / 3;
        for (int i = 1; i < n; i++) {
            xp = xp.multiply(x);
            xp = xp.setScale(prec, BigDecimal.ROUND_DOWN);
            y = y.add(xp.divide(fact[i], BigDecimal.ROUND_DOWN));
        }
        y = y.setScale(prec, BigDecimal.ROUND_DOWN);
//        System.out.println("exp_series x=" + x);
//        System.out.println("           y=" + y);
        return y;
    }

    static private BigDecimal exp(BigDecimal x) {
        int prec = 301;
        BigDecimal natural_e = naturalE(prec);

        BigDecimal y = new BigDecimal("1.0");  // sum of series on xc
        BigDecimal xc;                         // x - j
        BigDecimal ep = natural_e;             // e^j
        BigDecimal j = new BigDecimal("1");
        BigDecimal one = new BigDecimal("1.0");
        BigDecimal half = new BigDecimal("0.5");
        BigDecimal xp;                         // positive, then invert

        if (x.abs().compareTo(half) < 0) {
            return exp_series(x);
        }
        if (x.compareTo(new BigDecimal("0")) > 0) // positive
        {
            while (x.compareTo(j.add(one)) > 0) {
                ep = ep.multiply(natural_e);
                j = j.add(one);
            }
            xc = x.subtract(j);
            y = ep.multiply(exp_series(xc));
            y = y.setScale(prec, BigDecimal.ROUND_DOWN);
            return y;
        } else // negative
        {
            xp = x.negate();
            while (xp.compareTo(j.add(one)) > 0) {
                ep = ep.multiply(natural_e);
                j = j.add(one);
            }
            xc = xp.subtract(j);
            y = ep.multiply(exp_series(xc));
            y = y.setScale(prec, BigDecimal.ROUND_DOWN);

            BigDecimal epsilon = new BigDecimal("1");
            int bits = 1000; // precision in bits = about 3.32 precision in digits
            for (int i = 0; i < bits; i++) {
                epsilon = epsilon.multiply(new BigDecimal("0.5"));
            }
            epsilon = epsilon.setScale(prec, BigDecimal.ROUND_DOWN);

            return (one.add(epsilon)).divide(y, BigDecimal.ROUND_DOWN);
        }
    } // end exp

    static private BigDecimal exp(BigDecimal x, int prec) {

        BigDecimal y = new BigDecimal("1.0");  // sum of series on xc
        BigDecimal xc;                         // x - j
        BigDecimal ep = natural_e;             // e^j
        BigDecimal j = new BigDecimal("1");
        BigDecimal one = new BigDecimal("1.0");
        BigDecimal half = new BigDecimal("0.5");
        BigDecimal xp;                         // positive, then invert

        if (x.abs().compareTo(half) < 0) {
            return exp_series(x);
        }
        if (x.compareTo(new BigDecimal("0")) > 0) // positive
        {
            while (x.compareTo(j.add(one)) > 0) {
                ep = ep.multiply(natural_e);
                j = j.add(one);
            }
            xc = x.subtract(j);
            y = ep.multiply(exp_series(xc));
            y = y.setScale(prec, BigDecimal.ROUND_DOWN);
            return y;
        } else // negative
        {
            xp = x.negate();
            while (xp.compareTo(j.add(one)) > 0) {
                ep = ep.multiply(natural_e);
                j = j.add(one);
            }
            xc = xp.subtract(j);
            y = ep.multiply(exp_series(xc));
            y = y.setScale(prec, BigDecimal.ROUND_DOWN);

            return (one.add(epsilon)).divide(y, BigDecimal.ROUND_DOWN);
        }
    } // end exp

    static private void allPaths(int[] path, int[] states, int currentLength, ArrayList paths) {

        if (currentLength == path.length) {
            paths.add(path);
        } else if (currentLength == 0) {
            for (int i = 0; i < states.length; i++) {
                path[0] = states[i];
                allPaths(path, states, currentLength+1, paths);
            }
        } else {
            for (int i = 0; i < states.length; i++) {
                if (path[currentLength - 1] != states[i]) {
                    int[] newPath = (int[]) path.clone();
                    newPath[currentLength] = states[i];
                    allPaths(newPath, states, currentLength + 1, paths);
                }
            }
        }
    }

    static private double calculateProbability(ArrayList lstPaths, double[][] rm) {

        double prob = 0.0;
        Iterator it = lstPaths.iterator();
        while (it.hasNext()) {
            int[] path = (int[])it.next();
            System.out.print("Path: ");
            for (int i = 0; i < path.length; i++) {
                System.out.print(path[i]);
                System.out.print(" ");
            }
            System.out.println();

            double[] jumpRates = getJumpRates(path, rm);
            double[] exitRates = getExitRates(path, rm);
            prob += trajectoryLikelihood(jumpRates, exitRates,1.0);
        }

        return prob;

    }

    static private double[] getJumpRates(int[] path, double[][] rm) {
        double[] jumpRates = new double[path.length-1];
        for (int i = 0; i < path.length-1; i++) {
            jumpRates[i] = rm[path[i]][path[i+1]];

        }
        return jumpRates;
    }

    static private double[] getExitRates(int[] path, double[][] rm) {
        double[] exitRates = new double[path.length];
        for (int i = 0; i < path.length; i++) {
            exitRates[i] = -rm[path[i]][path[i]];

        }
        return exitRates;
    }
}

package lineageTree.substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Node;
import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.DoubleStream;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

public class GSC_Test {

    GeneralScarringLoss scarringModel;

    @Before
    public void setup(){

        RealParameter lossRate = new RealParameter("1.0");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);
    }

    @Test
    public void testRateMatrix(){

        //System.out.println("relative rates :\n" +
        //        Arrays.toString(scarringModel.getRelativeRates()) + "\n");
        System.out.println("rate matrix :");
        double[][] rateM = scarringModel.getRateMatrix();
        for(int i = 0; i < rateM.length; i++)
            System.out.println(Arrays.toString(rateM[i]));


        double[][] correctMatrix = new double[rateM.length][rateM.length];
        correctMatrix[0][0] = 0;
        for (int i=1; i<rateM.length-1; i++){
            // insert scar rates into first row
            correctMatrix[0][i] = scarringModel.sRate.getValues()[i-1];

            // insert loss rates into last column
            correctMatrix[i-1][rateM.length-1] = scarringModel.lRate.getValue();
        }
        correctMatrix[rateM.length-2][rateM.length-1] = scarringModel.lRate.getValue();

        for (int i=0; i<rateM.length; i++){
            assertArrayEquals("Assert matrix entries", correctMatrix[i], rateM[i], 1e-15);
        }

    }

    @Test
    public void testTransitionProbabilities1(){
        /* Test transition probabilities for scarring and loss window
        */
        double startTime = 25;
        double endTime = 23;
        double rate = 1.0;

        System.out.println("freqs = \n" + Arrays.toString(scarringModel.getFrequencies()) + "\n");

        int len = scarringModel.getStateCount();
        double[] prob = new double[len*len];
        scarringModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob);

        double[] correctMatrix = new double[]{
                0.002478752176666, 0.066428265529973, 0.066428265529973, 0.864664716763387,
                0, 0.135335283236613, 0, 0.864664716763387,
                0, 0, 0.135335283236613, 0.864664716763387,
                0, 0, 0, 1.0
        };

        System.out.println("\ntransition prob :\n" + Arrays.toString(prob));
        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals("Assert row sum == 1",1, sum, 1e-15);
        }

        assertArrayEquals("Assert correct transition probabilities",
                correctMatrix, prob, 1e-15);

    }

    @Test
    public void testTransitionProbabilities2(){
        /* Test transition probabilities for loss window
         */
        double startTime = 2;
        double endTime = 0;
        double rate = 1.0;

        System.out.println("freqs = \n" + Arrays.toString(scarringModel.getFrequencies()) + "\n");

        int len = scarringModel.getStateCount();
        double[] prob = new double[len*len];
        scarringModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob);

        double[] correctMatrix = new double[]{
                0.135335283236613, 0, 0, 0.864664716763387,
                0, 0.135335283236613, 0, 0.864664716763387,
                0, 0, 0.135335283236613, 0.864664716763387,
                0, 0, 0, 1.0
        };

        System.out.println("\ntransition prob :\n" + Arrays.toString(prob));
        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals("Assert row sum == 1",1, sum, 1e-15);
        }

        assertArrayEquals("Assert correct transition probabilities for loss",
                correctMatrix, prob, 0.01);

    }

    @Test
    public void test_different_parameters(){

        //set up scarring model
        RealParameter lossRate = new RealParameter("0.04");
        RealParameter scarRates = new RealParameter("10.0 20.0");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // set up time window
        double startTime = 25;
        double endTime = 23;
        double rate = 1.0;

        //set up correct transition probabilities
        double[] correctMatrix = new double[]{
                0.000000000000000, 0.307705448795545, 0.615410897591091, 0.076883653613364,
                0.000000000000000, 0.923116346386636, 0, 0.076883653613364,
                0, 0, 0.923116346386636, 0.076883653613364,
                0, 0, 0, 1
        };

        // retrieve transition probabilities from scarring model
        int len = scarringModel.getStateCount();;
        double[] prob = new double[len*len];
        scarringModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob);

        // test correctness of transition probabilities
        System.out.println("\ntransition prob :\n" + Arrays.toString(prob));
        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals("Assert row sum == 1",1, sum, 1e-15);
        }

        assertArrayEquals("Assert correct transition probabilities for loss",
                correctMatrix, prob, 0.01);


    }

    /*
    Test for single scarring rate
     */
    @Test
    public void test_different_parameters2(){

        //set up scarring model
        RealParameter lossRate = new RealParameter("0.0");
        RealParameter scarRates = new RealParameter("10");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // set up time window
        double startTime = 25;
        double endTime = 23;
        double rate = 1;

        //set up correct transition probabilities
        double[] correctMatrix = new double[]{
                0.000000002061154, 0.999999997938846, 0.0,
                0.000000000000000, 1.0, 0,
                0, 0, 1
        };

        // retrieve transition probabilities from scarring model
        int len = scarringModel.getStateCount();;
        double[] prob = new double[len*len];
        scarringModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob);

        // test correctness of transition probabilities
        System.out.println("\ntransition prob :\n" + Arrays.toString(prob));
        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals("Assert row sum == 1",1, sum, 1e-15);
        }

        assertArrayEquals("Assert correct transition probabilities for one scarring rate",
                correctMatrix, prob, 1e-15);
    }

    /*
    Test for single scarring rate
     */
    @Test
    public void test_clock_rate(){

        //set up scarring model
        RealParameter lossRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("10");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // set up time window
        double startTime = 25;
        double endTime = 23;
        double rate = 0.01;

        //set up correct transition probabilities
        double[] correctMatrix = new double[]{
                0.815462371187293, 0.180545618156699, 0.003992010656009,
                0.0, 0.996007989343991, 0.003992010656009,
                0, 0, 1
        };

        // retrieve transition probabilities from scarring model
        int len = scarringModel.getStateCount();;
        double[] prob = new double[len*len];
        scarringModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob);

        // test correctness of transition probabilities
        System.out.println("\ntransition prob :\n" + Arrays.toString(prob));
        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals("Assert row sum == 1",1, sum, 1e-15);
        }

        assertArrayEquals("Assert correct transition probabilities for one scarring rate",
                correctMatrix, prob, 1e-15);
    }


}

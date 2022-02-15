package lineageTree.substitutionmodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class GeneralScarringLoss extends SubstitutionModel.Base {

        final public Input<List<RealParameter>> scarringRatesInput = new Input<>("scarringRates",
                "Rates, at which scars types are introduced into the scarring region",
                new ArrayList<>(), Input.Validate.REQUIRED);
        final public Input<RealParameter> lossRateInput = new Input<>("lossRate",
                "Rate, at which scarring regions are lost", Input.Validate.REQUIRED);
        public Input<RealParameter> scarringHeightInput = new Input<>("scarringHeight",
                "Duration between the onset of scarring and sampling of the cells", Input.Validate.REQUIRED);
        public Input<RealParameter> scarringDurationInput = new Input<>("scarringDuration",
                "Duration of the scarring process", Input.Validate.REQUIRED);


        /**
         * flag to indicate matrix is up to date *
         */
        protected boolean updateMatrixScar = true;
        protected boolean updateMatrixLoss = true;


        RealParameter scarringHeightP;
        RealParameter scarringDurationP;
        double[] frequencies;
        double[][] rateMatrix;
        RealParameter sRate;
        RealParameter lRate;

    @Override
    public void initAndValidate() {

        // one state for each scar type + unedited + lost
        nrOfStates = scarringRatesInput.get().get(0).getDimension() + 2;
        rateMatrix = new double[nrOfStates][nrOfStates];

        sRate = scarringRatesInput.get().get(0);
        lRate = lossRateInput.get();

        // assert positive rates
        double sumScarringRates = 0;

        // add scar rates to rate matrix
        for (int i=0; i<sRate.getDimension(); i++){

            double scarringRate = sRate.getValues()[i];
            if (scarringRate < 0) {
                throw new RuntimeException("All scarring rates must be positive!");
            }
            sumScarringRates += scarringRate;
            rateMatrix[0][i+1] = scarringRate;
        }

        lRate = lossRateInput.get();
        double lossRate = lRate.getValue();

        if (lossRate < 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }
        for (int i = 0; i<nrOfStates-1; i++){
            rateMatrix[i][nrOfStates-1] = lossRate;
        }

        // center root frequency on unedited state
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;

        scarringHeightP = scarringHeightInput.get();
        scarringDurationP = scarringDurationInput.get();
    }

    /*
    Calculate transition probability matrix for loss without scarring
    The transition probability matrix for loss with scarring differs only in the first row.
     */
    public void getLossProbabilities(double[] matrix, double expOfDeltaLoss){

        // fill diagonal and final column
        for (int i=0; i<nrOfStates; i++){
            for (int j=0; j<nrOfStates; j++){

                if ( i==j ){
                    matrix[i*nrOfStates + j] = expOfDeltaLoss;
                }else if(j == nrOfStates-1){
                    matrix[i*nrOfStates + j] = 1 - expOfDeltaLoss;
                }else{
                    matrix[i*nrOfStates + j] = 0;
                }
            }
        }
        // set final diagonal element to 1
        matrix[nrOfStates * nrOfStates - 1] = 1;
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {

        double lossRate = lRate.getValue();
        Double[] scarRates = sRate.getValues();

        double scarringHeight = scarringHeightP.getValue();
        double scarringDuration = scarringDurationP.getValue();

        //multiply by joint branch rate from site model
        lossRate *= rate;
        for (int i=0; i < scarRates.length; i++){
            scarRates[i] *= rate;
        }

        double delta = startTime - endTime;
        double expOfDeltaLoss = Math.exp(-delta * lossRate);

        // calculate transition probabilities for loss process
        getLossProbabilities(matrix, expOfDeltaLoss);

        Stream<Double> scarSum = Stream.of(scarRates);
        Double scarRateSum = scarSum.reduce(0.0, (subtotal, element) -> subtotal + element);

        // for loss & scarring, add the scarring transition probabilities
        if ( (endTime >= (scarringHeight - scarringDuration)) & (endTime < scarringHeight) & (scarRateSum>0) ) {

                // fill first row
                matrix[0] = Math.exp(-delta * (lossRate + scarRateSum));
                for (int i = 0; i < nrOfStates - 2; i++) {
                    matrix[i + 1] = (scarRates[i] * expOfDeltaLoss - scarRates[i] * Math.exp(-delta * (lossRate + scarRateSum))) / scarRateSum;
                }
                matrix[nrOfStates - 1] = 1 - expOfDeltaLoss;

        }
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {return null;}

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof IntegerData;
    }

    @Override
    public double[] getFrequencies() {
        return frequencies;
    }

    public double getScarringHeight(){return scarringHeightP.getValue();}

    public double getScarringDuration(){return scarringDurationP.getValue();}

    public double[][] getRateMatrix(){return rateMatrix;}
}


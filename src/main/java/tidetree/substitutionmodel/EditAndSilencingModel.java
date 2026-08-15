package tidetree.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.type.RealVector;
import beast.base.spec.type.RealScalar;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.NonNegativeReal;

import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.spec.evolution.substitutionmodel.Base;
import beast.base.evolution.tree.Node;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

/**
 * @author Sophie Seidel
 **/
@Description("Computes the transition probabilities from an unedited state into edited states or a silenced state.")
public class EditAndSilencingModel extends Base {

        final public Input<RealVector<NonNegativeReal>> editRatesInput = new Input<>("editRates",
                "Rates at which edits are introduced into " +
                        "the genomic barcode during the editing window",

                Input.Validate.REQUIRED);

        final public Input<RealScalar<NonNegativeReal>> silencingRateInput = new Input<>("silencingRate",
                "Rate at which barcodes are silenced " +
                        "throughout the entire experiment", Input.Validate.REQUIRED);


        public Input<RealScalar<PositiveReal>> editHeightInput = new Input<>("editHeight",
                "Duration between the onset of edit and sampling of the cells", Input.Validate.REQUIRED);

        public Input<RealScalar<NonNegativeReal>> editDurationInput = new Input<>("editDuration",
                "Duration of the edit process", Input.Validate.REQUIRED);


        /**
         * flag to indicate matrix is up to date *
         */
        protected boolean updateMatrixScar = true;
        protected boolean updateMatrixLoss = true;


        RealScalar editHeightP;
        RealScalar editDurationP;
        double[] frequencies;
        double[][] rateMatrix;
        RealVector editRate_;
        RealScalar silencingRate_;

    @Override
    public void initAndValidate() {

        // one state for each edit type + unedited + lost
        nrOfStates = editRatesInput.get().size() + 2;
        rateMatrix = new double[nrOfStates][nrOfStates];

        editRate_ = editRatesInput.get();
        silencingRate_ = silencingRateInput.get();

        // assert positive rates
        //double sumEditRates = 0;

        // add edit rates to rate matrix
        for (int i=0; i<editRate_.size(); i++){

            double editRate = editRate_.get(i);
            if (editRate < 0) {
                throw new RuntimeException("All edit rates must be positive!");
            }
            //sumEditRates += editRate;
            rateMatrix[0][i+1] = editRate;
        }

        silencingRate_ = silencingRateInput.get();
        double silencingRate = silencingRate_.get();

        if (silencingRate < 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }
        for (int i = 0; i<nrOfStates-1; i++){
            rateMatrix[i][nrOfStates-1] = silencingRate;
        }

        // center root frequency on unedited state
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;

        editHeightP = editHeightInput.get();
        editDurationP = editDurationInput.get();
    }

    /*
    Calculate transition probability matrix for loss without edit
    The transition probability matrix for loss with edit differs only in the first row.
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

    //Rename time to height
    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {

        double silencingRate = silencingRate_.get();
        double[] editRates = new double[editRate_.size()];
	for (int i = 0; i < editRates.length; i++) {
	    editRates[i] = editRate_.get(i);
	}

        double editHeight = editHeightP.get();
        double editDuration = editDurationP.get();

        //multiply by joint branch rate from site model
        silencingRate *= rate;
        for (int i=0; i < editRates.length; i++){
            editRates[i] *= rate;
        }

        double delta = startTime - endTime;
        double expOfDeltaLoss = Math.exp(-delta * silencingRate);

        // calculate transition probabilities for loss process
        getLossProbabilities(matrix, expOfDeltaLoss);

        double editRateSum = Arrays.stream(editRates).sum();

        // for loss & edit, add the edit transition probabilities
        if ( (endTime >= (editHeight - editDuration)) & (endTime < editHeight) & (editRateSum>0) ) {

                // fill first row
                matrix[0] = Math.exp(-delta * (silencingRate + editRateSum));
                for (int i = 0; i < nrOfStates - 2; i++) {
                    matrix[i + 1] = (editRates[i] * expOfDeltaLoss - editRates[i] *
                            Math.exp(-delta * (silencingRate + editRateSum))) / editRateSum;
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

    @Override
    public double[] getRateMatrix(Node node) {
	double[] flat = new double[nrOfStates * nrOfStates];
      for (int i = 0; i < nrOfStates; i++) {
          System.arraycopy(rateMatrix[i], 0, flat, nrOfStates * i, nrOfStates);
      }
      return flat;
  }

    public double getEditHeight(){return editHeightP.get();}

    public double getEditDuration(){return editDurationP.get();}

    public double[][] getRateMatrix(){return rateMatrix;}
}


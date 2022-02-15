package lineageTree.distributions;

import beast.core.Input;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import lineageTree.substitutionmodel.GeneralScarringLoss;


public class CustomCore extends BeerLikelihoodCore {


    private double scarringStart;
    private double scarringStop;
    
    public CustomCore(int nrOfStates) {
        super(nrOfStates);
    }

    /**
     * initializes partial likelihood arrays.
     *  @param nodeCount           the number of nodes in the tree
     * @param patternCount        the number of patterns
     * @param matrixCount         the number of matrices (i.e., number of categories)
     * @param integrateCategories whether sites are being integrated over all matrices
     * @param m_siteModel
     */

    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities,
                           GeneralScarringLoss substitutionModel, SiteModelInterface.Base m_siteModel) {

        super.initialize(nodeCount, patternCount, matrixCount, integrateCategories, useAmbiguities);
        
        
        this.scarringStart = substitutionModel.getScarringHeight();
        this.scarringStop = scarringStart - substitutionModel.getScarringDuration();
    }

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param child1 the 'child 1' node
     * @param child2 the 'child 2' node
     * @param node the 'parent' node
     */
    public void calculatePartials(Node child1, Node child2, Node node) {

        int nodeIndex1 = child1.getNr();
        int nodeIndex2 = child2.getNr();
        int nodeIndex3 = node.getNr();

            if (states[nodeIndex1] != null) {
                if (states[nodeIndex2] != null) {
                    calculateStatesStatesPruning(
                            states[nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            states[nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                } else {
                    calculateStatesPartialsPruning(states[nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            partials[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                }
            } else {
                if (states[nodeIndex2] != null) {
                    calculateStatesPartialsPruning(states[nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                } else {
                    calculatePartialsPartialsPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                            partials[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                            partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
                }
            }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }

//
//        int k =0;
//        for (int i = 0; i < patternCount; i++) {
//            double f = 0.0;
//
//            for (int j = 0; j < stateCount; j++) {
//                f += partials[currentPartialsIndices[nodeIndex3]][nodeIndex3][k];
//                k++;
//            }
//            if (f == 0.0) {
//                Logger.getLogger("error").severe("A partial likelihood (node index = " + nodeIndex3 + ", pattern = "+ i +") is zero for all states.");
//            }
//        }
    }

    public void calculatePartials(Node child1, Node node) {

        int nodeIndex1 = child1.getNr();
        int nodeIndex3 = node.getNr();

        if (states[nodeIndex1] != null) {
            calculateStatesPruning(states[nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
        } else {
            calculatePartialsPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],
            matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1], partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
        }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }


    void calculatePartialsForCrossBranches(Node parent, Node child1, Node child2, boolean bool1, boolean bool2,
                                           SiteModelInterface.Base m_siteModel, GeneralScarringLoss substitutionModel,
                                           double branchRate){

        double[] nodePartials = new double[partialsSize];
        double[][] helperNodePartials = new double[4][partialsSize];

        double[] jointBranchRates0 = new double[nrOfMatrices];
        double[] jointBranchRates1 = new double [nrOfMatrices];

        for (int i = 0; i < nrOfMatrices; i++){
            jointBranchRates0[i] = m_siteModel.getRateForCategory(i, child1) * branchRate;
            jointBranchRates1[i] = m_siteModel.getRateForCategory(i, child2) * branchRate;
        }

        Node[] children = new Node[]{child1, child2};
        int[] childIndices = new int[]{child1.getNr(), child2.getNr()};
        double[][] jointbranchRates = new double[][]{jointBranchRates0, jointBranchRates1};

        double[] heightsBeforeParent = new double[2];
        boolean[] needsIntermediates = new boolean[2];

        double[] probs0 = new double[(nrOfStates) * (nrOfStates)];
        double[] probs1 = new double[(nrOfStates) * (nrOfStates)];
        double [] matrix0 = new double[nrOfMatrices * matrixSize];
        double [] matrix1 = new double[nrOfMatrices * matrixSize];

        // determine the number of helper nodes to calculate the partials over
        int nIntermediateNodes;
        double[] intNodeTimes;

        boolean [] parentBeforeChildAfterScarringStart = new boolean[] {
                (parent.getHeight() > scarringStart) & (child1.getHeight() < scarringStart),
                (parent.getHeight() > scarringStart) & (child2.getHeight() < scarringStart)
        };
        boolean [] parentBeforeChildAfterScarringStop = new boolean[]{
                (parent.getHeight() > scarringStop) & (child1.getHeight() < scarringStop),
                (parent.getHeight() > scarringStop) & (child2.getHeight() < scarringStop)
        };

        // for each child calculate partials upwards until the helper node below the parent
        // that's only necessary if the rate matrix changes on the branch!
        for (int i=0; i< children.length; i++){

            // if parent above scarring window (sw), child in sw
            if (parentBeforeChildAfterScarringStart[i] & !parentBeforeChildAfterScarringStop[i]){
                nIntermediateNodes = 1;
                intNodeTimes = new double[]{children[i].getHeight(), scarringStart, Double.NEGATIVE_INFINITY};
                heightsBeforeParent[i] = scarringStart;
                needsIntermediates[i] = true;

                helperNodePartials[i * 2 + 1] = calculatePartialsBeforeParent(parent, children[i], i, childIndices[i],
                        intNodeTimes, jointbranchRates[i], nIntermediateNodes, helperNodePartials, substitutionModel);

                // if parent in sw and child below sw
            }else if(!parentBeforeChildAfterScarringStart[i] & parentBeforeChildAfterScarringStop[i]){
                nIntermediateNodes = 1;
                intNodeTimes = new double[]{children[i].getHeight(), scarringStop, Double.NEGATIVE_INFINITY};
                heightsBeforeParent[i] = scarringStop;
                needsIntermediates[i] = true;

                helperNodePartials[i * 2 + 1] = calculatePartialsBeforeParent(parent, children[i], i, childIndices[i],
                        intNodeTimes, jointbranchRates[i], nIntermediateNodes, helperNodePartials, substitutionModel);

                // if parent above sw and child below sw
            }else if (parentBeforeChildAfterScarringStart[i] & parentBeforeChildAfterScarringStop[i]){
                nIntermediateNodes = 2;
                intNodeTimes = new double[]{children[i].getHeight(), scarringStop, scarringStart};
                heightsBeforeParent[i] = scarringStart;
                needsIntermediates[i] = true;

                helperNodePartials[i * 2 + 1] = calculatePartialsBeforeParent(parent, children[i], i, childIndices[i],
                        intNodeTimes, jointbranchRates[i], nIntermediateNodes, helperNodePartials, substitutionModel);

                // parent and child in the same rate matrix regime
            }else{
                heightsBeforeParent[i] = children[i].getHeight();
                needsIntermediates[i] = false;
            }
        }


        //calculate partials at parent
        // if intermediates are necessary their *partials* were computed above -> no leaf check necessary
        if (needsIntermediates[0]) {
            for (int k = 0; k < nrOfMatrices; k++) {
                substitutionModel.getTransitionProbabilities(null, parent.getHeight(), heightsBeforeParent[0],
                        jointBranchRates0[k], probs0);
                System.arraycopy(probs0, 0, matrix0, k * matrixSize, matrixSize);
            }

            if (needsIntermediates[1]) {
                for (int k = 0; k < nrOfMatrices; k++) {
                    substitutionModel.getTransitionProbabilities(null, parent.getHeight(), heightsBeforeParent[1],
                            jointBranchRates1[k], probs1);
                    System.arraycopy(probs1, 0, matrix1, k * matrixSize, matrixSize);
                }
                calculatePartialsPartialsPruning(helperNodePartials[1], matrix0, helperNodePartials[3], matrix1, nodePartials);

            } else {
                for (int k = 0; k < nrOfMatrices; k++) {
                    substitutionModel.getTransitionProbabilities(null, parent.getHeight(), child2.getHeight(),
                            jointBranchRates1[k], probs1);
                    System.arraycopy(probs1, 0, matrix1, k * matrixSize, matrixSize);
                }

                if (child2.isLeaf()) {
                    int[] states = new int[nrOfPatterns];
                    getNodeStates(childIndices[1], states);

                    calculateStatesPartialsPruning(states, matrix1, helperNodePartials[1], matrix0, nodePartials);

                } else {
                    double[] partials = new double[partialsSize];
                    getNodePartials(childIndices[1], partials);
                    calculatePartialsPartialsPruning(partials, matrix1, helperNodePartials[1], matrix0, nodePartials);
                }
            }
        } else {
            for (int k = 0; k < nrOfMatrices; k++) {
                substitutionModel.getTransitionProbabilities(null, parent.getHeight(), child1.getHeight(),
                        jointBranchRates0[k], probs0);
                System.arraycopy(probs0, 0, matrix0, k * matrixSize, matrixSize);
            }

            if (child1.isLeaf()){
                int[] states1 = new int[nrOfPatterns];
                getNodeStates(childIndices[0], states1);

                if (needsIntermediates[1]){
                    for (int k = 0; k < nrOfMatrices; k++) {
                        substitutionModel.getTransitionProbabilities(null, parent.getHeight(), heightsBeforeParent[1],
                                jointBranchRates1[k], probs1);
                        System.arraycopy(probs1, 0, matrix1, k * matrixSize, matrixSize);
                    }

                    // because helper node partials were computed, we can directly take them as input
                    calculateStatesPartialsPruning(states1, probs0, helperNodePartials[3], probs1, nodePartials);
                }else{
                    // because no helper nodes were needed, we have to check child2 being a leaf or not
                    for (int k = 0; k < nrOfMatrices; k++) {
                        substitutionModel.getTransitionProbabilities(null, parent.getHeight(), child2.getHeight(),
                                jointBranchRates1[k], probs1);
                        System.arraycopy(probs1, 0, matrix1, k * matrixSize, matrixSize);
                    }

                    if (child2.isLeaf()){
                        int[] states2 = new int [nrOfPatterns];
                        getNodeStates(childIndices[1], states2);
                        calculateStatesStatesPruning(states1, matrix0, states2, matrix1, nodePartials);
                    }else{
                        double[] partials2 = new double[partialsSize];
                        getNodePartials(childIndices[1], partials2);
                        calculateStatesPartialsPruning(states1, matrix0, partials2, matrix1, nodePartials);
                    }
                }
            }else {
                double[] partials1 = new double[partialsSize];
                getNodePartials(childIndices[0], partials1);

                if (needsIntermediates[1]){
                    for (int k = 0; k < nrOfMatrices; k++) {
                        substitutionModel.getTransitionProbabilities(null, parent.getHeight(), heightsBeforeParent[1],
                                jointBranchRates1[k], probs1);
                        System.arraycopy(probs1, 0, matrix1, k * matrixSize, matrixSize);
                    }
                    calculatePartialsPartialsPruning(partials1, matrix0, helperNodePartials[3], matrix1, nodePartials);
                }else{
                    for (int k = 0; k < nrOfMatrices; k++) {
                        substitutionModel.getTransitionProbabilities(null, parent.getHeight(), child2.getHeight(),
                                jointBranchRates1[k], probs1);
                        System.arraycopy(probs1, 0, matrix1, k * matrixSize, matrixSize);
                    }

                    if (child2.isLeaf()){
                        int[] states2 = new int [nrOfPatterns];
                        getNodeStates(childIndices[1], states2);
                        calculateStatesPartialsPruning(states2, matrix1, partials1, matrix0, nodePartials);
                    }else{
                        double[] partials2 = new double[partialsSize];
                        getNodePartials(childIndices[1], partials2);
                        calculatePartialsPartialsPruning(partials2, matrix1, partials1, matrix0, nodePartials);
                    }
                }
            }
        }
        setCurrentNodePartials(parent.getNr(), nodePartials);
    }

    void calculatePartialsForCrossBranches(Node parent, Node child, boolean bool1, boolean bool2,
                                           SiteModelInterface.Base m_siteModel, GeneralScarringLoss substitutionModel,
                                           double branchRate){

        double[] nodePartials = new double[partialsSize];
        double[][] helperNodePartials = new double[4][partialsSize];

        double[] jointBranchRates = new double[nrOfMatrices];


        for (int i = 0; i < nrOfMatrices; i++){
            jointBranchRates[i] = m_siteModel.getRateForCategory(i, child) * branchRate;
        }

        int childIndx = child.getNr();

        double heightsBeforeParent = -1.0;
        boolean needsIntermediates = false;

        double[] probs = new double[(nrOfStates) * (nrOfStates)];
        double [] matrix = new double[nrOfMatrices * matrixSize];

        // determine the number of helper nodes to calculate the partials over
        int nIntermediateNodes;
        double[] intNodeTimes;

        boolean  parentBeforeChildAfterScarringStart =
                (parent.getHeight() > scarringStart) & (child.getHeight() < scarringStart);
        boolean parentBeforeChildAfterScarringStop =
                (parent.getHeight() > scarringStop) & (child.getHeight() < scarringStop);

        // calculate partials upwards until the helper node below the parent
        // that's only necessary if the rate matrix changes on the branch!


        // if parent above scarring window (sw), child in sw
        if (parentBeforeChildAfterScarringStart & !parentBeforeChildAfterScarringStop){
                nIntermediateNodes = 1;
                intNodeTimes = new double[]{child.getHeight(), scarringStart, Double.NEGATIVE_INFINITY};
                heightsBeforeParent = scarringStart;
                needsIntermediates = true;

                helperNodePartials[1] = calculatePartialsBeforeParent(parent, child, 0, childIndx,
                        intNodeTimes, jointBranchRates, nIntermediateNodes, helperNodePartials, substitutionModel);

                // if parent in sw and child below sw
        }else if(!parentBeforeChildAfterScarringStart & parentBeforeChildAfterScarringStop){
                nIntermediateNodes = 1;
                intNodeTimes = new double[]{child.getHeight(), scarringStop, Double.NEGATIVE_INFINITY};
                heightsBeforeParent = scarringStop;
                needsIntermediates = true;

                helperNodePartials[1] = calculatePartialsBeforeParent(parent, child, 0, childIndx,
                        intNodeTimes, jointBranchRates, nIntermediateNodes, helperNodePartials, substitutionModel);

                // if parent above sw and child below sw
        }else if (parentBeforeChildAfterScarringStart & parentBeforeChildAfterScarringStop){
                nIntermediateNodes = 2;
                intNodeTimes = new double[]{child.getHeight(), scarringStop, scarringStart};
                heightsBeforeParent = scarringStart;
                needsIntermediates = true;

                helperNodePartials[1] = calculatePartialsBeforeParent(parent, child, 0, childIndx,
                        intNodeTimes, jointBranchRates, nIntermediateNodes, helperNodePartials, substitutionModel);

                // parent and child in the same rate matrix regime
        }else{
                heightsBeforeParent = child.getHeight();
                needsIntermediates = false;
        }



        //calculate partials at parent
        // if intermediates are necessary their *partials* were computed above -> no leaf check necessary
        if (needsIntermediates) {
            for (int k = 0; k < nrOfMatrices; k++) {
                substitutionModel.getTransitionProbabilities(null, parent.getHeight(), heightsBeforeParent,
                        jointBranchRates[k], probs);
                System.arraycopy(probs, 0, matrix, k * matrixSize, matrixSize);
            }
            calculatePartialsPruning(helperNodePartials[1], matrix, nodePartials);

        } else {
            for (int k = 0; k < nrOfMatrices; k++) {
                substitutionModel.getTransitionProbabilities(null, parent.getHeight(), child.getHeight(),
                        jointBranchRates[k], probs);
                System.arraycopy(probs, 0, matrix, k * matrixSize, matrixSize);
            }

            if (child.isLeaf()){
                int[] states = new int[nrOfPatterns];
                getNodeStates(childIndx, states);

                calculateStatesPruning(states, matrix, nodePartials);

            }else {
                double[] partials = new double[partialsSize];
                getNodePartials(childIndx, partials);

                calculatePartialsPruning(partials, matrix, nodePartials);
            }
        }
        setCurrentNodePartials(parent.getNr(), nodePartials);
    }

    double[] calculatePartialsBeforeParent(Node parent, Node child, int i, int childIndex, double[] intNodeTimes, double[] jointBranchRates, int nIntermediateNodes, double[][] helperNodePartials,
                                           GeneralScarringLoss substitutionModel) {

        double [] matrix = new double[nrOfMatrices * matrixSize];
        double[] probs = new double[matrixSize]; //TODO!! why is this +1??

        if (child.isLeaf()) {
            int[] states = new int[nrOfPatterns];
           getNodeStates(childIndex, states);

            if (nIntermediateNodes == 0) {
                for (int k = 0; k < nrOfMatrices; k++){
                    substitutionModel.getTransitionProbabilities(null, parent.getHeight(), child.getHeight(), jointBranchRates[k], probs);
                    System.arraycopy(probs, 0, matrix, k * matrixSize, matrixSize);
                }
                helperNodePartials[i * 2 + 1] = calculateStatesPruning(states, matrix, helperNodePartials[i * 2 + 1]);



            } else {
                for (int k = 0; k < nrOfMatrices; k++) {
                    substitutionModel.getTransitionProbabilities(null, intNodeTimes[1], intNodeTimes[0], jointBranchRates[k], probs);
                    System.arraycopy(probs, 0, matrix, k * matrixSize, matrixSize);
                }
                helperNodePartials[i * 2] = calculateStatesPruning(states, matrix, helperNodePartials[i * 2]);
                helperNodePartials[i * 2 + 1] = helperNodePartials[i * 2];


                if (nIntermediateNodes > 1) {
                    for (int k = 0; k < nrOfMatrices; k++) {
                        substitutionModel.getTransitionProbabilities(null, intNodeTimes[2], intNodeTimes[1], jointBranchRates[k], probs);
                        System.arraycopy(probs, 0, matrix, k * matrixSize, matrixSize);
                    }
                    helperNodePartials[i * 2 + 1] = calculatePartialsPruning(helperNodePartials[i * 2], matrix, helperNodePartials[i * 2 + 1]);

                }
            }
        } else {
            double[] partials = new double[partialsSize];
            getNodePartials(childIndex, partials);

            if (nIntermediateNodes == 0) {
                for (int k = 0; k < nrOfMatrices; k++) {
                    substitutionModel.getTransitionProbabilities(null, parent.getHeight(), child.getHeight(), jointBranchRates[k], probs);
                    System.arraycopy(probs, 0, matrix, k * matrixSize, matrixSize);
                }
                    helperNodePartials[i * 2 + 1] = calculatePartialsPruning(partials, matrix, helperNodePartials[i * 2 + 1]);
            } else {

                helperNodePartials[i * 2] = partials;
                for (int j = 0; j < nIntermediateNodes; j++) {
                    for (int k = 0; k < nrOfMatrices; k++) {
                        substitutionModel.getTransitionProbabilities(null, intNodeTimes[j + 1], intNodeTimes[j], jointBranchRates[k], probs);
                        System.arraycopy(probs, 0, matrix, k * matrixSize, matrixSize);
                    }
                    helperNodePartials[i * 2 + 1] = calculatePartialsPruning(helperNodePartials[i * 2], matrix, helperNodePartials[i * 2 + 1]);
                    helperNodePartials[i * 2] = helperNodePartials[i * 2 + 1];

                }
            }
        }

        return helperNodePartials[i * 2 + 1];
    }


    protected double[] calculateStatesPruning(int[] stateIndex1, double[] matrices1,
                                              double[] partials3) {
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int state1 = stateIndex1[k];

                int w = l * matrixSize;

                if (state1 < nrOfStates) {

                    for (int i = 0; i < nrOfStates; i++) {

                        partials3[v] = matrices1[w + state1];

                        v++;
                        w += nrOfStates;
                    }

                } else {
                    // single child has a gap or unknown state so set partials to 1
                    for (int j = 0; j < nrOfStates; j++) {
                        partials3[v] = 1.0;
                        v++;
                    }
                }
            }
        }
        return partials3;
    }

    // calculate partials for single child nodes
    protected double[] calculatePartialsPruning(double[] partials1, double[] matrices1,
                                                double[] partials3) {
        double sum1;

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int w = l * matrixSize;

                for (int i = 0; i < nrOfStates; i++) {

                    sum1 = 0.0;

                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += matrices1[w] * partials1[v + j];

                        w++;
                    }

                    partials3[u] = sum1;
                    u++;
                }
                v += nrOfStates;
            }
        }
        return partials3;
    }
}

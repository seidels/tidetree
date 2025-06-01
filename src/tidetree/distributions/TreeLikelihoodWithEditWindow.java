package tidetree.distributions;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import tidetree.substitutionmodel.EditAndSilencingModel;

import java.util.Arrays;

/**
 * @author Sophie Seidel
 **/
@Description("This class calculates the tree likelihood, i.e. the probability of seeing an alignment given a tree under" +
        "an editing model")
public class TreeLikelihoodWithEditWindow extends GenericTreeLikelihood {

    public Input<RealParameter> originInput = new Input<>("origin",
            "Start of the cell division process, usually start of the experiment.",
            Input.Validate.OPTIONAL);
    final public Input<Frequencies> rootFrequenciesInput =
            new Input<>("rootFrequencies", "prior state frequencies at root, optional", Input.Validate.OPTIONAL);

    final public Input<Boolean> useScalingInput = new Input<Boolean>("useScaling", "Whether or not to scale the log likelihood", false,
            Input.Validate.OPTIONAL);



    /**
     * INPUT objects.
     * Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected EditAndSilencingModel substitutionModel;
    protected double editStart;
    protected double editStop;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;
    protected RealParameter origin;
    protected boolean useOrigin = false;
    protected Node originNode;

    /**
     * Parameters derived from input objects.
     */
    // number of categories across site -> number of transition rate matrices
    protected int nrOfMatrices;
    // in the alignment matrix a pattern is a unique column of alignment entries; columns == sites with the same pattern
    // have the same likelihood. thus we save computation by not recomputing the site likelihood for identical patterns
    protected int nrOfPatterns;
    // number of states in substitution rate matrix
    //protected int nrOfStates;
    //protected int matrixSize;

    /**
     * MEMORY ALLOCATION
     */

    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;
    protected int matrixSize;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;


    /**
     * memory allocation for likelihoods for each of the patterns *
     */
   protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials
     */
   protected double[] m_fRootPartials;

    /**
     * CALCULATION ENGINE, containing methods that calculate the subtree likelihoods
     */
    protected CustomCore likelihoodCore;
    


    /**
     * Getters
     */
    public LikelihoodCore getLikelihoodCore() {
        return likelihoodCore;
    }

    public EditAndSilencingModel getSubstitutionModel() {return substitutionModel;}

    @Override
    public void initAndValidate() {
        if(logP == 0){
            int a = 1;
        }
        // Init input objects
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        int nodeCount = treeInput.get().getNodeCount();
        if (originInput.get() != null){
            origin = originInput.get();
            useOrigin = true;
            nodeCount = nodeCount + 1;
        }
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        Alignment alignment = dataInput.get();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = (EditAndSilencingModel) m_siteModel.substModelInput.get();
        editStart = substitutionModel.getEditHeight();
        editStop = editStart - substitutionModel.getEditDuration();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }

        // init parameters from input objects
        int nrOfStates = substitutionModel.getStateCount();
        /*if (nrOfStates != dataInput.get().getMaxStateCount()){
            throw new IllegalArgumentException("The number of sequence states provided in userDataType \n" +
                    " and the number of edit rates +2 have to be equal! \n" +
                    "Number of sequence states: " + dataInput.get().getMaxStateCount() + "\n"+
                    "Number of edit rates +2: " + nrOfStates);
        }*/
        nrOfPatterns = dataInput.get().getPatternCount();
        nrOfMatrices = m_siteModel.getCategoryCount();

        //Allocate memory
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];
        likelihoodCore = new CustomCore(nrOfStates);
        patternLogLikelihoods = new double[nrOfPatterns];
        m_fRootPartials = new double[nrOfPatterns * nrOfStates];

        matrixSize = (nrOfStates + 1) * (nrOfStates + 1);
        probabilities = new double[(nrOfStates + 1) * (nrOfStates + 1)];
        Arrays.fill(probabilities, 1.0);

        if (useScalingInput.get()) {
            boolean useScaling = true;
        }

        // init calculation engine
        initCore();

        // print startup message
        String className = getClass().getSimpleName();
        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));

    }

    protected void initCore() {
        // currently do not handle ambiguities
        boolean m_useAmbiguities = false;
        boolean m_useTipLikelihoods = false;

        int nodeCount = treeInput.get().getNodeCount();
        if (useOrigin){
            nodeCount = nodeCount + 1;
        }
        likelihoodCore.initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities,
                substitutionModel, m_siteModel);

        int extNodeCount;
        int intNodeCount;
        // origin node should not influence int and extNode calculations
        if (useOrigin){
            nodeCount = nodeCount -1;
            extNodeCount = nodeCount / 2 + 1;
            intNodeCount = nodeCount / 2 + 1;
        }else{
            extNodeCount = nodeCount / 2 + 1;
            intNodeCount = nodeCount / 2;
        }



        setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }

        if (useScalingInput.get()) {
            likelihoodCore.setUseScaling(1.01);
        }

    }

    /**
     * This method performs a recursion down the tree and computes the partial likelihood at each node
     * and stores them in the likelihood core
     * @param node
     * @return integer indicating whether the likelihood for the parent node has to be recalculated
     */
    protected int traverse(Node node) {

        int update = 2; //currently always recalculate for every node! (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();
        final double nodeHeight = node.getHeight();
        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        if ( (!node.isRoot() || (node.isRoot() && useOrigin)) &&  (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex]) ) {
            m_branchLengths[nodeIndex] = branchTime;
            Node parent = node.getParent();
            if(node.isRoot()){
                parent= new Node();
                parent.setHeight(origin.getValue());
            }

            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);

            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            int update1 = traverse(child1);
            final double childHeight1 = child1.getHeight();

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);
            final double childHeight2 = child2.getHeight();

            update1 = 2;

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {

                    boolean parentBeforeChildrenAfterEditHeight = (nodeHeight > editStart) &
                            ((childHeight1 < editStart) | (childHeight2 < editStart));
                    boolean parentBeforeChildrenAfterEditStop = (nodeHeight > editStop) &
                            ((childHeight1 < editStop) | (childHeight2 < editStop));

                    if (parentBeforeChildrenAfterEditHeight | parentBeforeChildrenAfterEditStop) {

                        likelihoodCore.calculatePartialsForCrossBranches(node, child1, child2,
                                parentBeforeChildrenAfterEditHeight, parentBeforeChildrenAfterEditStop,
                                m_siteModel, substitutionModel, branchRate);

                    }else {
                        likelihoodCore.calculatePartials(child1, child2, node);
                    }
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }


                if (node.isRoot()) {
                    // calculate the partials until the origin
                    if (useOrigin){
                        Double originHeight = origin.getValue();
                        originNode = new Node();
                        originNode.setHeight(originHeight);
                        originNode.setNr(node.getNr() + 1);

                        // calculate likelihood until the origin
                        //Node rootNode = node.copy(); //TODO copy here only the things necessary for proper calculation

                        Double rootHeight = node.getHeight();

                        //swap root with origin
                        //node.setHeight(originHeight);

                        if (m_siteModel.integrateAcrossCategories()) {

                            boolean parentBeforeChildrenAfterEditHeight = (originHeight > editStart) &
                                    (rootHeight < editStart);
                            boolean parentBeforeChildrenAfterEditStop = (originHeight > editStop) &
                                    (rootHeight < editStop);

                            if (parentBeforeChildrenAfterEditHeight | parentBeforeChildrenAfterEditStop) {

                                likelihoodCore.calculatePartialsForCrossBranches(originNode, node,
                                        parentBeforeChildrenAfterEditHeight, parentBeforeChildrenAfterEditStop,
                                        m_siteModel, substitutionModel, branchRate);

                            }else {
                                likelihoodCore.calculatePartials(node, originNode);
                            }
                        }
                        //node = rootNode;
                        // set height back to true root height
                        //  calculate the pattern likelihoods
                        final double[] proportions = m_siteModel.getCategoryProportions(originNode);
                        likelihoodCore.integratePartials(originNode.getNr(), proportions, m_fRootPartials);

                        double[] rootFrequencies = substitutionModel.getFrequencies();
                        if (rootFrequenciesInput.get() != null) {
                            rootFrequencies = rootFrequenciesInput.get().getFreqs();
                        }
                        likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);

                    }else{

                        //  calculate the pattern likelihoods at the root
                        final double[] proportions = m_siteModel.getCategoryProportions(node);
                        likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                        //likelihoodCore.getNodePartials(nodeIndex, m_fRootPartials);

                        double[] rootFrequencies = substitutionModel.getFrequencies();
                        if (rootFrequenciesInput.get() != null) {
                            rootFrequencies = rootFrequenciesInput.get().getFreqs();
                        }
                        likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                    }
                }
            }

        }
        return 2; //update;
    }



    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int i;
            int[] states = new int[patternCount];
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length == 1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }

    /**
     * @param taxon the taxon name as a string
     * @param data  the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     * or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
            if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {

        final TreeInterface tree = treeInput.get();

        if(tree.getNodeCount() == 1){
            return -10.0;
        }
        if(useOrigin) {
            Double originHeight = origin.getValue();
            if (tree.getRoot().getHeight() >= originHeight) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        try {
            if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
                calcLogP();
        } catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
            //TODO deleted possibility to turn off scaling
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    void calcLogP() {
        logP = 0.0;
        //if (useAscertainedSitePatterns) {
        //    final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
        //   for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
        //       logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
        //  }
        //} else {
        for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
            logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
        }
    }

    @Override
    public void store() {
        storedLogP = logP;
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }


    @Override
    public void restore() {
        logP = storedLogP;
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
    }


    @Override
    protected boolean requiresRecalculation() {

        return true;
    }
        /*hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();*/

}

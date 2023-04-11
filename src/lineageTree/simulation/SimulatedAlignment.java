package lineageTree.simulation;/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */



import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.BEASTClassLoader;
import beast.util.PackageManager;
import beast.util.Randomizer;
import feast.nexus.CharactersBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import lineageTree.substitutionmodel.EditAndSilencingModel;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Sophie Seidel
 * Adapted from feast's SimulatedAlignment class by Tim Vaughan
 * Main change is to introduce the editing window into the simulation
 */
@Description("A more flexible alignment simulator. Doesn't require " +
        "pre-specification of number of taxa.")
public class SimulatedAlignment extends Alignment {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree down which to simulate sequence evolution.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin", "Start of the process, usually the experiment",
            Input.Validate.OPTIONAL);

    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "Site model to use in simulation.",
            Input.Validate.REQUIRED);

    public Input<BranchRateModel> branchRateModelInput = new Input<>(
            "branchRateModel", "Clock model",
            Input.Validate.REQUIRED);

    public Input<Integer> sequenceLengthInput = new Input<>(
            "sequenceLength",
            "Length of sequence to simulate.",
            Input.Validate.REQUIRED);

    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "Name of file (if any) simulated alignment should be saved to.");

    private Tree tree;
    private Double originHeight;
    private Boolean useOrigin=Boolean.FALSE;
    private SiteModel siteModel;
    private BranchRateModel branchRateModel;
    private int seqLength;
    private DataType dataType;
    private String ancestralSeqStr;
    public int[] helperNodeSeq;
    public boolean[] branchType;

    public SimulatedAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        if (originInput.get() != null){
            useOrigin = Boolean.TRUE;
            originHeight = originInput.get().getValue();
        }

        siteModel = siteModelInput.get();
        branchRateModel = branchRateModelInput.get();
        seqLength = sequenceLengthInput.get();

        siteModel.getSubstitutionModel().getStateCount();

        sequences.clear();

        grabDataType();

        simulate();

        super.initAndValidate();

        // Write simulated alignment to disk if required
        if (outputFileNameInput.get() != null) {
            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                NexusBuilder nb = new NexusBuilder();
                nb.append(new TaxaBlock(new TaxonSet(this)));
                nb.append(new CharactersBlock(this));
                nb.write(pstream);
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to file "
                        + outputFileNameInput.get() + ".");
            }
        }
    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate() {
        int nTaxa = tree.getLeafNodeCount();

        double[] categoryProbs = siteModel.getCategoryProportions(tree.getRoot());

        int nCategories = siteModel.getCategoryCount();
        int nStates = siteModel.getSubstitutionModel().getStateCount();
        double[][] transitionProbs = new double[nCategories][nStates*nStates];

        int[][] alignment = new int[nTaxa][seqLength];

        int[] categories = new int[seqLength];
        for (int i=0; i<seqLength; i++)
            categories[i] = Randomizer.randomChoicePDF(categoryProbs);

        Node root = tree.getRoot();
        int[] parentSequence = new int[seqLength];
        helperNodeSeq = new int[seqLength];
        helperNodeSeq[0] = -2; // for testing purposes
        double[] frequencies = siteModel.getSubstitutionModel().getFrequencies();
        for (int i=0; i<parentSequence.length; i++)
            parentSequence[i] = Randomizer.randomChoicePDF(frequencies);

        ancestralSeqStr = dataType.encodingToString(parentSequence);

        traverse(root, parentSequence,
                categories, transitionProbs,
                alignment);

        for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
            String seqString = dataType.encodingToString(alignment[leafIdx]);

            String taxonName;
            if (tree.getNode(leafIdx).getID() != null)
                taxonName = tree.getNode(leafIdx).getID();
            else
                taxonName = "t" + leafIdx;

            sequenceInput.setValue(new Sequence(taxonName, seqString), this);
        }
    }

    /**
     * Traverse a tree, simulating a sequence alignment down it.
     *
     * @param node Node of the tree
     * @param parentSequence Sequence at the parent node in the tree
     * @param categories Mapping from sites to categories
     * @param transitionProbs transition probabilities
     * @param regionAlignment alignment for particular region
     */
    private void traverse(Node node,
            int[] parentSequence,
            int[] categories, double[][] transitionProbs,
            int[][] regionAlignment) {


        // get scarring time boundaries
        EditAndSilencingModel scarringModel = (EditAndSilencingModel)siteModel.getSubstitutionModel();
        double scarringStart = scarringModel.getEditHeight();
        double scarringStop = scarringStart -  scarringModel.getEditDuration();

        // if the process started at the origin and root is after the onset of scarring,
        // evolve parentSeq from origin to root
        if (node.isRoot() & useOrigin & node.getHeight() < scarringStart){

            boolean parentBeforeChildAfterScarringStart = (originHeight > scarringStart) &
                    (node.getHeight() < scarringStart);
            boolean parentBeforeChildAfterScarringStop = (originHeight > scarringStop) &
                    (node.getHeight() < scarringStop);

            // set up variables for child sequence
            int[] rootSequence = new int[parentSequence.length];
            int nStates = dataType.getStateCount();
            double[] charProb = new double[nStates];

            if(parentBeforeChildAfterScarringStart | parentBeforeChildAfterScarringStop){

                Node originNode = new Node();
                originNode.setHeight(originHeight);

                // change parent node and sequence to helper node sequence and continue as usual
                Node helperNode = evolveToHelperNode(parentBeforeChildAfterScarringStart, parentBeforeChildAfterScarringStop,
                        scarringStart, scarringStop, originNode, node, transitionProbs, parentSequence, categories);

                // Calculate transition probabilities
                for (int i=0; i<siteModel.getCategoryCount(); i++) {
                    final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRateModel.getRateForBranch(node);
                    siteModel.getSubstitutionModel().getTransitionProbabilities(
                            node, helperNode.getHeight(), node.getHeight(),
                            jointBranchRate,
                            transitionProbs[i]);
                }
                //System.out.println("Transition probabilities for final part" + Arrays.toString(transitionProbs[0]));

                // Draw characters on child sequence
                for (int i=0; i<rootSequence.length; i++) {
                    int category = categories[i];
                    System.arraycopy(transitionProbs[category],
                            helperNodeSeq[i]*nStates, charProb, 0, nStates);
                    rootSequence[i] = Randomizer.randomChoicePDF(charProb);
                }
            }else {

                // Calculate transition probabilities
                for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                    final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRateModel.getRateForBranch(node);
                    siteModel.getSubstitutionModel().getTransitionProbabilities(
                            node, originHeight, node.getHeight(),
                            jointBranchRate,
                            transitionProbs[i]);
                }

                // Draw characters on child sequence
                for (int i = 0; i < rootSequence.length; i++) {
                    int category = categories[i];
                    System.arraycopy(transitionProbs[category],
                            parentSequence[i] * nStates, charProb, 0, nStates);
                    rootSequence[i] = Randomizer.randomChoicePDF(charProb);
                }
            }
            parentSequence = rootSequence;
        }

        for (Node child : node.getChildren()) {

            boolean parentBeforeChildAfterScarringStart = (node.getHeight() > scarringStart) &
                    (child.getHeight() < scarringStart);
            boolean parentBeforeChildAfterScarringStop = (node.getHeight() > scarringStop) &
                    (child.getHeight() < scarringStop);

            // set up variables for child sequence
            int[] childSequence = new int[parentSequence.length];
            int nStates = dataType.getStateCount();
            double[] charProb = new double[nStates];

            if(parentBeforeChildAfterScarringStart | parentBeforeChildAfterScarringStop){

                // change parent node and sequence to helper node sequence and continue as usual
                Node helperNode = evolveToHelperNode(parentBeforeChildAfterScarringStart, parentBeforeChildAfterScarringStop,
                        scarringStart, scarringStop, node, child, transitionProbs, parentSequence, categories);

                // Calculate transition probabilities
                for (int i=0; i<siteModel.getCategoryCount(); i++) {
                    final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRateModel.getRateForBranch(node);
                    siteModel.getSubstitutionModel().getTransitionProbabilities(
                            child, helperNode.getHeight(), child.getHeight(),
                            jointBranchRate,
                            transitionProbs[i]);
                }
                 //System.out.println("Transition probabilities for final part" + Arrays.toString(transitionProbs[0]));

                // Draw characters on child sequence
                for (int i=0; i<childSequence.length; i++) {
                    int category = categories[i];
                    System.arraycopy(transitionProbs[category],
                            helperNodeSeq[i]*nStates, charProb, 0, nStates);
                    childSequence[i] = Randomizer.randomChoicePDF(charProb);
                }
            }else {

                // Calculate transition probabilities
                for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                    final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRateModel.getRateForBranch(node);
                    siteModel.getSubstitutionModel().getTransitionProbabilities(
                            child, node.getHeight(), child.getHeight(),
                            jointBranchRate,
                            transitionProbs[i]);
                }

                // Draw characters on child sequence
                for (int i = 0; i < childSequence.length; i++) {
                    int category = categories[i];
                    System.arraycopy(transitionProbs[category],
                            parentSequence[i] * nStates, charProb, 0, nStates);
                    childSequence[i] = Randomizer.randomChoicePDF(charProb);
                }
            }

            if (child.isLeaf()) {
                System.arraycopy(childSequence, 0,
                        regionAlignment[child.getNr()], 0, childSequence.length);
            } else {
                traverse(child, childSequence,
                        categories, transitionProbs,
                        regionAlignment);
            }
        }
    }

    /**
     * evolves the sequence to the helper node
     */
    private Node evolveToHelperNode(boolean boolScarringStart, boolean boolScarringStop, double scarringStart,
                                    double scarringStop, Node node, Node child, double[][] transitionProbs,
                                    int[] parentSequence, int[] categories){

        branchType = new boolean[]{boolScarringStart, boolScarringStop}; // for test
        Node helperNode = new Node();
        double helperNodeHeight;
        Node pNode = new Node();
        pNode.setHeight(node.getHeight());
        int nHelperNodes;

        if (boolScarringStart & boolScarringStop){
            nHelperNodes = 2;
            helperNodeHeight = scarringStart;

        }else if(boolScarringStart){
            nHelperNodes = 1;
            helperNodeHeight = scarringStart;

        }else if(boolScarringStop){
            nHelperNodes = 1;
            helperNodeHeight = scarringStop;

        }else {
            throw new RuntimeException("This should never evaluate!");
        }

        for (int h=0; h<nHelperNodes; h++){

            helperNode.setHeight(helperNodeHeight);

            // Calculate transition probabilities until helper node
            for (int i=0; i<siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRateModel.getRateForBranch(node);
                siteModel.getSubstitutionModel().getTransitionProbabilities(
                        helperNode, pNode.getHeight(), helperNode.getHeight(),
                        jointBranchRate,
                        transitionProbs[i]);
            }

            //System.out.println("Transition probabilities for " + Arrays.toString(transitionProbs[0]));

            // Draw characters on helper node sequence
            int nStates = dataType.getStateCount();
            double[] charProb = new double[nStates];
            for (int i=0; i<helperNodeSeq.length; i++) {
                int category = categories[i];
                System.arraycopy(transitionProbs[category],
                        parentSequence[i]*nStates, charProb, 0, nStates);
                helperNodeSeq[i] = Randomizer.randomChoicePDF(charProb);
            }
            pNode.setHeight(helperNodeHeight);
            helperNodeHeight = scarringStop;
            parentSequence = helperNodeSeq;
        }

        return helperNode;
    }

    /**
     * HORRIBLE function to identify data type from given description.
     */
    private void grabDataType() {
        if (userDataTypeInput.get() != null) {
            dataType = userDataTypeInput.get();
        } else {

            List<String> dataTypeDescList = new ArrayList<>();
            List<String> classNames = PackageManager.find(beast.evolution.datatype.DataType.class, "beast.evolution.datatype");
            for (String className : classNames) {
                try {
                    DataType thisDataType = (DataType) BEASTClassLoader.forName(className).newInstance();
                    if (dataTypeInput.get().equals(thisDataType.getTypeDescription())) {
                        dataType = thisDataType;
                        break;
                    }
                    dataTypeDescList.add(thisDataType.getTypeDescription());
                } catch (ClassNotFoundException
                    | InstantiationException
                    | IllegalAccessException e) {
                }
            }
            if (dataType == null) {
                throw new IllegalArgumentException("Data type + '"
                        + dataTypeInput.get()
                        + "' cannot be found.  Choose one of "
                        + Arrays.toString(dataTypeDescList.toArray(new String[0])));
            }
        }
    }

}

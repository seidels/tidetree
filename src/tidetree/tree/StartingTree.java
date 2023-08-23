package tidetree.tree;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.util.List;
import java.util.stream.IntStream;

/**
 * @author Sophie Seidel
 */
@Description("This class generates a valid starting tree given an alignment and information on the editing experiment.")
public class StartingTree extends Tree implements StateNodeInitialiser {
    final public Input<Alignment> taxaInput = new Input<>("taxa", "set of taxa to initialise tree specified by alignment", Input.Validate.REQUIRED);

    final public Input<RealParameter> rootHeightInput = new Input<>("rootHeight", "Time from beginning of the experiment until sequencing", Input.Validate.REQUIRED);
    final public Input<RealParameter> editDurationInput = new Input<>("editDuration", "Time duration from edit start to edit stop", Input.Validate.REQUIRED);
    final public Input<RealParameter> editHeightInput = new Input<>("editHeight", "Time from the onset of edit until sequencing", Input.Validate.REQUIRED);
    final public Input<Boolean> sequencesAreClusteredInput = new Input<>("sequencesAreClustered", "Is true, if identical sequences appear in clusters within the alignment", true, Input.Validate.REQUIRED);
    final public Input<Integer> nClustersInput =  new Input<>("nClusters", "Number of clusters, where each cluster consists of identical sequences. " +
            "This input is only used when sequences are clustered.");

    // set up useful parameters
    int[][] matchMatrix; // matrix indicating if taxa have identical sequences
    int[] taxaInTree; // vector indicating if taxa have (1) or have not (0) already been included in the tree
    int nClusters; // cluster is a set of identical sequences
    double editHeight;
    double editDuration;
    double rootHeight;
    int nTaxa;
    Alignment taxa;
    boolean sequencesAreClustered;
    int iIntNode; //running index counting the number of internal nodes included in the tree

    @Override
    public void initAndValidate() {

        // init taxa
        taxa = taxaInput.get();
        List<String> taxanames = taxa.getTaxaNames();

        nTaxa = taxanames.size();
        iIntNode = nTaxa;

        rootHeight = rootHeightInput.get().getValue();
        editDuration = editDurationInput.get().getValue();
        editHeight = editHeightInput.get().getValue();

        if (editHeight > rootHeight){
            throw new RuntimeException("editHeight has to be smaller or equal than rootHeight");
        }

        // number of clusters is only known a priory if sequences are clustered in the alignment
        sequencesAreClustered = sequencesAreClusteredInput.get();
        if (sequencesAreClustered){
            nClusters = nClustersInput.get();
            if (nClusters < 1){
                throw new RuntimeException("Cluster size has to be positive");
            }
        }else{
            nClusters = nTaxa;
        }

        matchMatrix = set_match_matrix(taxa);
        taxaInTree = new int[nTaxa];

        initStateNodes();
        super.initAndValidate();
    }

    @Override
    public void initStateNodes() {

        root = get_tree(rootHeight, editHeight, taxa, nClusters, matchMatrix);

        leafNodeCount = nTaxa;
        nodeCount = leafNodeCount *2 -1 ;
        internalNodeCount = leafNodeCount -1;

        initArrays();

        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(m_initial.get());
    }

    /**
     * Sets up the match matrix. The matrix indicates whether two taxa have the same (1)
     * sequence or not (0).
     * @param taxa alignment with sequence info
     * @return match matrix
     */
    private int[][] set_match_matrix(Alignment taxa){
        matchMatrix = new int[nTaxa][nTaxa];
        List<String> taxanames = taxa.getTaxaNames();

        for (int i=0; i<nTaxa; i++){
            for (int j=i+1; j<nTaxa; j++){

                // get sequence
                String seq_i = taxa.getSequenceAsString(taxanames.get(i));
                String seq_j = taxa.getSequenceAsString(taxanames.get(j));

                matchMatrix[i][j] = (seq_i.equals(seq_j)) ? 1 : 0;
                matchMatrix[j][i] = matchMatrix[i][j];
            }
        }
        return matchMatrix;
    }

    // reset the internal node counter
    private void reset_iIntNode (){
        iIntNode = nTaxa;
    }
    private void reset_taxa_in_tree () {
        for (int i=0; i<taxaInTree.length; i++){
        taxaInTree[i] = 0;
    }}

    /**
     * Collects the taxa that have the same sequence as the input taxon.
     * @param nMatches number of matches of the input taxon.
     * @param iTaxon input taxon
     * @return vector with matching taxa.
     */
    private int[] get_matching_taxa (int nMatches, int iTaxon){
        //match collector
        int[] matches = new int[nMatches];

        // prevent unnecessary calculation
        if(nMatches == 0){
            return matches;
        }

        int iMatch = 0;
        // collect taxa with matching sequences
        if(sequencesAreClustered){
            // given that taxa with matching sequences are clustered together and only forward
            // search is done i.e. taxon with smallest id is used to find matching taxa
            // with larger id
            for(int i=(iTaxon+1); i<=(iTaxon+nMatches); i++){
                matches[iMatch] = i;
                iMatch++;
            }
        }else{
            //given that taxa with matching sequences are distributed throughout the match matrix
            for(int i=(iTaxon+1); i<matchMatrix.length; i++){
                if (matchMatrix[iTaxon][i] == 1){
                    matches[iMatch] = i;
                    iMatch++;
                }
            }
        }
        return matches;
    }


    /**
     * Builds subtree for identical sequences where the parent of the subtree is below the
     * edit stop. The internal nodes between the present and the edit stop
     * are placed at equal intervals (divTimes) and then scaled by a random number (0,1] to
     * to prevent 2 internal nodes in different cluster trees to be places at the same height.
     * @param iTaxon i-th taxon
     * @param editStop time after which no edit events are allowed to occur.
     * @param taxa alignment
     * @param iCluster i-th cluster
     * @return subtree of identical sequences
     */
    private Node get_cluster_tree(int iTaxon, double editStop, Alignment taxa, int iCluster){
        int nMatches = IntStream.of(matchMatrix[iTaxon]).sum();
        int[] matches = get_matching_taxa(nMatches, iTaxon);
        double divTimes;
        if(nMatches != 0){
            divTimes =  editStop / nMatches;
        }else{
            divTimes = editStop;
        }
        List<String> taxaNames = taxa.getTaxaNames();

        // generate left node
        Node nodeLeft = new Node();
        nodeLeft.setHeight(0);
        nodeLeft.setNr(iTaxon);
        nodeLeft.setID(taxaNames.get(iTaxon));
        nodeLeft.setMetaData("cluster", iCluster);
        nodeLeft.metaDataString = ("cluster=" + iCluster);

        //mark node as integrated into the tree
        taxaInTree[iTaxon] = 1;

        for (int iMatch=0; iMatch < nMatches; iMatch++){
            //set up right node
            Node nodeRight = new Node();
            nodeRight.setHeight(0.0);
            nodeRight.setID(taxaNames.get(matches[iMatch]));
            nodeRight.setNr(matches[iMatch]);
            nodeRight.setMetaData("cluster", iCluster);
            nodeRight.metaDataString = ("cluster=" + iCluster);
            taxaInTree[matches[iMatch]] = 1; //mark matched taxon as integrated in the tree

            //set up parent node
            Node parent = new Node();
            // space internal nodes at unequal intervals until editHeight
            double height = divTimes * (iMatch+1) - Randomizer.uniform(0.01, 0.09);
            parent.setHeight(height);// * Randomizer.nextDouble());
            parent.setNr(iIntNode);

            parent.addChild(nodeLeft);
            parent.addChild(nodeRight);
            parent.setMetaData("cluster", iCluster);
            parent.metaDataString = ("cluster=" + iCluster);
            iIntNode++; // increase internal node counter

            nodeLeft = parent;
        }
        return nodeLeft;
    }

    /**
     * Moves from the current taxon in the alignment to the next taxon that has not yet been included in the tree.
     *
     * @param iTaxon
     * @return next taxon for which a subtree will be built.
     */
    private int update_iTaxon(int iTaxon){
        int nMatches;
        if (sequencesAreClustered) {
            // for clustered sequences, move iTaxon to the next cluster start by adding the #(matching sequences)
            nMatches = IntStream.of(matchMatrix[iTaxon]).sum();
            iTaxon += (nMatches + 1);

        }else{
            // for unclustered sequences, increase iTaxon, until a sequence not yet present in the cluster trees is found.
            while(iTaxon < taxaInTree.length && taxaInTree[iTaxon] == 1){
                iTaxon++;
            }
        }
        return iTaxon;
    }

    /**
     * Build starting tree from sequences. The root height is given by the length of the experiment and fixed. The
     * edit height indicates when edit of the sequences starts. The edit duration marks the edit stop,
     * the time when no more scars are allowed to be acccumulated. (For more, see @General In this tree, all taxa with identical sequences are clustered in a
     * subtree (also called clusterTree) below the edit stop. This is done to ensure that the tree log likelihood
     * does not evaluate to -Infinity.
     *
     * The codes takes a starting taxon and generates a subtree of matching taxa (i.e. taxa with the same sequence) for
     * it. Then, the next taxon (i-th taxon) with a different sequence is taken and, again, a subtree of its matching
     * taxa is created. These subtrees are joined by a parent node above the edit stop. Then, the procedure
     * continues for the next taxon with a different sequence.
     *
     * Hence, this methods relies on finding taxa with identical sequences and the next taxon with a different sequence.
     *
     * Both tasks are easier if the taxa are clustered, i.e. ordered by sequence identity in the alignment. Then,
     * matching taxa come in a contiguous block just behind the i-th taxon in the alignment. Also, the next different
     * sequence is found at position i-th taxon + #matches of that taxon.
     *
     * If the taxa are not clustered, the matching sequences can be in any position of the alignment. To identify them,
     * we construct the matchMatrix. To keep track of which taxa have already been looked up in the matchMatrix and
     * included in the tree we use the variable taxaInTree.
     *
     * @param rootHeight
     * @param editHeight
     * @param taxa
     * @param nClusters
     * @param matchMatrix
     * @return starting tree
     */
    public Node get_tree(double rootHeight, double editHeight, Alignment taxa, int nClusters, int[][] matchMatrix){

        if(nTaxa == 1){
            int onlyTaxon = 0;
            int onlyCluster = 0;

            List<String> taxaNames = taxa.getTaxaNames();
            Node onlyLeaf = new Node();
            onlyLeaf.setHeight(0);
            onlyLeaf.setNr(onlyTaxon);
            onlyLeaf.setID(taxaNames.get(onlyTaxon));
            onlyLeaf.setMetaData("cluster", onlyCluster);
            onlyLeaf.metaDataString = ("cluster=" + onlyCluster);
            /*Node parent = new Node();
            parent.setHeight(rootHeight);
            parent.setNr(1);
            parent.addChild(onlyLeaf);*/

            return onlyLeaf;
        }

        //define taxon and cluster counters
        int iTaxon = 0;
        int iCluster = 0;

        // time between the internal nodes connecting the cluster trees
        double divTime = (rootHeight - editHeight) / nClusters;

        double editStop = editHeight - editDuration;
        // get left subtree
        Node subtreeLeft = get_cluster_tree(iTaxon, editStop, taxa, iCluster);

        // update iTaxon to next cluster start
        iCluster++;
        iTaxon = update_iTaxon(iTaxon);

        //for (int iCluster=1; iCluster < nClusters; iCluster++)
        while (iTaxon < nTaxa){

            // get right subtree
            Node subtreeRight = get_cluster_tree(iTaxon, editStop, taxa, iCluster);
            iCluster++;
            if(iCluster > nClusters){
                throw new RuntimeException("The number of clusters is larger than specified!");
            }
            iTaxon = update_iTaxon(iTaxon);

            // get parent
            Node parent = new Node();
            parent.setHeight(editHeight + iCluster * divTime);
            // Number of nodes in the final tree - number of nodes that have yet to be created
            parent.setNr(iIntNode);
            parent.addChild(subtreeLeft);
            parent.addChild(subtreeRight);
            parent.setMetaData("cluster", 0);
            parent.metaDataString = "cluster=0";
            iIntNode++;

            // setup as left tree to end recursion
            subtreeLeft = parent;

        }

        if (iTaxon > nTaxa){
            //throw new RuntimeException("iTaxon reached " + iTaxon + " but it should never exceed the number of sequences: " + nTaxa + "!");
        }
        // MRCA height for starting tree should be above edit window in case scars are present in alignment
        subtreeLeft.setHeight(rootHeight);
        return subtreeLeft;
    }
}

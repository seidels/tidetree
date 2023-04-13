package lineageTree.tree;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.evolution.datatype.EditData;
import beast.base.util.Randomizer;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class StartingTreeTest {

    private Alignment singleLeaf, identicalSequenceAlignment, nonIdenticalSequenceAlignment, nonIdenticalNonClusteredSequenceAlignment;
    private StartingTree identicalSequenceTree, nonIdenticalSequenceTree, tooFewClustersTree, singleLeafTree;
    private StartingTree identicalSequenceTreeClusteringFalse, nonIdenticalSequenceTreeClusteringFalse,
            nonIdenticalNonClusteredSequenceTreeClusteringFalse;


    @Before
    public void setup(){

        Sequence a = new Sequence("1", "0,0");
        Sequence b = new Sequence("2", "0,0");
        Sequence c = new Sequence("3", "1,0");
        Sequence d = new Sequence("4", "1,0");
        Sequence e = new Sequence("5", "0,1");

        EditData scarDat = new EditData();
        scarDat.initByName("nrOfStates", 3);

        singleLeaf = new Alignment();
        singleLeaf.initByName("sequence", a, "userDataType", scarDat, "stateCount", 3);

        identicalSequenceAlignment = new Alignment();
        identicalSequenceAlignment.initByName("sequence", a, "sequence", b, "userDataType", scarDat, "stateCount", 3);

        nonIdenticalSequenceAlignment = new Alignment();
        nonIdenticalSequenceAlignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "sequence", e,
                "userDataType", scarDat, "stateCount", 3);

        nonIdenticalNonClusteredSequenceAlignment = new Alignment();
        nonIdenticalNonClusteredSequenceAlignment.initByName("sequence", a, "sequence", c, "sequence", b, "sequence", d, "sequence", e,
                "userDataType", scarDat, "stateCount", 3);

        // add sampling times
        TaxonSet singleTaxon = new TaxonSet();
        TaxonSet taxaWithIdenticalSequences = new TaxonSet();
        TaxonSet taxaWithNonIdenticalSequences = new TaxonSet();
        TaxonSet taxaWithNonIdenticalNonClusteredSequences = new TaxonSet();

        singleTaxon.initByName("alignment", singleLeaf);
        taxaWithIdenticalSequences.initByName("alignment", identicalSequenceAlignment);
        taxaWithNonIdenticalSequences.initByName("alignment", nonIdenticalSequenceAlignment);
        taxaWithNonIdenticalNonClusteredSequences.initByName("alignment", nonIdenticalNonClusteredSequenceAlignment);
        //TraitSet dateTraits = new TraitSet();
        //dateTraitsIdenticalSequences.initByName("taxa", taxaWithIdenticalSequences, "traitname", "date-forward", "value", "1=32");


        singleLeafTree = new StartingTree();
        singleLeafTree.initByName("rootHeight", 32.0, "editHeight", 25.0, "editDuration", 2.0,
                "sequencesAreClustered", true, "nClusters", 1, "taxa", singleLeaf);


        Randomizer.setSeed(6);
        identicalSequenceTree = new StartingTree();
        identicalSequenceTree.initByName("rootHeight", 32.0, "editHeight", 25.0, "editDuration", 2.0,
                "sequencesAreClustered", true, "nClusters", 1, "taxa", identicalSequenceAlignment);

        Randomizer.setSeed(6);
        identicalSequenceTreeClusteringFalse = new StartingTree();
        identicalSequenceTreeClusteringFalse.initByName("rootHeight", 32.0, "editHeight", 25.0, "editDuration", 2.0,
                "sequencesAreClustered", false, "taxa", identicalSequenceAlignment);

        Randomizer.setSeed(6);
        nonIdenticalSequenceTree = new StartingTree();
        nonIdenticalSequenceTree.initByName("rootHeight", 32.0, "editHeight", 25.0, "editDuration", 2.0,
                "sequencesAreClustered", true, "nClusters", 3, "taxa", nonIdenticalSequenceAlignment);



        tooFewClustersTree = new StartingTree();


        Randomizer.setSeed(6);
        nonIdenticalSequenceTreeClusteringFalse = new StartingTree();
        nonIdenticalSequenceTreeClusteringFalse.initByName("rootHeight", 32.0, "editHeight", 25.0, "editDuration", 2.0,
                "sequencesAreClustered", false, "taxa", nonIdenticalSequenceAlignment);

        Randomizer.setSeed(6);
        nonIdenticalNonClusteredSequenceTreeClusteringFalse = new StartingTree();
        nonIdenticalNonClusteredSequenceTreeClusteringFalse.initByName("rootHeight", 32.0, "editHeight", 25.0, "editDuration", 2.0,
                "sequencesAreClustered", false, "taxa", nonIdenticalNonClusteredSequenceAlignment);

    }

    @Test
    public void testSingleTaxonTree(){
        String trueTree = "0[&cluster=0]:0.0";
        assertEquals("singleLeafTree", singleLeafTree.toString(), trueTree);

    }

    @Test
   public void testTreeFromIdenticalSequences(){
        String trueTree = "(0[&cluster=0]:32.0,1[&cluster=0]:32.0)2[&cluster=0]:0.0";
        assertEquals("Correct starting tree from identical sequences, clustering=true",
                trueTree, identicalSequenceTree.toString());
    }

    @Test
    public void testTreeFromIdenticalSequencesClusteringFalse(){
        String trueTree = "(0[&cluster=0]:32.0,1[&cluster=0]:32.0)2[&cluster=0]:0.0";
        assertEquals("Correct starting tree from identical sequences, clustering=false",
                trueTree, identicalSequenceTreeClusteringFalse.toString());
    }

    @Test
    public void testTreeFromNonIdenticalSequence(){
        String trueTree = "(((0[&cluster=0]:22.937122117900326,1[&cluster=0]:22.937122117900326)5[&cluster=0]:6.7295445487663414,(2[&cluster=1]:22.91972105121713,3[&cluster=1]:22.91972105121713)6[&cluster=1]:6.746945615449537)7[&cluster=0]:2.333333333333332,4[&cluster=2]:32.0)8[&cluster=0]:0.0";
        assertEquals("Correct starting tree from identical sequences, clustering=true",
                trueTree, nonIdenticalSequenceTree.toString());
    }

    @Test
    public void testTreeFromNonIdenticalSequenceClusteringFalse(){
        String trueTree = "(((0[&cluster=0]:22.937122117900326,1[&cluster=0]:22.937122117900326)5[&cluster=0]:4.862877882099674,(2[&cluster=1]:22.91972105121713,3[&cluster=1]:22.91972105121713)6[&cluster=1]:4.88027894878287)7[&cluster=0]:4.199999999999999,4[&cluster=2]:32.0)8[&cluster=0]:0.0";
        assertEquals("Correct starting tree from identical sequences, clustering=false",
                trueTree, nonIdenticalSequenceTreeClusteringFalse.toString());
    }

    @Test
    public void testTreeWithWrongClusterInput() {
        String trueTree = "(((0[&cluster=0]:22.937122117900326,1[&cluster=0]:22.937122117900326)5[&cluster=0]:4.862877882099674,(2[&cluster=1]:22.91972105121713,3[&cluster=1]:22.91972105121713)6[&cluster=1]:4.88027894878287)7[&cluster=0]:4.199999999999999,4[&cluster=2]:32.0)8[&cluster=0]:0.0";
        try{
            tooFewClustersTree.initByName("rootHeight", 32.0, "editHeight", 25.0, "editDuration", 2.0,
                    "sequencesAreClustered", true, "nClusters", 2, "taxa", nonIdenticalSequenceAlignment);
        }catch (RuntimeException exception){
            assertEquals("initAndValidate() failed! The number of clusters is larger than specified!", exception.getMessage());
        };
    }

    @Test
    public void testTreeFromNonIdenticalNonClusteredSequenceClusteringFalse() {
        String trueTree = "(((0[&cluster=0]:22.937122117900326,2[&cluster=0]:22.937122117900326)5[&cluster=0]:4.862877882099674,(1[&cluster=1]:22.91972105121713,3[&cluster=1]:22.91972105121713)6[&cluster=1]:4.88027894878287)7[&cluster=0]:4.199999999999999,4[&cluster=2]:32.0)8[&cluster=0]:0.0";
        assertEquals("Correct starting tree from identical sequences, clustering=false",
                trueTree, nonIdenticalNonClusteredSequenceTreeClusteringFalse.toString());
    }



}

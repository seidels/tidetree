package lineageTree.simulation;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Node;
import lineageTree.simulation.SimulatedAlignment;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import lineageTree.substitutionmodel.GeneralScarringLoss;
import beast.evolution.datatype.ScarData;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

public class SimulatedAlignmentTest{

    @Test
    public void testBranchType00 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "scarData", "stateCount", 4);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        GeneralScarringLoss scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);


        SimulatedAlignment simAlignment = new SimulatedAlignment();
        int seqLength = 1;
        simAlignment.initByName("tree", tree1,
                "siteModel", siteM, "sequenceLength", seqLength,
                "branchRateModel", branchModel

        );

        int[] helperNode = simAlignment.helperNodeSeq;

        assertEquals(helperNode[0], -2);
    }



    @Test
    /*
    Test branch with parent above and child within scarring window
    Adjust the scarring model to do so.
     */
    public void testBranchType10 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "scarData", "stateCount", 4);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter scarringHeight = new RealParameter("3.0");
        RealParameter scarringDuration = new RealParameter("3.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        GeneralScarringLoss scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // clock model
        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);


        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree1, "siteModel", siteM, "sequenceLength", 1, "branchRateModel", branchModel);

        int[] helperNode = simAlignment.helperNodeSeq;
        assertTrue(simAlignment.branchType[0]);
        assertFalse(simAlignment.branchType[1]);

        //tested transition probabilities at runtime - works
    }

    @Test
    public void testBranchType01 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "scarData", "stateCount", 4);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter scarringHeight = new RealParameter("5.0");
        RealParameter scarringDuration = new RealParameter("3.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        GeneralScarringLoss scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // clock model
        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);

        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree1, "siteModel", siteM, "sequenceLength", 1, "branchRateModel", branchModel);

        assertFalse(simAlignment.branchType[0]);
        assertTrue(simAlignment.branchType[1]);

        //tested transition probabilities at runtime - works
    }

    @Test
    public void testBranchType11 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "scarData", "stateCount", 3);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        ScarData scarDat = new ScarData();
        scarDat.initByName("nrOfStates", 3);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.0");
        RealParameter scarRates = new RealParameter("10.0");
        RealParameter scarringHeight = new RealParameter("4.0");
        RealParameter scarringDuration = new RealParameter("3.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        GeneralScarringLoss scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // clock model
        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);

        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree1, "siteModel", siteM, "sequenceLength", 1, "userDataType", scarDat, "branchRateModel", branchModel);

        int[] helperNode = simAlignment.helperNodeSeq;
        assertTrue(simAlignment.branchType[0]);
        assertTrue(simAlignment.branchType[1]);
        //tested transition probabilities at runtime - works
    }
}
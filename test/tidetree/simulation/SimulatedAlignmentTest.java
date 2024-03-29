package tidetree.simulation;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import tidetree.substitutionmodel.EditAndSilencingModel;
import tidetree.evolution.datatype.EditData;
import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

public class SimulatedAlignmentTest{

    @Test
    public void testBranchType00 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        EditData editData = new EditData();
        editData.initByName("nrOfStates", 4);
        alignment.initByName("sequence", a, "userDataType", editData);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter silencingRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter editHeight = new RealParameter("25.0");
        RealParameter editDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel editAndSilencingModel = new EditAndSilencingModel();
        editAndSilencingModel.initByName("editRates", scarRates,
                "silencingRate", silencingRate,
                "editHeight", editHeight,
                "editDuration", editDuration, "frequencies", frequencies);

        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", editAndSilencingModel);


        SimulatedAlignment simAlignment = new SimulatedAlignment();
        int seqLength = 1;
        simAlignment.initByName("tree", tree1,
                "siteModel", siteM, "sequenceLength", seqLength,
                "branchRateModel", branchModel, "userDataType", editData

        );

        int[] helperNode = simAlignment.helperNodeSeq;

        assertEquals(helperNode[0], -2);
    }



    @Test
    /*
    Test branch with parent above and child within edit window
    Adjust the edit model to do so.
     */
    public void testBranchType10 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        EditData editData = new EditData();
        editData.initByName("nrOfStates", 4);
        alignment.initByName("sequence", a, "userDataType", editData, "stateCount", 4);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init edit model
        RealParameter silencingRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter editHeight = new RealParameter("3.0");
        RealParameter editDuration = new RealParameter("3.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel editAndSilencingModel = new EditAndSilencingModel();
        editAndSilencingModel.initByName("editRates", scarRates,
                "silencingRate", silencingRate,
                "editHeight", editHeight,
                "editDuration", editDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", editAndSilencingModel);

        // clock model
        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);


        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree1, "siteModel", siteM, "sequenceLength", 1,
                "branchRateModel", branchModel, "userDataType", editData);

        assertTrue(simAlignment.branchType[0]);
        assertFalse(simAlignment.branchType[1]);

        //tested transition probabilities at runtime - works
    }

    @Test
    public void testBranchType01 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        EditData editData = new EditData();
        editData.initByName("nrOfStates", 4);
        alignment.initByName("sequence", a, "userDataType", editData, "stateCount", 4);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init edit model
        RealParameter silencingRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter editHeight = new RealParameter("5.0");
        RealParameter editDuration = new RealParameter("3.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel editAndSilencingModel = new EditAndSilencingModel();
        editAndSilencingModel.initByName("editRates", scarRates,
                "silencingRate", silencingRate,
                "editHeight", editHeight,
                "editDuration", editDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", editAndSilencingModel);

        // clock model
        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);

        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree1, "siteModel", siteM,
                "sequenceLength", 1, "branchRateModel", branchModel, "userDataType", editData);

        assertFalse(simAlignment.branchType[0]);
        assertTrue(simAlignment.branchType[1]);

        //tested transition probabilities at runtime - works
    }

    @Test
    public void testBranchType11 (){

        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer", "stateCount", 3);

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:5.0)1[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        EditData scarDat = new EditData();
        scarDat.initByName("nrOfStates", 3);

        //init edit model
        RealParameter silencingRate = new RealParameter("0.0");
        RealParameter scarRates = new RealParameter("10.0");
        RealParameter editHeight = new RealParameter("4.0");
        RealParameter editDuration = new RealParameter("3.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel editAndSilencingModel = new EditAndSilencingModel();
        editAndSilencingModel.initByName("editRates", scarRates,
                "silencingRate", silencingRate,
                "editHeight", editHeight,
                "editDuration", editDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", editAndSilencingModel);

        // clock model
        RealParameter clock = new RealParameter("1.0");
        StrictClockModel branchModel = new StrictClockModel();
        branchModel.initByName("clock.rate", clock);

        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree1, "siteModel", siteM, "sequenceLength", 1,
                "userDataType", scarDat, "branchRateModel", branchModel);

        int[] helperNode = simAlignment.helperNodeSeq;
        assertTrue(simAlignment.branchType[0]);
        assertTrue(simAlignment.branchType[1]);
        //tested transition probabilities at runtime - works
    }
}
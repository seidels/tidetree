package tidetree.distributions;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import tidetree.substitutionmodel.EditAndSilencingModel;
import org.junit.Test;
import tidetree.util.AlignmentFromNexus;
import tidetree.util.NexusParser;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class TreeLikelihoodWithEditWindowTest {

    TreeLikelihoodWithEditWindow likelihood1, likelihood2, likelihood3;
    Tree tree1, tree2, tree3, treeImpossible;
    private TreeLikelihoodWithEditWindow likelihood4, likelihoodNegInf;

    @Test
    public void testSingleTaxon() {
        // init alignment
        Sequence a = new Sequence("0", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer", "stateCount", 4);

        //init tree
        tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "0[&cluster=0]:0.0",
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

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihood1 = new TreeLikelihoodWithEditWindow();
        likelihood1.initByName("data", alignment, "tree", tree1,
                "siteModel", siteM, "branchRateModel", clockModel);

        double logP = likelihood1.calculateLogP();
        assertEquals(-10.0, logP, 1e-13);
    }

    @Test
    public void testLikelihood1() {
        // parent above scarring window, children below scarring window

        // init alignment
        Sequence a = new Sequence("0", "0,");
        Sequence b = new Sequence("1", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "integer", "stateCount", 4);

        //init tree
        tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:28.5,1[&cluster=1]:28.5)2[&cluster=0]:0.0",
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

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihood1 = new TreeLikelihoodWithEditWindow();
        likelihood1.initByName("data", alignment, "tree", tree1,
                "siteModel", siteM, "branchRateModel", clockModel);

        Node parent = tree1.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);

        // test correct likelihood calculation for branches that need helper nodes ( that is branches, that
        // cross points where the transition rate matrix changes)

        /*double[] partials = likelihood1.calculatePartialsBeforeParent(parent, child1, 0, 0,
                new double[]{0, 23.0, 25.0}, 1.0, 2);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.0001234098040866, 0, 0, 0}, 1e-15);

        partials = likelihood1.calculatePartialsForCrossBranches(partials, parent, child1, child2,
                true, false);
        assertArrayEquals("Assert correct likelihood at parent", partials,
                new double[]{3.75566676593833e-9, 0, 0, 0}, 1e-15);
*/
        double logP = likelihood1.calculateLogP();
        assertEquals(Math.log(3.75566676593833e-9), logP, 1e-13);
    }

    @Test
    public void testLikelihood1b() {
        // parent above scarring window, children below scarring window, gammacategoryCount

        // init alignment
        Sequence a = new Sequence("0", "0,");
        Sequence b = new Sequence("1", "0,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "integer", "stateCount", 4);

        //init tree
        tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:28.5,1[&cluster=1]:28.5)2[&cluster=0]:0.0",
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

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        RealParameter shape = new RealParameter("1.0");
        SiteModel siteM = new SiteModel();
        siteM.initByName( "substModel", scarringModel, "gammaCategoryCount", 2,
                "shape", shape);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihood1 = new TreeLikelihoodWithEditWindow();
        likelihood1.initByName("data", alignment, "tree", tree1,
                "siteModel", siteM, "branchRateModel", clockModel);

        Node parent = tree1.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);

        // test correct likelihood calculation for branches that need helper nodes ( that is branches, that
        // cross points where the transition rate matrix changes)

        /*double[] partials = likelihood1.calculatePartialsBeforeParent(parent, child1, 0, 0,
                new double[]{0, 23.0, 25.0}, 1.0, 2);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.0001234098040866, 0, 0, 0}, 1e-15);

        partials = likelihood1.calculatePartialsForCrossBranches(partials, parent, child1, child2,
                true, false);
        assertArrayEquals("Assert correct likelihood at parent", partials,
                new double[]{3.75566676593833e-9, 0, 0, 0}, 1e-15);
*/
        double logP = likelihood1.calculateLogP();
        assertEquals(Math.log(0.0006354730093870854), logP, 1e-13);
    }

    @Test
    public void testLikelihood2(){
        // parent within scarring window, children below scarring window

        // init alignment
        Sequence a = new Sequence("0", "0,");
        Sequence b2 = new Sequence("1", "1,");

        Alignment alignment2 = new Alignment();
        alignment2.initByName("sequence", a, "sequence", b2, "dataType", "integer", "stateCount", 3);

        tree2 = new TreeParser();
        tree2.initByName("IsLabelledNewick", true, "taxa", alignment2, "newick",
                "(0[&cluster=0]:25,1[&cluster=1]:25)2[&cluster=0]:0.0",
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

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);
        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();


        likelihood2 = new TreeLikelihoodWithEditWindow();
        likelihood2.initByName("data", alignment2, "tree", tree2,
                "siteModel", siteM, "branchRateModel", clockModel);

        Node parent = tree2.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);

        /*double[] partials = likelihood2.calculatePartialsBeforeParent(parent, child1, 0, 0,
                new double[]{0, 23.0, Double.NEGATIVE_INFINITY}, 1.0, 1);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.01005183574463, 0, 0, 0}, 1e-12);

        partials= likelihood2.calculatePartialsForCrossBranches(partials, parent, child1, child2,
                false, true);
        assertArrayEquals("Assert correct likelihood at parent", partials,
                new double[]{4.081493696794465e-7, 0, 0, 0}, 1e-15);
*/
        double logP = likelihood2.calculateLogP();
        assertEquals(Math.log(4.081493696794465e-7), logP, 1e-13);
    }

    @Test
     public void testLikelihood3(){
        // parent above scarring window, children within scarring window

        // init alignment
        Sequence a = new Sequence("0", "0,");
        Sequence b3 = new Sequence("1", "3,");

        Alignment alignment3 = new Alignment();
        alignment3.initByName("sequence", a, "sequence", b3, "dataType", "integer", "stateCount", 3);

        //init tree
        tree3 = new TreeParser();
        tree3.initByName("IsLabelledNewick", true, "taxa", alignment3, "newick",
                "(0[&cluster=0]:25,1[&cluster=1]:25)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter scarringHeight = new RealParameter("2.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel scarringModel3 = new EditAndSilencingModel();
        scarringModel3.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM3 = new SiteModel();
        siteM3.initByName( "gammaCategoryCount", 0, "substModel", scarringModel3);
        StrictClockModel clockModel = new StrictClockModel();

        likelihood3 = new TreeLikelihoodWithEditWindow();
        likelihood3.initByName("data", alignment3, "tree", tree3,
                "siteModel", siteM3, "branchRateModel", clockModel);


        Node parent = tree3.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);

//        double[] partials = likelihood3.calculatePartialsBeforeParent(parent, child2, 1, 1,
//                new double[]{0, 2.0, Double.NEGATIVE_INFINITY}, 1.0, 1);
//        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
//                new double[]{0.329679953964361, 0.329679953964361, 0.329679953964361, 1.0}, 1e-15);
//
//        partials = likelihood3.calculatePartialsBeforeParent(parent, child1, 0, 0,
//                new double[]{0, 2.0, Double.NEGATIVE_INFINITY}, 1.0, 1);
//        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
//                new double[]{0.012277339903068, 0, 0, 0}, 1e-15);
//
//        partials= likelihood3.calculatePartialsForCrossBranches(partials, parent, child1, child2,
//                true, false);
//        assertArrayEquals("Assert correct likelihood at parent", partials,
//                new double[]{1.225782753675767e-04, 0, 0, 0}, 1e-15);

        double logP = likelihood3.calculateLogP();
        assertEquals(Math.log(1.225782753675767e-04), logP, 1e-14);
    }

    @Test
    public void testLikelihood3b(){
        // parent above scarring window, children within scarring window, 2 sites!

        // init alignment
        Sequence a = new Sequence("0", "0,0");
        Sequence b3 = new Sequence("1", "3,3");

        Alignment alignment3 = new Alignment();
        alignment3.initByName("sequence", a, "sequence", b3, "dataType", "integer", "stateCount", 4);

        //init tree
        tree3 = new TreeParser();
        tree3.initByName("IsLabelledNewick", true, "taxa", alignment3, "newick",
                "(0[&cluster=0]:25,1[&cluster=1]:25)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");
        RealParameter scarringHeight = new RealParameter("2.0");
        RealParameter scarringDuration = new RealParameter("2.0");
        //RealParameter scarRates = new RealParameter("0.0 0.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel scarringModel3 = new EditAndSilencingModel();
        scarringModel3.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM3 = new SiteModel();
        siteM3.initByName( "gammaCategoryCount", 0, "substModel", scarringModel3);
        StrictClockModel clockModel = new StrictClockModel();

        likelihood3 = new TreeLikelihoodWithEditWindow();
        likelihood3.initByName("data", alignment3, "tree", tree3,
                "siteModel", siteM3, "branchRateModel", clockModel);


        /*double[] partials = likelihood3.calculatePartialsBeforeParent(parent, child2, 1, 1,
                new double[]{0, 2.0, Double.NEGATIVE_INFINITY}, 1.0, 1);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.329679953964361, 0.329679953964361, 0.329679953964361, 1.0}, 1e-15);

        partials = likelihood3.calculatePartialsBeforeParent(parent, child1, 0, 0,
                new double[]{0, 2.0, Double.NEGATIVE_INFINITY}, 1.0, 1);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.012277339903068, 0, 0, 0}, 1e-15);

        partials= likelihood3.calculatePartialsForCrossBranches(partials, parent, child1, child2,
                true, false);
        assertArrayEquals("Assert correct likelihood at parent", partials,
                new double[]{1.225782753675767e-04, 0, 0, 0}, 1e-15);
*/
        double logP = likelihood3.calculateLogP();
        assertEquals(Math.log(Math.pow(1.225782753675767e-04,2)), logP, 1e-13);
    }

    @Test
    public void testLikelihood4(){
        //test likelihood calculation, when scarring window extends over the entire tree height
        // -> calculation does not require a changing rate matrix nor helper nodes

        // init alignment
        Sequence a = new Sequence("0", "0,");
        Sequence b2 = new Sequence("1", "1,");

        Alignment alignment2 = new Alignment();
        alignment2.initByName("sequence", a, "sequence", b2, "dataType", "integer", "stateCount", 3);

        //init tree
        tree2 = new TreeParser();
        tree2.initByName("IsLabelledNewick", true, "taxa", alignment2, "newick",
                "(0[&cluster=0]:25,1[&cluster=1]:25)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);


        // init scarring model
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        RealParameter scarringHeight = new RealParameter("100.0");
        RealParameter scarringDuration = new RealParameter("100.0");
        EditAndSilencingModel scarringModel4 = new EditAndSilencingModel();
        scarringModel4.initByName("editRates", new RealParameter("0.01 0.01"),
                "silencingRate", new RealParameter("0.01"),
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        SiteModel siteM4 = new SiteModel();
        siteM4.initByName("gammaCategoryCount", 0, "substModel", scarringModel4);
        StrictClockModel clockModel = new StrictClockModel();

        //init likelihood
        likelihood4 = new TreeLikelihoodWithEditWindow();
        likelihood4.initByName("data", alignment2, "tree", tree2,
                "siteModel", siteM4, "branchRateModel", clockModel);


        double logP = likelihood4.calculateLogP();
        assertEquals(Math.log(0.072374640511506), logP, 1e-14);

    }



    @Test
    public void testLikelihood4b(){
        //test likelihood calculation, when scarring window is of size 0
        // -> any tree should have 0 likelihood if it had scars

        // init alignment
        Sequence a = new Sequence("0", "1,");
        Sequence b2 = new Sequence("1", "2,");

        Alignment alignment2 = new Alignment();
        alignment2.initByName("sequence", a, "sequence", b2, "dataType", "integer", "stateCount", 3);

        //init tree
        tree2 = new TreeParser();
        tree2.initByName("IsLabelledNewick", true, "taxa", alignment2, "newick",
                "(0[&cluster=0]:25,1[&cluster=1]:25)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);


        // init scarring model
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("0.0");

        EditAndSilencingModel scarringModel4 = new EditAndSilencingModel();
        scarringModel4.initByName("editRates", new RealParameter("0.01 0.01"),
                "silencingRate", new RealParameter("0.01"),
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        SiteModel siteM4 = new SiteModel();
        siteM4.initByName("gammaCategoryCount", 0, "substModel", scarringModel4);
        StrictClockModel clockModel = new StrictClockModel();

        //init likelihood
        likelihood4 = new TreeLikelihoodWithEditWindow();
        likelihood4.initByName("data", alignment2, "tree", tree2,
                "siteModel", siteM4, "branchRateModel", clockModel);


        double logP = likelihood4.calculateLogP();
        assertEquals("Test 0 likelihood for tree with scars and scarring model without scarring window",
                Double.NEGATIVE_INFINITY, logP, 1e-14);

        //-------------------------------------------------------------------------------------------------------//

        scarringHeight.setValue(30.0);
        scarringDuration.setValue(2.0);
        EditAndSilencingModel scarringModel4b = new EditAndSilencingModel();
        scarringModel4b.initByName("editRates", new RealParameter("0.01 0.01"),
                "silencingRate", new RealParameter("0.01"),
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        siteM4.initByName("gammaCategoryCount", 0, "substModel", scarringModel4b);

        //init likelihood
        likelihood4.initByName("data", alignment2, "tree", tree2,
                "siteModel", siteM4, "branchRateModel", clockModel);


        logP = likelihood4.calculateLogP();
        assertEquals("Test 0 likelihood if scars segregate at parent and scarring window was before parent",
                Double.NEGATIVE_INFINITY, logP, 1e-14);

    }


    @Test
    public void testLikelihood5(){
        //test that impossible trees under the scarring model get a neg inf likelihood

        // init alignment
        Sequence a = new Sequence("0", "0,");
        Sequence b2 = new Sequence("1", "1,");
        Sequence c = new Sequence("2", "2,");

        Alignment alignment4 = new Alignment();
        alignment4.initByName("sequence", a, "sequence", b2, "sequence", c, "dataType", "integer", "stateCount", 3);


        //init trees
        treeImpossible = new TreeParser();
        treeImpossible.initByName("IsLabelledNewick", true, "taxa", alignment4, "newick",
                "((2[&cluster=0]:10,1[&cluster=1]:10)3[&cluster=0]:15.0,0[&cluster=2]:25)4[&cluster=0]:0.0",
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

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihoodNegInf = new TreeLikelihoodWithEditWindow();
        likelihoodNegInf.initByName("data", alignment4, "tree", treeImpossible,
                "siteModel", siteM, "branchRateModel", clockModel);

        double logP = likelihoodNegInf.calculateLogP();

        assertEquals(Double.NEGATIVE_INFINITY, logP, 1e-14);
    }

    @Test
    public void testLikelihood6(){
        //test tree that was found in tree posterior but should have been rejected

        // init alignment
        Sequence a = new Sequence("0", "1,");
        Sequence b = new Sequence("1", "1,");
        Sequence c = new Sequence("2", "2,");
        Sequence d = new Sequence("3", "2,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "dataType", "integer", "stateCount", 4);

        String wrongAcceptedTree = "((0[&cluster=0]:18.987792267733195,2[&cluster=1]:18.987792267733195)5[&cluster=1]:13.012207732266805,(1[&cluster=0]:30.196436618358874,3[&cluster=1]:30.196436618358874)4[&cluster=0]:1.8035633816411263)6[&cluster=0]:0.0;";
        //init trees
        treeImpossible = new TreeParser();
        treeImpossible.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                wrongAcceptedTree,
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.02");
        RealParameter scarRates = new RealParameter("1 1");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihoodNegInf = new TreeLikelihoodWithEditWindow();
        likelihoodNegInf.initByName("data", alignment, "tree", treeImpossible,
                "siteModel", siteM, "branchRateModel", clockModel);

        double logP = likelihoodNegInf.calculateLogP();

        assertEquals(Double.NEGATIVE_INFINITY, logP, 1e-14);
    }

    @Test
    public void testLikelihood7(){
        //test tree that contains all "branch classes"

        // init alignment
        // note: the second site is just reversing 1 and 2. Since the scar types have the same
        // scarring rate here, their partials will be the same - just the position in the partial
        // array will be swapped (1<->2). This also means that the tree likelihood is the tree
        // likelihood of one site x 2. To compare to the matlab script, just use the first site.
        Sequence a = new Sequence("0", "0,0");
        Sequence b = new Sequence("1", "1,2");
        Sequence c = new Sequence("2", "1,2");
        Sequence d = new Sequence("3", "2,1");
        Sequence e = new Sequence("4", "3,3");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "sequence", e, "dataType", "integer", "stateCount", 4);

        String newickTree = "((0[&cluster=0]:26.0,(1[&cluster=1]:13.0,2[&cluster=1]:13.0)5[&cluster=1]:13.0)7[&cluster=0]:6.0,(3[&cluster=2]:24.0,4[&cluster=3]:24.0)6[&cluster=1]:8.0)8:0.0;";
        //init trees
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newickTree,
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.02");
        RealParameter scarRates = new RealParameter("1 1");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihoodNegInf = new TreeLikelihoodWithEditWindow();
        likelihoodNegInf.initByName("data", alignment, "tree", tree,
                "siteModel", siteM, "branchRateModel", clockModel, "useScaling", true);

        Node parent = tree.getRoot();
        Node node_7 = parent.getChild(0);
        Node node_6 = parent.getChild(1);
        Node node_5 = node_7.getChild(1);

        //test partials at node 5
        likelihoodNegInf.traverse(tree.getRoot());
        double [] partials = new double[8];
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_5.getNr(), partials);
       assertArrayEquals("Assert correct likelihood at internal node 5:", partials,
                new double[]{0, 0.594520547970194, 0, 0, 0,0, 0.594520547970194,0}, 1e-15);

       //test partials at node 6
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_6.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 6:", partials,
                new double[]{0.101983098705779, 0, 0.235890505831029, 0, 0.101983098705779, 0.235890505831029 ,0,0}, 1e-15);

        //test partials at node 7
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_7.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 7:", partials,
                new double[]{0.002450084837716, 0, 0, 0, 0.002450084837716, 0, 0, 0}, 1e-15);


        double logP = likelihoodNegInf.calculateLogP();

        assertEquals(-8.447652794733264*2, logP, 1e-14);
    }
    @Test
    public void testLikelihood7b(){
        //test tree that contains all "branch classes"

        // init alignment
        // note: the second site is just reversing 1 and 2. Since the scar types have the same
        // scarring rate here, their partials will be the same - just the position in the partial
        // array will be swapped (1<->2). This also means that the tree likelihood is the tree
        // likelihood of one site x 2. To compare to the matlab script, just use the first site.
        Sequence a = new Sequence("0", "0,0");
        Sequence b = new Sequence("1", "1,2");
        Sequence c = new Sequence("2", "1,2");
        Sequence d = new Sequence("3", "2,1");
        Sequence e = new Sequence("4", "3,3");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "sequence", e, "dataType", "integer", "stateCount", 4);

        String newickTree = "((0[&cluster=0]:26.0,(1[&cluster=1]:13.0,2[&cluster=1]:13.0)5[&cluster=1]:13.0)7[&cluster=0]:6.0,(3[&cluster=2]:24.0,4[&cluster=3]:24.0)6[&cluster=1]:8.0)8:0.0;";
        //init trees
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newickTree,
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.02");
        RealParameter scarRates = new RealParameter("1 1");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        RealParameter shape = new RealParameter("1.0");
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 2, "substModel", scarringModel,
                "shape", shape);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihoodNegInf = new TreeLikelihoodWithEditWindow();
        likelihoodNegInf.initByName("data", alignment, "tree", tree,
                "siteModel", siteM, "branchRateModel", clockModel);

        Node parent = tree.getRoot();
        Node node_7 = parent.getChild(0);
        Node node_6 = parent.getChild(1);
        Node node_5 = node_7.getChild(1);

        //test partials at node 5
        likelihoodNegInf.traverse(tree.getRoot());
        double [] partials = new double[16];
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_5.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 5:", partials,
                new double[]{0, 0.836331904590560, 0, 0, 0, 0, 0.836331904590560, 0, 0, 0.422624893321294, 0, 0, 0,0, 0.422624893321294,0}, 1e-15);

        //test partials at node 6
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_6.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 6:", partials,
                new double[]{0.032054626428896, 0, 0.128958931608947, 0, 0.032054626428896, 0.128958931608947 ,0,0,
                        0.119317341270894, 0, 0.247654804887413, 0, 0.119317341270894, 0.247654804887413 ,0,0}, 1e-15);

        //test partials at node 7
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_7.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 7:", partials,
                new double[]{0.060425520687175, 0, 0, 0, 0.060425520687175, 0, 0, 0,
                        0.769152619166653e-4, 0, 0, 0, 0.769152619166653e-4, 0, 0, 0}, 1e-15);

        //test partials at root node
        likelihoodNegInf.getLikelihoodCore().getNodePartials(parent.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at root:", partials,
                new double[]{0.002643849069724, 0, 0, 0, 0.002643849069724, 0, 0, 0,
                        0.598196646044564e-5, 0, 0, 0, 0.598196646044564e-5, 0, 0, 0}, 1e-15);


        double logP = likelihoodNegInf.calculateLogP();

        assertEquals(-13.252813163014444, logP, 1e-14);
    }

    @Test
    public void testOrigin(){
        //test tree that contains all "branch classes" and calculates likelihood until origin

        // init alignment
        // note: the second site is just reversing 1 and 2. Since the scar types have the same
        // scarring rate here, their partials will be the same - just the position in the partial
        // array will be swapped (1<->2). This also means that the tree likelihood is the tree
        // likelihood of one site x 2. To compare to the matlab script, just use the first site.
        Sequence a = new Sequence("0", "0,0");
        Sequence b = new Sequence("1", "1,2");
        Sequence c = new Sequence("2", "1,2");
        Sequence d = new Sequence("3", "2,1");
        Sequence e = new Sequence("4", "3,3");

        RealParameter origin = new RealParameter("42");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "sequence", e, "dataType", "integer", "stateCount", 4);

        String newickTree = "((0[&cluster=0]:26.0,(1[&cluster=1]:13.0,2[&cluster=1]:13.0)5[&cluster=1]:13.0)7[&cluster=0]:6.0,(3[&cluster=2]:24.0,4[&cluster=3]:24.0)6[&cluster=1]:8.0)8:0.0;";
        //init trees
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newickTree,
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.02");
        RealParameter scarRates = new RealParameter("1 1");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("2.0");

        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        RealParameter shape = new RealParameter("1.0");
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 2, "substModel", scarringModel,
                "shape", shape);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihoodNegInf = new TreeLikelihoodWithEditWindow();
        likelihoodNegInf.initByName("data", alignment, "tree", tree,
                "siteModel", siteM, "branchRateModel", clockModel, "origin", origin, "useScaling", true);

        Node parent = tree.getRoot();
        Node node_7 = parent.getChild(0);
        Node node_6 = parent.getChild(1);
        Node node_5 = node_7.getChild(1);

        //test partials at node 5
        double logP = likelihoodNegInf.calculateLogP();
        logP = likelihoodNegInf.calculateLogP();
        assertEquals("1st likelihood calculation", -13.391340286074429, logP, 1e-14);

        double [] partials = new double[16];
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_5.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 5:", partials,
                new double[]{0, 0.836331904590560, 0, 0, 0, 0, 0.836331904590560, 0, 0, 0.422624893321294, 0, 0, 0,0, 0.422624893321294,0}, 1e-15);

        //test partials at node 6
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_6.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 6:", partials,
                new double[]{0.032054626428896, 0, 0.128958931608947, 0, 0.032054626428896, 0.128958931608947 ,0,0,
                         0.119317341270894, 0, 0.247654804887413, 0, 0.119317341270894, 0.247654804887413 ,0,0}, 1e-15);

        //test partials at node 7
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_7.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at internal node 7:", partials,
                new double[]{0.060425520687175, 0, 0, 0, 0.060425520687175, 0, 0, 0,
                        0.769152619166653e-4, 0, 0, 0, 0.769152619166653e-4, 0, 0, 0}, 1e-15);

        //test partials at root node
        /*likelihoodNegInf.getLikelihoodCore().getNodePartials(parent.getNr(), partials);
        assertArrayEquals("Assert correct likelihood at root:",
                partials,
                new double[]{0.002468211088543, 0, 0, 0, 0.002468211088543, 0, 0, 0,
                        0.429517181079105e-5, 0, 0, 0, 0.429517181079105e-5, 0, 0, 0}, 1e-15);
*/
        logP = likelihoodNegInf.calculateLogP();
        assertEquals("second calculation", -13.391340286074429, logP, 1e-14);

    }

    @Test
    public void testLikelihoodScaling(){
        //test alignment that is very unlikely given the tree, and assess scaling

        //String alignmentText = '#NEXUS\n\nbegin taxa;\n\tdimensions ntax=556;\n\ttaxlabels 0 1 10 100 101 102 103 104 105 106 107 108 109 11 110 111 112 113 114 115 116 117 118 119 12 120 121 122 123 124 125 126 127 128 129 13 130 131 132 133 134 135 136 137 138 139 14 140 141 142 143 144 145 146 147 148 149 15 150 151 152 153 154 155 156 157 158 159 16 160 161 162 163 164 165 166 167 168 169 17 170 171 172 173 174 175 176 177 178 179 18 180 181 182 183 184 185 186 187 188 189 19 190 191 192 193 194 195 196 197 198 199 2 20 200 201 202 203 204 205 206 207 208 209 21 210 211 212 213 214 215 216 217 218 219 22 220 221 222 223 224 225 226 227 228 229 23 230 231 232 233 234 235 236 237 238 239 24 240 241 242 243 244 245 246 247 248 249 25 250 251 252 253 254 255 256 257 258 259 26 260 261 262 263 264 265 266 267 268 269 27 270 271 272 273 274 275 276 277 278 279 28 280 281 282 283 284 285 286 287 288 289 29 290 291 292 293 294 295 296 297 298 299 3 30 300 301 302 303 304 305 306 307 308 309 31 310 311 312 313 314 315 316 317 318 319 32 320 321 322 323 324 325 326 327 328 329 33 330 331 332 333 334 335 336 337 338 339 34 340 341 342 343 344 345 346 347 348 349 35 350 351 352 353 354 355 356 357 358 359 36 360 361 362 363 364 365 366 367 368 369 37 370 371 372 373 374 375 376 377 378 379 38 380 381 382 383 384 385 386 387 388 389 39 390 391 392 393 394 395 396 397 398 399 4 40 400 401 402 403 404 405 406 407 408 409 41 410 411 412 413 414 415 416 417 418 419 42 420 421 422 423 424 425 426 427 428 429 43 430 431 432 433 434 435 436 437 438 439 44 440 441 442 443 444 445 446 447 448 449 45 450 451 452 453 454 455 456 457 458 459 46 460 461 462 463 464 465 466 467 468 469 47 470 471 472 473 474 475 476 477 478 479 48 480 481 482 483 484 485 486 487 488 489 49 490 491 492 493 494 495 496 497 498 499 5 50 500 501 502 503 504 505 506 507 508 509 51 510 511 512 513 514 515 516 517 518 519 52 520 521 522 523 524 525 526 527 528 529 53 530 531 532 533 534 535 536 537 538 539 54 540 541 542 543 544 545 546 547 548 549 55 550 551 552 553 554 555 56 57 58 59 6 60 61 62 63 64 65 66 67 68 69 7 70 71 72 73 74 75 76 77 78 79 8 80 81 82 83 84 85 86 87 88 89 9 90 91 92 93 94 95 96 97 98 99;\nend;\n\nbegin characters;\n\tdimensions nchar=5;\n\tformat datatype=integer;\n\tmatrix \n\t\t0 8,0,0,0,0\n\t\t1 8,12,12,0,0\n\t\t10 9,3,11,0,0\n\t\t100 12,13,12,0,0\n\t\t101 0,0,0,0,0\n\t\t102 8,0,0,0,0\n\t\t103 4,12,12,0,0\n\t\t104 12,0,0,0,0\n\t\t105 8,1,13,0,0\n\t\t106 0,0,0,0,0\n\t\t107 12,0,0,0,0\n\t\t108 0,0,0,0,0\n\t\t109 13,12,0,0,0\n\t\t11 8,12,0,0,0\n\t\t110 8,5,9,0,0\n\t\t111 9,8,0,0,0\n\t\t112 0,0,0,0,0\n\t\t113 8,0,0,0,0\n\t\t114 9,0,0,0,0\n\t\t115 8,0,0,0,0\n\t\t116 8,0,0,0,0\n\t\t117 4,0,0,0,0\n\t\t118 8,13,0,0,0\n\t\t119 8,0,0,0,0\n\t\t12 8,6,4,0,0\n\t\t120 8,4,9,12,0\n\t\t121 4,0,0,0,0\n\t\t122 0,0,0,0,0\n\t\t123 8,12,0,0,0\n\t\t124 8,11,0,0,0\n\t\t125 8,1,0,0,0\n\t\t126 8,4,0,0,0\n\t\t127 8,6,0,0,0\n\t\t128 6,12,0,0,0\n\t\t129 8,0,0,0,0\n\t\t13 13,13,0,0,0\n\t\t130 8,12,0,0,0\n\t\t131 8,0,0,0,0\n\t\t132 12,13,0,0,0\n\t\t133 8,6,0,0,0\n\t\t134 8,5,0,0,0\n\t\t135 8,12,0,0,0\n\t\t136 5,0,0,0,0\n\t\t137 0,0,0,0,0\n\t\t138 8,0,0,0,0\n\t\t139 8,0,0,0,0\n\t\t14 8,6,12,0,0\n\t\t140 8,0,0,0,0\n\t\t141 8,0,0,0,0\n\t\t142 12,4,0,0,0\n\t\t143 4,0,0,0,0\n\t\t144 0,0,0,0,0\n\t\t145 13,0,0,0,0\n\t\t146 8,0,0,0,0\n\t\t147 0,0,0,0,0\n\t\t148 0,0,0,0,0\n\t\t149 8,13,4,13,0\n\t\t15 0,0,0,0,0\n\t\t150 8,12,0,0,0\n\t\t151 9,0,0,0,0\n\t\t152 8,0,0,0,0\n\t\t153 8,0,0,0,0\n\t\t154 8,12,9,5,0\n\t\t155 8,13,4,0,0\n\t\t156 8,12,0,0,0\n\t\t157 8,0,0,0,0\n\t\t158 0,0,0,0,0\n\t\t159 8,12,0,0,0\n\t\t16 8,6,8,0,0\n\t\t160 0,0,0,0,0\n\t\t161 8,12,0,0,0\n\t\t162 0,0,0,0,0\n\t\t163 12,1,0,0,0\n\t\t164 8,7,0,0,0\n\t\t165 8,12,0,0,0\n\t\t166 0,0,0,0,0\n\t\t167 8,0,0,0,0\n\t\t168 9,0,0,0,0\n\t\t169 0,0,0,0,0\n\t\t17 8,0,0,0,0\n\t\t170 5,4,13,6,0\n\t\t171 8,6,0,0,0\n\t\t172 13,0,0,0,0\n\t\t173 8,0,0,0,0\n\t\t174 11,0,0,0,0\n\t\t175 12,0,0,0,0\n\t\t176 13,0,0,0,0\n\t\t177 8,11,0,0,0\n\t\t178 13,0,0,0,0\n\t\t179 8,9,9,12,0\n\t\t18 8,13,0,0,0\n\t\t180 8,9,12,0,0\n\t\t181 8,4,9,0,0\n\t\t182 8,0,0,0,0\n\t\t183 8,0,0,0,0\n\t\t184 8,0,0,0,0\n\t\t185 13,0,0,0,0\n\t\t186 8,0,0,0,0\n\t\t187 8,6,0,0,0\n\t\t188 0,0,0,0,0\n\t\t189 12,0,0,0,0\n\t\t19 8,12,0,0,0\n\t\t190 8,0,0,0,0\n\t\t191 8,12,0,0,0\n\t\t192 13,12,6,0,0\n\t\t193 9,0,0,0,0\n\t\t194 8,12,0,0,0\n\t\t195 12,0,0,0,0\n\t\t196 13,2,0,0,0\n\t\t197 8,11,12,0,0\n\t\t198 8,0,0,0,0\n\t\t199 8,12,0,0,0\n\t\t2 8,8,0,0,0\n\t\t20 0,0,0,0,0\n\t\t200 13,9,0,0,0\n\t\t201 4,12,12,0,0\n\t\t202 0,0,0,0,0\n\t\t203 9,0,0,0,0\n\t\t204 8,12,0,0,0\n\t\t205 4,0,0,0,0\n\t\t206 12,0,0,0,0\n\t\t207 8,0,0,0,0\n\t\t208 0,0,0,0,0\n\t\t209 8,13,0,0,0\n\t\t21 8,0,0,0,0\n\t\t210 8,12,0,0,0\n\t\t211 8,0,0,0,0\n\t\t212 8,0,0,0,0\n\t\t213 8,0,0,0,0\n\t\t214 12,12,12,0,0\n\t\t215 4,6,0,0,0\n\t\t216 8,4,9,12,9\n\t\t217 8,13,0,0,0\n\t\t218 5,0,0,0,0\n\t\t219 8,12,12,0,0\n\t\t22 8,0,0,0,0\n\t\t220 9,0,0,0,0\n\t\t221 0,0,0,0,0\n\t\t222 0,0,0,0,0\n\t\t223 8,0,0,0,0\n\t\t224 8,1,9,0,0\n\t\t225 0,0,0,0,0\n\t\t226 8,0,0,0,0\n\t\t227 12,0,0,0,0\n\t\t228 8,0,0,0,0\n\t\t229 8,9,9,0,0\n\t\t23 8,13,4,0,0\n\t\t230 8,0,0,0,0\n\t\t231 0,0,0,0,0\n\t\t232 10,0,0,0,0\n\t\t233 12,0,0,0,0\n\t\t234 8,12,0,0,0\n\t\t235 3,13,0,0,0\n\t\t236 8,12,0,0,0\n\t\t237 9,0,0,0,0\n\t\t238 8,8,0,0,0\n\t\t239 8,0,0,0,0\n\t\t24 8,0,0,0,0\n\t\t240 8,0,0,0,0\n\t\t241 8,0,0,0,0\n\t\t242 8,13,0,0,0\n\t\t243 12,0,0,0,0\n\t\t244 0,0,0,0,0\n\t\t245 0,0,0,0,0\n\t\t246 8,8,0,0,0\n\t\t247 8,12,0,0,0\n\t\t248 10,0,0,0,0\n\t\t249 12,0,0,0,0\n\t\t25 0,0,0,0,0\n\t\t250 0,0,0,0,0\n\t\t251 12,0,0,0,0\n\t\t252 8,0,0,0,0\n\t\t253 8,0,0,0,0\n\t\t254 9,0,0,0,0\n\t\t255 12,0,0,0,0\n\t\t256 12,0,0,0,0\n\t\t257 0,0,0,0,0\n\t\t258 12,0,0,0,0\n\t\t259 0,0,0,0,0\n\t\t26 4,0,0,0,0\n\t\t260 0,0,0,0,0\n\t\t261 12,0,0,0,0\n\t\t262 8,13,0,0,0\n\t\t263 0,0,0,0,0\n\t\t264 9,0,0,0,0\n\t\t265 0,0,0,0,0\n\t\t266 8,0,0,0,0\n\t\t267 9,0,0,0,0\n\t\t268 8,0,0,0,0\n\t\t269 8,13,8,0,0\n\t\t27 8,9,0,0,0\n\t\t270 8,12,0,0,0\n\t\t271 8,0,0,0,0\n\t\t272 0,0,0,0,0\n\t\t273 8,1,0,0,0\n\t\t274 13,0,0,0,0\n\t\t275 8,13,0,0,0\n\t\t276 11,12,12,12,0\n\t\t277 5,0,0,0,0\n\t\t278 8,12,13,0,0\n\t\t279 9,0,0,0,0\n\t\t28 0,0,0,0,0\n\t\t280 9,0,0,0,0\n\t\t281 8,0,0,0,0\n\t\t282 0,0,0,0,0\n\t\t283 8,13,0,0,0\n\t\t284 8,0,0,0,0\n\t\t285 8,12,0,0,0\n\t\t286 0,0,0,0,0\n\t\t287 0,0,0,0,0\n\t\t288 12,0,0,0,0\n\t\t289 0,0,0,0,0\n\t\t29 8,0,0,0,0\n\t\t290 9,0,0,0,0\n\t\t291 8,12,0,0,0\n\t\t292 8,12,0,0,0\n\t\t293 9,4,13,0,0\n\t\t294 8,0,0,0,0\n\t\t295 8,0,0,0,0\n\t\t296 8,0,0,0,0\n\t\t297 8,0,0,0,0\n\t\t298 8,0,0,0,0\n\t\t299 8,6,0,0,0\n\t\t3 0,0,0,0,0\n\t\t30 8,12,0,0,0\n\t\t300 8,4,0,0,0\n\t\t301 8,0,0,0,0\n\t\t302 8,9,0,0,0\n\t\t303 8,0,0,0,0\n\t\t304 8,0,0,0,0\n\t\t305 8,1,0,0,0\n\t\t306 8,12,0,0,0\n\t\t307 8,4,9,0,0\n\t\t308 10,0,0,0,0\n\t\t309 12,13,5,12,0\n\t\t31 8,0,0,0,0\n\t\t310 8,0,0,0,0\n\t\t311 11,9,0,0,0\n\t\t312 8,13,12,8,0\n\t\t313 8,0,0,0,0\n\t\t314 0,0,0,0,0\n\t\t315 12,0,0,0,0\n\t\t316 8,8,0,0,0\n\t\t317 0,0,0,0,0\n\t\t318 8,0,0,0,0\n\t\t319 12,0,0,0,0\n\t\t32 0,0,0,0,0\n\t\t320 0,0,0,0,0\n\t\t321 0,0,0,0,0\n\t\t322 0,0,0,0,0\n\t\t323 0,0,0,0,0\n\t\t324 7,0,0,0,0\n\t\t325 6,0,0,0,0\n\t\t326 8,0,0,0,0\n\t\t327 8,0,0,0,0\n\t\t328 0,0,0,0,0\n\t\t329 8,0,0,0,0\n\t\t33 12,0,0,0,0\n\t\t330 12,0,0,0,0\n\t\t331 8,0,0,0,0\n\t\t332 8,10,0,0,0\n\t\t333 8,0,0,0,0\n\t\t334 13,13,9,0,0\n\t\t335 8,0,0,0,0\n\t\t336 8,0,0,0,0\n\t\t337 8,9,9,9,0\n\t\t338 8,12,0,0,0\n\t\t339 0,0,0,0,0\n\t\t34 6,6,0,0,0\n\t\t340 8,0,0,0,0\n\t\t341 8,1,0,0,0\n\t\t342 0,0,0,0,0\n\t\t343 8,0,0,0,0\n\t\t344 8,6,12,6,0\n\t\t345 12,12,0,0,0\n\t\t346 9,0,0,0,0\n\t\t347 8,8,0,0,0\n\t\t348 8,0,0,0,0\n\t\t349 11,0,0,0,0\n\t\t35 8,0,0,0,0\n\t\t350 8,0,0,0,0\n\t\t351 8,12,0,0,0\n\t\t352 0,0,0,0,0\n\t\t353 0,0,0,0,0\n\t\t354 8,0,0,0,0\n\t\t355 8,0,0,0,0\n\t\t356 8,13,0,0,0\n\t\t357 8,3,2,0,0\n\t\t358 8,12,0,0,0\n\t\t359 0,0,0,0,0\n\t\t36 8,12,0,0,0\n\t\t360 0,0,0,0,0\n\t\t361 0,0,0,0,0\n\t\t362 8,0,0,0,0\n\t\t363 0,0,0,0,0\n\t\t364 8,0,0,0,0\n\t\t365 12,0,0,0,0\n\t\t366 10,0,0,0,0\n\t\t367 6,0,0,0,0\n\t\t368 13,0,0,0,0\n\t\t369 12,0,0,0,0\n\t\t37 13,0,0,0,0\n\t\t370 8,12,2,0,0\n\t\t371 8,4,9,12,0\n\t\t372 9,0,0,0,0\n\t\t373 9,0,0,0,0\n\t\t374 8,12,0,0,0\n\t\t375 4,12,0,0,0\n\t\t376 8,12,0,0,0\n\t\t377 0,0,0,0,0\n\t\t378 8,0,0,0,0\n\t\t379 0,0,0,0,0\n\t\t38 8,12,0,0,0\n\t\t380 12,2,0,0,0\n\t\t381 8,0,0,0,0\n\t\t382 12,8,0,0,0\n\t\t383 0,0,0,0,0\n\t\t384 4,0,0,0,0\n\t\t385 12,0,0,0,0\n\t\t386 8,9,0,0,0\n\t\t387 8,0,0,0,0\n\t\t388 8,13,0,0,0\n\t\t389 8,0,0,0,0\n\t\t39 8,9,0,0,0\n\t\t390 13,0,0,0,0\n\t\t391 8,0,0,0,0\n\t\t392 5,9,0,0,0\n\t\t393 8,0,0,0,0\n\t\t394 8,12,12,0,0\n\t\t395 12,0,0,0,0\n\t\t396 8,0,0,0,0\n\t\t397 12,6,0,0,0\n\t\t398 8,12,0,0,0\n\t\t399 8,6,0,0,0\n\t\t4 6,0,0,0,0\n\t\t40 8,0,0,0,0\n\t\t400 8,8,0,0,0\n\t\t401 8,8,8,0,0\n\t\t402 11,9,12,0,0\n\t\t403 8,6,0,0,0\n\t\t404 8,0,0,0,0\n\t\t405 9,0,0,0,0\n\t\t406 9,0,0,0,0\n\t\t407 0,0,0,0,0\n\t\t408 0,0,0,0,0\n\t\t409 13,4,0,0,0\n\t\t41 8,5,12,0,0\n\t\t410 4,0,0,0,0\n\t\t411 8,13,0,0,0\n\t\t412 12,8,12,0,0\n\t\t413 12,0,0,0,0\n\t\t414 0,0,0,0,0\n\t\t415 3,8,0,0,0\n\t\t416 8,0,0,0,0\n\t\t417 8,0,0,0,0\n\t\t418 12,0,0,0,0\n\t\t419 13,0,0,0,0\n\t\t42 8,0,0,0,0\n\t\t420 8,5,5,0,0\n\t\t421 4,0,0,0,0\n\t\t422 8,0,0,0,0\n\t\t423 12,0,0,0,0\n\t\t424 12,0,0,0,0\n\t\t425 0,0,0,0,0\n\t\t426 0,0,0,0,0\n\t\t427 0,0,0,0,0\n\t\t428 8,12,0,0,0\n\t\t429 8,0,0,0,0\n\t\t43 8,0,0,0,0\n\t\t430 8,0,0,0,0\n\t\t431 8,0,0,0,0\n\t\t432 8,12,12,0,0\n\t\t433 8,8,12,0,0\n\t\t434 0,0,0,0,0\n\t\t435 8,9,11,0,0\n\t\t436 8,12,0,0,0\n\t\t437 8,1,0,0,0\n\t\t438 9,0,0,0,0\n\t\t439 8,12,0,0,0\n\t\t44 12,0,0,0,0\n\t\t440 12,0,0,0,0\n\t\t441 8,4,0,0,0\n\t\t442 11,13,10,0,0\n\t\t443 6,9,0,0,0\n\t\t444 9,0,0,0,0\n\t\t445 12,0,0,0,0\n\t\t446 8,0,0,0,0\n\t\t447 8,12,9,0,0\n\t\t448 0,0,0,0,0\n\t\t449 8,4,0,0,0\n\t\t45 12,0,0,0,0\n\t\t450 12,0,0,0,0\n\t\t451 0,0,0,0,0\n\t\t452 10,12,0,0,0\n\t\t453 7,0,0,0,0\n\t\t454 9,0,0,0,0\n\t\t455 13,0,0,0,0\n\t\t456 12,12,0,0,0\n\t\t457 9,0,0,0,0\n\t\t458 0,0,0,0,0\n\t\t459 9,12,0,0,0\n\t\t46 8,0,0,0,0\n\t\t460 8,12,5,7,6\n\t\t461 12,12,13,0,0\n\t\t462 0,0,0,0,0\n\t\t463 0,0,0,0,0\n\t\t464 12,8,0,0,0\n\t\t465 12,0,0,0,0\n\t\t466 8,9,0,0,0\n\t\t467 8,1,0,0,0\n\t\t468 8,0,0,0,0\n\t\t469 8,11,12,8,0\n\t\t47 5,0,0,0,0\n\t\t470 12,4,0,0,0\n\t\t471 4,0,0,0,0\n\t\t472 8,12,8,0,0\n\t\t473 8,13,13,0,0\n\t\t474 8,0,0,0,0\n\t\t475 8,0,0,0,0\n\t\t476 8,12,13,0,0\n\t\t477 12,0,0,0,0\n\t\t478 12,0,0,0,0\n\t\t479 0,0,0,0,0\n\t\t48 0,0,0,0,0\n\t\t480 8,0,0,0,0\n\t\t481 8,13,0,0,0\n\t\t482 8,8,0,0,0\n\t\t483 8,8,0,0,0\n\t\t484 8,4,8,0,0\n\t\t485 0,0,0,0,0\n\t\t486 8,0,0,0,0\n\t\t487 8,0,0,0,0\n\t\t488 8,4,9,0,0\n\t\t489 8,8,9,0,0\n\t\t49 8,12,0,0,0\n\t\t490 12,0,0,0,0\n\t\t491 12,0,0,0,0\n\t\t492 7,0,0,0,0\n\t\t493 9,9,0,0,0\n\t\t494 8,0,0,0,0\n\t\t495 11,12,0,0,0\n\t\t496 4,0,0,0,0\n\t\t497 11,0,0,0,0\n\t\t498 12,0,0,0,0\n\t\t499 11,13,0,0,0\n\t\t5 8,0,0,0,0\n\t\t50 8,12,0,0,0\n\t\t500 2,0,0,0,0\n\t\t501 12,9,0,0,0\n\t\t502 8,0,0,0,0\n\t\t503 8,3,0,0,0\n\t\t504 12,0,0,0,0\n\t\t505 6,0,0,0,0\n\t\t506 8,0,0,0,0\n\t\t507 8,13,13,0,0\n\t\t508 8,11,12,0,0\n\t\t509 9,0,0,0,0\n\t\t51 8,0,0,0,0\n\t\t510 0,0,0,0,0\n\t\t511 6,12,0,0,0\n\t\t512 13,12,0,0,0\n\t\t513 8,0,0,0,0\n\t\t514 8,12,0,0,0\n\t\t515 5,0,0,0,0\n\t\t516 6,12,0,0,0\n\t\t517 8,0,0,0,0\n\t\t518 9,0,0,0,0\n\t\t519 8,0,0,0,0\n\t\t52 5,12,0,0,0\n\t\t520 8,0,0,0,0\n\t\t521 8,12,12,0,0\n\t\t522 8,8,4,12,1\n\t\t523 8,9,0,0,0\n\t\t524 8,0,0,0,0\n\t\t525 0,0,0,0,0\n\t\t526 7,0,0,0,0\n\t\t527 8,0,0,0,0\n\t\t528 8,0,0,0,0\n\t\t529 8,12,0,0,0\n\t\t53 13,8,6,0,0\n\t\t530 0,0,0,0,0\n\t\t531 8,0,0,0,0\n\t\t532 8,13,4,9,9\n\t\t533 8,9,0,0,0\n\t\t534 8,0,0,0,0\n\t\t535 12,0,0,0,0\n\t\t536 8,0,0,0,0\n\t\t537 8,0,0,0,0\n\t\t538 9,0,0,0,0\n\t\t539 8,12,0,0,0\n\t\t54 12,12,0,0,0\n\t\t540 0,0,0,0,0\n\t\t541 0,0,0,0,0\n\t\t542 8,12,0,0,0\n\t\t543 8,8,0,0,0\n\t\t544 8,0,0,0,0\n\t\t545 8,0,0,0,0\n\t\t546 0,0,0,0,0\n\t\t547 0,0,0,0,0\n\t\t548 12,0,0,0,0\n\t\t549 0,0,0,0,0\n\t\t55 9,0,0,0,0\n\t\t550 8,12,0,0,0\n\t\t551 12,0,0,0,0\n\t\t552 0,0,0,0,0\n\t\t553 8,12,12,0,0\n\t\t554 0,0,0,0,0\n\t\t555 8,12,12,0,0\n\t\t56 5,13,0,0,0\n\t\t57 0,0,0,0,0\n\t\t58 8,0,0,0,0\n\t\t59 8,0,0,0,0\n\t\t6 8,0,0,0,0\n\t\t60 8,12,9,0,0\n\t\t61 8,4,9,0,0\n\t\t62 8,12,12,12,11\n\t\t63 8,10,9,12,0\n\t\t64 8,0,0,0,0\n\t\t65 12,10,12,10,0\n\t\t66 8,4,9,12,0\n\t\t67 8,13,8,0,0\n\t\t68 8,0,0,0,0\n\t\t69 8,0,0,0,0\n\t\t7 0,0,0,0,0\n\t\t70 8,0,0,0,0\n\t\t71 8,4,9,2,0\n\t\t72 12,9,12,0,0\n\t\t73 4,0,0,0,0\n\t\t74 0,0,0,0,0\n\t\t75 0,0,0,0,0\n\t\t76 0,0,0,0,0\n\t\t77 0,0,0,0,0\n\t\t78 8,0,0,0,0\n\t\t79 8,6,0,0,0\n\t\t8 4,0,0,0,0\n\t\t80 0,0,0,0,0\n\t\t81 7,0,0,0,0\n\t\t82 8,9,0,0,0\n\t\t83 8,13,0,0,0\n\t\t84 12,12,0,0,0\n\t\t85 8,1,0,0,0\n\t\t86 8,0,0,0,0\n\t\t87 9,0,0,0,0\n\t\t88 8,9,0,0,0\n\t\t89 8,0,0,0,0\n\t\t9 12,13,0,0,0\n\t\t90 13,6,8,0,0\n\t\t91 0,0,0,0,0\n\t\t92 0,0,0,0,0\n\t\t93 8,0,0,0,0\n\t\t94 13,0,0,0,0\n\t\t95 8,7,0,0,0\n\t\t96 12,12,0,0,0\n\t\t97 0,0,0,0,0\n\t\t98 8,12,0,0,0\n\t\t99 0,0,0,0,0;\nend;';

        String fileName = "test/tidetree/distributions/simulate_alignment_and_tree.seed=2.1.alignment.nexus";

        AlignmentFromNexus alignment = new AlignmentFromNexus();
        alignment.initByName("fileName", fileName, "dataType", "integer");

        //Alignment alignment = NexusParser(alignmentText);
        // init alignment
        //Alignment alignment = new Alignment();
        //Alignment alignment = new AlignmentFromNexus(fileName);
        /*Sequence a = new Sequence("0", "0,0");
        Sequence b = new Sequence("1", "1,2");
        Sequence c = new Sequence("2", "1,2");
        Sequence d = new Sequence("3", "2,1");
        Sequence e = new Sequence("4", "3,3");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "sequence", e, "dataType", "integer", "stateCount", 4);
        */

        String newickTree = "(((((((275[&type=0]:19.91106435989237,((130[&type=0]:19.483630674546674,(150[&type=0]:18.007446644500877,(((60[&type=0]:14.828711130835476,273[&type=0]:14.828711130835476)689:1.844745369661748,(116[&type=0]:16.363707024544613,(119[&type=0]:15.077211738025817,534[&type=0]:15.077211738025817)703:1.2864952865187966)788:0.3097494759526107)813:0.20841147887067635,68[&type=0]:16.8818679793679)838:1.125578665132977)925:1.476184030045797)1009:0.31709911975748994,234[&type=0]:19.800729794304164)1025:0.11033456558820731)1031:1.290991750414726,((((431[&type=0]:15.864484161003832,(24[&type=0]:14.720572563429824,(210[&type=0]:13.54482195163459,508[&type=0]:13.54482195163459)631:1.1757506117952339)675:1.1439115975740073)747:1.9726566574452047,(187[&type=0]:16.012558156028675,347[&type=0]:16.012558156028675)761:1.8245826624203616)910:0.5087966037371991,(40[&type=0]:13.315347238791137,355[&type=0]:13.315347238791137)624:5.030590183395098)941:1.0976790797227842,((331[&type=0]:17.680407170040933,469[&type=0]:17.680407170040933)901:0.9746159292736323,468[&type=0]:18.655023099314565)959:0.7885934025944543)1006:1.7584396083980778)1068:2.2100324862133824,(23[&type=0]:23.224057679454056,(((((318[&type=0]:17.43694092975312,(93[&type=0]:16.843789250513097,542[&type=0]:16.843789250513097)836:0.5931516792400231)879:1.499955630336192,((439[&type=0]:14.726746755627884,((503[&type=0]:11.430681968430711,157[&type=0]:11.430681968430711)581:2.0574159028912025,388[&type=0]:13.488097871321914)629:1.2386488843059702)676:3.2331947272574606,338[&type=0]:17.959941482885345)922:0.9769550772039679)979:0.06863150772834103,((517[&type=0]:17.72159716216438,207[&type=0]:17.72159716216438)904:0.8638431950479948,(341[&type=0]:11.364009123515032,467[&type=0]:11.364009123515032)577:7.221431233697341)953:0.42008771060528005)983:3.0575188600402434,(((((304[&type=0]:14.731425992592838,252[&type=0]:14.731425992592838)677:1.8498840874566245,167[&type=0]:16.581310080049462)803:0.05150346623635116,(519[&type=0]:15.277799208153796,(326[&type=0]:14.907673371914067,124[&type=0]:14.907673371914067)694:0.3701258362397297)713:1.3550143381320172)808:1.261796328649897,(229[&type=0]:15.871822349101521,337[&type=0]:15.871822349101521)748:2.0227875258341896)915:1.8942769476775752,(219[&type=0]:19.644143940407712,(182[&type=0]:19.413070132511685,((((95[&type=0]:15.366969698906212,523[&type=0]:15.366969698906212)718:2.6890722422529674,((242[&type=0]:14.49313813269861,198[&type=0]:14.49313813269861)663:1.394950732160403,((126[&type=0]:12.856829327742382,449[&type=0]:12.856829327742382)610:3.030560555776752,476[&type=0]:15.887389883519134)749:6.989813398785572E-4)750:2.167953076300167)929:0.15512034216082782,(301[&type=0]:17.10638359642612,(354[&type=0]:11.172711462290458,335[&type=0]:11.172711462290458)573:5.933672134135662)851:1.1047786868938871)935:0.8740056353566246,(159[&type=0]:18.594786594219215,(446[&type=0]:15.969603463122098,(83[&type=0]:15.844279047410101,362[&type=0]:15.844279047410101)745:0.12532441571199726)755:2.625183131097117)954:0.49038132445741667)989:0.3279022138350527)1002:0.23107380789602772)1020:0.1447428822055734)1024:2.2741601052446114)1090:0.4285220129408813,(((537[&type=0]:15.765104488517842,482[&type=0]:15.765104488517842)738:1.0604179904079913,(507[&type=0]:13.997900448522097,153[&type=0]:13.997900448522097)642:2.8276220304037363)835:1.6280747997900633,(18[&type=0]:17.53410171536371,((487[&type=0]:13.182993991626756,152[&type=0]:13.182993991626756)621:3.7147047993578077,((6[&type=0]:14.771496583012741,127[&type=0]:14.771496583012741)681:1.2387845267415187,89[&type=0]:16.01028110975426)760:0.8874176812303034)839:0.6364029243791478)885:0.9194955633521857)947:4.037971662082882)1095:0.7324887386552774)1101:0.18803091706642405)1103:0.1990146331543876,(((((416[&type=0]:14.800007852194769,64[&type=0]:14.800007852194769)684:0.8340260335095326,125[&type=0]:15.634033885704302)735:4.02350543829488,((((262[&type=0]:15.250475252545812,(514[&type=0]:15.024301182939048,435[&type=0]:15.024301182939048)700:0.22617406960676334)711:1.0521623553348718,(283[&type=0]:12.31496999020384,544[&type=0]:12.31496999020384)600:3.987667617676843)784:1.9222099989588273,(((70[&type=0]:6.923365805115488,184[&type=0]:6.923365805115488)558:4.690137060440318,281[&type=0]:11.613502865555805)585:2.7554217481447836,370[&type=0]:14.368924613700589)659:3.855922993138922)938:0.351799867889536,((305[&type=0]:16.594534646255426,(433[&type=0]:14.699563852495977,(524[&type=0]:9.48853755431614,531[&type=0]:9.48853755431614)561:5.211026298179837)672:1.8949707937594482)805:1.8079931581703619,((299[&type=0]:16.022662793532383,(79[&type=0]:14.802321010752914,14[&type=0]:14.802321010752914)685:1.2203417827794691)762:0.7772997788094962,((133[&type=0]:15.830325471447008,(16[&type=0]:15.111896943074468,399[&type=0]:15.111896943074468)707:0.7184285283725398)743:0.28152845676801874,(12[&type=0]:13.480251862727366,171[&type=0]:13.480251862727366)628:2.6316020654876606)767:0.6881086441268529)831:1.6025652320839079)945:0.17411967030325926)952:1.0808918492701345)1021:0.8193740376777328,((417[&type=0]:17.178638020277532,186[&type=0]:17.178638020277532)860:2.045511067818424,((532[&type=0]:16.557284376495915,149[&type=0]:16.557284376495915)802:1.07168899053989,(356[&type=0]:12.615143997025248,67[&type=0]:12.615143997025248)605:5.013829370010557)896:1.5951757210601514)994:1.2527642735809579)1052:0.8299049111193746,((((((63[&type=0]:16.13348188875136,(520[&type=0]:14.97719203765199,(336[&type=0]:14.20445222058657,42[&type=0]:14.20445222058657)650:0.7727398170654212)697:1.156289851099368)769:0.6774973588611495,(502[&type=0]:12.943759250100316,(374[&type=0]:12.341900747757201,285[&type=0]:12.341900747757201)602:0.6018585023431147)611:3.8672199975121924)833:0.12586936318606945,(553[&type=0]:15.350067734647782,(155[&type=0]:14.82152635639175,(11[&type=0]:14.703860603486909,129[&type=0]:14.703860603486909)673:0.11766575290484127)688:0.5285413782560315)717:1.5867808761507955)840:3.55106490240944,(110[&type=0]:16.54784316394746,41[&type=0]:16.54784316394746)800:3.9400703492605587)1053:0.6872020461263872,(246[&type=0]:18.67255510150064,((50[&type=0]:15.574124777888109,29[&type=0]:15.574124777888109)732:3.097964884129464,(105[&type=0]:16.757744149543335,19[&type=0]:16.757744149543335)825:1.9143455124742381)961:4.654394830687636E-4)962:2.502560457833763)1067:0.09910294133269915,(((543[&type=0]:14.864752213196155,294[&type=0]:14.864752213196155)691:2.9677419604876363,(224[&type=0]:17.06094400703343,348[&type=0]:17.06094400703343)848:0.7715501666503606)909:1.5007675082157874,(513[&type=0]:17.394158294070557,266[&type=0]:17.394158294070557)874:1.9391033878290216)998:1.9409568187675248)1070:0.032599772129184856)1071:2.304284956878579)1105:0.5064073369656832,((((550[&type=0]:15.904148832523067,135[&type=0]:15.904148832523067)752:5.7024382544961085,((((240[&type=0]:16.58561432432156,474[&type=0]:16.58561432432156)804:1.0196393942872746,340[&type=0]:17.605253718608836)894:1.2704895763360398,(278[&type=0]:14.535791107714397,483[&type=0]:14.535791107714397)664:4.3399521872304785)973:1.6348845323086607,59[&type=0]:20.510627827253536)1056:1.095959259765639)1076:0.38217075556220337,(((472[&type=0]:18.006288931282484,529[&type=0]:18.006288931282484)924:1.0768553465708735,((447[&type=0]:15.69220657237755,((539[&type=0]:14.362875033482483,(394[&type=0]:13.133549799832226,292[&type=0]:13.133549799832226)617:1.2293252336502576)658:1.0992715971948712,521[&type=0]:15.462146630677355)724:0.23005994170019584)736:0.7648779146544342,(460[&type=0]:15.827851474123829,306[&type=0]:15.827851474123829)741:0.6292330129081556)793:2.6260597908213725)988:0.9351731850878551,(199[&type=0]:17.446792270034102,545[&type=0]:17.446792270034102)880:2.57152519290711)1040:1.9704403796401664)1088:1.4926125330578799,((((39[&type=0]:21.775244594863118,((((((386[&type=0]:11.127151022544542,36[&type=0]:11.127151022544542)570:5.104215754437536,(269[&type=0]:12.971888573206822,358[&type=0]:12.971888573206822)613:3.259478203775256)777:1.0009823393552786,(140[&type=0]:16.005851215641336,134[&type=0]:16.005851215641336)759:1.226497900696021)864:2.262881798828495,(376[&type=0]:17.046286171665656,78[&type=0]:17.046286171665656)846:2.4489447435001956)1011:0.6253852858734348,((436[&type=0]:19.94381655896472,69[&type=0]:19.94381655896472)1035:0.01616564762333894,((191[&type=0]:16.131079102962644,22[&type=0]:16.131079102962644)768:3.605657945454748,(378[&type=0]:17.54680186053027,236[&type=0]:17.54680186053027)887:2.189935187887123)1023:0.22324515817066626)1036:0.16063399445122784)1042:1.4491213301464434,(31[&type=0]:18.682696633948755,(((2[&type=0]:14.016223058618326,46[&type=0]:14.016223058618326)643:1.299984735153945,(333[&type=0]:13.362262206208543,295[&type=0]:13.362262206208543)626:1.9539455875637284)715:1.8532012716517556,313[&type=0]:17.169409065424027)855:1.5132875685247278)963:2.887040897236975)1075:0.2055070636773877)1080:0.4982005953559394,((((((113[&type=0]:10.29163247765084,329[&type=0]:10.29163247765084)564:5.5534396001471915,428[&type=0]:15.845072077798031)746:1.3294978333557648,238[&type=0]:17.174569911153796)858:1.7190732023769115,(389[&type=0]:17.57944561895723,343[&type=0]:17.57944561895723)889:1.3141974945734773)974:1.0355816561659843,(((209[&type=0]:14.88035134428572,(310[&type=0]:14.861154382804454,98[&type=0]:14.861154382804454)690:0.019196961481265262)693:2.920940718686177,161[&type=0]:17.801292062971896)907:1.0078631625378378,(271[&type=0]:15.32985592225559,141[&type=0]:15.32985592225559)716:3.4792993032541446)969:1.1200695441869577)1034:1.8216585222068815,((((((((165[&type=0]:15.630556040108786,302[&type=0]:15.630556040108786)734:1.0822200232159105,401[&type=0]:16.712776063324696)816:1.1541538310204835,(118[&type=0]:16.550753660506587,429[&type=0]:16.550753660506587)801:1.3161762338385934)912:0.8818180979141879,((180[&type=0]:17.192557630177248,555[&type=0]:17.192557630177248)861:0.06937301039395294,486[&type=0]:17.2619306405712)866:1.4868173516881669)967:0.43482819534603934,((156[&type=0]:15.483531452132919,432[&type=0]:15.483531452132919)726:0.7145171399875814,49[&type=0]:16.1980485921205)776:2.985527595484907)992:0.24824469930974757,(430[&type=0]:17.17358679002178,(241[&type=0]:16.796862700432733,(228[&type=0]:16.08526949542913,211[&type=0]:16.08526949542913)766:0.7115932050036022)830:0.3767240895890467)856:2.258234096893375)1003:0.40077704605889153,(183[&type=0]:19.561354895235613,((30[&type=0]:15.766101351652749,506[&type=0]:15.766101351652749)739:1.6567005407979298,437[&type=0]:17.42280189245068)877:2.138553002784935)1017:0.27124303773843295)1027:1.4097887860661231,(21[&type=0]:14.216536627234365,(88[&type=0]:11.944341618959998,420[&type=0]:11.944341618959998)592:2.272195008274368)651:7.025850091805804)1069:0.5084965728634039)1078:0.5225618983154838)1092:0.41097069426879784,217[&type=0]:22.684415884487855)1097:0.07207482481518213,(((((197[&type=0]:16.63761587073067,((194[&type=0]:16.136282956731712,441[&type=0]:16.136282956731712)770:0.1079487055088002,(82[&type=0]:11.998680674931139,154[&type=0]:11.998680674931139)593:4.2455509873093735)778:0.39338420849015776)809:1.4016290870597743,179[&type=0]:18.039244957790444)927:2.0228323392589758,((((131[&type=0]:15.4906414729165,291[&type=0]:15.4906414729165)727:2.1349278895290738,(85[&type=0]:13.788494608150977,387[&type=0]:13.788494608150977)637:3.8370747542945978)895:0.4140021260790192,(284[&type=0]:17.398769224579866,((344[&type=0]:16.598084532010695,0[&type=0]:16.598084532010695)806:0.13824904804273075,(((411[&type=0]:11.25843150050642,400[&type=0]:11.25843150050642)575:3.0813658172764313,357[&type=0]:14.339797317782851)657:2.0783147178454566,5[&type=0]:16.418112035628308)790:0.31822154442511774)822:0.6624356445264397)876:0.6408022639447282)928:1.554360317557915,((38[&type=0]:12.273318117902761,247[&type=0]:12.273318117902761)598:4.168702028092056,(43[&type=0]:15.28630461218008,489[&type=0]:15.28630461218008)714:1.1557155338147371)792:3.1519116600876913)1019:0.4681454909669114)1041:1.6040577843344863,(((204[&type=0]:17.63223523334537,475[&type=0]:17.63223523334537)897:0.10222840896025787,(((58[&type=0]:12.404497744648138,239[&type=0]:12.404497744648138)604:2.677083216531697,(123[&type=0]:13.762236354537231,(404[&type=0]:13.305413617401825,312[&type=0]:13.305413617401825)623:0.4568227371354059)635:1.3193446066426038)705:2.0445280829925387,(522[&type=0]:16.674704261549074,(270[&type=0]:14.332331354882436,139[&type=0]:14.332331354882436)655:2.342372906666638)814:0.45140478262329964)853:0.6083545981332534)905:2.1773241602614775,327[&type=0]:19.911787802567105)1032:1.7543472788168017)1077:0.23140125636593112,((((268[&type=0]:13.067388882242042,298[&type=0]:13.067388882242042)614:3.6633727326250174,(1[&type=0]:14.559837928317421,(297[&type=0]:13.96497607495714,391[&type=0]:13.96497607495714)640:0.5948618533602819)666:2.1709236865496386)820:0.7302321499696305,403[&type=0]:17.46099376483669)882:3.380930212050135,(253[&type=0]:20.45085739058461,((((138[&type=0]:19.17009699461107,(((190[&type=0]:11.352141041730093,422[&type=0]:11.352141041730093)576:5.682291998990475,((351[&type=0]:14.56662089123583,35[&type=0]:14.56662089123583)667:2.0890953129415077,(398[&type=0]:13.590454649414033,303[&type=0]:13.590454649414033)632:3.065261554763305)810:0.3787168365432301)844:1.9689166961120605,300[&type=0]:19.00334973683263)981:0.16674725777844301)991:0.21521149120972893,(86[&type=0]:17.897699658384433,((396[&type=0]:15.557523618805838,213[&type=0]:15.557523618805838)731:1.1683733449929363,(480[&type=0]:12.774443262065624,536[&type=0]:12.774443262065624)606:3.9514537017331506)819:1.171802694585658)916:1.4876088274363681)1000:0.9499246109789468,((481[&type=0]:18.353959613463466,((212[&type=0]:11.533204372074898,146[&type=0]:11.533204372074898)584:4.64968982401237,(223[&type=0]:16.005092044754683,316[&type=0]:16.005092044754683)758:0.17780215133258537)773:2.171065417376198)942:0.29841710915374264,(17[&type=0]:17.712251809199667,164[&type=0]:17.712251809199667)903:0.9401249134175416)958:1.6828563741825384)1048:0.05368672257383622,(484[&type=0]:18.80303390886273,((71[&type=0]:16.800021545701142,(488[&type=0]:16.503124699022884,(216[&type=0]:15.258310517371084,371[&type=0]:15.258310517371084)712:1.2448141816517992)796:0.29689684667825844)832:0.012817746426005527,(((181[&type=0]:14.01859516335786,66[&type=0]:14.01859516335786)644:0.2853181606416353,(120[&type=0]:13.958734974221247,307[&type=0]:13.958734974221247)639:0.3451783497782479)654:0.8943434272604396,61[&type=0]:15.198256751259935)709:1.6145825408672128)834:1.9901946167355824)968:1.5858859105108536)1050:0.0619375712110255)1051:0.391066586302216)1061:1.0556123608630124)1083:0.8589543715531995)1100:0.7248796663362214)1104:0.6361401910012923)1108:0.4819218763584239,((((418[&type=0]:17.696636525026644,450[&type=0]:17.696636525026644)902:0.690714516578371,(13[&type=0]:17.882093449819124,((51[&type=0]:16.32400305580682,76[&type=0]:16.32400305580682)785:0.6369848978538286,74[&type=0]:16.96098795366065)841:0.921105496158475)914:0.5052575917858917)944:4.005416799726127,(((((((226[&type=0]:10.299525301365168,296[&type=0]:10.299525301365168)565:3.4172681391565263,177[&type=0]:13.716793440521695)634:4.066642788435679,(364[&type=0]:11.40712972507203,381[&type=0]:11.40712972507203)579:6.376306503885344)906:1.6574609979245132,((((188[&type=0]:14.682822654717475,(527[&type=0]:12.779400997203512,92[&type=0]:12.779400997203512)607:1.903421657513963)671:1.2602133322519808,259[&type=0]:15.943035986969456)754:1.6543606100971893,451[&type=0]:17.597396597066645)893:0.3229708488939771,(250[&type=0]:14.176418052120676,227[&type=0]:14.176418052120676)649:3.743949393839946)920:1.5205297809212652)1005:0.8802664249767957,496[&type=0]:20.321163651858683)1046:1.2246252749794522,((528[&type=0]:20.79185220955404,(202[&type=0]:16.731453313668062,(101[&type=0]:12.250650673876343,473[&type=0]:12.250650673876343)596:4.480802639791719)821:4.060398895885978)1060:0.09546117500081408,((((168[&type=0]:12.318260056436163,193[&type=0]:12.318260056436163)601:3.6235685779302713,(425[&type=0]:14.29985438476521,34[&type=0]:14.29985438476521)653:1.6419742496012244)753:0.7222839093009128,(471[&type=0]:16.606008867237197,(258[&type=0]:14.789123919698312,121[&type=0]:14.789123919698312)682:1.8168849475388846)807:0.05810367643015013)811:0.6695287354855708,((99[&type=0]:14.771076001969165,463[&type=0]:14.771076001969165)680:2.4073099878144344,45[&type=0]:17.1783859897836)859:0.1552552893693182)870:3.5536721054019367)1062:0.6584755422832806)1074:0.3483269098940731,(((136[&type=0]:11.381445699707486,170[&type=0]:11.381445699707486)578:4.446805061388071,(218[&type=0]:13.881417311513822,392[&type=0]:13.881417311513822)638:1.9468334495817352)742:3.2185868726842664,(109[&type=0]:10.857508455135132,192[&type=0]:10.857508455135132)568:8.189329178644691)986:2.847278202952385)1082:0.49865200459893444)1093:1.5728601821233674,(((((208[&type=0]:18.962848441888518,(174[&type=0]:15.779272107297793,(442[&type=0]:8.337003709926137,499[&type=0]:8.337003709926137)559:7.442268397371656)740:3.1835763345907253)980:0.435013862916918,(551[&type=0]:18.602947067393544,(((470[&type=0]:12.086400951142181,233[&type=0]:12.086400951142181)594:4.429683627967105,426[&type=0]:16.516084579109286)798:0.16224088675339132,215[&type=0]:16.678325465862677)815:1.9246216015308661)955:0.7949152374118924)1001:0.4282039261278605,((261[&type=0]:15.501643076335812,345[&type=0]:15.501643076335812)729:2.3997081366600934,(525[&type=0]:12.403865826497775,15[&type=0]:12.403865826497775)603:5.497485386498131)917:1.924715017937391)1026:2.08454950703511,(((37[&type=0]:19.531184221009944,((220[&type=0]:17.12095172290147,((325[&type=0]:16.515608108602784,163[&type=0]:16.515608108602784)797:0.20376228267188878,232[&type=0]:16.719370391274673)818:0.4015813316267973)852:0.11898748611473664,(144[&type=0]:16.853503884020448,(77[&type=0]:12.111839861181394,375[&type=0]:12.111839861181394)595:4.741664022839053)837:0.38643532499575883)865:2.2912450119937375)1015:1.2400344525615097,(369[&type=0]:18.913060149695546,((505[&type=0]:11.704401532292222,530[&type=0]:11.704401532292222)589:6.917362843684808,145[&type=0]:18.62176437597703)956:0.29129577371851667)976:1.8581585238759075)1058:1.0328839035025439,(((((((377[&type=0]:13.086977995654143,103[&type=0]:13.086977995654143)615:3.6528332445324754,((361[&type=0]:16.028357934446188,453[&type=0]:16.028357934446188)764:0.6418428567777603,((540[&type=0]:13.352613528750881,526[&type=0]:13.352613528750881)625:2.058879584384817,457[&type=0]:15.411493113135698)720:1.2587076780882498)812:0.06961044896267055)823:0.43461104273378837,518[&type=0]:17.174422282920407)857:1.4021339254620777,(454[&type=0]:17.07174879302443,286[&type=0]:17.07174879302443)849:1.504807415358055)951:0.9685361796292682,(((408[&type=0]:11.471966239197897,479[&type=0]:11.471966239197897)582:5.286121643937513,((495[&type=0]:10.538486953199763,276[&type=0]:10.538486953199763)566:5.186971973574028,320[&type=0]:15.725458926773792)737:1.0326289563616182)826:0.8984128157141669,3[&type=0]:17.656500698849577)900:1.8885916891621761)1016:0.9877830219422279,((((498[&type=0]:17.904266244504164,((265[&type=0]:12.810872458293971,(405[&type=0]:11.68106853926804,137[&type=0]:11.68106853926804)587:1.1298039190259317)608:2.063120391773662,87[&type=0]:14.873992850067633)692:3.0302733944365308)918:0.04315541069703954,57[&type=0]:17.947421655201204)921:0.5818201829146226,424[&type=0]:18.529241838115826)949:1.800398804831513,((359[&type=0]:17.863511028088734,(384[&type=0]:13.113944545651409,75[&type=0]:13.113944545651409)616:4.749566482437325)911:1.642304177949562,(((368[&type=0]:12.284775166774303,(409[&type=0]:11.156362786784866,274[&type=0]:11.156362786784866)571:1.1284123799894363)599:0.8555916625485072,385[&type=0]:13.14036682932281)618:5.497961916231526,(249[&type=0]:14.08712142680258,84[&type=0]:14.08712142680258)647:4.551207318751755)957:0.8674864604839598)1012:0.8238254369090434)1047:0.20323476700664145)1057:0.9289556948318385,(148[&type=0]:15.60252036930112,243[&type=0]:15.60252036930112)733:5.8593107354847)1072:0.34227147228817856)1081:0.10651316089440854)1084:1.4317695995426583,((((334[&type=0]:19.484010776833227,455[&type=0]:19.484010776833227)1010:0.7718827376575916,(((379[&type=0]:19.436361808494823,373[&type=0]:19.436361808494823)1004:0.22461614781415662,(62[&type=0]:18.93134748827004,(263[&type=0]:17.546663452775825,492[&type=0]:17.546663452775825)886:1.3846840354942138)978:0.7296304680389412)1022:0.3008071469786948,(((339[&type=0]:17.137656313237656,554[&type=0]:17.137656313237656)854:1.5957142170490002,443[&type=0]:18.733370530286656)966:0.73112038382083,(((8[&type=0]:17.650779811049162,(497[&type=0]:16.054502418533353,393[&type=0]:16.054502418533353)765:1.596277392515809)899:0.49781592796968965,(501[&type=0]:15.0978004835091,478[&type=0]:15.0978004835091)706:3.050795255509751)933:1.1616211368359828,248[&type=0]:19.310216875854834)996:0.15427403825265174)1008:0.4972941891801881)1037:0.294108411203144)1045:1.5004337864235602,(((421[&type=0]:16.145091503011308,201[&type=0]:16.145091503011308)771:3.777119022017814,((122[&type=0]:18.714607685697413,(189[&type=0]:15.979812520298495,235[&type=0]:15.979812520298495)756:2.7347951653989178)965:0.18168744204223586,((4[&type=0]:18.274768586173998,((383[&type=0]:11.424213159090474,225[&type=0]:11.424213159090474)580:1.4074705617987053,254[&type=0]:12.83168372088918)609:5.443084865284819)939:0.04372864096163909,111[&type=0]:18.318497227135637)940:0.5777979006040113)975:1.0259153972894737)1033:0.2929998708730608,((((459[&type=0]:16.367121378514724,(10[&type=0]:16.286633678493143,(509[&type=0]:13.170385423126246,(493[&type=0]:12.269174606895751,151[&type=0]:12.269174606895751)597:0.9012108162304955)620:3.116248255366896)782:0.08048770002158179)789:2.150671391865366,((((231[&type=0]:14.742186277783128,510[&type=0]:14.742186277783128)678:0.31107431078451775,282[&type=0]:15.053260588567646)702:1.4073079283414884,(((466[&type=0]:13.63628517409881,169[&type=0]:13.63628517409881)633:1.7883024398796667,(173[&type=0]:14.335114550811824,47[&type=0]:14.335114550811824)656:1.0894730631666523)721:0.7342868236887039,535[&type=0]:16.15887443766718)772:0.3016940792419547)794:0.8492968736175541,363[&type=0]:17.30986539052669)868:1.2079273798534018)948:0.5275727869995137,(52[&type=0]:17.456839937680215,(128[&type=0]:16.260585995762135,328[&type=0]:16.260585995762135)779:1.19625394191808)881:1.5885256196993893)985:0.40042562294816975,(323[&type=0]:16.76147094939798,438[&type=0]:16.76147094939798)827:2.684320230929792)1007:0.7694192155744091)1044:1.5411169050121956)1079:0.2489692094960141,(((319[&type=0]:16.71408420467101,108[&type=0]:16.71408420467101)817:0.7492552845357388,((195[&type=0]:9.82281142054511,205[&type=0]:9.82281142054511)562:4.149549726706844,458[&type=0]:13.972361147251954)641:3.4909783419547953)883:2.5162784864232215,(293[&type=0]:17.394610832966066,((279[&type=0]:14.433322940944237,538[&type=0]:14.433322940944237)661:0.10761547571045149,237[&type=0]:14.540938416654688)665:2.853672416311378)875:2.5850071426639047)1039:2.025678534780422)1089:1.337088827100672)1102:0.6232426859434455)1107:0.6338044195444645)1109:0.24644614949909283,(((((((((456[&type=0]:16.420884156143444,515[&type=0]:16.420884156143444)791:2.247762739389639,221[&type=0]:18.668646895533083)960:0.39971598341578485,244[&type=0]:19.068362878948868)987:1.4389395567777328,((((26[&type=0]:13.366297035598688,504[&type=0]:13.366297035598688)627:0.3995537221213574,102[&type=0]:13.765850757720045)636:1.3502190725361185,80[&type=0]:15.116069830256164)708:5.251481230057069,((548[&type=0]:10.60422086551132,54[&type=0]:10.60422086551132)567:8.58548133961434,((7[&type=0]:12.963761247184797,106[&type=0]:12.963761247184797)612:5.471791282243645,330[&type=0]:18.435552529428442)946:0.7541496756972172)993:1.1778488551875732)1049:0.13975137541336835)1055:0.5869855314795309,(311[&type=0]:20.176346455214713,185[&type=0]:20.176346455214713)1043:0.9179415119914189)1066:1.383298841600336,((((114[&type=0]:17.311678100901236,((547[&type=0]:15.013978402043822,115[&type=0]:15.013978402043822)699:0.8181708467007613,380[&type=0]:15.832149248744583)744:1.479528852156653)869:0.27291455761033845,(512[&type=0]:17.202502668816976,((444[&type=0]:16.357646384036606,((56[&type=0]:14.572878102428158,104[&type=0]:14.572878102428158)668:0.9381474078473353,203[&type=0]:15.511025510275493)730:0.846620873761113)787:0.3960297081410822,402[&type=0]:16.753676092177688)824:0.4488265766392878)862:0.38208998969459884)891:0.587095083908725,(395[&type=0]:17.587597234459643,((494[&type=0]:15.401623605058301,317[&type=0]:15.401623605058301)719:0.8759645360011739,((461[&type=0]:13.293615517490968,280[&type=0]:13.293615517490968)622:1.374819519457498,44[&type=0]:14.668435036948466)670:1.6091531041110088)780:1.3100090934001685)892:0.5840905079606564)934:1.7066595650167287,((445[&type=0]:17.037236810428734,((((132[&type=0]:13.499368584870659,477[&type=0]:13.499368584870659)630:1.2676430920019417,96[&type=0]:14.7670116768726)679:0.24518124462852064,423[&type=0]:15.012192921501121)698:1.468868946273254,289[&type=0]:16.481061867774375)795:0.5561749426543585)845:1.8792207038606108,390[&type=0]:18.916457514289345)977:0.9618897931476837)1030:2.599239501369439)1094:0.14496045110080047,((((367[&type=0]:17.919449904067456,(((176[&type=0]:2.75640486664868,94[&type=0]:2.75640486664868)556:11.281450797391766,33[&type=0]:14.037855664040446)646:1.3870242612176078,((413[&type=0]:11.083521722711591,32[&type=0]:11.083521722711591)569:0.08808567376478749,546[&type=0]:11.171607396476379)572:4.253272528781675)722:2.494569978809402)919:0.9494670381934398,346[&type=0]:18.868916942260896)971:0.6609268393149961,(200[&type=0]:15.425718257668002,(255[&type=0]:11.633371355388979,251[&type=0]:11.633371355388979)586:3.7923469022790233)723:4.104125523907889)1014:2.40267912013703,((((260[&type=0]:15.468042429809527,500[&type=0]:15.468042429809527)725:1.9166636061995792,427[&type=0]:17.384706036009106)872:0.19626342272690778,434[&type=0]:17.580969458736014)890:1.2304602555505433,(((214[&type=0]:11.249388006564828,65[&type=0]:11.249388006564828)574:1.9013994172344226,160[&type=0]:13.150787423799251)619:3.0406213194448455,314[&type=0]:16.191408743244097)775:2.620020971042461)970:3.1210931874263643)1085:0.6900243581943464)1096:0.11825888603138779,(((91[&type=0]:19.024276411280347,365[&type=0]:19.024276411280347)984:1.924815407605731,((((410[&type=0]:16.98349644182141,73[&type=0]:16.98349644182141)843:2.580668860122362,((222[&type=0]:16.79392374282758,415[&type=0]:16.79392374282758)829:2.3434749440833897,((549[&type=0]:17.28511037579958,315[&type=0]:17.28511037579958)867:1.7203967826939213,((350[&type=0]:14.973618244035213,(533[&type=0]:14.707878283227192,230[&type=0]:14.707878283227192)674:0.26573996080802154)696:2.0083851171715583,414[&type=0]:16.982003361206772)842:2.023503797286729)982:0.13189152841746932)990:0.42676661503280044)1018:0.29355813947372766,(448[&type=0]:15.98411387892468,(324[&type=0]:9.83733179669723,81[&type=0]:9.83733179669723)563:6.14678208222745)757:3.8736095624928186)1029:0.9168783111044796,((196[&type=0]:15.039307094441492,147[&type=0]:15.039307094441492)701:1.2465773950863532,(162[&type=0]:11.797878178499774,112[&type=0]:11.797878178499774)590:4.488006311028071)781:4.4887172629941325)1059:0.17449006636410047)1063:1.7829042429814699,(((((((27[&type=0]:14.027142223941134,107[&type=0]:14.027142223941134)645:1.866448352647053,490[&type=0]:15.893590576588187)751:4.596177400614648,(9[&type=0]:18.143147955310983,((((322[&type=0]:6.302652471485192,267[&type=0]:6.302652471485192)557:5.535011764308619,257[&type=0]:11.837664235793811)591:5.985575244086823,287[&type=0]:17.823239479880634)908:0.21287321488665967,290[&type=0]:18.036112694767294)926:0.1070352605436895)932:2.3466200218918516)1054:0.5884809622556553,(((172[&type=0]:19.30493011381224,((452[&type=0]:16.027646605418326,245[&type=0]:16.027646605418326)763:1.3991678505943952,(440[&type=0]:14.810600498092386,142[&type=0]:14.810600498092386)686:2.616213957920335)878:1.878115657799519)995:0.20875287146798627,((25[&type=0]:15.20303230673458,349[&type=0]:15.20303230673458)710:2.90773170774378,397[&type=0]:18.11076401447836)931:1.402918970801867)1013:1.499210968438291,(((541[&type=0]:14.122737703825418,20[&type=0]:14.122737703825418)648:5.198639376404966,256[&type=0]:19.321377080230384)997:0.642084864583353,(((552[&type=0]:14.814452567813703,332[&type=0]:14.814452567813703)687:3.7395624713355424,(((360[&type=0]:14.79291612067981,272[&type=0]:14.79291612067981)683:2.2647192523209974,288[&type=0]:17.057635373000807)847:1.159845286490654,(166[&type=0]:17.481303402516822,(158[&type=0]:11.690039742426864,491[&type=0]:11.690039742426864)588:5.7912636600899585)884:0.7361772569746385)937:0.33653437965778465)950:0.31494975654343094,(((((117[&type=0]:15.080310754047124,((308[&type=0]:9.08311699406263,366[&type=0]:9.08311699406263)560:5.888907498974106,412[&type=0]:14.972024493036736)695:0.10828626101038807)704:1.1094351945269327,90[&type=0]:16.189745948574057)774:0.57431275504225,264[&type=0]:16.764058703616307)828:0.877817379130228,(97[&type=0]:16.35117429543616,485[&type=0]:16.35117429543616)786:1.2907017873103754)898:0.5696425374780283,(407[&type=0]:14.445179112029859,353[&type=0]:14.445179112029859)662:3.7663395081947044)936:0.6574461754681131)972:1.0944971491210609)1038:1.0494320089047804)1064:0.06535498573997245)1065:0.8604643796258884,(((72[&type=0]:17.088522790081434,(206[&type=0]:14.640200649237538,464[&type=0]:14.640200649237538)669:2.4483221408438958)850:1.6132956451491651,175[&type=0]:18.7018184352306)964:0.6641724304567411,(406[&type=0]:17.86971578915073,143[&type=0]:17.86971578915073)913:1.4962750765366089)999:2.572722453397038)1086:0.03566397894627116,(((((53[&type=0]:16.516346670119603,419[&type=0]:16.516346670119603)799:1.033400476789076,352[&type=0]:17.54974714690868)888:0.41905017318998006,(((55[&type=0]:11.509719455013212,372[&type=0]:11.509719455013212)583:5.712585969726,48[&type=0]:17.22230542473921)863:0.1407709337103995,28[&type=0]:17.36307635844961)871:0.6057209616490482)923:1.864690644404913,(((100[&type=0]:14.379272652028652,309[&type=0]:14.379272652028652)660:3.0093826490924744,465[&type=0]:17.388655301121126)873:0.7012536502788791,321[&type=0]:18.089908951400005)930:1.7435790131035667)1028:1.6922037020841998,(511[&type=0]:14.252422705354169,516[&type=0]:14.252422705354169)652:7.273268961233603)1073:0.44868563144287776)1087:0.2638765777568892,((178[&type=0]:15.495046610870862,462[&type=0]:15.495046610870862)728:2.8590280272520374,277[&type=0]:18.3540746381229)943:3.8841792376646396)1091:0.4937421860800093)1098:0.008810084071107838)1099:0.9227091499599176,(342[&type=0]:16.290322529364126,382[&type=0]:16.290322529364126)783:7.373192766534448)1106:1.182363296599494)1110:0.0";
        //init trees
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newickTree,
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.0");
        RealParameter scarRates = new RealParameter(" 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.025 0.025 0.025 0.025");
        RealParameter scarringHeight = new RealParameter("25.0");
        RealParameter scarringDuration = new RealParameter("25.0");

        RealParameter freqs = new RealParameter("1.0 0 0 0 0" +
                "0 0 0 0 0 " +
                "0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        EditAndSilencingModel scarringModel = new EditAndSilencingModel();
        scarringModel.initByName("editRates", scarRates,
                "silencingRate", lossRate,
                "editHeight", scarringHeight,
                "editDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        RealParameter shape = new RealParameter("1.0");
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 2, "substModel", scarringModel,
                "shape", shape);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihoodNegInf = new TreeLikelihoodWithEditWindow();
        likelihoodNegInf.initByName("data", alignment, "tree", tree,
                "siteModel", siteM, "branchRateModel", clockModel, "useScaling", true);

        Node parent = tree.getRoot();
        Node node_7 = parent.getChild(0);
        Node node_6 = parent.getChild(1);
        Node node_5 = node_7.getChild(1);

        //test partials at node 5
        likelihoodNegInf.traverse(tree.getRoot());
        double [] partials = new double[16];
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_5.getNr(), partials);
        //assertArrayEquals("Assert correct likelihood at internal node 5:", partials,
        //        new double[]{0, 0.836331904590560, 0, 0, 0, 0, 0.836331904590560, 0, 0, 0.422624893321294, 0, 0, 0,0, 0.422624893321294,0}, 1e-15);

        //test partials at node 6
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_6.getNr(), partials);
        //assertArrayEquals("Assert correct likelihood at internal node 6:", partials,
        //        new double[]{0.032054626428896, 0, 0.128958931608947, 0, 0.032054626428896, 0.128958931608947 ,0,0,
        //                0.119317341270894, 0, 0.247654804887413, 0, 0.119317341270894, 0.247654804887413 ,0,0}, 1e-15);

        //test partials at node 7
        likelihoodNegInf.getLikelihoodCore().getNodePartials(node_7.getNr(), partials);
        //assertArrayEquals("Assert correct likelihood at internal node 7:", partials,
        //        new double[]{0.060425520687175, 0, 0, 0, 0.060425520687175, 0, 0, 0,
        //                0.769152619166653e-4, 0, 0, 0, 0.769152619166653e-4, 0, 0, 0}, 1e-15);

        //test partials at root node
        likelihoodNegInf.getLikelihoodCore().getNodePartials(parent.getNr(), partials);
        //assertArrayEquals("Assert correct likelihood at root:", partials,
        //        new double[]{0.002643849069724, 0, 0, 0, 0.002643849069724, 0, 0, 0,
        //                0.598196646044564e-5, 0, 0, 0, 0.598196646044564e-5, 0, 0, 0}, 1e-15);


        double logP = likelihoodNegInf.calculateLogP();

        assertEquals(-13.252813163014444, logP, 1e-14);
    }
}

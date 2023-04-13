package lineageTree.distributions;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import lineageTree.substitutionmodel.EditAndSilencingModel;
import org.junit.Test;

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
        alignment.initByName("sequence", a, "dataType", "scarData", "stateCount", 4);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel3.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel3.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel4.initByName("scarringRates", new RealParameter("0.01 0.01"),
                "lossRate", new RealParameter("0.01"),
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel4.initByName("scarringRates", new RealParameter("0.01 0.01"),
                "lossRate", new RealParameter("0.01"),
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel4b.initByName("scarringRates", new RealParameter("0.01 0.01"),
                "lossRate", new RealParameter("0.01"),
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

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
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", scarringHeight,
                "scarringDuration", scarringDuration, "frequencies", frequencies);

        // init site model
        RealParameter shape = new RealParameter("1.0");
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 2, "substModel", scarringModel,
                "shape", shape);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihoodNegInf = new TreeLikelihoodWithEditWindow();
        likelihoodNegInf.initByName("data", alignment, "tree", tree,
                "siteModel", siteM, "branchRateModel", clockModel, "origin", origin);

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
}

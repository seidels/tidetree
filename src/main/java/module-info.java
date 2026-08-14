open module tidetree {
    requires beast.pkgmgmt;
    requires beast.base;
    requires feast;
    requires static beast.fx;
    requires static javafx.controls;

    exports tidetree.distributions;
    exports tidetree.evolution.datatype;
    exports tidetree.simulation;
    exports tidetree.substitutionmodel;
    exports tidetree.tree;
    exports tidetree.util;
    exports tidetree.app.beauti;


    provides beast.base.core.BEASTInterface with
	tidetree.evolution.datatype.EditData,
	tidetree.substitutionmodel.EditAndSilencingModel,
	tidetree.tree.StartingTree,
	tidetree.distributions.TreeLikelihoodWithEditWindow,
	tidetree.simulation.SimulatedAlignment,
	tidetree.util.AlignmentFromNexus;

      provides beastfx.app.inputeditor.AlignmentImporter with
          tidetree.app.beauti.NexusImporter;
}

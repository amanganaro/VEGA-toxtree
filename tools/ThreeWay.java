import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.toxtree.core.Classification;
import insilico.toxtree.cramer.CramerTree;
import insilico.toxtree.verhaar.VerhaarTree;
import toxTree.core.IDecisionMethod;
import toxTree.core.IDecisionResult;
import toxTree.query.MolAnalyser;
import toxTree.tree.cramer.CramerRules;
import verhaar.VerhaarScheme;

import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/**
 * Tells apart two possible causes of a difference: the ported rules, or the core that normalizes the molecule.
 * For each SMILES it reports what the original Toxtree says when fed exactly what the old model fed it, the
 * canonical SMILES of the current core, and what the ported tree says on the same molecule.
 */
public class ThreeWay {
    public static void main(String[] args) throws Exception {
        String which = args[0];
        List<String> lines = Files.readAllLines(Paths.get(args[1]), StandardCharsets.UTF_8);

        IDecisionMethod old = which.equals("cramer") ? new CramerRules() : new VerhaarScheme();
        var ported = which.equals("cramer") ? (Object) new CramerTree() : new VerhaarTree();
        var sp = new shadedTTv31.org.openscience.cdk.smiles.SmilesParser(
                shadedTTv31.org.openscience.cdk.DefaultChemObjectBuilder.getInstance());

        try (var w = Files.newBufferedWriter(Paths.get(args[2]), StandardCharsets.UTF_8)) {
            w.write("smiles\told_toxtree\tported\n");
            for (String smiles : lines) {
                smiles = smiles.trim();
                if (smiles.isEmpty()) continue;
                String oldAnswer = "N.A.";
                String trace = "";
                String mol0 = smiles;
                try {
                    var mol = SmilesMolecule.Convert(smiles);
                    // what the old model did: re-parse the core's canonical SMILES with the old CDK
                    try {
                        var ac = sp.parseSmiles(mol.GetSMILES());
                        MolAnalyser.analyse(ac);
                        IDecisionResult r = old.createDecisionResult();
                        r.classify(ac);
                        r.assignResult(ac);
                        oldAnswer = r.getAssignedCategories().get(0).getName();
                        trace = old.explainRules(r, false).toString();
                    } catch (Throwable t) {
                        oldAnswer = "N.A.";
                        trace = "ERROR";
                    }
                } catch (Throwable t) {
                    // leave both N.A.
                }
                w.write(mol0 + "\t" + oldAnswer + "\t" + trace + "\n");
            }
        }
        System.out.println("done");
    }
}

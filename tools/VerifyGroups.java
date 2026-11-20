import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.toxtree.chem.Smarts;
import shadedTTv31.org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import shadedTTv31.org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import toxTree.query.FunctionalGroups;
import toxTree.query.MolAnalyser;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Checks each functional group of the port against the one Toxtree builds, atom for atom.
 * <p>
 * Both sides are given the same canonical SMILES, so the heavy atoms carry the same positions on both. For every
 * molecule and every group the two sets of matched heavy atoms are compared. This catches a SMARTS that reads
 * the name of a group rather than its shape, even on molecules where the difference changes no classification.
 */
public class VerifyGroups {

    private record Group(String name, QueryAtomContainer query, String smarts) {
    }

    public static void main(String[] args) throws Exception {
        List<String> smilesList = Files.readAllLines(Paths.get(args[0]), StandardCharsets.UTF_8);
        int limit = args.length > 1 ? Integer.parseInt(args[1]) : smilesList.size();

        List<Group> groups = groups();
        Map<String, int[]> stats = new LinkedHashMap<>();   // checked, differing
        Map<String, List<String>> examples = new LinkedHashMap<>();
        for (Group g : groups) {
            stats.put(g.name(), new int[5]);
            examples.put(g.name(), new ArrayList<>());
        }

        var sp = new shadedTTv31.org.openscience.cdk.smiles.SmilesParser(
                shadedTTv31.org.openscience.cdk.DefaultChemObjectBuilder.getInstance());
        UniversalIsomorphismTester uit = new UniversalIsomorphismTester();

        int done = 0;
        for (String raw : smilesList) {
            if (done >= limit) {
                break;
            }
            String smiles = raw.trim();
            if (smiles.isEmpty()) {
                continue;
            }
            String canonical;
            org.openscience.cdk.interfaces.IAtomContainer mine;
            org.openscience.cdk.interfaces.IAtomContainer kekule;
            shadedTTv31.org.openscience.cdk.interfaces.IAtomContainer theirs;
            try {
                var v = SmilesMolecule.Convert(smiles);
                canonical = v.GetSMILES();
                mine = SmilesMolecule.CreateCDKMolecule(canonical, null);
                org.openscience.cdk.tools.manipulator.AtomContainerManipulator
                        .convertImplicitToExplicitHydrogens(mine);
                theirs = sp.parseSmiles(canonical);
                MolAnalyser.analyse(theirs);
                // Toxtree marks the CHn carbons before matching any group that asks for one.
                FunctionalGroups.preProcess(theirs);
                kekule = mine.clone();
                for (var a : kekule.atoms()) {
                    a.setIsAromatic(false);
                }
                for (var b : kekule.bonds()) {
                    b.setIsAromatic(false);
                }
            } catch (Throwable t) {
                continue;
            }
            done++;

            for (Group g : groups) {
                Set<String> aromaticView;
                Set<String> kekuleView;
                Set<String> toxtree;
                try {
                    aromaticView = mySets(mine, g.smarts());
                    kekuleView = mySets(kekule, g.smarts());
                    toxtree = theirSets(uit, theirs, g.query());
                } catch (Throwable t) {
                    continue;
                }
                int[] s = stats.get(g.name());
                s[0]++;
                if (aromaticView.equals(toxtree)) {
                    s[1]++;
                }
                if (kekuleView.equals(toxtree)) {
                    s[2]++;
                }
                boolean presentToxtree;
                try {
                    presentToxtree = uit.isSubgraph(theirs, g.query());
                } catch (Throwable t) {
                    presentToxtree = !toxtree.isEmpty();
                }
                if (presentToxtree == !aromaticView.isEmpty() || presentToxtree == !kekuleView.isEmpty()) {
                    s[4]++;
                }
                if (!aromaticView.equals(toxtree) && !kekuleView.equals(toxtree)) {
                    s[3]++;
                    List<String> ex = examples.get(g.name());
                    if (ex.size() < 3) {
                        ex.add(String.format("      %-50s toxtree=%s mine=%s",
                                canonical.length() > 48 ? canonical.substring(0, 48) + ".." : canonical,
                                toxtree, aromaticView));
                    }
                }
            }
        }

        System.out.printf("%n=== functional groups checked atom for atom on %d molecules%n%n", done);
        System.out.printf("%-42s %7s %9s %9s %8s %9s%n", "group", "checked", "aromatic", "kekule", "neither", "presence");
        for (Group g : groups) {
            int[] s = stats.get(g.name());
            String view = s[3] > 0 ? "   <-- differs"
                    : (s[1] == s[0] ? "   aromatic view" : "   KEKULE view only");
            System.out.printf("%-42s %7d %9d %9d %8d %9d%s%n", g.name(), s[0], s[1], s[2], s[3], s[4], view);
            examples.get(g.name()).forEach(System.out::println);
        }
    }

    /** The heavy atoms of each match, as sorted position lists. */
    private static Set<String> mySets(org.openscience.cdk.interfaces.IAtomContainer mol, String smarts) {
        Set<String> sets = new TreeSet<>();
        for (int[] match : Smarts.of(smarts).matchAll(mol).uniqueAtoms()) {
            Set<Integer> heavy = new TreeSet<>();
            for (int i : match) {
                if (!"H".equals(mol.getAtom(i).getSymbol())) {
                    heavy.add(i);
                }
            }
            if (!heavy.isEmpty()) {
                sets.add(heavy.toString());
            }
        }
        return sets;
    }

    private static Set<String> theirSets(UniversalIsomorphismTester uit,
                                         shadedTTv31.org.openscience.cdk.interfaces.IAtomContainer mol,
                                         QueryAtomContainer query) throws Exception {
        Set<String> sets = new TreeSet<>();
        List<?> maps = uit.getSubgraphAtomsMaps(mol, query);
        for (Object o : maps) {
            Set<Integer> heavy = new TreeSet<>();
            for (Object r : (List<?>) o) {
                int index = ((shadedTTv31.org.openscience.cdk.isomorphism.mcss.RMap) r).getId1();
                if (!"H".equals(mol.getAtom(index).getSymbol())) {
                    heavy.add(index);
                }
            }
            if (!heavy.isEmpty()) {
                sets.add(heavy.toString());
            }
        }
        return sets;
    }

    private static List<Group> groups() throws Exception {
        String chn = Smarts.CHN;
        String metals = "Na,K,Ca,Mg,Al";
        List<Group> g = new ArrayList<>();
        g.add(new Group("alcohol", FunctionalGroups.alcohol(false), chn + "-[#8]-[#1]"));
        g.add(new Group("aldehyde", FunctionalGroups.aldehyde(), "[#6,#1]-[#6]([#1])=[#8]"));
        g.add(new Group("acetal", FunctionalGroups.acetal(), "[#6]-[#8]-[#6]-[#8]-[#6]"));
        g.add(new Group("acyclic_acetal", FunctionalGroups.acyclic_acetal(),
                "[#6]-!@[#8]-!@[#6]-!@[#8]-!@[#6]"));
        g.add(new Group("ketone", FunctionalGroups.ketone(),
                "[!#1;!#8;!#16;!#7]-[#6](=[#8])-[!#1;!#8;!#16;!#7]"));
        g.add(new Group("ketone_ring", FunctionalGroups.ketone_ring(), "[#6;R]-[#6;R](=[#8])-[#6;R]"));
        g.add(new Group("sidechain_ketone", FunctionalGroups.sidechain_ketone(),
                "[#6]-[#6;!R](=[#8;!R])-[#6]"));
        g.add(new Group("carboxylicAcid", FunctionalGroups.carboxylicAcid(), "[*]-[#6](-[#8]-[#1])=[#8]"));
        g.add(new Group("ester", FunctionalGroups.ester(), "[*]-[#6](-[#8]-[!#1])=[#8]"));
        g.add(new Group("mercaptan", FunctionalGroups.mercaptan(), "[#6]-[#16]-[#1]"));
        g.add(new Group("sulphide", FunctionalGroups.sulphide(), "[#6]-[#16]-[#6]"));
        g.add(new Group("methylether", FunctionalGroups.methylether(), "[#8](-[#6]([#1])([#1])[#1])-[#6]"));
        g.add(new Group("hydroxy_ring", FunctionalGroups.hydroxy_ring(), "[#6;R]-[#8]-[#1]"));
        g.add(new Group("hydroxy1", FunctionalGroups.hydroxy1(), "[#6;!a]-[#8]-[#1]"));
        g.add(new Group("ether", FunctionalGroups.ether(), chn + "-[#8]-" + chn));
        g.add(new Group("alkoxy", FunctionalGroups.alkoxy(), "[#6]-[#8]-[#6]"));
        g.add(new Group("carbonate", FunctionalGroups.carbonate(), "[*]-[#8]-[#6](-[#8]-[!#1])=[#8]"));
        g.add(new Group("thioester", FunctionalGroups.thioester(),
                "[!#1;!#16;!#8]-[#6](-[#16]-[!#1])=[#8,#16]"));
        g.add(new Group("lactone(false)", FunctionalGroups.lactone(false), "[#6]@[#6](-@[#8]@[#6])=!@[#8]"));
        g.add(new Group("anhydride", FunctionalGroups.anhydride(), "[#6]@[#6](-@[#8]@[#6]=!@[#8])=!@[#8]"));
        g.add(new Group("isopreneUnit", FunctionalGroups.isopreneUnit(), "[#6]-[#6](-[#6])-[#6]-[#6]"));
        g.add(new Group("primaryAmine(false)", FunctionalGroups.primaryAmine(false), "[#6]-[#7](-[#1])-[#1]"));
        g.add(new Group("secondaryAmine(false)", FunctionalGroups.secondaryAmine(false),
                "[#6]-[#7+0](-[#6])-[#1]"));
        g.add(new Group("secondaryAmine(true)", FunctionalGroups.secondaryAmine(true),
                chn + "-!@[#7+0](-!@" + chn + ")-[#1]"));
        g.add(new Group("tertiaryAmine", FunctionalGroups.tertiaryAmine(), "[!#1]-[#7](-[!#1])-[!#1]"));
        g.add(new Group("nitro1double", FunctionalGroups.nitro1double(), "[!#1]-[#7+](=[#8])-[#8-]"));
        g.add(new Group("nitro2double", FunctionalGroups.nitro2double(), "[!#1]-[#7](=[#8])=[#8]"));
        g.add(new Group("cyano", FunctionalGroups.cyano(), "[#7,#6,#1]-[#6]#[#7]"));
        g.add(new Group("Nnitroso", FunctionalGroups.Nnitroso(), "[!#1]-[#7](-[!#1])-[#7]=[#8]"));
        g.add(new Group("diAzo", FunctionalGroups.diAzo(), "[#7]#[#7]"));
        g.add(new Group("triAzeno", FunctionalGroups.triAzeno(), "[#6]-[#7]=[#7]-[#7](-[#1])-[#1]"));
        g.add(new Group("quaternaryNitrogen1(true)", FunctionalGroups.quaternaryNitrogen1(true),
                "[*]-[#7+](-[*])(-[*])-[*]"));
        g.add(new Group("quaternaryNitrogen2(true)", FunctionalGroups.quaternaryNitrogen2(true),
                "[*]-[#7+](-[*])=[*]"));
        g.add(new Group("quarternaryNitrogenException", FunctionalGroups.quarternaryNitrogenException(),
                "[*]-[#7+](-[*])=[#6](-[!#1])-[!#1]"));
        g.add(new Group("vicinalDiKetone", FunctionalGroups.vicinalDiKetone(),
                "[#8]=[#6](-[!#8;!#16])-[#6](=[#8])-[!#8;!#16]"));
        g.add(new Group("ketoneAttachedToTerminalVinyl", FunctionalGroups.ketoneAttachedToTerminalVinyl(),
                "[!#8;!#16;!#7]-[#6](=[#8])-[#6]([#1])=[#6]([#1])[#1]"));
        g.add(new Group("ketalAttachedToTerminalVinyl", FunctionalGroups.ketalAttachedToTerminalVinyl(),
                "[!#8;!#1]-[#6](-[#8]-[!#8;!#1])(-[#8]-[!#8;!#1])-[#6]=[#6]([#1])[#1]"));
        g.add(new Group("alcoholSecondary...Vinyl", FunctionalGroups.alcoholSecondaryAttachedToTerminalVinyl(),
                "[#1]-[#8]-[#6](-[#1])(-[#6])-[#6](-[#1])-[#6](-[#1])-[#1]"));
        g.add(new Group("esterOfAlcoholSecondary...", FunctionalGroups.esterOfalcoholSecondaryAttachedToTerminalVinyl(),
                "[#6]-[#6](-[#6](=[#8])-[#8]-[#6])-[#6]=[#6]"));
        g.add(new Group("allylAlcoholAcetal", FunctionalGroups.allylAlcoholAcetal(),
                "[#6](-[#8]-[#6]-[#6]=[#6])(-[#8]-[#6]-[#6]=[#6])-[#1]"));
        g.add(new Group("allylAlcoholEsterDerivative", FunctionalGroups.allylAlcoholEsterDerivative(),
                "[#8]=[#6](-[#6])-[#8]-[#6]([#1])([#1])-[#6]([#1])=[#6]([#1])[#1]"));
        g.add(new Group("stericallyHindered", FunctionalGroups.stericallyHindered(),
                "[#6]~[#6]~[#6](-[#6](-[!#1])(-[!#1])-[!#1])~[#6](~[!#1])~[#6](-[#6](-[!#1])(-[!#1])-[!#1])~[#6]"));
        g.add(new Group("phenol", verhaar.query.FunctionalGroups.phenol(),
                "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-,:[#8]-[#1]"));
        g.add(new Group("aniline", verhaar.query.FunctionalGroups.aniline(),
                "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-,:[#7]-[#1]"));
        g.add(new Group("pyridine_character", verhaar.query.FunctionalGroups.pyridine_character(),
                "[#7]1:[#6,#7]:[#6,#7]:[#6,#7]:[#6,#7]:[#6,#7]:1"));
        g.add(new Group("halogen(Cl,Br,F)", verhaar.query.FunctionalGroups.halogen(new String[]{"Cl", "Br", "F"}),
                "[*]~[Cl,Br,F]"));
        g.add(new Group("halogenAtBetaFromUnsaturation",
                verhaar.query.FunctionalGroups.halogenAtBetaFromUnsaturation(new String[]{"Cl", "Br", "F", "I"}),
                "[*]=,#,:[*]-[*]-[F,Cl,Br,I]"));
        g.add(new Group("ringSubstituted(Cl)", verhaar.query.FunctionalGroups.ringSubstituted("Cl"),
                "[*]@[#6](@[*])!@[Cl]"));
        g.add(new Group("saltOfCarboxylicAcid", FunctionalGroups.saltOfCarboxylicAcid(metals.split(",")),
                "[#6,#1]-[#6](-[#8]-[" + metals + "])=[#8]"));
        String me = "Na,K,Ca";
        String[] meArr = {"Na", "K", "Ca"};
        g.add(new Group("alcohol(true)", FunctionalGroups.alcohol(true),
                "[$([#6;!a]~[#1]),$([#6;!a;!$([#6]~[!#6])])]-[#8]-[#1]"));
        g.add(new Group("lactone(true)", FunctionalGroups.lactone(true),
                "[#6]=@[#6]@[#6](-@[#8]@[#6])=!@[#8]"));
        g.add(new Group("methoxy_ring", FunctionalGroups.methoxy_ring(),
                "[*]-[#8]-[#6]([#1])([#1])[#1]"));
        g.add(new Group("sulphonate(me)", FunctionalGroups.sulphonate(meArr),
                "[!#8;!#1;!#7]-[#16](=[#8])(=[#8])-[#8]-[" + me + "]"));
        g.add(new Group("sulphonate(me,false)", FunctionalGroups.sulphonate(meArr, false),
                "[!#8;!#1;!#7]-[#16](=[#8])(=[#8])-[#8]"));
        g.add(new Group("sulphamate(me)", FunctionalGroups.sulphamate(meArr),
                "[#7]-[#16](=[#8])(=[#8])-[#8]-[" + me + "]"));
        g.add(new Group("sulphamate(null)", FunctionalGroups.sulphamate(null),
                "[#7]-[#16](=[#8])(=[#8])-[#8]"));
        g.add(new Group("sulphate(me)", FunctionalGroups.sulphate(meArr),
                "[" + me + "]-[#8]-[#16](=[#8])(=[#8])-[#8]"));
        g.add(new Group("sulphateOfAmine(0)", FunctionalGroups.sulphateOfAmine(0),
                "[#16](-[#8-]~[#7+])(-[#8])(=[#8])=[#8]"));
        g.add(new Group("hydrochlorideOfAmine(1)", FunctionalGroups.hydrochlorideOfAmine(1),
                "[#6]-[#7+](-[#1])(-[#1])(-[#1])~[#17-]"));
        g.add(new Group("hydrochlorideOfAmine(3)", FunctionalGroups.hydrochlorideOfAmine(3),
                "[#6]-[#7+](-[#6])(-[#6])(-[#1])~[#17-]"));
        g.add(new Group("polyoxyethylene(1)", FunctionalGroups.polyoxyethylene(1),
                "[*]-[#8]-[#6]([#1])([#1])-[#6]([#1])([#1])-[*]"));
        g.add(new Group("polyoxyethylene(2)", FunctionalGroups.polyoxyethylene(2),
                "[*]-[#8]-[#6]([#1])([#1])-[#6]([#1])([#1])-[#8]-[#6]([#1])([#1])-[#6]([#1])([#1])-[*]"));
        g.add(new Group("quaternaryNitrogen1(false)", FunctionalGroups.quaternaryNitrogen1(false),
                "[*]-[#7](-[*])(-[*])-[*]"));
        g.add(new Group("quaternaryNitrogen2(false)", FunctionalGroups.quaternaryNitrogen2(false),
                "[*]-[#7](-[*])=[*]"));
        g.add(new Group("ringSubstituted(null)", verhaar.query.FunctionalGroups.ringSubstituted(null),
                "[*]@[#6]@[*]"));
        g.add(new Group("ionicGroup", verhaar.query.FunctionalGroups.ionicGroup(), "[*]~[*]"));
        return g;
    }
}

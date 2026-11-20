package insilico.toxtree;

import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.toxtree.core.Classification;
import insilico.toxtree.cramer.CramerTree;
import org.junit.jupiter.api.Test;

import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Scanner;

/**
 * Lists the molecules the port puts in a different class from Toxtree, so that each can be looked at on its own.
 * It asserts nothing: the counts are asserted by {@link TraceComparisonTest}, and this is for reading.
 */
class CategoryDiffTest {
    private static String shorten(String category) {
        int i = category.indexOf(" (Class");
        return i < 0 ? category : category.substring(0, i);
    }

    @Test
    void list() throws Exception {
        CramerTree tree = new CramerTree();
        for (String file : new String[]{"reference/cramer-core.tsv", "reference/cramer-logp.tsv"}) {
            System.out.println("### " + file);
            try (InputStream in = getClass().getClassLoader().getResourceAsStream(file);
                 Scanner sc = new Scanner(in, StandardCharsets.UTF_8)) {
                if (sc.hasNextLine()) sc.nextLine();
                while (sc.hasNextLine()) {
                    String[] p = sc.nextLine().split("\t", -1);
                    if (p.length < 3 || p[2].startsWith("ERROR")) continue;
                    String got;
                    try {
                        Classification r = tree.classify(SmilesMolecule.CreateCDKMolecule(p[0], null));
                        got = r.isClassified() ? r.categoryName() : "N.A.";
                    } catch (Exception x) {
                        got = "N.A.";
                    }
                    if (!got.equals(p[1])) {
                        System.out.printf("%-50s %-14s -> %-14s %s%n",
                                p[0].length() > 48 ? p[0].substring(0, 48) + ".." : p[0],
                                shorten(p[1]), shorten(got),
                                p[2].contains("(") ? "residui" : "");
                    }
                }
            }
        }
    }
}

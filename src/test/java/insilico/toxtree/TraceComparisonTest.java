/*
 * Copyright (C) 2026 ProtoQSAR / Alberto Manganaro
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 */
package insilico.toxtree;

import insilico.core.molecule.conversion.SmilesMolecule;
import insilico.toxtree.core.Classification;
import insilico.toxtree.core.DecisionTree;
import insilico.toxtree.cramer.CramerTree;
import insilico.toxtree.verhaar.VerhaarTree;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Compares this port against decision traces captured from the original Toxtree.
 * <p>
 * The point is not a pass/fail on the final category but a measurement: for every rule, how often this port
 * answers what Toxtree answered. While rules are still missing the numbers simply show how far the work has got;
 * once they are all in, any rule below 100% is a defect with a known address.
 * <p>
 * The reference file is produced by the harness described in the README and lives in
 * {@code src/test/resources/reference}.
 */
class TraceComparisonTest {

    private record Reference(String smiles, String category, String trace) {
    }

    @Test
    void verhaarReproducesToxtree() throws Exception {
        Agreement agreement = report(new VerhaarTree(), "reference/verhaar-core.tsv");
        assertEquals(agreement.compared(), agreement.sameTrace(),
                "Every decision trace must match the one Toxtree produced");
    }

    @Test
    void verhaarReproducesToxtreeOnTheVerhaar2000Set() throws Exception {
        // An independent set: the 121 molecules of Verhaar et al. 2000, which Toxtree ships as its own
        // validation data. Nothing in the port was tuned on them.
        Agreement agreement = report(new VerhaarTree(), "reference/verhaar-2000.tsv");
        assertEquals(agreement.compared(), agreement.sameTrace(),
                "Every decision trace must match the one Toxtree produced");
    }

    @Test
    void cramerProgress() throws Exception {
        // Not an assertion yet: the port is still being written, and the report says how far it has got.
        report(new CramerTree(), "reference/cramer-core.tsv");
    }

    /** @param compared molecules classified, @param sameTrace how many took the very same path */
    private record Agreement(int compared, int sameCategory, int sameTrace) {
    }

    private Agreement report(DecisionTree tree, String resource) throws IOException {
        List<Reference> references = load(resource);
        assertTrue(references.size() > 100, "The reference file looks empty: " + resource);

        Map<String, int[]> perRule = new TreeMap<>(TraceComparisonTest::compareRuleIds);
        Map<String, List<String>> disagreements = new TreeMap<>(TraceComparisonTest::compareRuleIds);
        int sameCategory = 0;
        int sameTrace = 0;
        int usable = 0;
        int reachedUnimplemented = 0;

        for (Reference reference : references) {
            IAtomContainer mol;
            try {
                mol = SmilesMolecule.CreateCDKMolecule(reference.smiles(), null);
            } catch (Exception x) {
                continue;
            }
            usable++;

            Classification result = tree.classify(mol);
            if (result.trace().contains("-")) {
                reachedUnimplemented++;
            }
            if (result.categoryName().equals(reference.category())) {
                sameCategory++;
            }
            if (result.trace().equals(reference.trace())) {
                sameTrace++;
            }

            Map<String, String> expected = answersOf(reference.trace());
            Map<String, String> actual = answersOf(result.trace());
            for (Map.Entry<String, String> entry : expected.entrySet()) {
                String got = actual.get(entry.getKey());
                if (got == null || "-".equals(got)) {
                    continue;
                }
                int[] counts = perRule.computeIfAbsent(entry.getKey(), k -> new int[2]);
                counts[1]++;
                if (got.equals(entry.getValue())) {
                    counts[0]++;
                } else {
                    disagreements.computeIfAbsent(entry.getKey(), k -> new ArrayList<>())
                            .add(reference.smiles() + "  (expected " + entry.getValue() + ", got " + got + ")");
                }
            }
        }

        StringBuilder out = new StringBuilder();
        out.append(String.format("%n=== %s against %s%n", tree.getTitle(), resource));
        out.append(String.format("molecules compared : %d%n", usable));
        out.append(String.format("identical category : %d (%.1f%%)%n", sameCategory, percent(sameCategory, usable)));
        out.append(String.format("identical trace    : %d (%.1f%%)%n", sameTrace, percent(sameTrace, usable)));
        out.append(String.format("paths that hit a missing rule: %d%n%n", reachedUnimplemented));
        out.append(String.format("%-10s %8s %8s   %s%n", "rule", "agree", "compared", "agreement"));
        for (Map.Entry<String, int[]> entry : perRule.entrySet()) {
            int[] counts = entry.getValue();
            out.append(String.format("%-10s %8d %8d   %6.1f%%%s%n", entry.getKey(), counts[0], counts[1],
                    percent(counts[0], counts[1]), counts[0] == counts[1] ? "" : "   <-- differs"));
        }
        if (!disagreements.isEmpty()) {
            out.append(String.format("%nmolecules where this port and Toxtree disagree:%n"));
            disagreements.forEach((rule, list) -> {
                out.append(String.format("  rule %s (%d)%n", rule, list.size()));
                list.stream().limit(8).forEach(line -> out.append("    ").append(line).append(System.lineSeparator()));
            });
        }
        System.out.println(out);
        return new Agreement(usable, sameCategory, sameTrace);
    }

    private static double percent(int part, int total) {
        return total == 0 ? 0 : 100.0 * part / total;
    }

    /** Splits a trace such as {@code "0.1Y,0.2N"} into rule id to answer. */
    private static Map<String, String> answersOf(String trace) {
        Map<String, String> answers = new LinkedHashMap<>();
        for (String step : trace.split(",")) {
            step = step.trim();
            if (step.length() < 2) {
                continue;
            }
            answers.put(step.substring(0, step.length() - 1), step.substring(step.length() - 1));
        }
        return answers;
    }

    /** Orders rule ids the way a chemist reads them: 1.4 before 1.4.1, and 1.5 before 1.10. */
    private static int compareRuleIds(String a, String b) {
        String[] left = a.split("\\.");
        String[] right = b.split("\\.");
        for (int i = 0; i < Math.min(left.length, right.length); i++) {
            int c = Integer.compare(parse(left[i]), parse(right[i]));
            if (c != 0) {
                return c;
            }
        }
        return Integer.compare(left.length, right.length);
    }

    private static int parse(String value) {
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException x) {
            return Integer.MAX_VALUE;
        }
    }

    private List<Reference> load(String resource) throws IOException {
        List<Reference> references = new ArrayList<>();
        try (InputStream in = getClass().getClassLoader().getResourceAsStream(resource)) {
            if (in == null) {
                throw new IOException("Missing reference file: " + resource);
            }
            try (Scanner scanner = new Scanner(in, StandardCharsets.UTF_8)) {
                if (scanner.hasNextLine()) {
                    scanner.nextLine();
                }
                while (scanner.hasNextLine()) {
                    String[] parts = scanner.nextLine().split("\t", -1);
                    if (parts.length >= 3 && !parts[2].startsWith("ERROR")) {
                        references.add(new Reference(parts[0], parts[1], parts[2]));
                    }
                }
            }
        }
        return references;
    }
}

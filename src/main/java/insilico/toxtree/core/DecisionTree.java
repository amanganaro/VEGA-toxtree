/*
 * Copyright (C) 2026 ProtoQSAR / Alberto Manganaro
 *
 * Derived from Toxtree, Copyright Ideaconsult Ltd. (C) 2005-2021.
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 */
package insilico.toxtree.core;

import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.StringJoiner;

/**
 * Walks a decision tree described by a transition table, the way Toxtree's {@code UserDefinedTree} does.
 * <p>
 * The table has one row per rule, each row holding four entries:
 * <pre>
 *     { nextIfNo, nextIfYes, categoryIfNo, categoryIfYes }
 * </pre>
 * All four are 1-based, and 0 means "nothing here". Walking starts at rule 1 and stops as soon as a row assigns
 * a category; a row that neither assigns a category nor names a next rule leaves the molecule unclassified.
 * <p>
 * Subclasses supply the rules, the table and the category names, and nothing else: keeping the tree declarative
 * is what makes it cheap to check against the original.
 */
public abstract class DecisionTree {

    /** Guards against a malformed table that loops; a real path is far shorter than this. */
    private static final int MAX_STEPS = 200;

    private final DecisionRule[] rules;
    private final int[][] transitions;
    private final String[] categories;

    protected DecisionTree(DecisionRule[] rules, int[][] transitions, String[] categories) {
        if (rules.length != transitions.length) {
            throw new IllegalArgumentException(
                    "The tree has " + rules.length + " rules but " + transitions.length + " transition rows");
        }
        for (int i = 0; i < transitions.length; i++) {
            if (transitions[i].length != 4) {
                throw new IllegalArgumentException("Transition row " + (i + 1) + " must have 4 entries");
            }
        }
        this.rules = rules;
        this.transitions = transitions;
        this.categories = categories;
    }

    /** The name of the tree, as it appears in model output. */
    public abstract String getTitle();

    /** The rules, in the order of the transition table. */
    public DecisionRule[] getRules() {
        return rules.clone();
    }

    public String[] getCategories() {
        return categories.clone();
    }

    /**
     * Classifies one molecule.
     * <p>
     * The molecule must already be normalized: this library does not parse SMILES, add hydrogens or perceive
     * aromaticity, it uses what it is given. A failure inside a rule stops the walk and is reported in the
     * result rather than thrown, so that one bad molecule in a batch does not stop the batch.
     *
     * @param mol a normalized molecule
     */
    public Classification classify(IAtomContainer mol) {
        RuleContext ctx = new RuleContext(mol);
        StringJoiner trace = new StringJoiner(",");

        int node = 1;
        for (int step = 0; step < MAX_STEPS; step++) {
            DecisionRule rule = rules[node - 1];

            boolean answer;
            try {
                answer = rule.isImplemented() && rule.verify(mol, ctx);
            } catch (Exception | StackOverflowError x) {
                trace.add(rule.getId() + "?");
                return Classification.failed(trace.toString(), x);
            }

            trace.add(rule.getId() + (rule.isImplemented() ? (answer ? "Y" : "N") : "-"));

            int[] row = transitions[node - 1];
            int category = answer ? row[3] : row[2];
            if (category != 0) {
                return Classification.of(category, categories[category - 1], trace.toString());
            }

            int next = answer ? row[1] : row[0];
            if (next == 0) {
                return Classification.unclassified(trace.toString());
            }
            node = next;
        }

        return Classification.failed(trace.toString(),
                new IllegalStateException("The decision path did not terminate after " + MAX_STEPS + " steps"));
    }
}

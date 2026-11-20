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

import java.util.List;
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

    /** How often a molecule may be replaced by pieces of itself before the walk gives up. */
    private static final int MAX_RESIDUE_DEPTH = 4;

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
        return walk(mol, 1, 0, null);
    }

    /**
     * Walks the tree from one question, on one molecule.
     *
     * @param from  the question to start at, 1-based
     * @param depth how many times the molecule has already been replaced by pieces of itself
     */
    private Classification walk(IAtomContainer mol, int from, int depth, RuleContext parent) {
        RuleContext ctx = new RuleContext(mol, parent);
        StringJoiner trace = new StringJoiner(",");

        int node = from;
        for (int step = 0; step < MAX_STEPS; step++) {
            DecisionRule rule = rules[node - 1];

            boolean answer;
            try {
                answer = rule.isImplemented() && rule.verify(mol, ctx);
            } catch (Exception | StackOverflowError x) {
                if (rule.appearsInTrace()) {
                    trace.add(rule.getId() + "?");
                }
                return Classification.failed(trace.toString(), x);
            }

            if (rule.appearsInTrace()) {
                trace.add(rule.getId() + (rule.isImplemented() ? (answer ? "Y" : "N") : "-"));
            }

            List<IAtomContainer> residues = ctx.takeResidues();
            boolean showResidues = ctx.residuesShowInTrace();

            int[] row = transitions[node - 1];
            int category = answer ? row[3] : row[2];
            if (category != 0) {
                return Classification.of(category, categories[category - 1], trace.toString());
            }

            int next = answer ? row[1] : row[0];
            if (next == 0) {
                return Classification.unclassified(trace.toString());
            }

            if (residues != null && !residues.isEmpty() && depth < MAX_RESIDUE_DEPTH) {
                // The rule has replaced the molecule with pieces of itself. Every later question is asked of
                // each piece in turn, and the last piece to reach a category decides, as in the original.
                Classification last = null;
                StringBuilder tail = new StringBuilder();
                for (IAtomContainer residue : residues) {
                    Classification sub = walk(residue, next, depth + 1, ctx.parentOrSelf());
                    if (showResidues) {
                        // The original sets these off in brackets, with no comma before them.
                        tail.append('(').append(sub.trace()).append(')');
                    } else {
                        if (tail.length() > 0) {
                            tail.append(',');
                        }
                        tail.append(sub.trace());
                    }
                    last = sub;
                }
                String whole = showResidues ? trace + tail.toString() : trace + "," + tail;
                if (last == null) {
                    return Classification.unclassified(trace.toString());
                }
                return new Classification(last.categoryIndex(), last.categoryName(), whole, last.error());
            }
            node = next;
        }

        return Classification.failed(trace.toString(),
                new IllegalStateException("The decision path did not terminate after " + MAX_STEPS + " steps"));
    }
}

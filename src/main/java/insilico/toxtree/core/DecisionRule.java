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

/**
 * A single yes/no question of a decision tree.
 * <p>
 * Implementations must be stateless and safe to share between molecules: everything that depends on the molecule
 * being classified belongs in the {@link RuleContext}. Patterns and other immutable data should be built once, in
 * the constructor or in a static initializer.
 */
public interface DecisionRule {

    /**
     * The identifier used in the decision trace, as in the original Toxtree: {@code "2"} for Cramer,
     * {@code "1.4.1"} for Verhaar. It is part of the contract, because traces are compared against the reference
     * implementation.
     */
    String getId();

    /** Short human-readable description of the question. */
    String getTitle();

    /**
     * Answers the question for one molecule.
     *
     * @param mol the molecule, already normalized by the caller
     * @param ctx per-molecule state shared by the rules of one classification
     * @return true when the answer is "yes"
     */
    boolean verify(IAtomContainer mol, RuleContext ctx);

    /**
     * False for rules that are declared by the original tree but have no implementation. An unimplemented rule
     * answers "no" and the trace marks it, so that it can be told apart from a genuine negative answer.
     */
    default boolean isImplemented() {
        return true;
    }

    /**
     * False for a rule that is answered but not recorded in the decision trace.
     * <p>
     * The Cramer tree reuses a handful of generic questions, "is it an open chain", "is it aromatic", as
     * shortcuts from several branches, and the original leaves them out of the trace it prints. Keeping them out
     * here too is what lets a trace be compared with the original's character for character.
     */
    default boolean appearsInTrace() {
        return true;
    }

    /** The same rule, answered as before but left out of the decision trace. */
    static DecisionRule silent(DecisionRule rule) {
        return new DecisionRule() {
            @Override
            public String getId() {
                return rule.getId();
            }

            @Override
            public String getTitle() {
                return rule.getTitle();
            }

            @Override
            public boolean verify(org.openscience.cdk.interfaces.IAtomContainer mol, RuleContext ctx) {
                return rule.verify(mol, ctx);
            }

            @Override
            public boolean isImplemented() {
                return rule.isImplemented();
            }

            @Override
            public boolean appearsInTrace() {
                return false;
            }
        };
    }
}

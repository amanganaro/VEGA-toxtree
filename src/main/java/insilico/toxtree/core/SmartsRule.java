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

import insilico.toxtree.chem.Smarts;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.List;

/**
 * A rule that answers by matching SMARTS substructures, which is how most of the Verhaar rules work.
 * <p>
 * Each pattern can be required to match or required <em>not</em> to match, and the patterns are combined with
 * either "any" or "all", reproducing Toxtree's {@code AbstractRuleSmartSubstructure}. The patterns are compiled
 * once, when the rule is built, so a bad SMARTS fails immediately rather than on the first molecule.
 */
public class SmartsRule extends AbstractRule {

    /** How the individual pattern outcomes are combined. */
    public enum Combine {
        /** The rule answers yes when at least one pattern is satisfied. */
        ANY,
        /** The rule answers yes only when every pattern is satisfied. */
        ALL
    }

    private record Pattern(String smarts, boolean negated) {
    }

    private final List<Pattern> patterns = new ArrayList<>();
    private final Combine combine;

    protected SmartsRule(String id, String title, Combine combine) {
        super(id, title);
        this.combine = combine;
    }

    /** Adds a pattern that must match. */
    protected final SmartsRule require(String... smarts) {
        for (String s : smarts) {
            Smarts.of(s);
            patterns.add(new Pattern(s, false));
        }
        return this;
    }

    /** Adds a pattern that must not match. */
    protected final SmartsRule reject(String... smarts) {
        for (String s : smarts) {
            Smarts.of(s);
            patterns.add(new Pattern(s, true));
        }
        return this;
    }

    /** The SMARTS this rule was built with, for reporting and tests. */
    public final List<String> getSmarts() {
        return patterns.stream().map(Pattern::smarts).toList();
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (patterns.isEmpty()) {
            return false;
        }
        IAtomContainer target = ctx.withExplicitHydrogens();
        for (Pattern pattern : patterns) {
            boolean satisfied = Smarts.matches(target, pattern.smarts()) != pattern.negated();
            if (combine == Combine.ANY && satisfied) {
                return true;
            }
            if (combine == Combine.ALL && !satisfied) {
                return false;
            }
        }
        return combine == Combine.ALL;
    }
}

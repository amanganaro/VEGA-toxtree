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
package insilico.toxtree.verhaar.rules;

import insilico.toxtree.core.DecisionRule;
import insilico.toxtree.core.RuleContext;
import insilico.toxtree.core.SmartsRule;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Verhaar rule 1.7.1: halogenated compounds that satisfy one of the 1.5 rules, but are not alpha- or
 * beta-halogen substituted.
 * <p>
 * The rule delegates to the four rules of branch 1.5, so it is given them rather than building its own copies:
 * the tree holds a single instance of each.
 */
public final class Rule171 extends SmartsRule {

    private static final String HALOGEN = "[Cl,Br,I,F][!$(C=*);!$(CC=*)]";

    private final DecisionRule[] branch15;

    public Rule171(DecisionRule... branch15) {
        super("1.7.1",
                "Are halogenated compounds that comply with rule Q.1.5 but not alpha-, beta- halogen substituted compounds",
                Combine.ALL);
        require(HALOGEN);
        this.branch15 = branch15;
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (!super.verify(mol, ctx)) {
            return false;
        }
        for (DecisionRule rule : branch15) {
            if (rule.verify(mol, ctx)) {
                return true;
            }
        }
        return false;
    }
}

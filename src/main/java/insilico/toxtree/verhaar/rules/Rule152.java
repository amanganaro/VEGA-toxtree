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

import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.RuleContext;
import insilico.toxtree.core.SmartsRule;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Verhaar rule 1.5.2: aliphatic alcohols, but not allylic or propargylic ones.
 * <p>
 * Note that "aliphatic" carries Toxtree's meaning: the molecule must have no ring at all. A cyclohexanol is
 * therefore not an aliphatic alcohol here.
 */
public final class Rule152 extends SmartsRule {

    private static final String ALCOHOL = "[OH1]C[!$(C=C);!$(C#C)]";
    private static final String DOUBLE_BOND = "C=C";

    public Rule152() {
        super("1.5.2", "Be aliphatic alcohols but not allylic/propargylic alcohols", Combine.ALL);
        require(ALCOHOL);
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (!ctx.flags().isAliphatic()) {
            return false;
        }
        if (!super.verify(mol, ctx)) {
            return false;
        }
        // A triple bond anywhere makes it propargylic, a double bond allylic; both are excluded.
        if (ctx.flags().isAcetylenic()) {
            return false;
        }
        return !Smarts.matches(ctx.withExplicitHydrogens(), DOUBLE_BOND);
    }
}

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

import insilico.toxtree.core.RuleContext;
import insilico.toxtree.core.SmartsRule;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Verhaar rule 2.4: primary alkyl amines, made of carbon, hydrogen and nitrogen only.
 * <p>
 * The element test looks at the heavy atoms alone, so hydrogen plays no part, and "aliphatic" again means
 * "without rings".
 */
public final class Rule24 extends SmartsRule {

    private static final String PRIMARY_AMINE = "A[NX3;H2]";

    public Rule24() {
        super("2.4", "Be primary alkyl amines (containing only C,H,N)", Combine.ALL);
        require(PRIMARY_AMINE);
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        for (IAtom atom : mol.atoms()) {
            String symbol = atom.getSymbol();
            if (!"C".equals(symbol) && !"N".equals(symbol) && !"H".equals(symbol)) {
                return false;
            }
        }
        if (!ctx.flags().isAliphatic()) {
            return false;
        }
        return super.verify(mol, ctx);
    }
}

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
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Verhaar rule 1.4.1: acyclic compounds that carry no halogen two bonds away from an unsaturation, which is to
 * say no allylic or propargylic halogen.
 */
public final class Rule141 extends AbstractRule {

    /**
     * An unsaturated bond, then two single bonds, then a halogen. The first bond may also be aromatic: Toxtree's
     * "unsaturated" bond matches double, triple and aromatic alike.
     */
    static final String HALOGEN_BETA_FROM_UNSATURATION = "[*]=,#,:[*]-[*]-[F,Cl,Br,I]";

    public Rule141() {
        super("1.4.1", "Be acyclic compound NOT containing halogen at beta-positions from unsaturations "
                + "(e.g. allylic / propargilyc halogens)");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (!ctx.flags().isOpenChain()) {
            return false;
        }
        return !Smarts.matches(ctx.withExplicitHydrogens(), HALOGEN_BETA_FROM_UNSATURATION);
    }
}

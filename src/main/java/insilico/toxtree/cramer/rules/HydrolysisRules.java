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
package insilico.toxtree.cramer.rules;

import insilico.toxtree.chem.Hydrolysis;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.DecisionRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.List;

/** The Cramer questions that ask what water does to the substance. */
public final class HydrolysisRules {

    /** The key under which the fragments are kept, so that the questions after the first reuse them. */
    private static final String PRODUCTS = "hydrolysisProducts";

    private HydrolysisRules() {
    }

    /** The fragments of the molecule after hydrolysis, worked out once per molecule. */
    static List<IAtomContainer> productsOf(RuleContext ctx) {
        return ctx.cached(PRODUCTS, () -> Hydrolysis.products(ctx.kekule()));
    }

    /**
     * Questions 15 and 29: readily hydrolysed, and to pieces that are at most mononuclear.
     * <p>
     * "Mononuclear" is read off the smallest set of smallest rings: no fragment may carry more than one.
     */
    public static DecisionRule toMononuclearResidues(String id) {
        return new AbstractRule(id, "Readily hydrolysed to mononuclear residues") {
            @Override
            public boolean verify(IAtomContainer mol, RuleContext ctx) {
                List<IAtomContainer> products = productsOf(ctx);
                if (products == null) {
                    return false;
                }
                for (IAtomContainer product : products) {
                    if (Cycles.sssr(product).toRingSet().getAtomContainerCount() > 1) {
                        return false;
                    }
                }
                return true;
            }
        };
    }

    /** Question 17: readily hydrolysed, and at least one of the pieces is a common terpene. */
    public static DecisionRule toCommonTerpene(String id) {
        return new AbstractRule(id, "Readily hydrolysed to a common terpene") {
            @Override
            public boolean verify(IAtomContainer mol, RuleContext ctx) {
                List<IAtomContainer> products = productsOf(ctx);
                if (products == null) {
                    return false;
                }
                for (IAtomContainer product : products) {
                    int rings = Cycles.relevant(product).toRingSet().getAtomContainerCount();
                    if (CommonTerpeneRule.isCommonTerpene(product, rings)) {
                        return true;
                    }
                }
                return false;
            }
        };
    }
}

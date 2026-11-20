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

import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.DecisionRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;

/** The two questions Cramer asks about rings themselves rather than about what they are made of. */
public final class RingRules {

    private RingRules() {
    }

    /**
     * "Is there more than one aromatic ring?", asked as questions 14 and 28.
     * <p>
     * The count comes from the essential cycles, not from the complete set the rest of the tree works with. That
     * is the original's choice and it matters for fused systems: naphthalene has two essential cycles but three
     * cycles in all.
     */
    public static DecisionRule manyAromaticRings(String id) {
        return new AbstractRule(id, "More than one aromatic ring") {
            @Override
            public boolean verify(IAtomContainer mol, RuleContext ctx) {
                IRingSet rings = Cycles.essential(mol).toRingSet();
                int aromatic = 0;
                for (IAtomContainer ring : rings.atomContainers()) {
                    boolean allAromatic = true;
                    for (IAtom atom : ring.atoms()) {
                        if (!atom.isAromatic()) {
                            allAromatic = false;
                            break;
                        }
                    }
                    if (allAromatic) {
                        aromatic++;
                    }
                }
                return aromatic > 1;
            }
        };
    }

    /**
     * "Does any ring carry a substituent?", asked as questions 13 and 27.
     * <p>
     * A substituent is any bond from a ring atom that is not itself a ring bond and does not lead to a hydrogen.
     */
    public static DecisionRule ringsHaveSubstituents(String id, String title) {
        return new AbstractRule(id, title) {
            @Override
            public boolean verify(IAtomContainer mol, RuleContext ctx) {
                IAtomContainer target = ctx.withExplicitHydrogens();
                IRingSet rings = ctx.explicitFlags().getRingSet();
                for (int i = 0; i < rings.getAtomContainerCount(); i++) {
                    IAtomContainer ring = rings.getAtomContainer(i);
                    for (IAtom ringAtom : ring.atoms()) {
                        int index = target.indexOf(ringAtom);
                        if (index < 0) {
                            continue;
                        }
                        for (IBond bond : target.getConnectedBondsList(target.getAtom(index))) {
                            if ("H".equals(bond.getBegin().getSymbol()) || "H".equals(bond.getEnd().getSymbol())) {
                                continue;
                            }
                            if (!bond.isInRing()) {
                                return true;
                            }
                        }
                    }
                }
                return false;
            }
        };
    }
}

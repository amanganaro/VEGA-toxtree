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

import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.HashSet;
import java.util.Set;

/**
 * Cramer question 9: is the lactone fused to another ring, or is it a five- or six-membered one with an
 * unsaturation next to the carbonyl?
 * <p>
 * A lone lactone ring is judged on its own unsaturation; a lactone that shares atoms with another ring is caught
 * whatever its unsaturation.
 */
public final class FusedOrUnsaturatedLactoneRule extends AbstractRule {

    private static final String LACTONE = "[#6]@[#6](-@[#8]@[#6])=!@[#8]";
    /** The same, with a double bond alpha to the carbonyl. */
    private static final String UNSATURATED_LACTONE = "[#6]=@[#6]@[#6](-@[#8]@[#6])=!@[#8]";

    private static final int CARBONYL_CARBON = 1;
    private static final int RING_OXYGEN = 2;

    public FusedOrUnsaturatedLactoneRule() {
        super("9", "Lactone fused to another ring, or a 5- or 6-membered unsaturated lactone");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.kekule();
        if (!Smarts.matches(target, LACTONE)) {
            return false;
        }

        IRingSet rings = ctx.explicitFlags().getRingSet();
        if (rings.getAtomContainerCount() == 1) {
            IAtomContainer ring = rings.getAtomContainer(0);
            int size = ring.getAtomCount();
            if (size < 5 || size > 6) {
                return false;
            }
            for (IBond bond : ring.bonds()) {
                if (bond.getOrder() == IBond.Order.DOUBLE || bond.getOrder() == IBond.Order.TRIPLE) {
                    return Smarts.matches(target, UNSATURATED_LACTONE);
                }
            }
            return false;
        }

        IAtomContainer lactoneRing = ringWithLactoneBond(ctx, target, rings);
        return lactoneRing != null && sharesAtomsWithAnotherRing(ctx, rings, lactoneRing);
    }

    private static IAtomContainer ringWithLactoneBond(RuleContext ctx, IAtomContainer target, IRingSet rings) {
        Set<Long> lactoneBonds = new HashSet<>();
        for (int[] match : Smarts.of(LACTONE).matchAll(target).uniqueAtoms()) {
            lactoneBonds.add(key(match[CARBONYL_CARBON], match[RING_OXYGEN]));
        }
        for (int i = 0; i < rings.getAtomContainerCount(); i++) {
            IAtomContainer ring = rings.getAtomContainer(i);
            for (IBond bond : ring.bonds()) {
                if (lactoneBonds.contains(key(ctx.indexOf(bond.getBegin()), ctx.indexOf(bond.getEnd())))) {
                    return ring;
                }
            }
        }
        return null;
    }

    private static boolean sharesAtomsWithAnotherRing(RuleContext ctx, IRingSet rings, IAtomContainer lactoneRing) {
        Set<Integer> lactoneAtoms = new HashSet<>();
        for (IAtom atom : lactoneRing.atoms()) {
            lactoneAtoms.add(ctx.indexOf(atom));
        }
        for (int i = 0; i < rings.getAtomContainerCount(); i++) {
            IAtomContainer ring = rings.getAtomContainer(i);
            if (ring == lactoneRing) {
                continue;
            }
            for (IAtom atom : ring.atoms()) {
                if (lactoneAtoms.contains(ctx.indexOf(atom))) {
                    return true;
                }
            }
        }
        return false;
    }

    private static long key(int one, int other) {
        return (long) Math.min(one, other) << 32 | Math.max(one, other);
    }
}

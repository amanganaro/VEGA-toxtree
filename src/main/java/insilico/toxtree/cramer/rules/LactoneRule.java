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
import org.openscience.cdk.interfaces.IRingSet;

import java.util.HashSet;
import java.util.Set;

/**
 * Cramer question 8: is it a lactone or a cyclic diester?
 * <p>
 * A lactone here is an ester whose single-bonded oxygen and whose carbonyl carbon both lie in a ring, while the
 * carbonyl oxygen does not. Having found one, the question is whether its ring is made of carbon and that oxygen
 * alone: a ring carrying another heteroatom is a heterocycle and goes elsewhere in the tree. An anhydride
 * disqualifies the molecule outright.
 */
public final class LactoneRule extends AbstractRule {

    /** Ring carbon, ring carbonyl carbon, ring oxygen, ring carbon; the carbonyl oxygen is outside the ring. */
    private static final String LACTONE = "[#6]@[#6](-@[#8]@[#6])=!@[#8]";
    private static final String ANHYDRIDE = "[#6]@[#6](-@[#8]@[#6]=!@[#8])=!@[#8]";

    private static final int CARBONYL_CARBON = 1;
    private static final int RING_OXYGEN = 2;
    private static final int CARBONYL_OXYGEN = 4;

    public LactoneRule() {
        super("8", "Lactone or cyclic diester");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.kekule();

        Set<Long> lactoneBonds = new HashSet<>();
        Set<Integer> lactoneOxygens = new HashSet<>();
        for (int[] match : Smarts.of(LACTONE).matchAll(target).uniqueAtoms()) {
            lactoneBonds.add(key(match[CARBONYL_CARBON], match[RING_OXYGEN]));
            lactoneOxygens.add(match[RING_OXYGEN]);
            lactoneOxygens.add(match[CARBONYL_OXYGEN]);
        }
        if (lactoneBonds.isEmpty()) {
            return false;
        }
        if (Smarts.matches(target, ANHYDRIDE)) {
            return false;
        }

        IRingSet rings = ctx.explicitFlags().getRingSet();
        for (int i = 0; i < rings.getAtomContainerCount(); i++) {
            IAtomContainer ring = rings.getAtomContainer(i);
            if (ring.getAtomCount() <= 3 || !containsLactoneBond(ctx, ring, lactoneBonds)) {
                continue;
            }
            for (IAtom atom : ring.atoms()) {
                String symbol = atom.getSymbol();
                if ("C".equals(symbol)) {
                    continue;
                }
                if ("O".equals(symbol) && lactoneOxygens.contains(ctx.indexOf(atom))) {
                    continue;
                }
                // Another heteroatom in the same ring: a heterocycle, not a lactone.
                return false;
            }
            return true;
        }
        return false;
    }

    private static boolean containsLactoneBond(RuleContext ctx, IAtomContainer ring, Set<Long> bonds) {
        for (var bond : ring.bonds()) {
            if (bonds.contains(key(ctx.indexOf(bond.getBegin()), ctx.indexOf(bond.getEnd())))) {
                return true;
            }
        }
        return false;
    }

    private static long key(int one, int other) {
        return (long) Math.min(one, other) << 32 | Math.max(one, other);
    }
}

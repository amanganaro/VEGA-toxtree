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

import insilico.toxtree.chem.FunctionalGroup;
import insilico.toxtree.chem.GroupAllocation;
import insilico.toxtree.chem.RingSubstituents;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;

/**
 * Cramer question 30: setting aside a hydroxyl or a methoxyl on the ring itself, does an aromatic ring carry
 * anything beyond short aliphatic groups?
 * <p>
 * A substituent is acceptable when it carries no ring of its own, when its longest carbon chain is at most five,
 * and when its heteroatoms all belong to one of the simple groups. Anything else was meant to make the answer
 * yes.
 * <p>
 * <b>In practice it never does, and that is reproduced on purpose.</b> Like question 6, this question asks
 * Toxtree to consider only the aromatic rings, and the selection that follows never matches any: the aromatic
 * flag is stored as if it were a bit mask and then tested against ring flag positions. No ring is selected, the
 * loop over substituents never runs, and the answer is no whatever the molecule. In the golden dataset the
 * question is asked 25 times and answered no every time, including for benzamide and phenyl isocyanate, whose
 * substituents are plainly not on the list.
 * <p>
 * The intended logic is kept below and can be switched on by {@link #SELECT_AROMATIC_RINGS}.
 */
public final class AromaticRingComplexSubstituentsRule extends AbstractRule {

    private static final int MAX_CHAIN = 5;
    private static final String METHOXY_RING = "methoxyRing";

    private static final List<FunctionalGroup> ALLOWED_GROUPS = List.of(
            CramerGroups.hydroxyOnRing(),
            FunctionalGroup.of(METHOXY_RING, "[#6;R]-[#8]-[#6]([#1])([#1])[#1]", 0),
            CramerGroups.alcohol(),
            CramerGroups.ketone(),
            CramerGroups.aldehyde(),
            CramerGroups.carboxylicAcid(),
            CramerGroups.ester());

    private static final List<String> ALLOWED = ALLOWED_GROUPS.stream()
            .map(FunctionalGroup::id)
            .distinct()
            .toList();

    public AromaticRingComplexSubstituentsRule() {
        super("30", "Aromatic ring with complex substituents");
    }

    /** Whether the aromatic rings are selected at all. False reproduces Toxtree; see the class comment. */
    private static final boolean SELECT_AROMATIC_RINGS = false;

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (!SELECT_AROMATIC_RINGS) {
            return false;
        }
        List<IRing> aromatic = ctx.explicitFlags().getAromaticRings();
        if (aromatic.isEmpty()) {
            return false;
        }

        IAtomContainer target = ctx.withExplicitHydrogens();
        GroupAllocation allocation = GroupAllocation.of(target, ALLOWED_GROUPS);

        for (IRing ring : aromatic) {
            Set<Integer> ringAtoms = new HashSet<>();
            for (IAtom atom : ring.atoms()) {
                ringAtoms.add(ctx.indexOf(atom));
            }
            for (IAtomContainer substituent : RingSubstituents.of(target, ringAtoms)) {
                if (Cycles.sssr(substituent).toRingSet().getAtomContainerCount() > 0) {
                    return true;
                }
                if (longestCarbonChain(substituent) > MAX_CHAIN) {
                    return true;
                }
                if (!allocation.onlyTheseGroups(substituent, ALLOWED)) {
                    return true;
                }
            }
        }
        return false;
    }

    /** The longest run of carbons in the fragment, counted in atoms. */
    private static int longestCarbonChain(IAtomContainer fragment) {
        int longest = 0;
        for (IAtom start : fragment.atoms()) {
            if ("C".equals(start.getSymbol())) {
                longest = Math.max(longest, walk(fragment, start, Collections.newSetFromMap(new IdentityHashMap<>())));
            }
        }
        return longest;
    }

    private static int walk(IAtomContainer fragment, IAtom atom, Set<IAtom> visited) {
        if (!visited.add(atom)) {
            return 0;
        }
        int best = 0;
        for (IAtom neighbour : neighbours(fragment, atom)) {
            if ("C".equals(neighbour.getSymbol())) {
                best = Math.max(best, walk(fragment, neighbour, visited));
            }
        }
        visited.remove(atom);
        return best + 1;
    }

    private static List<IAtom> neighbours(IAtomContainer fragment, IAtom atom) {
        List<IAtom> list = new ArrayList<>();
        fragment.getConnectedBondsList(atom).forEach(bond -> list.add(bond.getOther(atom)));
        return list;
    }
}

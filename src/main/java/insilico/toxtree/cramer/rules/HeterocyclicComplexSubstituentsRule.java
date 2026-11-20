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
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRing;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Cramer question 11: does a heterocyclic ring bear substituents beyond the simple ones?
 * <p>
 * The heteroatoms of the ring itself are set aside; what is judged is everything hanging off it. Allowed are
 * plain hydrocarbons and a short list of simple groups: alcohols, aldehydes, acetals, ketones, acids, esters,
 * mercaptans, sulphides, methyl ethers and a hydroxyl on a ring. Anything else makes the answer yes.
 * <p>
 * Where there is more than one heterocyclic ring, each is considered in turn with the rest of the structure as
 * its substituents.
 */
public final class HeterocyclicComplexSubstituentsRule extends AbstractRule {

    private static final List<FunctionalGroup> ALLOWED_GROUPS = List.of(
            CramerGroups.alcohol(),
            CramerGroups.aldehyde(),
            CramerGroups.acetal(),
            CramerGroups.ketone(),
            CramerGroups.carboxylicAcid(),
            CramerGroups.ester(),
            CramerGroups.mercaptan(),
            CramerGroups.sulphide(),
            CramerGroups.methylEther(),
            CramerGroups.hydroxyOnRing());

    private static final List<String> ALLOWED = ALLOWED_GROUPS.stream()
            .map(FunctionalGroup::id)
            .distinct()
            .toList();

    public HeterocyclicComplexSubstituentsRule() {
        super("11", "Has a heterocyclic ring with complex substituents");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        List<IRing> rings = ctx.explicitFlags().getHeterocyclicRings();
        if (rings.isEmpty()) {
            return false;
        }

        // The Kekule view: Toxtree reads bond orders without regard to aromaticity, so the sulphur of a
        // thiophene counts as an ordinary sulphide and is accounted for.
        IAtomContainer target = ctx.kekule();
        GroupAllocation allocation = GroupAllocation.of(target, ALLOWED_GROUPS);

        for (IRing ring : rings) {
            Set<Integer> ringAtoms = new HashSet<>();
            for (IAtom atom : ring.atoms()) {
                ringAtoms.add(ctx.indexOf(atom));
            }
            for (IAtomContainer substituent : RingSubstituents.of(target, ringAtoms)) {
                if (!allocation.onlyTheseGroups(substituent, ALLOWED)) {
                    return true;
                }
            }
        }
        return false;
    }
}

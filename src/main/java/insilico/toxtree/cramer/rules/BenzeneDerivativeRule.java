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
import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Cramer question 6: a benzene derivative bearing only hydrocarbon chains, hydroxy or hydroxy-ester chains, and
 * alkoxy groups, with at least one alkoxy para to a chain. It is the question that was meant to put safrole and
 * myristicin in the high class.
 * <p>
<b>This question can never be answered yes, and that is reproduced on purpose.</b>
 * <p>
 * Both this question and question 30 ask Toxtree to consider only the aromatic rings, through
 * {@code setAnalyzeAromatic(true)} together with {@code setAnalyzeOnlyRingsWithFlagSet(true)}. The selection
 * that follows never matches: the rule stores the aromatic flag as if it were a bit mask, then tests it against
 * loop indices that are ring flag positions, and the two never line up. No ring is ever selected, the loop over
 * substituents never runs, and the answer is no whatever the molecule.
 * <p>
 * It is demonstrable rather than inferred. Anethole and estragole, an alkoxy group para to a hydrocarbon chain,
 * are precisely what this question describes, and Toxtree answers no to both; across the 603 molecules of the
 * golden dataset the question is asked 243 times and answered no every time, as is question 30 for its 25.
 * <p>
 * The intended logic is kept below and can be switched on by {@link #SELECT_AROMATIC_RINGS}, so that whoever
 * decides to correct this can see what the question was meant to do and what it would change.
 */public final class BenzeneDerivativeRule extends AbstractRule {

    private static final String HYDROXY1 = "hydroxy1";
    private static final String ETHER = "ether";

    private static final List<FunctionalGroup> GROUPS = List.of(
            FunctionalGroup.of(HYDROXY1, "[#6;!a]-[#8]-[#1]", 0),
            CramerGroups.alcohol(),
            CramerGroups.ester(),
            FunctionalGroup.of(ETHER, "[#6]-[#8]-[#6]"));

    private static final List<String> CHAIN_GROUPS =
            List.of(HYDROXY1, CramerGroups.ALCOHOL, CramerGroups.ESTER);
    private static final List<String> ETHER_GROUP = List.of(ETHER);

    private enum Substituent {
        HYDROCARBON, HYDROXY_CHAIN, ESTER_CHAIN, ALKOXY
    }

    public BenzeneDerivativeRule() {
        super("6", "Benzene derivative with certain substituents");
    }

    /** Whether the aromatic rings are selected at all. False reproduces Toxtree; see the class comment. */
    private static final boolean SELECT_AROMATIC_RINGS = false;

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (!SELECT_AROMATIC_RINGS) {
            return false;
        }
        IAtomContainer ring = singleBenzeneRing(ctx);
        if (ring == null) {
            return false;
        }

        IAtomContainer target = ctx.withExplicitHydrogens();
        List<Integer> positions = cyclicOrder(ctx, ring);
        if (positions.size() != 6) {
            return false;
        }

        GroupAllocation allocation = GroupAllocation.of(target, GROUPS);
        Substituent[] places = new Substituent[6];
        int substituents = 0;

        for (IAtomContainer substituent : RingSubstituents.of(target, new HashSet<>(positions))) {
            substituents++;
            Substituent kind = classify(substituent, allocation);
            if (kind == null) {
                return false;
            }
            for (int place : attachmentPlaces(target, substituent, positions)) {
                places[place] = kind;
            }
        }
        if (substituents == 0) {
            return false;
        }

        for (int i = 0; i < 6; i++) {
            if (places[i] != Substituent.ALKOXY) {
                continue;
            }
            Substituent para = places[(i + 3) % 6];
            if (para == Substituent.HYDROCARBON || para == Substituent.ESTER_CHAIN) {
                return true;
            }
        }
        return false;
    }

    private static Substituent classify(IAtomContainer substituent, GroupAllocation allocation) {
        if (allocation.onlyTheseGroups(substituent, List.of())) {
            // Nothing but carbon and hydrogen.
            return Substituent.HYDROCARBON;
        }
        if (allocation.onlyTheseGroups(substituent, CHAIN_GROUPS)) {
            return Smarts.matches(substituent, GROUPS.get(0).smarts())
                    ? Substituent.HYDROXY_CHAIN
                    : Substituent.ESTER_CHAIN;
        }
        if (allocation.onlyTheseGroups(substituent, ETHER_GROUP)) {
            return Substituent.ALKOXY;
        }
        return null;
    }

    /** The one benzene ring, or null when there is not exactly one. */
    private static IAtomContainer singleBenzeneRing(RuleContext ctx) {
        IRingSet rings = ctx.explicitFlags().getRingSet();
        IAtomContainer benzene = null;
        int count = 0;
        for (int i = 0; i < rings.getAtomContainerCount(); i++) {
            IAtomContainer ring = rings.getAtomContainer(i);
            if (ring.getAtomCount() != 6) {
                continue;
            }
            boolean aromaticCarbons = true;
            for (IAtom atom : ring.atoms()) {
                if (!atom.isAromatic() || !"C".equals(atom.getSymbol())) {
                    aromaticCarbons = false;
                    break;
                }
            }
            if (aromaticCarbons) {
                count++;
                benzene = ring;
            }
        }
        return count == 1 ? benzene : null;
    }

    /** The ring atoms in the order they are joined, as positions in the molecule. */
    private static List<Integer> cyclicOrder(RuleContext ctx, IAtomContainer ring) {
        List<Integer> order = new ArrayList<>();
        Set<Integer> remaining = new HashSet<>();
        for (IAtom atom : ring.atoms()) {
            remaining.add(ctx.indexOf(atom));
        }
        Integer current = remaining.iterator().next();
        while (current != null) {
            order.add(current);
            remaining.remove(current);
            Integer next = null;
            for (IBond bond : ring.bonds()) {
                int begin = ctx.indexOf(bond.getBegin());
                int end = ctx.indexOf(bond.getEnd());
                if (begin == current && remaining.contains(end)) {
                    next = end;
                    break;
                }
                if (end == current && remaining.contains(begin)) {
                    next = begin;
                    break;
                }
            }
            current = next;
        }
        return order;
    }

    /** The ring positions a substituent is attached to; a bridge is attached to two. */
    private static List<Integer> attachmentPlaces(IAtomContainer target, IAtomContainer substituent,
                                                  List<Integer> positions) {
        Set<Integer> substituentAtoms = new HashSet<>();
        for (IAtom atom : substituent.atoms()) {
            substituentAtoms.add(target.indexOf(atom));
        }
        List<Integer> places = new ArrayList<>();
        for (int place = 0; place < positions.size(); place++) {
            IAtom ringAtom = target.getAtom(positions.get(place));
            for (IBond bond : target.getConnectedBondsList(ringAtom)) {
                IAtom other = bond.getOther(ringAtom);
                if (substituentAtoms.contains(target.indexOf(other))) {
                    places.add(place);
                    break;
                }
            }
        }
        return places;
    }
}

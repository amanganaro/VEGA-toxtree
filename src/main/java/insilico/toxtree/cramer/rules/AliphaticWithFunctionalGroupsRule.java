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
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.List;

/**
 * Cramer question 20: an open-chain substance bearing only common functional groups, and not too many of them.
 * <p>
 * Two allowances, as in the original. The oxygen groups a fat or a sugar would carry, and the polyoxyethylene
 * chains of a surfactant, are tolerated up to four times each; the remaining groups only once. Beyond that, and
 * for any heteroatom that belongs to none of them, the answer is no.
 */
public final class AliphaticWithFunctionalGroupsRule extends AbstractRule {

    private static final String CARBONATE = "carbonate";
    private static final String THIOESTER = "thioester";
    private static final String PRIMARY_AMINE = "primaryAmine";
    private static final String TERTIARY_AMINE = "tertiaryAmine";
    private static final String POLYOXYETHYLENE = "polyoxyethylene";

    /** Allowed up to four times each. */
    private static final List<FunctionalGroup> UP_TO_FOUR = List.of(
            FunctionalGroup.of(CARBONATE, "[*]-[#8]-[#6](-[#8]-[!#1])=[#8]", 0, 4),
            CramerGroups.carboxylicAcid(),
            CramerGroups.ester(),
            CramerGroups.aldehyde(),
            CramerGroups.alcohol());

    /** Allowed only once each. */
    private static final List<FunctionalGroup> AT_MOST_ONE = List.of(
            CramerGroups.ketone(),
            CramerGroups.acetal(),
            FunctionalGroup.of(THIOESTER, "[!#1;!#16;!#8]-[#6](-[#16]-[!#1])=[#8,#16]", 0, 3),
            CramerGroups.sulphide(),
            CramerGroups.mercaptan(),
            FunctionalGroup.of(PRIMARY_AMINE, "[#6]-[#7](-[#1])-[#1]", 0),
            FunctionalGroup.of(TERTIARY_AMINE, "[!#1]-[#7](-[!#1])-[!#1]"));

    public AliphaticWithFunctionalGroupsRule() {
        super("20", "Open chain with only simple functional groups");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (!ctx.flags().isAliphatic()) {
            return false;
        }
        IAtomContainer target = ctx.kekule();

        List<FunctionalGroup> groups = new ArrayList<>();
        groups.addAll(UP_TO_FOUR);
        groups.addAll(AT_MOST_ONE);
        // Longest first, so that a long chain is allocated whole rather than in pieces.
        for (int length = 4; length >= 1; length--) {
            groups.add(polyoxyethylene(length));
        }

        GroupAllocation allocation = GroupAllocation.of(target, groups);
        for (FunctionalGroup group : UP_TO_FOUR) {
            if (allocation.count(group.id()) > 4) {
                return false;
            }
        }
        for (FunctionalGroup group : AT_MOST_ONE) {
            if (allocation.count(group.id()) > 1) {
                return false;
            }
        }
        if (allocation.count(POLYOXYETHYLENE) > 4) {
            return false;
        }

        List<String> allowed = new ArrayList<>();
        groups.forEach(group -> allowed.add(group.id()));
        return allocation.onlyTheseGroups(target, allowed);
    }

    /** A chain of {@code length} ethylene oxide units, the atoms at either end being only the anchors. */
    private static FunctionalGroup polyoxyethylene(int length) {
        StringBuilder smarts = new StringBuilder("[*]");
        for (int i = 0; i < length; i++) {
            smarts.append("-[#8]-[#6]([#1])([#1])-[#6]([#1])([#1])");
        }
        smarts.append("-[*]");
        return FunctionalGroup.of(POLYOXYETHYLENE, smarts.toString(), 0, 1 + length * 7);
    }
}

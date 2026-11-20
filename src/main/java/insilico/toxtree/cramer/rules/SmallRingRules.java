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
import insilico.toxtree.core.DecisionRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.ArrayList;
import java.util.List;

/** Cramer questions 25 and 26, both about small rings carrying only the substituents question 24 allows. */
public final class SmallRingRules {

    private static final String ACYCLIC_ACETAL = "acyclicAcetal";
    private static final String SIDECHAIN_KETONE = "sidechainKetone";
    private static final String RING_KETONE = "ringKetone";
    private static final String SULPHONATE = "sulphonate";
    private static final String SULPHAMATE = "sulphamate";

    private static final String ALIPHATIC_CHN = "[$([#6;!a]~[#1]),$([#6;!a;!$([#6]~[!#6])])]";

    private static final FunctionalGroup ALIPHATIC_ALCOHOL =
            FunctionalGroup.of(CramerGroups.ALCOHOL, ALIPHATIC_CHN + "-[#8]-[#1]", 0);
    private static final FunctionalGroup SIDECHAIN_KETONE_GROUP =
            FunctionalGroup.of(SIDECHAIN_KETONE, "[#6]-[#6;!R](=[#8;!R])-[#6]", 0, 3);
    private static final FunctionalGroup RING_KETONE_GROUP =
            FunctionalGroup.of(RING_KETONE, "[#6;R]-[#6;R](=[#8])-[#6;R]", 0, 3);
    private static final FunctionalGroup SULPHONATE_GROUP =
            FunctionalGroup.of(SULPHONATE, "[!#8;!#1;!#7]-[#16](=[#8])(=[#8])-[#8]", 0);
    private static final FunctionalGroup SULPHAMATE_GROUP =
            FunctionalGroup.of(SULPHAMATE, "[#7]-[#16](=[#8])(=[#8])-[#8]");

    /** The ring ketone that question 26 insists a single ring must carry. */
    private static final String RING_KETONE_SMARTS = "[#6;R]-[#6;R](=[#8])-[#6;R]";

    private SmallRingRules() {
    }

    /**
     * Question 25: a cyclopropane or cyclobutane carrying only the substituents of question 24.
     */
    public static DecisionRule cyclopropaneOrCyclobutane(String id) {
        List<FunctionalGroup> groups = List.of(
                ALIPHATIC_ALCOHOL,
                CramerGroups.aldehyde(),
                CramerGroups.acetal(),
                SIDECHAIN_KETONE_GROUP,
                CramerGroups.carboxylicAcid(),
                CramerGroups.ester(),
                SULPHONATE_GROUP);
        List<String> allowed = ids(groups);

        return new AbstractRule(id, "Cyclopropane or cyclobutane with simple substituents") {
            @Override
            public boolean verify(IAtomContainer mol, RuleContext ctx) {
                IRingSet rings = ctx.explicitFlags().getRingSet();
                if (rings.getAtomContainerCount() != 1) {
                    return false;
                }
                IAtomContainer ring = rings.getAtomContainer(0);
                if (ring.getAtomCount() > 4) {
                    return false;
                }
                for (IAtom atom : ring.atoms()) {
                    if (!"C".equals(atom.getSymbol())) {
                        return false;
                    }
                }

                IAtomContainer target = ctx.withExplicitHydrogens();
                GroupAllocation allocation = GroupAllocation.of(target, groups);
                if (allocation.onlyTheseGroups(target, allowed)) {
                    return true;
                }
                for (IAtomContainer substituent : RingSubstituents.of(target, ring)) {
                    if (!allocation.onlyTheseGroups(substituent, allowed)) {
                        return false;
                    }
                }
                return true;
            }
        };
    }

    /**
     * Question 26: a monocycloalkanone, or any bicyclic compound, carrying nothing beyond the groups of
     * question 24. A single ring has to carry a ketone in the ring itself; two rings need not.
     */
    public static DecisionRule monocycloalkanoneOrBicyclic(String id) {
        List<FunctionalGroup> groups = List.of(
                ALIPHATIC_ALCOHOL,
                CramerGroups.aldehyde(),
                FunctionalGroup.of(ACYCLIC_ACETAL, "[#6]-!@[#8]-!@[#6]-!@[#8]-!@[#6]"),
                SIDECHAIN_KETONE_GROUP,
                RING_KETONE_GROUP,
                CramerGroups.carboxylicAcid(),
                CramerGroups.ester(),
                SULPHONATE_GROUP,
                SULPHAMATE_GROUP);
        List<String> allowed = ids(groups);

        return new AbstractRule(id, "Monocycloalkanone or a bicyclic compound") {
            @Override
            public boolean verify(IAtomContainer mol, RuleContext ctx) {
                IAtomContainer target = ctx.withExplicitHydrogens();
                if (!GroupAllocation.of(target, groups).onlyTheseGroups(target, allowed)) {
                    return false;
                }
                int rings = ctx.explicitFlags().getRingCount();
                if (rings == 1) {
                    return !ctx.flags().isHeterocyclic() && Smarts.matches(target, RING_KETONE_SMARTS);
                }
                return rings == 2;
            }
        };
    }

    private static List<String> ids(List<FunctionalGroup> groups) {
        List<String> ids = new ArrayList<>();
        groups.forEach(group -> {
            if (!ids.contains(group.id())) {
                ids.add(group.id());
            }
        });
        return ids;
    }
}

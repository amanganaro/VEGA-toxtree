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
import org.openscience.cdk.interfaces.IRingSet;

import java.util.List;

/**
 * Cramer question 24: a single carbocyclic ring with nothing but simple substituents.
 * <p>
 * The ring must be the only one, must have at least five atoms, which leaves out cyclopropane and cyclobutane,
 * and must be made of carbon. What may hang off it is an aliphatic alcohol, an aldehyde, an acyclic acetal, a
 * ketone outside the ring, an acid, an ester, or a sulphonate or sulphamate.
 */
public final class MonocarbocyclicRule extends AbstractRule {

    private static final String ACYCLIC_ACETAL = "acyclicAcetal";
    private static final String SIDECHAIN_KETONE = "sidechainKetone";
    private static final String SULPHONATE = "sulphonate";
    private static final String SULPHAMATE = "sulphamate";

    /** An alcohol on a carbon that is not aromatic. */
    private static final String ALIPHATIC_CHN = "[$([#6;!a]~[#1]),$([#6;!a;!$([#6]~[!#6])])]";

    private static final List<FunctionalGroup> ALLOWED_GROUPS = List.of(
            FunctionalGroup.of(CramerGroups.ALCOHOL, ALIPHATIC_CHN + "-[#8]-[#1]", 0),
            CramerGroups.aldehyde(),
            FunctionalGroup.of(ACYCLIC_ACETAL, "[#6]-!@[#8]-!@[#6]-!@[#8]-!@[#6]"),
            FunctionalGroup.of(SIDECHAIN_KETONE, "[#6]-[#6;!R](=[#8;!R])-[#6]", 0, 3),
            CramerGroups.carboxylicAcid(),
            CramerGroups.ester(),
            FunctionalGroup.of(SULPHONATE, "[!#8;!#1;!#7]-[#16](=[#8])(=[#8])-[#8]", 0),
            FunctionalGroup.of(SULPHAMATE, "[#7]-[#16](=[#8])(=[#8])-[#8]"));

    private static final List<String> ALLOWED = ALLOWED_GROUPS.stream()
            .map(FunctionalGroup::id)
            .distinct()
            .toList();

    public MonocarbocyclicRule() {
        super("24", "Monocarbocyclic with simple substituents");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IRingSet rings = ctx.explicitFlags().getRingSet();
        if (rings.getAtomContainerCount() != 1) {
            return false;
        }
        IAtomContainer ring = rings.getAtomContainer(0);
        if (ring.getAtomCount() < 5) {
            return false;
        }
        for (IAtom atom : ring.atoms()) {
            if (!"C".equals(atom.getSymbol())) {
                return false;
            }
        }

        IAtomContainer target = ctx.withExplicitHydrogens();
        GroupAllocation allocation = GroupAllocation.of(target, ALLOWED_GROUPS);
        if (allocation.onlyTheseGroups(target, ALLOWED)) {
            return true;
        }
        for (IAtomContainer substituent : RingSubstituents.of(target, ring)) {
            if (!allocation.onlyTheseGroups(substituent, ALLOWED)) {
                return false;
            }
        }
        return true;
    }

    static {
        // Fail early if any pattern is malformed.
        ALLOWED_GROUPS.forEach(group -> Smarts.of(group.smarts()));
    }
}

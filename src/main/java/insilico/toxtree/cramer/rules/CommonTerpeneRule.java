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
import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.List;

/**
 * Cramer questions 16 and 38: is the substance a common terpene?
 * <p>
 * The test is one of composition rather than of skeleton. The carbons must come in multiples of five and there
 * must be at least two such units; the hydrogens must come to eight per unit; nothing but carbon, hydrogen and
 * oxygen is allowed; at least one isoprene unit must be present; and any oxygen must belong to an alcohol, an
 * aldehyde or a carboxylic acid. Rings bring one more condition, on how many double bonds are left over.
 * <p>
 * The hydrogen test uses integer division, so it is satisfied by a range rather than by an exact count. That is
 * the original's arithmetic and it is kept.
 */
public final class CommonTerpeneRule extends AbstractRule {

    private static final String ISOPRENE = "isoprene";
    private static final String ALCOHOL = "alcohol";
    private static final String ALDEHYDE = "aldehyde";
    private static final String ACID = "acid";

    /** Five carbons: a chain of three whose middle atom carries a branch, the chain bonds single or double. */
    private static final FunctionalGroup ISOPRENE_UNIT =
            FunctionalGroup.overlapping(ISOPRENE, "[#6]-,=[#6](-[#6])-,=[#6]-,=[#6]");

    private static final List<FunctionalGroup> OXYGEN_GROUPS = List.of(
            FunctionalGroup.overlapping(ALCOHOL, Smarts.CHN + "-[#8]-[#1]", 0),
            FunctionalGroup.overlapping(ALDEHYDE, "[#6,#1]-[#6]([#1])=[#8]", 0),
            FunctionalGroup.overlapping(ACID, "[*]-[#6](-[#8]-[#1])=[#8]", 0));

    public CommonTerpeneRule(String id) {
        super(id, "Common terpene");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        return isCommonTerpene(ctx.kekule(), ctx.explicitFlags().getRingCount());
    }

    /**
     * The test itself, so that question 17 can put the fragments of a hydrolysis to it.
     *
     * @param target    a molecule with explicit hydrogens
     * @param ringCount how many rings it has
     */
    public static boolean isCommonTerpene(IAtomContainer target, int ringCount) {

        int carbon = 0;
        int hydrogen = 0;
        int oxygen = 0;
        for (IAtom atom : target.atoms()) {
            switch (atom.getSymbol()) {
                case "C" -> carbon++;
                case "H" -> hydrogen++;
                case "O" -> oxygen++;
                default -> {
                }
            }
        }

        if (carbon % 5 != 0) {
            return false;
        }
        int units = carbon / 5;
        if (units < 2 || hydrogen / 8 != units) {
            return false;
        }
        if (carbon + hydrogen + oxygen < target.getAtomCount()) {
            return false;
        }

        List<FunctionalGroup> groups = new ArrayList<>();
        groups.add(ISOPRENE_UNIT);
        List<String> allowed = new ArrayList<>();
        allowed.add(ISOPRENE);
        if (oxygen > 0) {
            groups.addAll(OXYGEN_GROUPS);
            allowed.add(ALCOHOL);
            allowed.add(ALDEHYDE);
            allowed.add(ACID);
        }

        GroupAllocation allocation = GroupAllocation.of(target, groups);
        if (allocation.count(ISOPRENE) == 0) {
            return false;
        }

        if (ringCount > 0 && carbon + hydrogen == target.getAtomCount()) {
            int doubleBonds = 0;
            for (IBond bond : target.bonds()) {
                if (bond.getOrder() == IBond.Order.DOUBLE) {
                    doubleBonds++;
                }
            }
            if (ringCount == 1 && doubleBonds != 2) {
                return false;
            }
            if (ringCount == 2 && doubleBonds != 1) {
                return false;
            }
        }

        return allocation.onlyTheseGroups(target, allowed);
    }
}

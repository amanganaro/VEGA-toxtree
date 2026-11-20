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
import org.openscience.cdk.interfaces.IRingSet;

import java.util.List;

/**
 * Cramer question 5: a simply branched aliphatic hydrocarbon, or a common carbohydrate.
 * <p>
 * Two quite different tests behind one question. A molecule of nothing but carbon and hydrogen is judged on its
 * shape alone: it passes when it has no ring and no triple bond. Anything else is put to the carbohydrate test,
 * which asks for carbon, hydrogen and oxygen only, at least three hydroxyls, and no functional group beyond the
 * ones a sugar is allowed: an ether if it is cyclic, at most one aldehyde and one ketone if it is not.
 */
public final class SimplyBranchedOrCarbohydrateRule extends AbstractRule {

    private static final String ALCOHOL = "alcohol";
    private static final String ETHER = "ether";
    private static final String ALDEHYDE = "aldehyde";
    private static final String KETONE = "ketone";

    /** The carbon is only the anchor and stays unallocated, so it remains free for another group. */
    private static final FunctionalGroup ALCOHOL_GROUP =
            FunctionalGroup.of(ALCOHOL, Smarts.CHN + "-[#8]-[#1]", 0);
    private static final FunctionalGroup ETHER_GROUP =
            FunctionalGroup.of(ETHER, Smarts.CHN + "-[#8]-" + Smarts.CHN);
    private static final FunctionalGroup ALDEHYDE_GROUP =
            FunctionalGroup.of(ALDEHYDE, "[#6,#1]-[#6]([#1])=[#8]", 0);
    private static final FunctionalGroup KETONE_GROUP =
            FunctionalGroup.of(KETONE, "[!#1;!#8;!#16;!#7]-[#6](=[#8])-[!#1;!#8;!#16;!#7]", 0, 3);

    private static final List<FunctionalGroup> CYCLIC = List.of(ALCOHOL_GROUP, ETHER_GROUP);
    private static final List<FunctionalGroup> ACYCLIC = List.of(ALCOHOL_GROUP, ALDEHYDE_GROUP, KETONE_GROUP);

    public SimplyBranchedOrCarbohydrateRule() {
        super("5", "Simply branched aliphatic hydrocarbon or a common carbohydrate");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.kekule();

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

        if (target.getAtomCount() == carbon + hydrogen) {
            // A hydrocarbon: "aliphatic" already means no ring, and it excludes the acetylenic ones.
            return ctx.flags().isAliphatic() && !ctx.flags().isAcetylenic();
        }
        if (target.getAtomCount() != carbon + hydrogen + oxygen || oxygen == 0) {
            return false;
        }

        IRingSet rings = ctx.explicitFlags().getRingSet();
        if (rings.getAtomContainerCount() > 0) {
            if (ctx.flags().isAromatic()) {
                return false;
            }
            for (int i = 0; i < rings.getAtomContainerCount(); i++) {
                int size = rings.getAtomContainer(i).getAtomCount();
                if (size < 5 || size > 6) {
                    return false;
                }
            }
            GroupAllocation allocation = GroupAllocation.of(target, CYCLIC);
            if (allocation.count(ALCOHOL) < 3 || allocation.count(ETHER) == 0) {
                return false;
            }
            return allocation.onlyTheseGroups(target, List.of(ALCOHOL, ETHER));
        }

        GroupAllocation allocation = GroupAllocation.of(target, ACYCLIC);
        if (allocation.count(ALCOHOL) < 3) {
            return false;
        }
        if (allocation.count(ALDEHYDE) > 1 || allocation.count(KETONE) > 1) {
            return false;
        }
        return allocation.onlyTheseGroups(target, List.of(ALCOHOL, ALDEHYDE, KETONE));
    }
}

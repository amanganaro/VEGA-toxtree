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

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Cramer question 21: does the substance carry three or more different functional groups?
 * <p>
 * The counting is done element by element, and it is generous in places and strict in others, both of which the
 * original intends. A nitrogen in an amine counts once however many amines there are, while every hydroxyl
 * counts separately, which is what puts a sugar over the threshold. An element that is on none of the lists
 * counts as a group all by itself, a shortcut the original admits to.
 */
public final class ManyFunctionalGroupsRule extends AbstractRule {

    private static final int THRESHOLD = 3;
    private static final String METALS = "Na,K,Ca,Mg,Al";

    private static final String AMINE_SECONDARY = "[#6]-[#7+0](-[#6])-[#1]";
    private static final String AMINE_TERTIARY = "[!#1]-[#7](-[!#1])-[!#1]";
    private static final String AMINE_PRIMARY = "[#6]-[#7](-[#1])-[#1]";
    private static final String ACETYLENIC = "[#6]#[#6]";

    private static final List<FunctionalGroup> NITROGEN_GROUPS = List.of(
            FunctionalGroup.of("cyano", "[#7,#6,#1]-[#6]#[#7]"),
            FunctionalGroup.of("nitroCharged", "[!#1]-[#7+](=[#8])-[#8-]"),
            FunctionalGroup.of("nitroNeutral", "[!#1]-[#7](=[#8])=[#8]"),
            FunctionalGroup.of("nNitroso", "[!#1]-[#7](-[!#1])-[#7]=[#8]"),
            FunctionalGroup.of("diazo", "[#7]#[#7]"),
            FunctionalGroup.of("triazeno", "[#6]-[#7]=[#7]-[#7](-[#1])-[#1]"));

    private static final List<FunctionalGroup> OXYGEN_GROUPS = List.of(
            FunctionalGroup.of("salt", "[#6,#1]-[#6](-[#8]-[" + METALS + "])=[#8]", 0),
            CramerGroups.acetal(),
            FunctionalGroup.of("ether", Smarts.CHN + "-[#8]-" + Smarts.CHN),
            CramerGroups.aldehyde(),
            CramerGroups.ketone(),
            CramerGroups.alcohol());

    private static final List<FunctionalGroup> ACID_OR_ESTER = List.of(
            CramerGroups.carboxylicAcid(),
            CramerGroups.ester());

    private static final List<FunctionalGroup> SULPHUR_GROUPS = List.of(
            CramerGroups.sulphide(),
            CramerGroups.mercaptan(),
            FunctionalGroup.of("thioester", "[!#1;!#16;!#8]-[#6](-[#16]-[!#1])=[#8,#16]", 0, 3),
            FunctionalGroup.of("sulphate", "[" + METALS + "]-[#8]-[#16](=[#8])(=[#8])-[#8]"),
            FunctionalGroup.of("sulphamate", "[#7]-[#16](=[#8])(=[#8])-[#8]-[" + METALS + "]"),
            FunctionalGroup.of("sulphonate", "[!#8;!#1;!#7]-[#16](=[#8])(=[#8])-[#8]-[" + METALS + "]", 0));

    public ManyFunctionalGroupsRule() {
        super("21", "Three or more different functional groups");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.kekule();

        Set<String> elements = new HashSet<>();
        for (IAtom atom : target.atoms()) {
            elements.add(atom.getSymbol());
        }

        int groups = 0;
        for (String element : elements) {
            switch (element) {
                case "H" -> {
                }
                case "C" -> {
                    if (Smarts.matches(target, ACETYLENIC)) {
                        groups++;
                    }
                }
                case "N" -> {
                    if (Smarts.matchesAny(target, AMINE_SECONDARY, AMINE_TERTIARY, AMINE_PRIMARY)) {
                        groups++;
                    }
                    groups += countAll(target, NITROGEN_GROUPS);
                }
                case "O" -> {
                    if (countAll(target, ACID_OR_ESTER) > 0) {
                        groups++;
                    }
                    groups += countAll(target, OXYGEN_GROUPS);
                }
                case "S" -> groups += countAll(target, SULPHUR_GROUPS);
                // Anything else counts as a group in its own right.
                default -> groups++;
            }
            if (groups >= THRESHOLD) {
                return true;
            }
        }
        return groups >= THRESHOLD;
    }

    /**
     * How many of these groups are present, each counting once however many times it occurs. An acetal with two
     * ether linkages is an acetal and an ether, which is two, not three.
     */
    private static int countAll(IAtomContainer target, List<FunctionalGroup> groups) {
        GroupAllocation allocation = GroupAllocation.of(target, groups);
        int total = 0;
        for (FunctionalGroup group : groups) {
            if (allocation.count(group.id()) > 0) {
                total++;
            }
        }
        return total;
    }
}

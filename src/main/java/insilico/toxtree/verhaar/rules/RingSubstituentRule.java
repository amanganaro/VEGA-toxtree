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
package insilico.toxtree.verhaar.rules;

import insilico.toxtree.chem.FunctionalGroup;
import insilico.toxtree.chem.GroupAllocation;
import insilico.toxtree.chem.RingSubstituents;
import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRing;

import java.util.ArrayList;
import java.util.List;

/**
 * The shape shared by Verhaar's class 2 rules: a named aromatic system whose substituents may only be a limited
 * number of nitro groups, a limited number of halogens, and otherwise plain alkyl.
 * <p>
 * Rules 2.1, 2.2, 2.3 and 2.5 are this one rule with different settings:
 * <table>
 *   <caption>Settings per rule</caption>
 *   <tr><th>rule</th><th>main structure</th><th>halogens</th><th>max halogens</th><th>max nitro</th></tr>
 *   <tr><td>2.1</td><td>phenol</td><td>Cl, Br, F</td><td>3</td><td>1</td></tr>
 *   <tr><td>2.2</td><td>aniline</td><td>Cl</td><td>3</td><td>1</td></tr>
 *   <tr><td>2.3</td><td>none, aromatic rings only</td><td>Cl, Br</td><td>2</td><td>1, and exactly 1 required</td></tr>
 *   <tr><td>2.5</td><td>pyridine</td><td>Cl</td><td>2</td><td>0</td></tr>
 * </table>
 * <p>
 * Each ring is cut out in turn and every substituent is judged twice: once by what kind of group it is, which is
 * where the counting happens, and once by whether all its heteroatoms belong to an allowed group at all. The
 * counters run across the whole molecule, not per ring, as in the original.
 * <p>
 * <b>The main structure is never required to be present.</b> That looks wrong, and chemically it is: Toxtree
 * answers "yes, a non- or weakly acidic phenol" for thiophene and for thiirane. It is a defect in the original,
 * and a subtle one. {@code RuleRingMainStrucSubstituents} assigns the main structure inside {@code initQuery()},
 * which the superclass constructor calls, and then declares the field as {@code mainStructure = null}: in Java a
 * field initializer runs after the superclass constructor, so the assignment is wiped before the rule is ever
 * used. The query itself survives in the allowed-group list, because it was registered while the field still
 * held it. The behaviour is reproduced here on purpose, since the point of the port is to answer as the
 * published models answer today; 18 of the 603 molecules in the golden dataset depend on it.
 */
public class RingSubstituentRule extends AbstractRule {

    private static final String NITRO = "nitro";
    private static final String HALOGEN = "halogen";
    private static final String MAIN = "main";

    private static final FunctionalGroup NITRO_CHARGED =
            FunctionalGroup.of(NITRO, "[!#1]-[N+](=[#8])-[#8-]");
    private static final FunctionalGroup NITRO_NEUTRAL =
            FunctionalGroup.of(NITRO, "[!#1]-[#7](=[#8])=[#8]");

    /** A benzene ring carrying a hydroxyl. */
    public static final String PHENOL = "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]-[#1]";
    /** A benzene ring carrying an amine. */
    public static final String ANILINE = "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#7]-[#1]";
    /** An aromatic six-membered ring that starts at a nitrogen and is otherwise carbon or nitrogen. */
    public static final String PYRIDINE = "[#7]1:[#6,#7]:[#6,#7]:[#6,#7]:[#6,#7]:[#6,#7]:1";

    private final String mainStructure;
    private final List<FunctionalGroup> groups = new ArrayList<>();
    private final List<String> allowedIds = new ArrayList<>();
    private final List<String> mainIds = new ArrayList<>();
    private final List<String> nitroIds = List.of(NITRO);
    private final List<String> halogenIds = List.of(HALOGEN);
    private final int maxHalogens;
    private final int maxNitroGroups;
    private final boolean aromaticRingsOnly;
    private final boolean exactlyOneNitroGroup;

    public RingSubstituentRule(String id, String title, String mainStructure, String halogens,
                                  int maxHalogens, int maxNitroGroups, boolean aromaticRingsOnly,
                                  boolean exactlyOneNitroGroup) {
        super(id, title);
        this.mainStructure = mainStructure;
        this.maxHalogens = maxHalogens;
        this.maxNitroGroups = maxNitroGroups;
        this.aromaticRingsOnly = aromaticRingsOnly;
        this.exactlyOneNitroGroup = exactlyOneNitroGroup;

        groups.add(NITRO_CHARGED);
        groups.add(NITRO_NEUTRAL);
        groups.add(FunctionalGroup.of(HALOGEN, "[*]~[" + halogens + "]"));
        allowedIds.add(NITRO);
        allowedIds.add(HALOGEN);
        if (mainStructure != null) {
            groups.add(FunctionalGroup.of(MAIN, mainStructure));
            allowedIds.add(MAIN);
        }
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.withExplicitHydrogens();

        // No ring, no answer: Toxtree only records a ring set when the molecule has at least one ring, and these
        // rules give up when they are handed none. An acyclic molecule therefore never satisfies them.
        if (ctx.explicitFlags().getRingCount() == 0) {
            return false;
        }

        List<IRing> rings = ringsOf(ctx);
        GroupAllocation allocation = GroupAllocation.of(target, groups);

        int[] counters = new int[2]; // nitro groups, halogens
        for (IRing ring : rings) {
            if (aromaticRingsOnly && !isAromatic(ring)) {
                continue;
            }
            for (IAtomContainer substituent : RingSubstituents.of(target, ring)) {
                if (!substituentIsAllowed(substituent, allocation, counters)) {
                    return false;
                }
                if (!allocation.onlyTheseGroups(substituent, allowedIds)) {
                    return false;
                }
            }
        }
        return !exactlyOneNitroGroup || counters[0] == 1;
    }

    private List<IRing> ringsOf(RuleContext ctx) {
        List<IRing> rings = new ArrayList<>();
        var ringSet = ctx.explicitFlags().getRingSet();
        for (int i = 0; i < ringSet.getAtomContainerCount(); i++) {
            rings.add((IRing) ringSet.getAtomContainer(i));
        }
        return rings;
    }

    private static boolean isAromatic(IRing ring) {
        for (IAtom atom : ring.atoms()) {
            if (!atom.isAromatic()) {
                return false;
            }
        }
        return true;
    }

    /**
     * Judges one substituent and keeps the running counts. The order of the tests is the original's, and it
     * matters: a substituent without heteroatoms satisfies the first test and is accepted as alkyl before any
     * counter is touched.
     */
    private boolean substituentIsAllowed(IAtomContainer substituent, GroupAllocation allocation, int[] counters) {
        if (allocation.onlyTheseGroups(substituent, mainIds)) {
            return true;
        }
        if (allocation.onlyTheseGroups(substituent, nitroIds)) {
            counters[0]++;
            return counters[0] <= maxNitroGroups;
        }
        if (allocation.onlyTheseGroups(substituent, halogenIds)) {
            counters[1]++;
            return counters[1] <= maxHalogens;
        }
        for (IAtom atom : substituent.atoms()) {
            String symbol = atom.getSymbol();
            if (!"C".equals(symbol) && !"H".equals(symbol)) {
                return false;
            }
        }
        return true;
    }
}

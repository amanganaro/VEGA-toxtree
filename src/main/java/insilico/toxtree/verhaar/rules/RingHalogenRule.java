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
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.List;

/**
 * Verhaar rules 1.4.3 and 1.4.4: cyclic compounds that are unsubstituted, or whose substituents are acyclic and
 * made of carbon, hydrogen and halogens only. The two differ in nothing but the number of rings they accept.
 * <p>
 * A benzylic halogen disqualifies the molecule outright, before the rings are looked at, because such compounds
 * do not satisfy rule 1.4.1.
 */
public class RingHalogenRule extends AbstractRule {

    private static final String HALOGEN = "halogen";
    private static final List<FunctionalGroup> GROUPS =
            List.of(new FunctionalGroup(HALOGEN, "[*]~[F,Cl,Br,I]"));
    private static final List<String> ALLOWED = List.of(HALOGEN);

    private final boolean polycyclic;

    public RingHalogenRule(String id, String title, boolean polycyclic) {
        super(id, title);
        this.polycyclic = polycyclic;
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.withExplicitHydrogens();
        if (Smarts.matches(target, Rule141.HALOGEN_BETA_FROM_UNSATURATION)) {
            return false;
        }

        int rings = ctx.explicitFlags().getRingCount();
        if (rings == 0 || (polycyclic ? rings <= 1 : rings != 1)) {
            return false;
        }

        GroupAllocation allocation = GroupAllocation.of(target, GROUPS);
        // When nothing but halogens stands out in the whole molecule there is no need to take it apart.
        if (allocation.onlyTheseGroups(target, ALLOWED)) {
            return true;
        }

        IRingSet ringSet = ctx.explicitFlags().getRingSet();
        for (int i = 0; i < ringSet.getAtomContainerCount(); i++) {
            for (IAtomContainer substituent : RingSubstituents.of(target, ringSet.getAtomContainer(i))) {
                if (!allocation.onlyTheseGroups(substituent, ALLOWED)) {
                    return false;
                }
            }
        }
        return true;
    }
}

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

import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import java.util.HashSet;
import java.util.Set;

/**
 * Verhaar rule 1.2: the molecule must not contain ionic groups.
 * <p>
 * Toxtree looks for a salt rather than for a charge: it splits the structure into connected components and tries
 * to pair oppositely charged atoms that sit in <em>different</em> components, drawing a virtual bond between the
 * ones it can pair. A zwitterion, whose charges balance inside a single component, is therefore not an ionic
 * group, while sodium cyclamate is.
 * <p>
 * The original walks the pairing atom by atom, consuming charge as it goes, because it also wants the bonds
 * themselves for highlighting. The answer to this question only depends on whether at least one pair exists, so
 * that is what is computed here.
 */
public final class RuleIonicGroups extends AbstractRule {

    public RuleIonicGroups() {
        super("1.2", "Not contain ionic groups");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        return !hasIonicGroup(mol);
    }

    private static boolean hasIonicGroup(IAtomContainer mol) {
        IAtomContainerSet components = ConnectivityChecker.partitionIntoMolecules(mol);
        int count = components.getAtomContainerCount();
        if (count < 2) {
            return false;
        }

        Set<Integer> withPositive = new HashSet<>();
        Set<Integer> withNegative = new HashSet<>();
        for (int i = 0; i < count; i++) {
            for (IAtom atom : components.getAtomContainer(i).atoms()) {
                Integer charge = atom.getFormalCharge();
                if (charge == null || charge == 0) {
                    continue;
                }
                (charge > 0 ? withPositive : withNegative).add(i);
            }
        }

        if (withPositive.isEmpty() || withNegative.isEmpty()) {
            return false;
        }
        // A pair is possible unless every charge of both signs sits in one and the same component.
        return withPositive.size() > 1 || !withPositive.equals(withNegative);
    }
}

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

import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.Set;

/**
 * Cramer question 3: does the structure contain elements other than carbon, hydrogen, oxygen, nitrogen or
 * <em>divalent</em> sulphur?
 * <p>
 * Sulphur is the interesting part, and the original decides it in two different ways depending on the molecule.
 * It first asks CDK for the atom's valency, which for an aromatic sulphur such as thiophene's is two. But the
 * moment it meets a sulphur whose valency is not two it tries to record the atom in a container that the caller
 * did not supply, throws a null pointer, swallows it, and starts again on an older path that adds up the bond
 * orders instead, counting an aromatic bond as one and a half and including the implicit hydrogens.
 * <p>
 * The two disagree exactly where it matters: by bond order an aromatic sulphur comes to three, not two, so
 * thiophene would be rejected. Which path runs is decided by the whole molecule, not by the atom being looked
 * at, so both are implemented and chosen the same way.
 */
public final class ElementsOtherThanRule extends AbstractRule {

    private static final Set<String> ALLOWED = Set.of("C", "H", "O", "N", "S");

    public ElementsOtherThanRule() {
        super("3", "Contains elements other than C, H, O, N, divalent S");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.withExplicitHydrogens();
        boolean byValency = true;
        for (IAtom atom : target.atoms()) {
            if ("S".equals(atom.getSymbol())) {
                Integer valency = atom.getValency();
                if (valency == null || valency != 2) {
                    byValency = false;
                    break;
                }
            }
        }

        for (IAtom atom : target.atoms()) {
            String symbol = atom.getSymbol();
            if ("H".equals(symbol)) {
                continue;
            }
            if (!ALLOWED.contains(symbol)) {
                return true;
            }
            if ("S".equals(symbol) && !byValency && !isDivalentByBondOrder(target, atom)) {
                return true;
            }
        }
        return false;
    }

    private static boolean isDivalentByBondOrder(IAtomContainer mol, IAtom atom) {
        double order = 0;
        for (IBond bond : mol.getConnectedBondsList(atom)) {
            if (bond.isAromatic()) {
                order += 1.5;
            } else {
                switch (bond.getOrder()) {
                    case SINGLE -> order += 1;
                    case DOUBLE -> order += 2;
                    case TRIPLE -> order += 3;
                    default -> {
                    }
                }
            }
        }
        Integer implicit = atom.getImplicitHydrogenCount();
        if (implicit != null) {
            order += implicit;
        }
        return order - 2 < 0.1;
    }
}

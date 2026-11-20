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
package insilico.toxtree.chem;

import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;

/**
 * Takes a group off a molecule and returns what is left.
 * <p>
 * The bond that joins the group to the rest is broken and both ends are filled with hydrogen, then the molecule
 * is split into whatever pieces that leaves. Question 4 uses it to judge a salt as the free acid, question 9 to
 * judge a lactone as the hydroxy acid it would open into.
 * <p>
 * The hydrogen is counted as the original counts it, by subtracting the neighbours an atom now has from the
 * valency it was given, and that is worth knowing because it does not always give a sensible molecule. The
 * carbon of a benzenesulphonate had two neighbours in the ring and the sulphur; losing the sulphur leaves it
 * with two, and four minus two is two hydrogens, so the benzene that comes back carries seven. Nothing will call
 * that ring aromatic afterwards, and the questions that follow are answered accordingly. It is reproduced
 * because it is what the published models do.
 */
public final class Detachment {

    private Detachment() {
    }

    /**
     * Breaks the molecule where the named group attaches and returns the pieces.
     *
     * @param mol      a molecule with explicit hydrogens
     * @param smarts   the group, written so that one of its atoms is the attachment point
     * @param anchor   position of the attachment atom in the pattern
     * @param attached position of the group atom it is bonded to
     * @return the pieces, or null when the group is not there
     */
    public static List<IAtomContainer> detach(IAtomContainer mol, String smarts, int anchor, int attached) {
        List<int[]> bonds = new ArrayList<>();
        Set<IBond> taken = Collections.newSetFromMap(new IdentityHashMap<>());
        for (int[] match : Smarts.of(smarts).matchAll(mol).uniqueAtoms()) {
            IBond bond = mol.getBond(mol.getAtom(match[anchor]), mol.getAtom(match[attached]));
            if (bond == null || !taken.add(bond)) {
                continue;
            }
            bonds.add(new int[]{match[anchor], match[attached]});
        }
        if (bonds.isEmpty()) {
            return null;
        }

        IAtomContainer broken;
        try {
            broken = mol.clone();
        } catch (CloneNotSupportedException x) {
            return null;
        }

        Set<Integer> anchors = new HashSet<>();
        Set<Integer> ends = new HashSet<>();
        for (int[] pair : bonds) {
            anchors.add(pair[0]);
            IBond bond = broken.getBond(broken.getAtom(pair[0]), broken.getAtom(pair[1]));
            if (bond == null) {
                continue;
            }
            broken.removeBond(bond);
            ends.add(pair[0]);
            ends.add(pair[1]);
        }
        for (int end : ends) {
            fillWithHydrogen(broken, broken.getAtom(end));
        }

        // Only the pieces the group hung off are kept. The group itself is not judged again: taking a
        // sulphonate off benzenesulphonic acid leaves benzene to carry on with, not the sulphonic acid.
        List<IAtomContainer> pieces = new ArrayList<>();
        IAtomContainerSet parts = ConnectivityChecker.partitionIntoMolecules(broken);
        for (int i = 0; i < parts.getAtomContainerCount(); i++) {
            IAtomContainer part = parts.getAtomContainer(i);
            for (IAtom atom : part.atoms()) {
                if (anchors.contains(broken.indexOf(atom))) {
                    pieces.add(part);
                    break;
                }
            }
        }
        return pieces;
    }

    /** Adds as many hydrogens as the atom's valency says are now missing. */
    private static void fillWithHydrogen(IAtomContainer mol, IAtom atom) {
        Integer valency = atom.getValency();
        if (valency == null) {
            return;
        }
        int missing = valency - mol.getConnectedBondsCount(atom);
        for (int i = 0; i < missing; i++) {
            IAtom hydrogen = mol.getBuilder().newAtom();
            hydrogen.setSymbol("H");
            hydrogen.setAtomicNumber(1);
            hydrogen.setImplicitHydrogenCount(0);
            mol.addAtom(hydrogen);
            mol.addBond(mol.indexOf(atom), mol.indexOf(hydrogen), IBond.Order.SINGLE);
        }
    }
}

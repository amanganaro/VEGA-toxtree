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
import java.util.List;

/**
 * Breaks a molecule at the bonds the body is expected to break: an azo linkage, an imine, and the bonds around a
 * tertiary amine. Question 33 uses it to decide whether every piece of a large sulphonated dye would still carry
 * a sulphonate of its own.
 * <p>
 * As in the original, the pieces are broken again in their turn until nothing more can be, and the ends are
 * capped with hydrogen.
 */
public final class Metabolism {

    /** @param smarts the linkage, @param first and @param second the positions of the bond to break */
    private record Reaction(String smarts, int first, int second) {
    }

    private static final List<Reaction> REACTIONS = List.of(
            new Reaction("[*]-[#7]=[#7]-[*]", 1, 2),                    // azo
            new Reaction("[*]-[#7](-[*])=[#6](-[*])-[*]", 1, 3),        // imine
            new Reaction("[#6]-[#7](-[#6])-[#6]", 0, 1));               // tertiary amine

    /** How deep the breaking may go, so that a pathological structure cannot run away. */
    private static final int MAX_DEPTH = 10;

    private Metabolism() {
    }

    /**
     * The pieces the molecule breaks into, or null when none of the linkages is present.
     *
     * @param mol a molecule with explicit hydrogens
     */
    public static List<IAtomContainer> products(IAtomContainer mol) {
        return products(mol, 0);
    }

    private static List<IAtomContainer> products(IAtomContainer mol, int depth) {
        if (depth >= MAX_DEPTH) {
            return null;
        }
        for (Reaction reaction : REACTIONS) {
            List<IAtomContainer> pieces = apply(mol, reaction);
            if (pieces == null) {
                continue;
            }
            List<IAtomContainer> all = new ArrayList<>();
            for (IAtomContainer piece : pieces) {
                List<IAtomContainer> further = products(piece, depth + 1);
                if (further == null) {
                    all.add(piece);
                } else {
                    all.addAll(further);
                }
            }
            return all;
        }
        return null;
    }

    private static List<IAtomContainer> apply(IAtomContainer mol, Reaction reaction) {
        List<int[]> bonds = new ArrayList<>();
        for (int[] match : Smarts.of(reaction.smarts()).matchAll(mol).uniqueAtoms()) {
            bonds.add(new int[]{match[reaction.first()], match[reaction.second()]});
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

        for (int[] pair : bonds) {
            IAtom one = broken.getAtom(pair[0]);
            IAtom other = broken.getAtom(pair[1]);
            IBond bond = broken.getBond(one, other);
            if (bond == null) {
                continue;
            }
            broken.removeBond(bond);
            addHydrogen(broken, one);
            addHydrogen(broken, other);
        }

        List<IAtomContainer> pieces = new ArrayList<>();
        IAtomContainerSet parts = ConnectivityChecker.partitionIntoMolecules(broken);
        for (int i = 0; i < parts.getAtomContainerCount(); i++) {
            pieces.add(parts.getAtomContainer(i));
        }
        return pieces.size() > 1 ? pieces : null;
    }

    private static void addHydrogen(IAtomContainer mol, IAtom to) {
        IAtom hydrogen = mol.getBuilder().newAtom();
        hydrogen.setSymbol("H");
        hydrogen.setAtomicNumber(1);
        hydrogen.setImplicitHydrogenCount(0);
        mol.addAtom(hydrogen);
        mol.addBond(mol.indexOf(to), mol.indexOf(hydrogen), IBond.Order.SINGLE);
    }
}

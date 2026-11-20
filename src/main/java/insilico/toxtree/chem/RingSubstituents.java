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
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Cuts a ring out of a molecule and returns what hangs off it.
 * <p>
 * Several rules are phrased as "a ring whose substituents are all of this kind". They are answered by removing
 * the ring's atoms, splitting what is left into connected pieces, and judging each piece on its own. A lone
 * hydrogen is not a substituent and is dropped, so an unsubstituted ring yields nothing to judge and passes.
 * <p>
 * The pieces keep the atom objects of the molecule they came from, so marks worked out on the whole molecule
 * still apply to them.
 */
public final class RingSubstituents {

    private RingSubstituents() {
    }

    /**
     * The substituents of one ring.
     *
     * @param mol  the molecule
     * @param ring the ring to remove, as perceived on that same molecule
     */
    public static List<IAtomContainer> of(IAtomContainer mol, IAtomContainer ring) {
        // Positions, not object identity: the rings CDK hands back wrap the molecule's atoms rather than being
        // the very same objects, so == would never match while indexOf does.
        Set<Integer> ringAtoms = new HashSet<>();
        for (IAtom atom : ring.atoms()) {
            int index = mol.indexOf(atom);
            if (index >= 0) {
                ringAtoms.add(index);
            }
        }
        return of(mol, ringAtoms);
    }

    /**
     * The substituents of a ring given by the positions of its atoms, for callers that work on a view of the
     * molecule the ring was not perceived on. The two must hold the same atoms in the same order.
     */
    public static List<IAtomContainer> of(IAtomContainer mol, Set<Integer> ringAtoms) {

        IAtomContainer rest = SilentChemObjectBuilder.getInstance().newAtomContainer();
        for (int i = 0; i < mol.getAtomCount(); i++) {
            if (!ringAtoms.contains(i)) {
                rest.addAtom(mol.getAtom(i));
            }
        }
        for (IBond bond : mol.bonds()) {
            if (!ringAtoms.contains(mol.indexOf(bond.getBegin()))
                    && !ringAtoms.contains(mol.indexOf(bond.getEnd()))) {
                rest.addBond(bond);
            }
        }

        List<IAtomContainer> substituents = new ArrayList<>();
        IAtomContainerSet parts = ConnectivityChecker.partitionIntoMolecules(rest);
        for (int i = 0; i < parts.getAtomContainerCount(); i++) {
            IAtomContainer part = parts.getAtomContainer(i);
            if (part.getAtomCount() == 1 && "H".equals(part.getAtom(0).getSymbol())) {
                continue;
            }
            substituents.add(part);
        }
        return substituents;
    }
}

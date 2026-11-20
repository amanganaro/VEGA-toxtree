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

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.ArrayList;
import java.util.List;

/**
 * Computes the {@link MolFlags} of a molecule and marks its ring atoms and bonds.
 * <p>
 * This replaces Toxtree's {@code MolAnalyser}, but only the part the rules actually read. Atom typing, hydrogen
 * handling and aromaticity perception are <em>not</em> done here: the molecule is expected to arrive already
 * normalized, as VEGA's {@code InsilicoMolecule} delivers it. That is the deliberate difference from Toxtree,
 * which made hydrogens explicit and ran its own Hückel detector.
 */
public final class Preprocessing {

    private Preprocessing() {
    }

    /**
     * Perceives the rings, marks ring atoms and bonds and derives the structural flags.
     * <p>
     * Rings come from {@link Cycles#all(IAtomContainer)}, with the SSSR set as a fallback when the complete set
     * is intractable: the same choice, and the same fallback, as the original.
     */
    public static MolFlags analyse(IAtomContainer mol) {
        IRingSet ringSet;
        try {
            ringSet = Cycles.all(mol).toRingSet();
        } catch (Intractable x) {
            ringSet = Cycles.sssr(mol).toRingSet();
        }

        int ringCount = ringSet.getAtomContainerCount();
        if (ringCount == 0) {
            boolean acetylenic = false;
            for (IBond bond : mol.bonds()) {
                if (bond.getOrder() == IBond.Order.TRIPLE) {
                    acetylenic = true;
                    break;
                }
            }
            // No ring: aliphatic unless a triple bond is present. Toxtree's own definition.
            return new MolFlags(true, !acetylenic, acetylenic, false, false, false, false, 0, ringSet,
                    List.of(), List.of());
        }

        boolean heterocyclic = false;
        boolean heterocyclic3 = false;
        boolean heteroaromatic = false;
        int aromaticRingCount = 0;
        List<IRing> heterocyclicRings = new ArrayList<>();
        List<IRing> aromaticRings = new ArrayList<>();

        for (int i = 0; i < ringCount; i++) {
            IRing ring = (IRing) ringSet.getAtomContainer(i);

            for (IBond bond : ring.bonds()) {
                bond.setIsInRing(true);
            }

            int aromaticAtoms = 0;
            boolean heterocyclicRing = false;
            for (IAtom atom : ring.atoms()) {
                atom.setIsInRing(true);
                if (atom.isAromatic()) {
                    aromaticAtoms++;
                }
                if (!"C".equals(atom.getSymbol())) {
                    heterocyclicRing = true;
                }
            }

            boolean aromaticRing = aromaticAtoms == ring.getAtomCount();
            ring.setFlag(CDKConstants.ISAROMATIC, aromaticRing);

            if (heterocyclicRing) {
                heterocyclic = true;
                heterocyclicRings.add(ring);
                if (ring.getAtomCount() == 3) {
                    heterocyclic3 = true;
                }
            }
            if (aromaticRing) {
                aromaticRingCount++;
                aromaticRings.add(ring);
                if (heterocyclicRing) {
                    heteroaromatic = true;
                }
            }
        }

        return new MolFlags(false, false, false, aromaticRingCount > 0, heterocyclic, heterocyclic3,
                heteroaromatic, aromaticRingCount, ringSet, heterocyclicRings, aromaticRings);
    }
}

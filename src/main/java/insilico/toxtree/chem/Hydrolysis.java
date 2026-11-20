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
 * Splits a molecule the way water would, for the Cramer questions that ask whether a substance is readily
 * hydrolysed and what it hydrolyses to.
 * <p>
 * Toxtree keeps six reactions as CML files and applies the first one that matches: the esters of carboxylic and
 * thiocarboxylic acids in their four combinations of oxygen and sulphur, phosphate esters and sulphate esters.
 * They are written here as patterns instead, each naming the bond to break. The rule for which bond that is comes
 * from the original: the first single bond of the pattern whose two atoms are both named elements rather than
 * anonymous attachment points.
 * <p>
 * Breaking it adds a hydroxyl to the side that carries a double bond and a hydrogen to the other, so an ester
 * gives back its acid and its alcohol.
 */
public final class Hydrolysis {

    /**
     * @param smarts the group that can be hydrolysed
     * @param first  position in the pattern of one atom of the bond to break
     * @param second position of the other
     */
    private record Reaction(String smarts, int first, int second) {
    }

    private static final List<Reaction> REACTIONS = List.of(
            new Reaction("[*]-[#6](-[#8]-[*])=[#8]", 1, 2),    // ester of a carboxylic acid
            new Reaction("[*]-[#6](-[#16]-[*])=[#8]", 1, 2),   // thioester
            new Reaction("[*]-[#6](-[#8]-[*])=[#16]", 1, 2),   // thiono ester
            new Reaction("[*]-[#6](-[#16]-[*])=[#16]", 1, 2),  // dithio ester
            new Reaction("[*]=[#15](-[#8]-[*])(-[#8]-[*])-[#8]-[*]", 1, 2), // phosphate ester
            new Reaction("[*]-[#16](-[#8]-[*])(=[#8])=[#8]", 1, 2));        // sulphate ester

    private Hydrolysis() {
    }

    /**
     * The fragments the molecule falls into, or null when none of the reactions applies.
     * <p>
     * Only the first reaction that matches is used, as in the original, but every occurrence of it is split.
     *
     * @param mol a molecule with explicit hydrogens
     */
    public static List<IAtomContainer> products(IAtomContainer mol) {
        for (Reaction reaction : REACTIONS) {
            List<IAtomContainer> products = apply(mol, reaction);
            if (products != null) {
                return products;
            }
        }
        return null;
    }

    private static List<IAtomContainer> apply(IAtomContainer mol, Reaction reaction) {
        List<int[]> bondsToBreak = new ArrayList<>();
        Set<IBond> taken = Collections.newSetFromMap(new IdentityHashMap<>());

        for (int[] match : Smarts.of(reaction.smarts()).matchAll(mol).uniqueAtoms()) {
            List<IBond> induced = inducedBonds(mol, match);
            boolean overlaps = false;
            for (IBond bond : induced) {
                if (taken.contains(bond)) {
                    overlaps = true;
                    break;
                }
            }
            if (overlaps) {
                continue;
            }
            taken.addAll(induced);
            bondsToBreak.add(new int[]{match[reaction.first()], match[reaction.second()]});
        }
        if (bondsToBreak.isEmpty()) {
            return null;
        }

        IAtomContainer broken;
        try {
            broken = mol.clone();
        } catch (CloneNotSupportedException x) {
            return null;
        }

        for (int[] pair : bondsToBreak) {
            IAtom one = broken.getAtom(pair[0]);
            IAtom other = broken.getAtom(pair[1]);
            IBond bond = broken.getBond(one, other);
            if (bond == null) {
                continue;
            }
            broken.removeBond(bond);

            // The hydroxyl goes to whichever side carries a double bond: the carbonyl carbon of an ester, the
            // phosphorus of a phosphate.
            IAtom hydroxylSide = hasDoubleBond(broken, one) ? one : other;
            IAtom hydrogenSide = hydroxylSide == one ? other : one;

            addHydrogen(broken, hydrogenSide);
            IAtom oxygen = broken.getBuilder().newAtom();
            oxygen.setSymbol("O");
            oxygen.setAtomicNumber(8);
            oxygen.setImplicitHydrogenCount(0);
            broken.addAtom(oxygen);
            broken.addBond(broken.indexOf(hydroxylSide), broken.indexOf(oxygen), IBond.Order.SINGLE);
            addHydrogen(broken, oxygen);
        }

        List<IAtomContainer> products = new ArrayList<>();
        IAtomContainerSet parts = ConnectivityChecker.partitionIntoMolecules(broken);
        for (int i = 0; i < parts.getAtomContainerCount(); i++) {
            products.add(parts.getAtomContainer(i));
        }
        return products;
    }

    private static void addHydrogen(IAtomContainer mol, IAtom to) {
        IAtom hydrogen = mol.getBuilder().newAtom();
        hydrogen.setSymbol("H");
        hydrogen.setAtomicNumber(1);
        hydrogen.setImplicitHydrogenCount(0);
        mol.addAtom(hydrogen);
        mol.addBond(mol.indexOf(to), mol.indexOf(hydrogen), IBond.Order.SINGLE);
    }

    private static boolean hasDoubleBond(IAtomContainer mol, IAtom atom) {
        for (IBond bond : mol.getConnectedBondsList(atom)) {
            if (bond.getOrder() == IBond.Order.DOUBLE) {
                return true;
            }
        }
        return false;
    }

    private static List<IBond> inducedBonds(IAtomContainer mol, int[] match) {
        Set<Integer> atoms = new HashSet<>();
        for (int atom : match) {
            atoms.add(atom);
        }
        List<IBond> bonds = new ArrayList<>();
        for (IBond bond : mol.bonds()) {
            if (atoms.contains(mol.indexOf(bond.getBegin())) && atoms.contains(mol.indexOf(bond.getEnd()))) {
                bonds.add(bond);
            }
        }
        return bonds;
    }
}

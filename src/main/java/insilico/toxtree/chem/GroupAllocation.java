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

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Works out which functional group each atom of a molecule belongs to, and answers the question several rules are
 * built on: is every heteroatom of this fragment accounted for by one of the groups we allow?
 * <p>
 * This is Toxtree's marking mechanism. The original writes a property on each atom as it matches a group and then
 * walks the fragment looking for an atom that carries none; the marks are kept here instead, which leaves the
 * caller's molecule untouched and lets two threads classify at once.
 * <p>
 * Carbon and hydrogen are never questioned: only heteroatoms have to belong to an allowed group. That is why a
 * purely alkyl substituent satisfies every such rule, and why the carbon marks the original also records make no
 * difference to the answer.
 */
public final class GroupAllocation {

    private final IAtomContainer molecule;
    private final List<Set<String>> marksByAtom;

    private GroupAllocation(IAtomContainer molecule, List<Set<String>> marksByAtom) {
        this.molecule = molecule;
        this.marksByAtom = marksByAtom;
    }

    /** Matches every group against the molecule and records what each atom was allocated to. */
    public static GroupAllocation of(IAtomContainer mol, Collection<FunctionalGroup> groups) {
        int count = mol.getAtomCount();
        List<Set<String>> marks = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            marks.add(null);
        }
        for (FunctionalGroup group : groups) {
            Set<IBond> taken = Collections.newSetFromMap(new IdentityHashMap<>());
            for (int[] match : Smarts.of(group.smarts()).matchAll(mol).uniqueAtoms()) {
                // Toxtree throws away matches of the same group that overlap, and it counts them by bond, not
                // by atom. Both halves of that matter. On a hexahalocyclohexane the matches of "halogen on a
                // ring carbon" share ring bonds, so only the first survives and the other halogens are left
                // unaccounted for, which is why Toxtree rejects that molecule. On a trichloromethyl group the
                // three matches of "halogen" share the carbon but no bond, so all three are kept.
                List<IBond> bonds = inducedBonds(mol, match);
                boolean overlaps = false;
                for (IBond bond : bonds) {
                    if (taken.contains(bond)) {
                        overlaps = true;
                        break;
                    }
                }
                if (overlaps) {
                    continue;
                }
                taken.addAll(bonds);
                for (int position = 0; position < match.length; position++) {
                    if (group.isUnmarked(position)) {
                        continue;
                    }
                    int atom = match[position];
                    Set<String> set = marks.get(atom);
                    if (set == null) {
                        set = new LinkedHashSet<>(2);
                        marks.set(atom, set);
                    }
                    set.add(group.id());
                }
            }
        }
        return new GroupAllocation(mol, marks);
    }

    /** The bonds of the molecule that run between two atoms of the same match. */
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

    /**
     * True when every heteroatom of the fragment belongs to at least one of the named groups.
     * <p>
     * A fragment without heteroatoms satisfies this trivially, which is deliberate: that is how an alkyl
     * substituent is accepted by rules that otherwise only allow, say, a nitro group.
     */
    public boolean onlyTheseGroups(IAtomContainer fragment, Collection<String> ids) {
        for (IAtom atom : fragment.atoms()) {
            String symbol = atom.getSymbol();
            if ("H".equals(symbol) || "C".equals(symbol)) {
                continue;
            }
            Set<String> marks = marksOf(atom);
            if (marks == null) {
                return false;
            }
            boolean allowed = false;
            for (String id : ids) {
                if (marks.contains(id)) {
                    allowed = true;
                    break;
                }
            }
            if (!allowed) {
                return false;
            }
        }
        return true;
    }

    /** The groups an atom was allocated to, or null when it was allocated to none. */
    public Set<String> marksOf(IAtom atom) {
        // By position: CDK sometimes hands out wrappers around an atom rather than the atom itself, and those
        // compare equal without being identical.
        int index = molecule.indexOf(atom);
        return index < 0 ? null : marksByAtom.get(index);
    }
}

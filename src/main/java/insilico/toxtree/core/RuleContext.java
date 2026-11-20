/*
 * Copyright (C) 2026 ProtoQSAR / Alberto Manganaro
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 */
package insilico.toxtree.core;

import insilico.toxtree.chem.MolFlags;
import insilico.toxtree.chem.Preprocessing;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Supplier;

/**
 * Everything the rules of one classification share: the structural flags and a scratch area for results that are
 * expensive to compute and used by more than one rule.
 * <p>
 * One instance per molecule. Rules are stateless; all mutable state lives here. This is what lets the library be
 * used from several threads, which the original could not do.
 */
public final class RuleContext {

    private final IAtomContainer molecule;
    private final MolFlags flags;
    private final Map<String, Object> cache = new HashMap<>();
    private IAtomContainer explicitHydrogens;
    private MolFlags explicitFlags;
    private IAtomContainer kekule;

    public RuleContext(IAtomContainer mol) {
        this.molecule = mol;
        this.flags = Preprocessing.analyse(mol);
    }

    /** The structural flags of the molecule. */
    public MolFlags flags() {
        return flags;
    }

    /**
     * The molecule to run SMARTS against: a copy in which hydrogens are explicit atoms.
     * <p>
     * This is needed for fidelity, not for convenience. Toxtree makes every hydrogen explicit before it starts,
     * and a good number of its patterns rely on it: {@code C=;!:CC(O[H])} names a hydrogen atom outright, and
     * {@code [OH1]C[!$(C=C);!$(C#C)]} only matches methanol because the third atom can be a hydrogen. Matching
     * those against a molecule with implicit hydrogens silently answers "no".
     * <p>
     * The copy is made once per molecule and keeps the aromaticity and ring flags of the original, which is what
     * lets the patterns run without re-preparation.
     */
    /**
     * The structural flags of {@link #withExplicitHydrogens()}.
     * <p>
     * Rules that take a ring apart need rings perceived on the very molecule they will cut, because the pieces
     * must share its atom objects. Hydrogens never sit in a ring, so these are the same rings as
     * {@link #flags()} describes, only expressed over the copy.
     */
    public MolFlags explicitFlags() {
        if (explicitFlags == null) {
            explicitFlags = Preprocessing.analyse(withExplicitHydrogens());
        }
        return explicitFlags;
    }

    /**
     * The molecule with explicit hydrogens and the aromatic flags cleared, so that bond orders are read the way
     * Toxtree reads them.
     * <p>
     * Many of Toxtree's functional groups are built from bonds that ask only for an order, single or double,
     * and take no notice of whether the bond is aromatic. SMARTS is stricter: {@code -} does not match an
     * aromatic bond. Pyridinium is the case that shows it, where Toxtree finds a quaternary nitrogen by reading
     * the Kekulé orders of the ring. Clearing the flags reproduces that reading exactly, and it must not be
     * used for patterns that ask about aromaticity.
     */
    public IAtomContainer kekule() {
        if (kekule == null) {
            try {
                IAtomContainer copy = withExplicitHydrogens().clone();
                for (IAtom atom : copy.atoms()) {
                    atom.setIsAromatic(false);
                }
                for (IBond bond : copy.bonds()) {
                    bond.setIsAromatic(false);
                }
                kekule = copy;
            } catch (CloneNotSupportedException x) {
                kekule = withExplicitHydrogens();
            }
        }
        return kekule;
    }

    public IAtomContainer withExplicitHydrogens() {
        if (explicitHydrogens == null) {
            try {
                IAtomContainer copy = molecule.clone();
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(copy);
                explicitHydrogens = copy;
            } catch (CloneNotSupportedException x) {
                explicitHydrogens = molecule;
            }
        }
        return explicitHydrogens;
    }

    /**
     * Returns a cached value, computing it on first request. Used for things such as the hydrolysis products,
     * which several rules ask for in turn.
     */
    @SuppressWarnings("unchecked")
    public <T> T cached(String key, Supplier<T> supplier) {
        if (cache.containsKey(key)) {
            return (T) cache.get(key);
        }
        T value = supplier.get();
        cache.put(key, value);
        return value;
    }
}

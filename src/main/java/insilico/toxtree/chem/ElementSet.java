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

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.TreeSet;

/**
 * The set of chemical elements of a molecule, or the set a rule allows, with the two conventions Toxtree's
 * {@code SmartElementsList} uses. Both are easy to miss and both change the answers:
 * <ul>
 *   <li><b>Hydrogen is ignored.</b> Adding {@code "H"} to the set does nothing, so hydrogen never takes part in
 *       any comparison. A happy side effect: it makes these rules indifferent to whether hydrogens are explicit
 *       or implicit, which is the main difference between Toxtree's preprocessing and ours.</li>
 *   <li><b>Halogens collapse to a single symbol {@code "X"}</b>, but only those listed as halogens for that
 *       rule. A rule that declares only {@code Cl, Br, F} as halogens keeps iodine as a distinct element, so a
 *       molecule containing iodine is no longer "made only of C, H, N, O, S and halogens". Verhaar's rule 0.1
 *       relies on exactly this.</li>
 * </ul>
 */
public final class ElementSet {

    /** The symbol every recognised halogen collapses to. */
    public static final String HALOGEN = "X";

    /** The halogens assumed when a rule does not say otherwise. */
    public static final Set<String> DEFAULT_HALOGENS = Set.of("F", "Cl", "Br", "I");

    private final Set<String> symbols;
    private final Set<String> halogens;

    private ElementSet(Set<String> symbols, Set<String> halogens) {
        this.symbols = symbols;
        this.halogens = halogens;
    }

    /** The set a rule allows, using the default halogens. */
    public static ElementSet of(String... elements) {
        return of(DEFAULT_HALOGENS, elements);
    }

    /** The set a rule allows, naming which elements count as halogens. */
    public static ElementSet of(Set<String> halogens, String... elements) {
        Set<String> symbols = new TreeSet<>();
        for (String element : elements) {
            add(symbols, halogens, element);
        }
        return new ElementSet(symbols, new TreeSet<>(halogens));
    }

    /** The elements a molecule is made of, under the same conventions. */
    public static ElementSet ofMolecule(IAtomContainer mol, Set<String> halogens) {
        Set<String> symbols = new TreeSet<>();
        for (IAtom atom : mol.atoms()) {
            add(symbols, halogens, atom.getSymbol());
        }
        return new ElementSet(symbols, new TreeSet<>(halogens));
    }

    private static void add(Set<String> symbols, Set<String> halogens, String element) {
        if ("H".equals(element)) {
            return;
        }
        symbols.add(halogens.contains(element) ? HALOGEN : element);
    }

    /** True when this set contains the element, treating any declared halogen as {@code "X"}. */
    public boolean contains(String element) {
        if ("H".equals(element)) {
            return true;
        }
        return symbols.contains(halogens.contains(element) ? HALOGEN : element);
    }

    /** True when every element of {@code other} is allowed by this set. */
    public boolean containsAll(ElementSet other) {
        return symbols.containsAll(other.symbols);
    }

    /** True when the two sets hold exactly the same elements. */
    public boolean isSameAs(ElementSet other) {
        return symbols.equals(other.symbols) && halogens.equals(other.halogens);
    }

    /** True when the two sets have at least one element in common. */
    public boolean intersects(ElementSet other) {
        for (String symbol : other.symbols) {
            if (symbols.contains(symbol)) {
                return true;
            }
        }
        return false;
    }

    public Set<String> getHalogens() {
        return new LinkedHashSet<>(halogens);
    }

    public Set<String> symbols() {
        return new LinkedHashSet<>(symbols);
    }

    @Override
    public String toString() {
        return Arrays.toString(symbols.toArray());
    }
}

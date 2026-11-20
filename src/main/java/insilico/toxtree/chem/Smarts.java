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
package insilico.toxtree.chem;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Compiles and caches SMARTS patterns.
 * <p>
 * Patterns are created with {@code setPrepare(false)}. This matters: with preparation on, CDK recomputes ring
 * membership and aromaticity on the molecule it is given, rewriting state that the caller owns. VEGA hit exactly
 * that problem once, where a model silently changed the molecule for every model that ran after it. The molecules
 * reaching this library are already normalized, so preparation is both unnecessary and harmful.
 */
public final class Smarts {

    private static final Map<String, SmartsPattern> CACHE = new ConcurrentHashMap<>();

    /**
     * A carbon as Toxtree's {@code markCHn} recognises one: either it bears at least one hydrogen, or it bears
     * none and every neighbour it has is a carbon.
     * <p>
     * Several of Toxtree's groups are built from this kind of atom rather than from a plain carbon, and the
     * difference decides real cases: the carbonyl carbon of an amide is not one, so a urea nitrogen is not an
     * aliphatic secondary amine.
     */
    public static final String CHN = "[$([#6]~[#1]),$([#6;!$([#6]~[!#6])])]";

    private Smarts() {
    }

    /**
     * Returns the compiled pattern for a SMARTS string, compiling it on first use.
     *
     * @throws IllegalArgumentException if the SMARTS cannot be parsed
     */
    public static SmartsPattern of(String smarts) {
        return CACHE.computeIfAbsent(smarts, s -> {
            try {
                return SmartsPattern.create(s, SilentChemObjectBuilder.getInstance()).setPrepare(false);
            } catch (Exception x) {
                throw new IllegalArgumentException("Invalid SMARTS: " + s, x);
            }
        });
    }

    /** True when the molecule contains the substructure at least once. */
    public static boolean matches(IAtomContainer mol, String smarts) {
        return of(smarts).matches(mol);
    }

    /** True when the molecule contains every one of the substructures. */
    public static boolean matchesAll(IAtomContainer mol, String... smarts) {
        for (String s : smarts) {
            if (!matches(mol, s)) {
                return false;
            }
        }
        return true;
    }

    /** True when the molecule contains at least one of the substructures. */
    public static boolean matchesAny(IAtomContainer mol, String... smarts) {
        for (String s : smarts) {
            if (matches(mol, s)) {
                return true;
            }
        }
        return false;
    }

    /** Number of distinct matches of the substructure, counting each set of atoms once. */
    public static int count(IAtomContainer mol, String smarts) {
        return of(smarts).matchAll(mol).uniqueAtoms().count();
    }
}

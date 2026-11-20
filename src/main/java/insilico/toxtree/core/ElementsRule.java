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
package insilico.toxtree.core;

import insilico.toxtree.chem.ElementSet;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Set;

/**
 * A rule about which chemical elements the molecule is made of. Seven of the Verhaar rules are of this kind.
 * <p>
 * See {@link ElementSet} for the two conventions that govern the comparison: hydrogen is ignored, and the
 * halogens declared for the rule collapse into one symbol.
 */
public class ElementsRule extends AbstractRule {

    /** How the molecule's elements are compared with the ones the rule names. */
    public enum Mode {
        /** Yes when the molecule contains nothing outside the named elements. */
        ONLY,
        /** Yes when the molecule's elements are exactly the named ones, no more and no fewer. */
        EXACTLY,
        /** Yes when the molecule contains at least one of the named elements. */
        ANY
    }

    private final Mode mode;
    private final ElementSet elements;
    private final Set<String> halogens;

    protected ElementsRule(String id, String title, Mode mode, String... elements) {
        this(id, title, mode, ElementSet.DEFAULT_HALOGENS, elements);
    }

    protected ElementsRule(String id, String title, Mode mode, Set<String> halogens, String... elements) {
        super(id, title);
        this.mode = mode;
        this.halogens = halogens;
        this.elements = ElementSet.of(halogens, elements);
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        ElementSet actual = ElementSet.ofMolecule(mol, halogens);
        return switch (mode) {
            case ONLY -> elements.containsAll(actual);
            case EXACTLY -> actual.isSameAs(elements);
            case ANY -> actual.intersects(elements);
        };
    }
}

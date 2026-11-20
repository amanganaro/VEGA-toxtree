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
package insilico.toxtree.verhaar.rules;

import insilico.toxtree.core.ElementsRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Set;

/**
 * Verhaar rule 1.1: the molecule must not contain iodine. The question is phrased negatively, so the element
 * test is inverted, exactly as in the original.
 */
public final class Rule11 extends ElementsRule {

    public Rule11() {
        super("1.1", "Not contain I", Mode.ANY, Set.of("Cl", "Br", "F"), "I");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        return !super.verify(mol, ctx);
    }
}

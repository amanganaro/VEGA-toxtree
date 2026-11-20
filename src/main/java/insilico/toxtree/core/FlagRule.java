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

import insilico.toxtree.chem.MolFlags;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.function.Predicate;

/**
 * A rule that is answered by one of the structural flags alone, such as "is it aromatic" or "is it
 * heterocyclic". Several of the Cramer questions are of this kind, and the tree reuses some of them under more
 * than one number.
 */
public final class FlagRule extends AbstractRule {

    private final Predicate<MolFlags> question;

    public FlagRule(String id, String title, Predicate<MolFlags> question) {
        super(id, title);
        this.question = question;
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        return question.test(ctx.flags());
    }
}

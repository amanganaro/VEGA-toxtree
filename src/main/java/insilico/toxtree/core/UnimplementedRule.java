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

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * A rule that is declared by the tree but not written yet.
 * <p>
 * It always answers "no", and the trace marks it with {@code -} instead of {@code N}, so that a port in progress
 * can be run end to end and compared with the reference: every molecule whose path reaches an unimplemented rule
 * is visibly incomplete rather than quietly wrong.
 */
public final class UnimplementedRule extends AbstractRule {

    public UnimplementedRule(String id, String title) {
        super(id, title);
    }

    @Override
    public boolean isImplemented() {
        return false;
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        return false;
    }
}

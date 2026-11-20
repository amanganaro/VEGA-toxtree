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

import insilico.toxtree.chem.FunctionalGroup;
import insilico.toxtree.chem.GroupAllocation;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.List;
import java.util.function.Predicate;

/**
 * A rule that asks whether the whole molecule is made of nothing but the groups it allows: every heteroatom has
 * to belong to one of them, carbon and hydrogen being free. This is Toxtree's
 * {@code RuleOnlyAllowedSubstructures}, used by rules 1.4.2 and 1.6.1.
 * <p>
 * A precondition is checked first, since both rules only apply to a certain shape of molecule.
 */
public class CoverageRule extends AbstractRule {

    private final Predicate<RuleContext> precondition;
    private final List<FunctionalGroup> groups;
    private final List<String> allowedIds;

    public CoverageRule(String id, String title, Predicate<RuleContext> precondition,
                           List<FunctionalGroup> groups) {
        super(id, title);
        this.precondition = precondition;
        this.groups = List.copyOf(groups);
        this.allowedIds = groups.stream().map(FunctionalGroup::id).distinct().toList();
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        if (!precondition.test(ctx)) {
            return false;
        }
        IAtomContainer target = ctx.withExplicitHydrogens();
        return GroupAllocation.of(target, groups).onlyTheseGroups(target, allowedIds);
    }
}

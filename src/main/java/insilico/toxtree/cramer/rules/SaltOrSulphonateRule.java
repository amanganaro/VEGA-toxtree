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
package insilico.toxtree.cramer.rules;

import insilico.toxtree.chem.FunctionalGroup;
import insilico.toxtree.chem.GroupAllocation;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.util.List;
import java.util.Set;

/**
 * Cramer question 4: do the elements that question 3 objected to occur only as a harmless salt?
 * <p>
 * Reached only when question 3 found something beyond carbon, hydrogen, oxygen, nitrogen and divalent sulphur,
 * it asks whether that something is a sodium, potassium, calcium, magnesium or ammonium salt of a carboxylic
 * acid, the sulphate or hydrochloride of an amine, or a sulphonate, sulphamate or sulphate. If so the substance
 * is judged as the free acid or amine would be.
 * <p>
 * The original marks every oxygen and nitrogen with its own symbol before it starts, and every divalent sulphur
 * besides, so those never need a group to vouch for them. What has to be accounted for is the rest: the metals,
 * the chlorine of a hydrochloride, and the sulphur of a sulphonate or sulphate, which is not divalent.
 */
public final class SaltOrSulphonateRule extends AbstractRule {

    private static final String METALS = "Na,K,Ca,Mg,N";
    private static final String SALT = "salt";

    private static final List<FunctionalGroup> GROUPS = List.of(
            // salts of a carboxylic acid, as a neutral metal ester, as a separated anion, and the variant the
            // original writes with a carbanion
            FunctionalGroup.of(SALT, "[#6,#1]-[#6](-[#8]-[" + METALS + "])=[#8]", 0),
            FunctionalGroup.of(SALT, "[#6]-[#8-]~[" + METALS + "]"),
            FunctionalGroup.of(SALT, "[#6]-[#6-]~[" + METALS + "]"),
            // hydrochloride of an amine, primary through tertiary
            FunctionalGroup.of(SALT, "[#7+]~[#17-]"),
            // sulphonate, sulphamate, sulphate, and the sulphate of an amine
            FunctionalGroup.of(SALT, "[!#8;!#1;!#7]-[#16](=[#8])(=[#8])-[#8]", 0),
            FunctionalGroup.of(SALT, "[#7]-[#16](=[#8])(=[#8])-[#8]"),
            // a sulphate only counts when it carries a metal: a dialkyl sulphate does not qualify
            FunctionalGroup.of(SALT, "[Na,K,Ca]-[#8]-[#16](=[#8])(=[#8])-[#8]"),
            FunctionalGroup.of(SALT, "[#16](-[#8-]~[#7+])(-[#8])(=[#8])=[#8]"));

    private static final Set<String> NEVER_QUESTIONED = Set.of("C", "H", "O", "N");

    public SaltOrSulphonateRule() {
        super("4", "Elements not listed in question 3 occur only as a salt, sulphonate, sulphamate or sulphate");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.kekule();
        GroupAllocation allocation = GroupAllocation.of(target, GROUPS);

        for (IAtom atom : target.atoms()) {
            String symbol = atom.getSymbol();
            if (NEVER_QUESTIONED.contains(symbol)) {
                continue;
            }
            if ("S".equals(symbol) && isDivalent(target, atom)) {
                continue;
            }
            Set<String> marks = allocation.marksOf(atom);
            if (marks == null || !marks.contains(SALT)) {
                return false;
            }
        }
        return true;
    }

    private static boolean isDivalent(IAtomContainer mol, IAtom atom) {
        double order = 0;
        for (IBond bond : mol.getConnectedBondsList(atom)) {
            switch (bond.getOrder()) {
                case SINGLE -> order += 1;
                case DOUBLE -> order += 2;
                case TRIPLE -> order += 3;
                case QUADRUPLE -> order += 4;
                default -> {
                }
            }
        }
        return order - 2 < 0.1;
    }
}

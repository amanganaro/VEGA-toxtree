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

import insilico.toxtree.chem.Metabolism;
import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.List;

/**
 * Cramer question 33: does the substance carry enough sulphonate or sulphamate groups to be excreted unchanged?
 * <p>
 * It must carry at least one to begin with. Then the molecule is broken at the bonds the body is expected to
 * break, the azo and imine linkages and the tertiary amines, and each piece has to carry a sulphur in a
 * sulphonate or sulphamate of its own and to be no larger than twenty carbons. A substance that cannot be broken
 * at all is judged by its own size.
 */
public final class SulphonateGroupsRule extends AbstractRule {

    private static final int MAX_CARBONS = 20;

    private static final String SULPHONATE = "[!#8;!#1;!#7]-[#16](=[#8])(=[#8])-[#8]";
    private static final String SULPHAMATE = "[#7]-[#16](=[#8])(=[#8])-[#8]";

    public SulphonateGroupsRule() {
        super("33", "Has a sufficient number of sulphonate or sulphamate groups");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.kekule();
        if (!Smarts.matchesAny(target, SULPHONATE, SULPHAMATE)) {
            return false;
        }

        List<IAtomContainer> residues = Metabolism.products(target);
        if (residues == null) {
            return carbons(target) <= MAX_CARBONS;
        }
        for (IAtomContainer residue : residues) {
            if (count(residue, "S") == 0) {
                return false;
            }
            if (!Smarts.matchesAny(residue, SULPHONATE, SULPHAMATE)) {
                return false;
            }
            if (carbons(residue) > MAX_CARBONS) {
                return false;
            }
        }
        return true;
    }

    private static int carbons(IAtomContainer mol) {
        return count(mol, "C");
    }

    private static int count(IAtomContainer mol, String symbol) {
        int n = 0;
        for (IAtom atom : mol.atoms()) {
            if (symbol.equals(atom.getSymbol())) {
                n++;
            }
        }
        return n;
    }
}

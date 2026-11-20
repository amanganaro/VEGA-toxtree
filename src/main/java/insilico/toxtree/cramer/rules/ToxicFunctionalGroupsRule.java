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

import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Cramer question 2: does the substance contain a functional group associated with enhanced toxicity?
 * <p>
 * An aliphatic secondary amine, a nitro group, a cyano, an N-nitroso, a diazo, a triazeno or a quaternary
 * nitrogen. The quaternary nitrogen is the only one with a way out: it does not count when it appears as
 * {@code >C=N+R2}, or as the hydrochloride or sulphate salt of a primary or tertiary amine.
 */
public final class ToxicFunctionalGroupsRule extends AbstractRule {

    private static final String SECONDARY_AMINE_ALIPHATIC =
            Smarts.CHN + "-!@[#7+0](-!@" + Smarts.CHN + ")-[#1]";
    private static final String NITRO_CHARGED = "[!#1]-[#7+](=[#8])-[#8-]";
    private static final String NITRO_NEUTRAL = "[!#1]-[#7](=[#8])=[#8]";
    private static final String CYANO = "[#7,#6,#1]-[#6]#[#7]";
    private static final String N_NITROSO = "[!#1]-[#7](-[!#1])-[#7]=[#8]";
    private static final String DIAZO = "[#7]#[#7]";
    private static final String TRIAZENO = "[#6]-[#7]=[#7]-[#7](-[#1])-[#1]";

    private static final String[] QUATERNARY_NITROGEN = {
            "[*]-[#7+](-[*])(-[*])-[*]",
            "[*]-[#7](-[*])(-[*])-[*]",
            "[*]-[#7+](-[*])=[*]",
            "[*]-[#7](-[*])=[*]",
    };

    /** The {@code >C=N+R2} form, which is not counted as a toxic quaternary nitrogen. */
    private static final String QUATERNARY_NITROGEN_EXCEPTION = "[*]-[#7+](-[*])=[#6](-[!#1])-[!#1]";

    /** Salts of a primary or a tertiary amine, which are not counted either. */
    private static final String[] AMINE_SALTS = {
            // hydrochloride, primary then tertiary
            "[#6]-[#7+](-[#1])(-[#1])(-[#1])~[#17-]",
            "[#6]-[#7+](-[#6])(-[#6])(-[#1])~[#17-]",
            // sulphate, primary then tertiary
            "[#16](-[#8-]~[#7+](-[#6])(-[#1])(-[#1])-[#1])(-[#8])(=[#8])=[#8]",
            "[#16](-[#8-]~[#7+](-[#6])(-[#6])(-[#6])-[#1])(-[#8])(=[#8])=[#8]",
    };

    public ToxicFunctionalGroupsRule() {
        super("2", "Contains functional groups associated with enhanced toxicity");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        // Order-only groups: read on the Kekule view, as Toxtree does.
        IAtomContainer target = ctx.kekule();

        if (ctx.flags().isAliphatic() && Smarts.matches(target, SECONDARY_AMINE_ALIPHATIC)) {
            return true;
        }
        if (Smarts.matchesAny(target, NITRO_CHARGED, NITRO_NEUTRAL, CYANO, N_NITROSO, DIAZO, TRIAZENO)) {
            return true;
        }
        if (!Smarts.matchesAny(target, QUATERNARY_NITROGEN)) {
            return false;
        }
        if (Smarts.matches(target, QUATERNARY_NITROGEN_EXCEPTION)) {
            return false;
        }
        return !Smarts.matchesAny(target, AMINE_SALTS);
    }
}

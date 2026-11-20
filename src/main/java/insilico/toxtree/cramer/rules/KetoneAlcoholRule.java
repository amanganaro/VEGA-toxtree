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
import insilico.toxtree.chem.Smarts;
import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Cramer question 18: is it a ketone, alcohol, acid or ester of the kinds the tree treats as more than simply
 * inert?
 * <p>
 * The question is answered by a list of named structures rather than by a principle. Some are whole substances
 * in their own right, allyl alcohol or acrolein, and the substance has to <em>be</em> one of them; the rest are
 * groups the substance has to contain. At the end come the acetylenic compounds, the open-chain ketones and
 * ketals with at least four carbons on one side, and the sterically hindered rings.
 */
public final class KetoneAlcoholRule extends AbstractRule {

    /** Substances the molecule has to be, not merely contain. */
    private static final String[] SUBSTANCES = {
            "C=CCO",        // allyl alcohol
            "SCC=C",        // allyl mercaptan
            "C=CCSCC=C",    // allyl sulphide
            "C=CCN",        // allyl amine
            "C=CC=O",       // acrolein
            "C=C(C)C=O",    // methacrolein
            "C=CC(O)=O",    // acrylic acid
            "C=C(C)C(O)=O", // methacrylic acid
    };

    private static final String[] GROUPS = {
            // vicinal diketone
            "[#8]=[#6](-[!#8;!#16])-[#6](=[#8])-[!#8;!#16]",
            // ketone, ketal, secondary alcohol or its ester, on a terminal vinyl
            "[!#8;!#16;!#7]-[#6](=[#8])-[#6]([#1])=[#6]([#1])[#1]",
            "[!#8;!#1]-[#8]-[#6](-[#8]-[!#8;!#1])(-[!#8;!#1])-[#6]=[#6]",
            "[#1]-[#8]-[#6]([#1])(-[#6])-[#6]([#1])=[#6]([#1])[#1]",
            "[#6]-[#6](-[#6](=[#8])-[#8]-[#6])-[#6]=[#6]",
            // allyl alcohol as an acetal or an ester
            "[#6](-[#8]-[#6]-[#6]=[#6])(-[#8]-[#6]-[#6]=[#6])[#1]",
            "[#8]=[#6](-[#6])-[#8]-[#6]([#1])([#1])-[#6]([#1])=[#6]([#1])[#1]",
            // allyl thioester
            "[#6](=[#8])-[#16]-[#6]-[#6]=[#6]",
            // methacrolein acetal
            "[#6]=[#6](-[#6])-[#6]([#1])(-[#8]-[#6])-[#8]-[#6]",
            // a ring whose every carbon carries something other than hydrogen
            "[#6]1(-[!#1])~[#6](-[!#1])~[#6](-[!#1])~[#6](-[!#1])~[#6](-[!#1])~[#6]1-[!#1]",
    };

    /** An open-chain ketone or ketal with at least four carbons on one side. */
    private static final String KETONE_4 = "[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#6]";
    private static final String KETAL_4 = "[#6]-[#6]-[#6]-[#6]-[#6](-[#8]-[#6])(-[#8]-[#6])-[#6]";
    private static final String ALIPHATIC_CHN = "[$([#6;!a]~[#1]),$([#6;!a;!$([#6]~[!#6])])]";

    private static List<String> substanceInchis;

    public KetoneAlcoholRule() {
        super("18", "Ketone, alcohol, acid or ester of a reactive kind");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        IAtomContainer target = ctx.kekule();

        if (Smarts.matchesAny(target, GROUPS)) {
            return true;
        }
        if (isOneOfTheNamedSubstances(mol)) {
            return true;
        }
        if (ctx.flags().isAcetylenic()) {
            return true;
        }
        if (!ctx.flags().isAliphatic()) {
            return false;
        }

        // Open chain: a long ketone or ketal is enough, on its own or with an alcohol.
        FunctionalGroup ketone = FunctionalGroup.of("ketone4", KETONE_4);
        FunctionalGroup ketal = FunctionalGroup.of("ketal4", KETAL_4);
        FunctionalGroup alcohol = FunctionalGroup.of(CramerGroups.ALCOHOL, ALIPHATIC_CHN + "-[#8]-[#1]", 0);

        if (Smarts.matches(target, KETONE_4)) {
            GroupAllocation alone = GroupAllocation.of(target, List.of(ketone));
            if (alone.onlyTheseGroups(target, List.of(ketone.id()))) {
                return true;
            }
            GroupAllocation withAlcohol = GroupAllocation.of(target, List.of(ketone, alcohol));
            if (Smarts.matches(target, alcohol.smarts())
                    && withAlcohol.onlyTheseGroups(target, List.of(ketone.id(), alcohol.id()))) {
                return true;
            }
        }
        if (Smarts.matches(target, KETAL_4)) {
            GroupAllocation alone = GroupAllocation.of(target, List.of(ketal));
            return alone.onlyTheseGroups(target, List.of(ketal.id()));
        }
        return false;
    }

    private boolean isOneOfTheNamedSubstances(IAtomContainer mol) {
        try {
            String inchi = InChIGeneratorFactory.getInstance().getInChIGenerator(mol).getInchi();
            return inchi != null && substanceInchis().contains(inchi);
        } catch (Exception x) {
            return false;
        }
    }

    private static synchronized List<String> substanceInchis() {
        if (substanceInchis == null) {
            List<String> inchis = new ArrayList<>();
            SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
            for (String smiles : SUBSTANCES) {
                try {
                    IAtomContainer substance = parser.parseSmiles(smiles);
                    for (IBond bond : substance.bonds()) {
                        bond.setIsAromatic(false);
                    }
                    inchis.add(InChIGeneratorFactory.getInstance().getInChIGenerator(substance).getInchi());
                } catch (Exception x) {
                    // a substance we cannot express is simply never matched
                }
            }
            substanceInchis = inchis;
        }
        return substanceInchis;
    }

    static {
        Set<String> seen = new HashSet<>();
        for (String smarts : GROUPS) {
            if (seen.add(smarts)) {
                Smarts.of(smarts);
            }
        }
    }
}

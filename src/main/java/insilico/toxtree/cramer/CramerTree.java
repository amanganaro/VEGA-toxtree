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
package insilico.toxtree.cramer;

import insilico.toxtree.chem.MolFlags;
import insilico.toxtree.core.DecisionRule;
import insilico.toxtree.core.DecisionTree;
import insilico.toxtree.core.FlagRule;
import insilico.toxtree.core.UnimplementedRule;
import insilico.toxtree.cramer.rules.AliphaticWithFunctionalGroupsRule;
import insilico.toxtree.cramer.rules.BenzeneDerivativeRule;
import insilico.toxtree.cramer.rules.CommonTerpeneRule;
import insilico.toxtree.cramer.rules.ElementsOtherThanRule;
import insilico.toxtree.cramer.rules.HeterocyclicComplexSubstituentsRule;
import insilico.toxtree.cramer.rules.HydrolysisRules;
import insilico.toxtree.cramer.rules.KetoneAlcoholRule;
import insilico.toxtree.cramer.rules.LactoneRule;
import insilico.toxtree.cramer.rules.ManyFunctionalGroupsRule;
import insilico.toxtree.cramer.rules.MonocarbocyclicRule;
import insilico.toxtree.cramer.rules.RingRules;
import insilico.toxtree.cramer.rules.SaltOrSulphonateRule;
import insilico.toxtree.cramer.rules.SmallRingRules;
import insilico.toxtree.cramer.rules.SimplyBranchedOrCarbohydrateRule;
import insilico.toxtree.cramer.rules.StructureListRule;
import insilico.toxtree.cramer.rules.SulphonateGroupsRule;
import insilico.toxtree.cramer.rules.ToxicFunctionalGroupsRule;

/**
 * The Cramer decision tree, which sorts substances into three classes of presumed toxic hazard.
 * <p>
 * Thirty-nine questions, of which 34 are the original ones and the last six are shortcuts that other branches
 * reuse: questions 28 and 29 repeat 14 and 15, and 34 to 39 are the generic "is it an open chain", "is it
 * heterocyclic", "is it aromatic" and "is it a common terpene" questions borrowed by several branches.
 */
public final class CramerTree extends DecisionTree {

    private static final String[] CATEGORIES = {
            "Low (Class I)",
            "Intermediate (Class II)",
            "High (Class III)"
    };

    /** {@code { ifNo, ifYes, assignIfNo, assignIfYes }}, 1-based, 0 meaning "nothing". */
    private static final int[][] TRANSITIONS = {
            {2, 0, 0, 1},    //  1  Q1
            {3, 0, 0, 3},    //  2  Q2
            {5, 4, 0, 0},    //  3  Q3
            {0, 7, 3, 0},    //  4  Q4
            {6, 0, 0, 1},    //  5  Q5
            {7, 0, 0, 3},    //  6  Q6
            {16, 8, 0, 0},   //  7  Q7
            {10, 9, 0, 0},   //  8  Q8
            {34, 0, 0, 3},   //  9  Q9
            {11, 0, 0, 3},   // 10  Q10
            {12, 33, 0, 0},  // 11  Q11
            {22, 13, 0, 0},  // 12  Q12
            {0, 14, 3, 0},   // 13  Q13
            {22, 15, 0, 0},  // 14  Q14
            {33, 37, 0, 0},  // 15  Q15
            {17, 0, 0, 1},   // 16  Q16
            {19, 38, 0, 0},  // 17  Q17
            {0, 0, 1, 2},    // 18  Q18
            {23, 20, 0, 0},  // 19  Q19
            {22, 21, 0, 0},  // 20  Q20
            {18, 0, 0, 3},   // 21  Q21
            {33, 0, 0, 2},   // 22  Q22
            {24, 27, 0, 0},  // 23  Q23
            {25, 18, 0, 0},  // 24  Q24
            {26, 0, 0, 2},   // 25  Q25
            {22, 0, 0, 2},   // 26  Q26
            {0, 28, 3, 0},   // 27  Q27
            {30, 29, 0, 0},  // 28  Q28
            {33, 36, 0, 0},  // 29  Q29
            {39, 31, 0, 0},  // 30  Q30
            {32, 39, 0, 0},  // 31  Q31
            {22, 0, 0, 2},   // 32  Q32
            {0, 0, 3, 1},    // 33  Q33, the last of the original questions
            {35, 20, 0, 0},  // 34  from Q9
            {23, 10, 0, 0},  // 35
            {19, 30, 0, 0},  // 36  from Q29
            {16, 22, 0, 0},  // 37  from Q15
            {19, 18, 0, 0},  // 38  from Q17
            {19, 18, 0, 0},  // 39  from Q30 or Q31
    };

    public CramerTree() {
        super(buildRules(), TRANSITIONS, CATEGORIES);
    }

    @Override
    public String getTitle() {
        return "Cramer rules";
    }

    private static DecisionRule[] buildRules() {
        return new DecisionRule[]{
                StructureListRule.normalBodyConstituent(),
                new ToxicFunctionalGroupsRule(),
                new ElementsOtherThanRule(),
                new SaltOrSulphonateRule(),
                new SimplyBranchedOrCarbohydrateRule(),
                new BenzeneDerivativeRule(),
                new FlagRule("7", "Heterocyclic", MolFlags::isHeterocyclic),
                new LactoneRule(),
                new UnimplementedRule("9", "Fused or unsaturated lactone"),
                new FlagRule("10", "Three-membered heterocycle", MolFlags::isHeterocyclic3),
                new HeterocyclicComplexSubstituentsRule(),
                new FlagRule("12", "Heteroaromatic", MolFlags::isHeteroaromatic),
                RingRules.ringsHaveSubstituents("13", "Ring with substituents"),
                RingRules.manyAromaticRings("14"),
                HydrolysisRules.toMononuclearResidues("15"),
                new CommonTerpeneRule("16"),
                HydrolysisRules.toCommonTerpene("17"),
                new KetoneAlcoholRule(),
                new FlagRule("19", "Open chain", MolFlags::isOpenChain),
                new AliphaticWithFunctionalGroupsRule(),
                new ManyFunctionalGroupsRule(),
                StructureListRule.commonComponentOfFood(),
                new FlagRule("23", "Aromatic", MolFlags::isAromatic),
                new MonocarbocyclicRule(),
                SmallRingRules.cyclopropaneOrCyclobutane("25"),
                SmallRingRules.monocycloalkanoneOrBicyclic("26"),
                RingRules.ringsHaveSubstituents("27", "Rings with substituents"),
                RingRules.manyAromaticRings("28"),
                HydrolysisRules.toMononuclearResidues("29"),
                new UnimplementedRule("30", "Ring with complex substituents"),
                new UnimplementedRule("31", "Acyclic acetal or ester of question 30"),
                new UnimplementedRule("32", "Only the functional groups of question 30"),
                new SulphonateGroupsRule(),
                new FlagRule("34", "Open chain", MolFlags::isOpenChain),
                new FlagRule("35", "Heterocyclic", MolFlags::isHeterocyclic),
                new FlagRule("36", "Aromatic", MolFlags::isAromatic),
                new FlagRule("37", "Heterocyclic", MolFlags::isHeterocyclic),
                new CommonTerpeneRule("38"),
                new FlagRule("39", "Aromatic", MolFlags::isAromatic),
        };
    }
}

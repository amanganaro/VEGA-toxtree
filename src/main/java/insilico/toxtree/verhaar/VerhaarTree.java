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
package insilico.toxtree.verhaar;

import insilico.toxtree.core.DecisionRule;
import insilico.toxtree.core.DecisionTree;
import insilico.toxtree.core.DescriptorRangeRule;
import insilico.toxtree.core.ElementsRule;
import insilico.toxtree.core.SmartsRule;
import insilico.toxtree.core.UnimplementedRule;
import insilico.toxtree.verhaar.rules.Rule11;
import insilico.toxtree.verhaar.rules.RuleIonicGroups;
import insilico.toxtree.verhaar.rules.Rule152;
import insilico.toxtree.verhaar.rules.Rule171;
import insilico.toxtree.verhaar.rules.Rule24;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;

import java.util.Set;

/**
 * The Verhaar scheme for predicting the mode of action in aquatic toxicity, in its original form.
 * <p>
 * This is the scheme VEGA uses. Toxtree also ships a modified order, from Enoch et al. 2008, which asks the
 * class 4 question first; it uses the very same rules, so it would only need a second transition table.
 * <p>
 * Rules are listed in the order of the transition table, and the comment on each row carries the question number
 * used in the decision trace.
 */
public final class VerhaarTree extends DecisionTree {

    private static final String[] CATEGORIES = {
            "Class 1 (narcosis or baseline toxicity)",
            "Class 2 (less inert compounds)",
            "Class 3 (unspecific reactivity)",
            "Class 4 (compounds and groups of compounds acting by a specific mechanism)",
            "Class 5 (Not possible to classify according to these rules)"
    };

    /** {@code { ifNo, ifYes, assignIfNo, assignIfYes }}, 1-based, 0 meaning "nothing". */
    private static final int[][] TRANSITIONS = {
            {34, 2, 0, 0},   //  1  Q0.1
            {34, 3, 0, 0},   //  2  Q0.2
            {34, 4, 0, 0},   //  3  Q0.3
            {34, 5, 0, 0},   //  4  Q1.1
            {34, 6, 0, 0},   //  5  Q1.2
            {7, 0, 0, 1},    //  6  Q1.3
            {12, 8, 0, 0},   //  7  Q1.4
            {9, 0, 0, 1},    //  8  Q1.4.1
            {10, 0, 0, 1},   //  9  Q1.4.2
            {11, 0, 0, 1},   // 10  Q1.4.3
            {21, 0, 0, 1},   // 11  Q1.4.4
            {17, 13, 0, 0},  // 12  Q1.5
            {14, 0, 0, 1},   // 13  Q1.5.1
            {15, 0, 0, 1},   // 14  Q1.5.2
            {16, 0, 0, 1},   // 15  Q1.5.3
            {21, 0, 0, 1},   // 16  Q1.5.4
            {19, 18, 0, 0},  // 17  Q1.6
            {21, 0, 0, 1},   // 18  Q1.6.1
            {21, 20, 0, 0},  // 19  Q1.7
            {21, 0, 0, 1},   // 20  Q1.7.1
            {22, 0, 0, 2},   // 21  Q2.1
            {23, 0, 0, 2},   // 22  Q2.2
            {24, 0, 0, 2},   // 23  Q2.3
            {25, 0, 0, 2},   // 24  Q2.4
            {26, 0, 0, 2},   // 25  Q2.5
            {27, 0, 0, 3},   // 26  Q3.1
            {28, 0, 0, 3},   // 27  Q3.2
            {29, 0, 0, 3},   // 28  Q3.3
            {30, 0, 0, 3},   // 29  Q3.4
            {31, 0, 0, 3},   // 30  Q3.5
            {32, 0, 0, 3},   // 31  Q3.6
            {33, 0, 0, 3},   // 32  Q3.7
            {34, 0, 0, 3},   // 33  Q3.8
            {0, 0, 5, 4},    // 34  Q4
    };

    public VerhaarTree() {
        super(buildRules(), TRANSITIONS, CATEGORIES);
    }

    @Override
    public String getTitle() {
        return "Verhaar scheme";
    }

    private static DecisionRule[] buildRules() {
        DecisionRule rule151 = smarts("1.5.1",
                "Be linear ethers or monocyclic mono ethers, but not epoxides or peroxides",
                "[$([#6]);!$(C=O)]-[$([O;!R]),$([O;R;!r3])]-[$([#6]);!$(C=O)]");
        DecisionRule rule152 = new Rule152();
        DecisionRule rule153 = new SmartsRule("1.5.3",
                "Be alcohols with aromatic moieties, but NOT phenols or benzylic alcohols", SmartsRule.Combine.ALL) {
            {
                reject("c[OH1]", "cC[OH1]");
                require("(C[OH1].a)");
            }
        };
        DecisionRule rule154 = smarts("1.5.4",
                "Be ketones, but not alpha-, beta unsaturated ketones (e.g. 1-butenone or acetophenone)",
                "O=C([C;!$(C=C)])([C;!$(C=C)])");

        return new DecisionRule[]{
                // Preconditions: outside these, the scheme does not apply and the molecule goes to Q4.
                elements("0.1", "Consists only of C,H,N,O,S,halogens (exluding I)",
                        ElementsRule.Mode.ONLY, Set.of("Cl", "Br", "F"), "C", "H", "N", "O", "S", "X"),
                new DescriptorRangeRule("0.2", "Have a logKow between 0 and 6",
                        new XLogPDescriptor(), -0.75, 6.75),
                new DescriptorRangeRule("0.3", "Have a molecular mass (MW) not more than 600 Daltons",
                        new WeightDescriptor(), 0, 600),

                // Class 1: narcosis or baseline toxicity.
                new Rule11(),
                new RuleIonicGroups(),
                elements("1.3", "Contain only C&H", ElementsRule.Mode.EXACTLY, "C", "H"),
                elements("1.4", "Contain only C,H and halogen", ElementsRule.Mode.EXACTLY,
                        "C", "H", "F", "Br", "I", "Cl"),
                new UnimplementedRule("1.4.1", "Be acyclic and not have halogens at beta from unsaturation"),
                new UnimplementedRule("1.4.2", "Be monocyclic with allowed ring substituents"),
                new UnimplementedRule("1.4.3", "Be monocyclic with allowed substituents"),
                new UnimplementedRule("1.4.4", "Be polycyclic with allowed substituents"),
                elements("1.5", "Contain C,H & O", ElementsRule.Mode.EXACTLY, "C", "H", "O"),
                rule151,
                rule152,
                rule153,
                rule154,
                elements("1.6", "Contain only C,H,N", ElementsRule.Mode.EXACTLY, "C", "H", "N"),
                new UnimplementedRule("1.6.1", "Be aliphatic secondary or tertiary amines"),
                elements("1.7", "Contain only C,H,O and halogen", ElementsRule.Mode.EXACTLY,
                        "C", "H", "O", "Cl", "Br", "I", "F"),
                new Rule171(rule151, rule152, rule153, rule154),

                // Class 2: less inert compounds.
                new UnimplementedRule("2.1", "Be phenols with allowed substituents"),
                new UnimplementedRule("2.2", "Be anilines with allowed substituents"),
                new UnimplementedRule("2.3", "Be mononitroaromatics with allowed substituents"),
                new Rule24(),
                new UnimplementedRule("2.5", "Be pyridines with allowed substituents"),

                // Class 3: unspecific reactivity.
                smartsAny("3.1", "Possess alylic/propargylic activation",
                        "C=;!:CC([Cl,Br,I])",
                        "C=;!:CC(C#N)",
                        "C=;!:CC(O[H])",
                        "C=;!:CC(C(C)=O)",
                        "C=;!:CC([CH1]=O)",
                        "C#C[C;v4;H1][Cl,Br,I]",
                        "C#C[C;v4;H1]C#N",
                        "C#C[C;v4;H1]O[H]",
                        "C#C[C;v4;H1]C(C)=O",
                        "CC#C[C;v4;H1](C)[CH1]=O"),
                smartsAny("3.2", "Possess benzylic activation",
                        "a[$([CX4H1]),$([CX4][C,N,O,S,Cl,Br,I])][Cl,Br,I]",
                        "a[$([CX4H1]),$([CX4][C,N,O,S,Cl,Br])]C#N",
                        "a[$([CX4H1]),$([CX4][C,N,O,S,Cl,Br])][OH1]",
                        "a[$([CX4H1]),$([CX4][C,N,O,S,Cl,Br])]C(C)=O",
                        "a[$([CX4H1]),$([CX4][C,N,O,S,Cl,Br])][CH1]=O"),
                smartsAny("3.3", "Be other compound with a (good) leaving group at an alpha position",
                        "A=,#AC([C;!R])[Cl,Br,I]",
                        "A=,#AC([C;!R])[C;!R]#N",
                        "A=,#AC([C;!R])[OH1;!R]",
                        "A=,#AC([C;!R])[C;!R]C(O)=O",
                        "A=,#AC([C;!R])[C;!R][CH1]=O"),
                smartsAny("3.4", "Possess a three-membered heterocyclic ring",
                        "C1OC1", "C1NC1", "[r3;!$([#6])]"),
                smartsAny("3.5", "Possess activated C-C double/triple bonds",
                        "[C!R]=,#[C!R]C(=O)O",
                        "[C!R]=,#[C!R]C(=O)",
                        "[C!R]=,#[C!R]C(=O)[N!R]([C!R])[C!R]",
                        "[C!R]=,#[C!R]C#N",
                        "[C!R]=,#[C!R](N(=O)=O)",
                        "[C!R]=,#[C!R]([N+](=O)[O-])",
                        "[C!R]=,#[C!R]S(=O)(=O)",
                        "O=C1C=CC(=O)C=C1"),
                smartsAny("3.6", "Be hydrazines or other compounds with a single, double or triple N-N linkage",
                        "[NH1;X3][NH1;X3]", "[NH0;X2]=[NH0;X2]", "[NX3]N#N"),
                smartsAny("3.7", "Be activated nitriles",
                        "[OH1][CX4]C#N", "[CX2]#CC#N"),
                smartsAny("3.8", "Contain one of the following structural entities",
                        "C(=O)OC(=O)",
                        "[C;R][O;R][C;R](=O)[O;R]",
                        "C(=O)Cl",
                        "NC(=O)Cl",
                        "C=C=O",
                        "[CX3H1](=O)",
                        "[CX3H2](=O)",
                        "N=C=O",
                        "SC#N",
                        "N=C=S",
                        "[S;X2][S;X2]",
                        "S(=O)(=O)[OX2H0]",
                        "[O;H1]S(=O)(=O)[OX2H]",
                        "[C;R][S;R]([O;R])(=O)(=O)",
                        "[O;R][S;R]([O;R])(=O)(=O)",
                        "[O;!R]C([Cl,Br,F,I])[!$([Cl,Br,F,I])]",
                        "[O;!R]CC([Cl,Br,F,I])[!$([Cl,Br,F,I])]",
                        "[N;!R][C;!R][C;!R]([Cl,F,Br,I])",
                        "[S;!R][C;!R][C;!R]([Cl,F,Br,I])"),

                // Class 4: compounds acting by a specific mechanism.
                smartsAny("4", "Compounds acting by a specific mechanism",
                        "C[Sn]",
                        "c1ccc(cc1)C(c2ccccc2)=,-C(Cl)Cl",
                        "[R](Cl)[R][R](Cl)[R][R](Cl)",
                        "Claa-aaCl",
                        "O=C(O)C1C(C=C)C1(C)(C)",
                        "OC(=O)[NX3]",
                        "[NX3]C(=O)S",
                        "[NX3]C(=S)S",
                        "O=P(Oc1ccccc1)(Oc2ccccc2)Oc3ccccc3",
                        "P(O)(O)(O)(=S)",
                        "P(O)(O)(S)(=S)",
                        "P(O)(C)(S)(=S)"),
        };
    }

    private static DecisionRule elements(String id, String title, ElementsRule.Mode mode, String... elements) {
        return new ElementsRule(id, title, mode, elements) {
        };
    }

    private static DecisionRule elements(String id, String title, ElementsRule.Mode mode, Set<String> halogens,
                                         String... elements) {
        return new ElementsRule(id, title, mode, halogens, elements) {
        };
    }

    private static DecisionRule smarts(String id, String title, String... smarts) {
        return new SmartsRule(id, title, SmartsRule.Combine.ALL) {
            {
                require(smarts);
            }
        };
    }

    private static DecisionRule smartsAny(String id, String title, String... smarts) {
        return new SmartsRule(id, title, SmartsRule.Combine.ANY) {
            {
                require(smarts);
            }
        };
    }
}

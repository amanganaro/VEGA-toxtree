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
package insilico.toxtree.chem;

import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;

import java.util.List;

/**
 * The structural flags the rules ask about, as computed by {@link Preprocessing}.
 * <p>
 * The semantics are those of Toxtree's {@code MolAnalyser}, which are not always the obvious ones. In particular
 * {@link #isAliphatic()} means "has no ring at all", not "is not aromatic", and it is switched off as soon as a
 * triple bond is present. Rules were written against those meanings, so they are kept verbatim.
 */
public final class MolFlags {

    private final boolean openChain;
    private final boolean aliphatic;
    private final boolean acetylenic;
    private final boolean aromatic;
    private final boolean heterocyclic;
    private final boolean heterocyclic3;
    private final boolean heteroaromatic;
    private final int aromaticRingCount;
    private final IRingSet ringSet;
    private final List<IRing> heterocyclicRings;
    private final List<IRing> aromaticRings;

    MolFlags(boolean openChain, boolean aliphatic, boolean acetylenic, boolean aromatic, boolean heterocyclic,
             boolean heterocyclic3, boolean heteroaromatic, int aromaticRingCount, IRingSet ringSet,
             List<IRing> heterocyclicRings, List<IRing> aromaticRings) {
        this.openChain = openChain;
        this.aliphatic = aliphatic;
        this.acetylenic = acetylenic;
        this.aromatic = aromatic;
        this.heterocyclic = heterocyclic;
        this.heterocyclic3 = heterocyclic3;
        this.heteroaromatic = heteroaromatic;
        this.aromaticRingCount = aromaticRingCount;
        this.ringSet = ringSet;
        this.heterocyclicRings = heterocyclicRings;
        this.aromaticRings = aromaticRings;
    }

    /** True when the molecule has no ring. */
    public boolean isOpenChain() {
        return openChain;
    }

    /**
     * True when the molecule has no ring <em>and</em> no triple bond. This is Toxtree's meaning, not the usual
     * chemical one: a cyclohexane is not "aliphatic" here.
     */
    public boolean isAliphatic() {
        return aliphatic;
    }

    /** True when the molecule has no ring and at least one triple bond. */
    public boolean isAcetylenic() {
        return acetylenic;
    }

    /** True when at least one ring has all its atoms flagged aromatic. */
    public boolean isAromatic() {
        return aromatic;
    }

    /** True when at least one ring contains an atom other than carbon. */
    public boolean isHeterocyclic() {
        return heterocyclic;
    }

    /** True when at least one three-membered ring contains an atom other than carbon. */
    public boolean isHeterocyclic3() {
        return heterocyclic3;
    }

    /** True when at least one ring is both aromatic and heterocyclic. */
    public boolean isHeteroaromatic() {
        return heteroaromatic;
    }

    public int getAromaticRingCount() {
        return aromaticRingCount;
    }

    /** The perceived rings; never null, possibly empty. */
    public IRingSet getRingSet() {
        return ringSet;
    }

    public int getRingCount() {
        return ringSet.getAtomContainerCount();
    }

    public List<IRing> getHeterocyclicRings() {
        return heterocyclicRings;
    }

    public List<IRing> getAromaticRings() {
        return aromaticRings;
    }
}

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
import insilico.toxtree.chem.Smarts;

/**
 * The functional groups the Cramer questions allocate atoms to.
 * <p>
 * Each is written as Toxtree builds it, including which of its atoms are only the anchor that says where the
 * group attaches. An anchor is left unallocated so that it stays available to another group, which is what lets
 * a single carbon belong to the chain and to nothing else.
 */
public final class CramerGroups {

    public static final String ALCOHOL = "alcohol";
    public static final String ALDEHYDE = "aldehyde";
    public static final String ACETAL = "acetal";
    public static final String KETONE = "ketone";
    public static final String ACID = "acid";
    public static final String ESTER = "ester";
    public static final String MERCAPTAN = "mercaptan";
    public static final String SULPHIDE = "sulphide";
    public static final String METHYL_ETHER = "methylether";
    public static final String HYDROXY_RING = "hydroxyRing";

    private CramerGroups() {
    }

    public static FunctionalGroup alcohol() {
        return FunctionalGroup.of(ALCOHOL, Smarts.CHN + "-[#8]-[#1]", 0);
    }

    public static FunctionalGroup aldehyde() {
        return FunctionalGroup.of(ALDEHYDE, "[#6,#1]-[#6]([#1])=[#8]", 0);
    }

    /** Two ethers on the same carbon. */
    public static FunctionalGroup acetal() {
        return FunctionalGroup.of(ACETAL, "[#6]-[#8]-[#6]-[#8]-[#6]");
    }

    public static FunctionalGroup ketone() {
        return FunctionalGroup.of(KETONE, "[!#1;!#8;!#16;!#7]-[#6](=[#8])-[!#1;!#8;!#16;!#7]", 0, 3);
    }

    public static FunctionalGroup carboxylicAcid() {
        return FunctionalGroup.of(ACID, "[*]-[#6](-[#8]-[#1])=[#8]", 0);
    }

    public static FunctionalGroup ester() {
        return FunctionalGroup.of(ESTER, "[*]-[#6](-[#8]-[!#1])=[#8]", 0, 3);
    }

    public static FunctionalGroup mercaptan() {
        return FunctionalGroup.of(MERCAPTAN, "[#6]-[#16]-[#1]");
    }

    public static FunctionalGroup sulphide() {
        return FunctionalGroup.of(SULPHIDE, "[#6]-[#16]-[#6]");
    }

    public static FunctionalGroup methylEther() {
        return FunctionalGroup.of(METHYL_ETHER, "[#8](-[#6]([#1])([#1])[#1])-[#6]");
    }

    /** A hydroxyl on a ring carbon; the carbon is only the anchor. */
    public static FunctionalGroup hydroxyOnRing() {
        return FunctionalGroup.of(HYDROXY_RING, "[#6;R]-[#8]-[#1]", 0);
    }
}

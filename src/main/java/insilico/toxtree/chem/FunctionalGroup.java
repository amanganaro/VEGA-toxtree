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

import java.util.Arrays;

/**
 * A functional group the rules can allocate atoms to.
 * <p>
 * Toxtree builds these as hand-made query containers; here they are SMARTS. What matters besides the pattern is
 * which of its atoms count as <em>belonging</em> to the group. A group is usually written with an anchor atom
 * that only says where it attaches, such as the {@code [*]} in {@code [*]C(=O)O[!#1]}: the original leaves such
 * atoms unmarked so that they stay available for another group, and so does this.
 *
 * @param id       the name the rules refer to the group by
 * @param smarts   the pattern
 * @param unmarked positions in the pattern, counted from zero, that do not belong to the group
 */
public record FunctionalGroup(String id, String smarts, int... unmarked) {

    public FunctionalGroup {
        Smarts.of(smarts);
        unmarked = unmarked.clone();
        Arrays.sort(unmarked);
    }

    /** True when the atom at this position in the pattern is only context and must not be marked. */
    public boolean isUnmarked(int position) {
        return Arrays.binarySearch(unmarked, position) >= 0;
    }

    @Override
    public String toString() {
        return id + " " + smarts;
    }
}

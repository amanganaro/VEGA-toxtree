/*
 * Copyright (C) 2026 ProtoQSAR / Alberto Manganaro
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 */
package insilico.toxtree.core;

/**
 * The outcome of running a decision tree on one molecule: the assigned category and the path taken to reach it.
 * <p>
 * The trace is the sequence of answers, as in the original Toxtree, for example
 * {@code "1N,2N,3N,5N,6N,7N,16N,17N,19Y,20Y,21N,18N"}. It is what makes a port testable: two implementations can
 * be compared rule by rule rather than only on the final answer.
 *
 * @param categoryIndex 1-based index of the assigned category, or 0 when none was reached
 * @param categoryName  the category name, or {@link #UNCLASSIFIED} when none was reached
 * @param trace         the decision path
 * @param error         the failure that interrupted the classification, or null
 */
public record Classification(int categoryIndex, String categoryName, String trace, Throwable error) {

    public static final String UNCLASSIFIED = "N.A.";

    /** True when a category was assigned. */
    public boolean isClassified() {
        return categoryIndex > 0;
    }

    /** True when the classification was interrupted by a failure. */
    public boolean isError() {
        return error != null;
    }

    static Classification of(int index, String name, String trace) {
        return new Classification(index, name, trace, null);
    }

    static Classification unclassified(String trace) {
        return new Classification(0, UNCLASSIFIED, trace, null);
    }

    static Classification failed(String trace, Throwable error) {
        return new Classification(0, UNCLASSIFIED, trace, error);
    }

    @Override
    public String toString() {
        return categoryName + " [" + trace + "]";
    }
}

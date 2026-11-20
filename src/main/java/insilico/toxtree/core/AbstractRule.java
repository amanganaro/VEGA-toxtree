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

/** Holds the identifier and the title, so that rules only have to supply their logic. */
public abstract class AbstractRule implements DecisionRule {

    private final String id;
    private final String title;

    protected AbstractRule(String id, String title) {
        this.id = id;
        this.title = title;
    }

    @Override
    public final String getId() {
        return id;
    }

    @Override
    public final String getTitle() {
        return title;
    }

    @Override
    public String toString() {
        return id + " " + title;
    }
}

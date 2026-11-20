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

import insilico.toxtree.core.AbstractRule;
import insilico.toxtree.core.RuleContext;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.HashSet;
import java.util.Set;

/**
 * Cramer questions 1 and 22, which are not about structure at all but about membership of a list: is the
 * substance a normal constituent of the body, or a common component of food?
 * <p>
 * The lists are expert judgements, 402 and 103 substances, and they are matched by exact standard InChI. The
 * aromatic flags are cleared from a copy of the molecule before the InChI is generated, as the original does,
 * because the generator is sensitive to them.
 * <p>
 * The lists come from Toxtree unchanged. The answers they give are therefore only as complete as the lists are,
 * which the original is candid about: a substance absent from them is not thereby absent from food.
 */
public final class StructureListRule extends AbstractRule {

    private final String resource;
    private volatile Set<String> inchis;

    private StructureListRule(String id, String title, String resource) {
        super(id, title);
        this.resource = resource;
    }

    /** Question 1: a normal constituent of the body, or an optical isomer of one. */
    public static StructureListRule normalBodyConstituent() {
        return new StructureListRule("1", "Normal constituent of the body", "/insilico/toxtree/cramer/bodymol.inchi");
    }

    /** Question 22: a common component of food, or structurally closely related to one. */
    public static StructureListRule commonComponentOfFood() {
        return new StructureListRule("22", "Common component of food", "/insilico/toxtree/cramer/foodmol.inchi");
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        Set<String> list = list();
        try {
            IAtomContainer copy = mol.clone();
            for (IBond bond : copy.bonds()) {
                bond.setIsAromatic(false);
            }
            String inchi = InChIGeneratorFactory.getInstance().getInChIGenerator(copy).getInchi();
            return inchi != null && list.contains(inchi);
        } catch (Exception x) {
            return false;
        }
    }

    private Set<String> list() {
        Set<String> local = inchis;
        if (local == null) {
            synchronized (this) {
                local = inchis;
                if (local == null) {
                    local = read();
                    inchis = local;
                }
            }
        }
        return local;
    }

    private Set<String> read() {
        Set<String> values = new HashSet<>();
        try (InputStream in = StructureListRule.class.getResourceAsStream(resource)) {
            if (in == null) {
                throw new IllegalStateException("Missing structure list: " + resource);
            }
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    line = line.trim();
                    if (!line.isEmpty()) {
                        values.add(line);
                    }
                }
            }
        } catch (IOException x) {
            throw new UncheckedIOException(x);
        }
        return values;
    }
}

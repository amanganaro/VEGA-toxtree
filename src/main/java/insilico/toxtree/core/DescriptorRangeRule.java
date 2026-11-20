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
package insilico.toxtree.core;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerResult;

/**
 * A rule that answers yes when a CDK molecular descriptor falls inside a closed range.
 * <p>
 * The descriptor is used with its default parameters. That is not an oversight: Toxtree's own logP rule declares
 * parameters and then never passes them to the descriptor, so the reference values are the ones produced by the
 * defaults. Changing them here would silently move the boundary for every molecule.
 */
public class DescriptorRangeRule extends AbstractRule {

    private final IMolecularDescriptor descriptor;
    private final double min;
    private final double max;

    public DescriptorRangeRule(String id, String title, IMolecularDescriptor descriptor, double min, double max) {
        super(id, title);
        this.descriptor = descriptor;
        this.min = min;
        this.max = max;
    }

    public double getMin() {
        return min;
    }

    public double getMax() {
        return max;
    }

    /** The descriptor value for a molecule, for tests and for reporting. */
    public double value(IAtomContainer mol) {
        DescriptorValue value;
        synchronized (descriptor) {
            value = descriptor.calculate(mol);
        }
        if (value.getException() != null) {
            throw new IllegalStateException("Descriptor " + getId() + " failed", value.getException());
        }
        IDescriptorResult result = value.getValue();
        if (result instanceof DoubleResult d) {
            return d.doubleValue();
        }
        if (result instanceof IntegerResult i) {
            return i.intValue();
        }
        throw new IllegalStateException("Unsupported descriptor result: " + result.getClass().getName());
    }

    @Override
    public boolean verify(IAtomContainer mol, RuleContext ctx) {
        double d = value(mol);
        return d >= min && d <= max;
    }
}

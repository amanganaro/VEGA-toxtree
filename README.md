# VEGA - ToxTree rules

The **Cramer** and **Verhaar** decision rules, ported to CDK 2.x.

They replace `ToxTreeWrapper` (72 MB) and `toxtree-cramer`, and with them ambit, CDK 1.5.13, beam 0.9.1 and the
IdeaConsult Nexus. The VEGA models `CramerToxTree` and `VerhaarToxtree` depend on this artifact and contain no logic
of their own: they only run the tree.

## Usage

```java
VerhaarTree tree = new VerhaarTree();
Classification result = tree.classify(molecule.GetStructure());

result.categoryName();   // "Class 1 (narcosis or baseline toxicity)"
result.trace();          // "0.1Y,0.2Y,0.3Y,1.1Y,1.2Y,1.3Y"
result.isClassified();
```

**The molecule must arrive already normalised.** The library does not read SMILES, does not add hydrogens and does
not perceive aromaticity: it uses what it receives. This is deliberate, and it is the difference from the original,
which did its own preprocessing. In the VEGA models the normalisation is the one of `InsilicoMolecule`.

The only dependency is CDK, in `provided` scope: its version is chosen by whoever uses the library, which in VEGA is
the core.

## Structure

| Package | Content |
|---|---|
| `insilico.toxtree.core` | shared scaffolding: `DecisionTree`, `DecisionRule`, `RuleContext`, `Classification` and the rule base classes |
| `insilico.toxtree.chem` | shared chemistry: structural flags, compiled SMARTS, element sets |
| `insilico.toxtree.verhaar` | the Verhaar scheme and its rules |
| `insilico.toxtree.cramer` | the Cramer scheme and its rules |

Each tree is declarative, as in the original: a list of rules and a transition table
`{ nextIfNo, nextIfYes, categoryIfNo, categoryIfYes }`. Rules are stateless; everything that depends on the molecule
lives in the `RuleContext`, so the library can be used from several threads, which the original did not allow.

## The decision trace

Besides the category, every classification returns the path followed through the tree, in the same format as
Toxtree:

```
1N,2N,3N,5N,6N,7N,16N,17N,19Y,20Y,21N,18N
```

This is what makes the port verifiable: two implementations are compared rule by rule, not only on the final
result. A rule not yet written shows up as `-` instead of `N`, so an incomplete tree is told apart at a glance from
one that answers wrongly.

## Status

| Tree | Questions implemented | Agreement with Toxtree |
|---|---|---|
| Verhaar | **34 of 34** | **100%**: 586 molecules out of 586, identical trace and category |
| Cramer | **39 of 39** | golden dataset: identical category on all 584; logP training set: 9927 out of 9937 (99.9%) |

Both trees are complete.

## Verification

`TraceComparisonTest` compares this port with the traces captured from the original Toxtree on the 603 molecules of
the VEGA golden dataset, and prints the agreement of every single rule plus the molecules on which the two
implementations diverge.

The reference files are in `src/test/resources/reference`:

| File | Molecules | Source |
|---|---|---|
| `cramer-core.tsv`, `verhaar-core.tsv` | 603 | VEGA golden dataset |
| `verhaar-2000.tsv` | 121 | Toxtree validation data, Verhaar et al. 2000 |
| `cramer-logp.tsv`, `verhaar-logp.tsv` | 9961 | logP training set, run through VEGA 1.2.4 |

They are regenerated with the tools in [`tools/`](tools/README.md), the only place that still needs the old wrapper
and the IdeaConsult Nexus.

### Group-by-group verification

The trace comparison tells whether a rule answers as it should, but it does not see a badly written functional
group until the molecule that exposes it comes along. That is why there is `tools/VerifyGroups`, which compares
**every group, atom by atom**, with what Toxtree builds: same molecule, same atoms matched.

Of 66 groups, 50 match exactly and 11 match on the Kekulé view, i.e. they are defined on bond order. The three
remaining warnings have no effect on any classification and are explained in `tools/README.md`.

## Quirks of the original, reproduced on purpose

The port reproduces the current behaviour, not the one you would expect. These points are worth writing down,
because they look like bugs and partly are:

- **Cramer questions 6 and 30 can never answer yes.** Both ask Toxtree to consider only aromatic rings, and the
  selection never finds one: the aromatic flag is stored as if it were a bit mask and then compared with the
  positions of the ring flags. Anethole and estragole, an alkoxy group para to a hydrocarbon chain, are exactly what
  question 6 describes, and Toxtree answers no to both; likewise benzamide and phenyl isocyanate at question 30. The
  intended logic is still in the code, behind a switch in each class.
- **Question 3 decides sulphur divalence in two different ways** depending on the molecule: it asks CDK for the
  valence, which is 2 for an aromatic sulphur, but falls back on a sum of bond orders, where aromatic sulphur counts
  3, as soon as it meets a sulphur that is not divalent. Which path runs is decided by the whole molecule, not by
  the atom.
- **Toxtree carbons are not just any carbons:** many groups require a carbon bearing at least one hydrogen, or none
  and only carbon neighbours. The carbonyl carbon of an amide does not qualify, so a urea nitrogen is not an
  aliphatic secondary amine.
- **Bond queries look only at the order and ignore aromaticity**, which SMARTS does not allow: that is why the groups
  defined on bond order run on a Kekulé view with the aromatic flags cleared. This is how pyridinium turns out to
  have a quaternary nitrogen and the sulphur of thiophene a sulphide.
- **The main structure is never required.** Toxtree answers "yes, it is a weakly acidic phenol" for thiophene and
  for thiirane. `RuleRingMainStrucSubstituents` assigns the main structure inside `initQuery()`, which the superclass
  constructor calls, and then declares the field as `mainStructure = null`: in Java, field initialisers run after the
  superclass constructor, so the assignment is wiped out before the rule is used. The query survives among the
  allowed groups, because it was registered while the field still held it. 18 of the 603 molecules of the golden
  dataset depend on this.
- **No ring, no answer.** Toxtree records the ring set only if there is at least one ring, and the class 2 rules give
  up when they receive none. An acyclic molecule therefore never satisfies them.
- **Overlapping matches of the same group are discarded, counting by bond.** On a hexahalogenated cyclohexane the
  matches of "halogen on a ring carbon" share ring bonds: only the first survives and the other halogens stay
  uncovered, which is why Toxtree rejects that molecule. On a trichloromethyl group, instead, the three matches share
  the carbon but no bond, so they all count.

## Known differences from Toxtree

- **Residues are implemented**, but not everything about them can be reproduced. When question 4 detaches a
  sulphonate, or question 9 opens a lactone, Toxtree fills the broken bonds by counting `valence − neighbours`, and
  the count is wrong: benzenesulphonate yields a benzene with **seven** hydrogens. Aromaticity perception fails on
  that benzene, and the following questions answer accordingly. The port reproduces the count and therefore almost
  always the verdict too, but on some fused systems the two aromaticity detectors disagree on an impossible
  molecule. This causes part of the 10 differences left out of 9937.
- **An impossible carbon.** `C(C)(C)(C)(C)C`, with five carbons attached, makes the atom typing of the old CDK fail
  and Toxtree miscounts its hydrogens; CDK 2.3 types it without complaint. It is the only difference not caused by
  residues.
- **Hydrogens.** Toxtree makes them all explicit; VEGA keeps them implicit. Several SMARTS depend on explicit ones
  (`C=;!:CC(O[H])` names a hydrogen, and `[OH1]C[!$(C=C);!$(C#C)]` finds methanol only because the third atom can be
  a hydrogen), so matching runs on a copy with explicit hydrogens, built once per molecule.
- **Aromaticity.** Toxtree uses ambit's old Hückel detector; here the aromaticity is the one VEGA decides. This is
  the remaining divergence, and the trace test measures it.

## License

GPL v2 or later. The code derives from [Toxtree](https://toxtree.sourceforge.net/) by Ideaconsult Ltd., released
under that license; this is why the library lives in a repository separate from the VEGA core, which is LGPL v3.

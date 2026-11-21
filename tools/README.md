# Verification tools against the original ToxTree

Two programs used to check the port, not to run it. They are not part of the library and are not compiled with it:
they need the old `ToxTreeWrapper`, which is exactly what the port removes.

To use them you need a separate Maven project with these dependencies:

```xml
<dependency><groupId>InSilico</groupId><artifactId>ToxTreeWrapper</artifactId>
           <version>1.0</version><classifier>shaded</classifier></dependency>
<dependency><groupId>toxtree</groupId><artifactId>toxtree-cramer</artifactId><version>3.1.0</version></dependency>
<dependency><groupId>uk.ac.ebi.beam</groupId><artifactId>beam-core</artifactId><version>1.3.3</version></dependency>
<dependency><groupId>InSilico</groupId><artifactId>Core</artifactId><version>1.3.23</version></dependency>
<dependency><groupId>InSilico</groupId><artifactId>ToxTreeRules</artifactId><version>0.1-SNAPSHOT</version></dependency>
```

plus the repository `https://nexus.ideaconsult.net/content/repositories/releases` for toxtree-cramer.

## ThreeWay

Generates the reference files in `src/test/resources/reference`. For each SMILES it normalises the molecule with
the VEGA core, takes its canonical SMILES, parses it again with CDK 1.5.13 — exactly what the old models did — and
saves the category and the decision trace.

```
java -cp ... ThreeWay cramer smiles.txt cramer-reference.tsv
java -cp ... ThreeWay verhaar smiles.txt verhaar-reference.tsv
```

It also tells which side a difference is on: if today's ToxTree agrees with a published VEGA output, the difference
belongs to the port and not to the core.

## VerifyGroups

Compares **every functional group**, atom by atom, with what ToxTree builds. Both sides receive the same canonical
SMILES, so the heavy atoms have the same positions, and for each molecule the two sets of matched atoms are
compared.

```
java -cp ... VerifyGroups smiles.txt [count]
```

For each group it prints on how many molecules the set matches on the aromatic view, on the Kekulé view, on
neither, and whether at least the presence of the group matches. It is needed because a test on the final
classifications does not see a wrong group until the molecule that exposes it comes along: this one sees it at once.

Two things to know when reading the output:

- **A single-bond group cannot be compared by atom sets.** The CDK method that returns the mappings does not produce
  them for such small queries, so the "neither" column fills up while the presence column stays correct. This is the
  case of `halogen`.
- **`ionicGroup` cannot be compared at all:** it relies on virtual bonds that ToxTree adds between charged
  fragments, and that do not exist in SMARTS.

## Last run

Of 66 groups, 50 match exactly on both views and 11 only on the Kekulé view — that is, groups defined on bond order,
which ToxTree reads regardless of aromaticity, and the rules that use them must go through `RuleContext.kekule()`.

Three warnings remain, none with any effect on a classification: `halogen` and `ionicGroup` for the reasons above,
and `esterOfalcoholSecondaryAttachedToTerminalVinyl`, where on 3 large molecules out of 9938 the port finds a match
that ToxTree does not.

# VEGA - ToxTree rules

Le regole decisionali di **Cramer** e **Verhaar**, portate su CDK 2.x.

Sostituiscono `ToxTreeWrapper` (72 MB) e `toxtree-cramer`, e con essi ambit, CDK 1.5.13, beam 0.9.1 e il Nexus di
IdeaConsult. I modelli VEGA `CramerToxTree` e `VerhaarToxtree` dipendono da questo artefatto e non contengono
nessuna logica: si limitano a far girare l'albero.

## Come si usa

```java
VerhaarTree tree = new VerhaarTree();
Classification result = tree.classify(molecule.GetStructure());

result.categoryName();   // "Class 1 (narcosis or baseline toxicity)"
result.trace();          // "0.1Y,0.2Y,0.3Y,1.1Y,1.2Y,1.3Y"
result.isClassified();
```

**La molecola deve arrivare già normalizzata.** La libreria non legge SMILES, non aggiunge idrogeni e non percepisce
l'aromaticità: usa quello che riceve. Questo è voluto, ed è la differenza con l'originale, che faceva il
preprocessing per conto suo. Nei modelli VEGA la normalizzazione è quella di `InsilicoMolecule`.

L'unica dipendenza è CDK, in scope `provided`: la versione la decide chi usa la libreria, che nel caso di VEGA è il
core.

## Struttura

| Package | Contenuto |
|---|---|
| `insilico.toxtree.core` | impalcatura comune: `DecisionTree`, `DecisionRule`, `RuleContext`, `Classification` e le classi base delle regole |
| `insilico.toxtree.chem` | chimica condivisa: flag strutturali, SMARTS compilati, insiemi di elementi |
| `insilico.toxtree.verhaar` | lo schema Verhaar e le sue regole |
| `insilico.toxtree.cramer` | lo schema Cramer e le sue regole |

Ogni albero è dichiarativo, come nell'originale: un elenco di regole e una tabella di transizioni
`{ seNo, seSi, assegnaSeNo, assegnaSeSi }`. Le regole sono senza stato; tutto ciò che dipende dalla molecola sta nel
`RuleContext`, così la libreria si può usare da più thread, cosa che l'originale non permetteva.

## La traccia di decisione

Oltre alla categoria, ogni classificazione restituisce il percorso seguito nell'albero, nello stesso formato di
Toxtree:

```
1N,2N,3N,5N,6N,7N,16N,17N,19Y,20Y,21N,18N
```

È il motivo per cui il porting è verificabile: due implementazioni si confrontano regola per regola, non solo sul
risultato finale. Una regola non ancora scritta compare come `-` invece di `N`, così un albero incompleto si
distingue a colpo d'occhio da uno che risponde male.

## Stato

| Albero | Regole fatte | Concordanza con Toxtree |
|---|---|---|
| Verhaar | **34 su 34** | **100%**: 586 molecole su 586, traccia e categoria identiche |
| Cramer | 0 su 34 | — |

Verhaar è finito. Manca Cramer, che riusa l'impalcatura e in più ha bisogno della simulazione di idrolisi e delle
liste di strutture in InChI.

## Verifica

`TraceComparisonTest` confronta questo porting con le tracce catturate dal Toxtree originale sulle 603 molecole del
dataset golden di VEGA, e stampa la concordanza di ogni singola regola più le molecole su cui le due implementazioni
divergono.

I file di riferimento stanno in `src/test/resources/reference` e si rigenerano con il profilo `toxtree-reference`,
che è l'unico punto in cui serve ancora il Nexus di IdeaConsult.

## Tre particolarità dell'originale, riprodotte apposta

Il porting riproduce il comportamento attuale, non quello che ci si aspetterebbe. Tre punti meritano di essere
scritti, perché sembrano errori e in parte lo sono:

- **La struttura principale non viene mai richiesta.** Toxtree risponde "sì, è un fenolo debolmente acido" per il
  tiofene e per la tiirana. `RuleRingMainStrucSubstituents` assegna la struttura principale dentro `initQuery()`,
  che il costruttore della superclasse chiama, e poi dichiara il campo come `mainStructure = null`: in Java gli
  inizializzatori di campo girano dopo il costruttore della superclasse, quindi l'assegnazione viene cancellata
  prima che la regola venga usata. La query sopravvive fra i gruppi ammessi, perché era stata registrata quando il
  campo la conteneva ancora. Da questo dipendono 18 delle 603 molecole del dataset golden.
- **Nessun anello, nessuna risposta.** Toxtree registra l'insieme degli anelli solo se ce n'è almeno uno, e le
  regole di classe 2 si arrendono quando non ne ricevono nessuno. Una molecola aciclica quindi non le soddisfa mai.
- **I match sovrapposti si scartano, contandoli per legame.** Su un cicloesano esalogenato i match di "alogeno su
  un carbonio dell'anello" condividono legami dell'anello: sopravvive solo il primo e gli altri alogeni restano
  scoperti, e per questo Toxtree rifiuta quella molecola. Su un gruppo triclorometile invece i tre match
  condividono il carbonio ma nessun legame, quindi valgono tutti.

## Due differenze note rispetto a Toxtree

- **Idrogeni.** Toxtree li rende tutti espliciti; VEGA li tiene impliciti. Parecchi SMARTS dipendono da quelli
  espliciti (`C=;!:CC(O[H])` nomina un idrogeno, e `[OH1]C[!$(C=C);!$(C#C)]` trova il metanolo solo perché il terzo
  atomo può essere un idrogeno), quindi il matching gira su una copia con gli idrogeni espliciti, costruita una volta
  per molecola.
- **Aromaticità.** Toxtree usa il vecchio rilevatore di Hückel di ambit; qui arriva quella decisa da VEGA. È la
  divergenza che resta, e si misura con il test delle tracce.

## Licenza

GPL v2 o successiva. Il codice deriva da [Toxtree](https://toxtree.sourceforge.net/) di Ideaconsult Ltd., rilasciato
con quella licenza; per questo la libreria sta in un repository separato dal core VEGA, che è LGPL v3.

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

| Albero | Domande fatte | Concordanza con Toxtree |
|---|---|---|
| Verhaar | **34 su 34** | **100%**: 586 molecole su 586, traccia e categoria identiche |
| Cramer | **39 su 39** | **categoria identica su tutte le 584 molecole**; traccia identica su 576 (98,6%) |

Entrambi gli alberi sono completi. Le 8 tracce che differiscono riguardano 4 molecole, per due motivi soli,
descritti sotto: i residui e un carbonio impossibile.

## Verifica

`TraceComparisonTest` confronta questo porting con le tracce catturate dal Toxtree originale sulle 603 molecole del
dataset golden di VEGA, e stampa la concordanza di ogni singola regola più le molecole su cui le due implementazioni
divergono.

I file di riferimento stanno in `src/test/resources/reference`:

| File | Molecole | Origine |
|---|---|---|
| `cramer-core.tsv`, `verhaar-core.tsv` | 603 | dataset golden di VEGA |
| `verhaar-2000.tsv` | 121 | dati di validazione di Toxtree, Verhaar et al. 2000 |
| `cramer-logp.tsv`, `verhaar-logp.tsv` | 9961 | training set di logP, fatto girare con VEGA 1.2.4 |

Si rigenerano con gli strumenti in [`tools/`](tools/README.md), l'unico punto in cui serve ancora il wrapper
vecchio e il Nexus di IdeaConsult.

### La verifica gruppo per gruppo

Il confronto delle tracce dice se una regola risponde come deve, ma non vede un gruppo funzionale scritto male
finché non capita la molecola che lo mette in crisi. Per questo c'è `tools/VerifyGroups`, che confronta **ogni
gruppo, atomo per atomo**, con quello che Toxtree costruisce: stessa molecola, stessi atomi toccati.

Su 66 gruppi, 50 coincidono esattamente e 11 coincidono sulla vista di Kekulé, cioè sono definiti sull'ordine dei
legami. Le tre segnalazioni che restano non hanno effetto su nessuna classificazione e sono spiegate in
`tools/README.md`.

## Particolarità dell'originale, riprodotte apposta

Il porting riproduce il comportamento attuale, non quello che ci si aspetterebbe. Questi punti meritano di essere
scritti, perché sembrano errori e in parte lo sono:

- **Le domande 6 e 30 di Cramer non possono mai rispondere sì.** Entrambe chiedono a Toxtree di considerare solo
  gli anelli aromatici, e la selezione non ne trova mai uno: il flag aromatico viene memorizzato come se fosse una
  maschera di bit e poi confrontato con le posizioni dei flag degli anelli. L'anetolo e l'estragolo, un alcossile
  in para a una catena idrocarburica, sono esattamente ciò che la domanda 6 descrive, e Toxtree risponde no a
  entrambi; così la benzammide e il fenil isocianato alla domanda 30. La logica prevista resta nel codice, dietro
  un interruttore in ciascuna classe.
- **La domanda 3 decide la divalenza dello zolfo in due modi diversi** a seconda della molecola: chiede la valenza
  a CDK, che per uno zolfo aromatico dà 2, ma ripiega su una somma degli ordini di legame, dove lo zolfo aromatico
  fa 3, appena incontra uno zolfo non divalente. Quale percorso giri lo decide la molecola intera, non l'atomo.
- **I carboni di Toxtree non sono carboni qualsiasi:** molti gruppi richiedono un carbonio che porti almeno un
  idrogeno, oppure nessuno e solo vicini carboni. Il carbonio carbonilico di un'ammide non qualifica, e quindi un
  azoto ureico non è un'ammina secondaria alifatica.
- **Le query sui legami guardano solo l'ordine e ignorano l'aromaticità**, cosa che SMARTS non permette: per questo
  i gruppi definiti sull'ordine girano su una vista di Kekulé con i flag aromatici azzerati. È così che il
  piridinio risulta avere un azoto quaternario e lo zolfo del tiofene un solfuro.
- **I match sovrapposti di uno stesso gruppo si scartano, contandoli per legame.**

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

- **Residui.** Quando la domanda 4 di Cramer risponde sì, Toxtree stacca il gruppo e pone tutte le domande
  successive al residuo, non alla molecola intera; la sua traccia lo segna con delle parentesi. Il porting
  continua a interrogare la molecola intera. Riguarda 7 tracce su 584 e **non cambia nessuna classificazione**,
  perché le categorie coincidono comunque.
- **Un carbonio impossibile.** `C(C)(C)(C)(C)C`, con cinque carboni attaccati, fa fallire l'atom typing del CDK
  vecchio e Toxtree ne conta male gli idrogeni; CDK 2.3 lo tipizza senza storcere il naso. È l'unica differenza
  non dovuta ai residui.
- **Idrogeni.** Toxtree li rende tutti espliciti; VEGA li tiene impliciti. Parecchi SMARTS dipendono da quelli
  espliciti (`C=;!:CC(O[H])` nomina un idrogeno, e `[OH1]C[!$(C=C);!$(C#C)]` trova il metanolo solo perché il terzo
  atomo può essere un idrogeno), quindi il matching gira su una copia con gli idrogeni espliciti, costruita una volta
  per molecola.
- **Aromaticità.** Toxtree usa il vecchio rilevatore di Hückel di ambit; qui arriva quella decisa da VEGA. È la
  divergenza che resta, e si misura con il test delle tracce.

## Licenza

GPL v2 o successiva. Il codice deriva da [Toxtree](https://toxtree.sourceforge.net/) di Ideaconsult Ltd., rilasciato
con quella licenza; per questo la libreria sta in un repository separato dal core VEGA, che è LGPL v3.

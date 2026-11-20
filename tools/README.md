# Strumenti di verifica contro il ToxTree originale

Due programmi che servono a controllare il porting, non a farlo girare. Non fanno parte della libreria e non
vengono compilati con essa: hanno bisogno del vecchio `ToxTreeWrapper`, cioè proprio di quello che il porting
elimina.

Per usarli serve un progetto Maven a parte con queste dipendenze:

```xml
<dependency><groupId>InSilico</groupId><artifactId>ToxTreeWrapper</artifactId>
           <version>1.0</version><classifier>shaded</classifier></dependency>
<dependency><groupId>toxtree</groupId><artifactId>toxtree-cramer</artifactId><version>3.1.0</version></dependency>
<dependency><groupId>uk.ac.ebi.beam</groupId><artifactId>beam-core</artifactId><version>1.3.3</version></dependency>
<dependency><groupId>InSilico</groupId><artifactId>Core</artifactId><version>1.3.23</version></dependency>
<dependency><groupId>InSilico</groupId><artifactId>ToxTreeRules</artifactId><version>0.1-SNAPSHOT</version></dependency>
```

più il repository `https://nexus.ideaconsult.net/content/repositories/releases` per toxtree-cramer.

## ThreeWay

Genera i file di riferimento che stanno in `src/test/resources/reference`. Per ogni SMILES normalizza la molecola
con il core di VEGA, ne prende lo SMILES canonico, lo riparsa con il CDK 1.5.13 — cioè fa esattamente quello che
facevano i modelli vecchi — e salva categoria e traccia di decisione.

```
java -cp ... ThreeWay cramer smiles.txt cramer-riferimento.tsv
java -cp ... ThreeWay verhaar smiles.txt verhaar-riferimento.tsv
```

Serve anche a stabilire da che parte sta una differenza: se ToxTree di oggi concorda con un'uscita di VEGA
pubblicata, allora la differenza è del porting e non del core.

## VerifyGroups

Confronta **ogni gruppo funzionale**, atomo per atomo, con quello che ToxTree costruisce. Entrambe le parti
ricevono lo stesso SMILES canonico, così gli atomi pesanti hanno le stesse posizioni, e per ogni molecola si
confrontano i due insiemi di atomi toccati.

```
java -cp ... VerifyGroups smiles.txt [quante]
```

Per ogni gruppo stampa su quante molecole l'insieme coincide con la vista aromatica, con quella di Kekulé, su
nessuna delle due, e se almeno la presenza del gruppo coincide. Serve perché un test sulle classificazioni finali
non vede un gruppo sbagliato finché non capita la molecola che lo mette in crisi: questo lo vede subito.

Due cose da sapere leggendo l'uscita:

- **Un gruppo con un solo legame non è confrontabile per insiemi di atomi.** Il metodo di CDK che restituisce le
  mappature non le produce per query così piccole, quindi la colonna "nessuna delle due" si riempie mentre la
  colonna della presenza resta a posto. È il caso di `halogen`.
- **`ionicGroup` non è confrontabile affatto:** si basa su legami virtuali che ToxTree aggiunge fra frammenti
  carichi, e che in SMARTS non esistono.

## Com'è andata l'ultima volta

Su 66 gruppi, 50 coincidono esattamente su entrambe le viste e 11 soltanto sulla vista di Kekulé — cioè sono
gruppi definiti sull'ordine dei legami, che ToxTree legge senza badare all'aromaticità, e le regole che li usano
devono passare da `RuleContext.kekule()`.

Restano tre segnalazioni, tutte senza effetto su nessuna classificazione: `halogen` e `ionicGroup` per i motivi
sopra, e `esterOfalcoholSecondaryAttachedToTerminalVinyl`, dove su 3 molecole grandi su 9938 il porting trova una
corrispondenza che ToxTree non trova.

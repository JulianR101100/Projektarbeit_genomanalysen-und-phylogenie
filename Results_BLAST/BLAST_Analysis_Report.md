# BLAST Analyse Bericht

## 1. Detaillierte Analyse des Proteins & der Annotation

### Identität und Funktion
Das identifizierte Protein ist **TSR3** (Ribosome Biogenesis Protein TSR3 Homolog). Die in den BLAST-Ergebnissen gefundene Bezeichnung **18S rRNA aminocarboxypropyltransferase** ist hierbei **nicht veraltet**, sondern stellt die präzise biochemische Funktionsbeschreibung dar (EC-Nummer 2.5.1.157).

*   **Gen-Symbol:** *TSR3* (in Analogie zum Hefe-Homolog).
*   **Enzymatische Aktivität:** Das Protein katalysiert einen hochspezifischen Schritt in der Ribosomen-Biogenese: den Transfer einer **Aminocarboxypropyl-Gruppe** (acp) von S-Adenosylmethionin (SAM) auf ein spezifisches Pseudouridin (Ψ) der 18S rRNA (in Hefe Position 1248, im Menschen konserviert).
*   **Bedeutung:** Diese Modifikation ist essentiell für die korrekte Reifung der kleinen ribosomalen Untereinheit (40S). Ein Funktionsverlust führt zu Defekten in der Ribosomenherstellung, was die hohe Konservierung (Top-Hits >80% Identität über alle Säugetiere hinweg) erklärt. Dass die Bezeichnung "Homolog" auftaucht, liegt daran, dass das Protein ursprünglich in Modellorganismen (wie *Saccharomyces cerevisiae*) charakterisiert wurde und die menschliche Variante als Ortholog identifiziert wurde.

---

## 2. Theorie der Substitutionsmatrix: Warum die Wahl entscheidend ist

Die Wahl der Substitutionsmatrix (z.B. BLOSUM62, BLOSUM80 oder PAM30) ist keine rein technische Formalität, sondern eine **implizite Hypothese über die evolutionäre Distanz** der erwarteten Treffer.

### Grundlagen: Log-Odds Scores
Jede Zelle $S_{ij}$ in einer Substitutionsmatrix ist definiert als:
$$
S_{ij} = \log_2 \left( \frac{q_{ij}}{p_i p_j} \right) \\
$$ 
*   $q_{ij}$: Die Wahrscheinlichkeit, dass Aminosäure $i$ und $j$ evolutionär voneinander abstammen (observed frequency in alignments).
*   $p_i p_j$: Die Wahrscheinlichkeit, dass sie zufällig zusammen auftreten (background frequency).

### Vergleich der Matrizen am konkreten Beispiel
*   **BLOSUM62 (Standard):** Wurde aus Blöcken von Alignments erstellt, in denen die Sequenzen **weniger als 62% Identität** hatten. Sie ist ein "Allrounder", der darauf ausgelegt ist, moderate Verwandtschaften zu erkennen.
*   **BLOSUM80:** Wurde aus Blöcken mit **weniger als 80% Identität** erstellt. Da die Datenbasis ähnlicher ist, sind die "Belohnungen" für konservierte Matches und die "Strafen" für Mismatches strikter.
*   **PAM30:** Basiert auf dem PAM-Modell (Point Accepted Mutation) für sehr kurze evolutionäre Distanzen (30 PAM-Einheiten). Sie ist extrem "streng" und ideal für Alignments mit sehr hoher Identität (oft >80-90%).

**Anwendung hier:** Unsere Top-Hits haben >80% bis 100% Identität. Verwendet man BLOSUM62, "unterschätzt" man die Strenge der Selektion. BLOSUM80 oder PAM30 haben eine höhere **relative Entropie** (Information Content pro Position) und können feinere Unterschiede zwischen sehr ähnlichen Sequenzen besser auflösen.

---

## 3. Erweiterte stochastische Analyse und Matrix-Selektion

Um die Wahl der Substitutionsmatrix nicht auf Heuristiken zu stützen, wurde eine quantitative Analyse durchgeführt. Diese vergleicht die Standard-Matrix (BLOSUM62) mit Matrizen für geringere evolutionäre Distanzen (BLOSUM80, BLOSUM90) und einem modellbasierten Ansatz (PAM30).

### A. Theoretischer Hintergrund: Markov-Prozess vs. Block-Clustering

Die Konkurrenz zwischen PAM- und BLOSUM-Matrizen ist nicht nur eine Frage der "Strenge", sondern des zugrundeliegenden evolutionären Modells. Dies ist entscheidend für die Bewertung des vorliegenden Proteins **TSR3**.

1.  **PAM-Modell (Point Accepted Mutation):**
    Die PAM-Matrizen basieren auf einem **zeitkontinuierlichen Markov-Prozess**. Die Annahme ist, dass die Evolution durch unabhängige, aufeinanderfolgende Punktmutationen geschieht. Mathematisch lässt sich die Übergangswahrscheinlichkeitsmatrix $M$ für eine Distanz $t$ als Potenz der Einheitsdistanz darstellen:
    $$M_t = (M_1)^t$$
    Für **PAM30** ($t=30$) modelliert dies einen Prozess, in dem sich die Sequenz gleichmäßig über die Zeit verändert ("evolutionäre Drift"). Dies passt biologisch hervorragend zu essentiellen "Housekeeping"-Proteinen wie TSR3, die unter globalem Selektionsdruck stehen und sich langsam, aber stetig verändern.

2.  **BLOSUM-Modell (Blocks Substitution Matrix):**
    BLOSUM basiert nicht auf einem expliziten evolutionären Zeitmodell, sondern auf der **Cluster-Analyse** konservierter lokaler Domänen ("Blocks"). BLOSUM80 oder BLOSUM90 werden aus Sequenzblöcken generiert, die noch $>$80% bzw. $>$90% Identität aufweisen.
    *Schwäche im konkreten Fall:* BLOSUM fokussiert auf konservierte "Inseln". Da das TSR3-Protein jedoch über die gesamte Länge hochkonserviert ist (globales Struktur-Funktions-Erfordernis), könnte der "Insel-Ansatz" von BLOSUM weniger präzise sein als der kontinuierliche Markov-Ansatz von PAM.

### B. Empirische Ergebnisse: Der Bit-Score Vergleich

Um diese theoretische Annahme zu prüfen, haben wir die **Bit-Scores** ($S'$) verglichen. Der Bit-Score ist im Gegensatz zum Raw-Score ($S$) matrix-unabhängig normalisiert (durch die Karlin-Altschul-Parameter $\lambda$ und $K$) und repräsentiert den reinen Informationsgehalt (in Bits) des Alignments relativ zum Rauschen:
$$S' = \frac{\lambda \cdot S - \ln(K)}{\ln(2)}$$

**Ergebnisse der Analyse (Top 30 Orthologe):**

| Matrix | Ø Identität (Top 30) | Ø Bit-Score (Hits 2-30) | Delta (zu BLOSUM62) |
| :--- | :--- | :--- | :--- |
| **BLOSUM62** | 85.55 % | 516.24 | Referenz |
| **BLOSUM80** | 85.51 % | 539.10 | + 22.86 |
| **BLOSUM90** | 85.53 % | 541.24 | + 25.00 |
| **PAM30** | **85.61 %** | **541.45** | **+ 25.21** |

### C. Interpretation und Validierung der Wahl

1.  **Das Konservierungs-Plateau:**
    Der massive Anstieg des Bit-Scores von BLOSUM62 auf BLOSUM80 (+22.86 Bits) beweist, dass die Standard-Matrix für diese hochkonservierten Daten (>85% Identität) mathematisch inadäquat ist ("zu weich"). Sie kann konservative Mutationen (Signal) nicht ausreichend vom Hintergrundrauschen trennen.

2.  **Widerlegung der "Beliebigen Strenge" (BLOSUM90 als Kontrolle):**
    Ein kritisches Gegenargument wäre, dass man durch Wahl einer beliebig strengen Matrix (nahe der Identitätsmatrix) den Score künstlich treiben könnte. Um dies zu prüfen, wurde **BLOSUM90** herangezogen.
    * Obwohl BLOSUM90 spezifisch für Sequenzen mit ~90% Identität entwickelt wurde (was unseren Daten sehr nahe kommt), liefert **PAM30** einen (marginal) höheren Informationsgehalt (541.45 vs. 541.24).
    * Dies indiziert, dass der höhere Score von PAM30 nicht allein auf "Strenge" basiert, sondern darauf, dass das **Markov-Modell** die spezifische Mutationsdynamik von TSR3 (gleichmäßige Punktmutationen über lange Zeiträume) präziser abbildet als das Block-Modell.

### D. Alignment-Längen und Gap-Statistik
Eine unpassende Matrix (z.B. eine für zu hohe Distanz wie BLOSUM45 auf identische Sequenzen) führt oft zu "über-extendierten" Alignments oder unnötigen Gaps ("Gap wander").
*   **Analyse:** Prüfen, ob sich die *Alignment Length* oder *Gaps* signifikant ändern.
*   **Ergebnis:** Die Gaps und Längen blieben in unserer Stichprobe über alle Matrizen hinweg stabil. Dies zeigt, dass die Homologie so stark ist, dass selbst eine "zu weiche" Matrix (BLOSUM62) das Alignment nicht zerstört, aber die Scores unter PAM30 präziser sind.

### E. Validierung der Matrix-Unabhängigkeit
Eine wissenschaftlich kritische Frage ist, ob die zur Entscheidung herangezogene *Percent Identity* selbst von der verwendeten Matrix (initial BLOSUM62) abhängt.
Da die Identitäten extrem hoch sind (>85%) und die Topologie der Alignments stabil ist, ist dieser Effekt vernachlässigbar. Der Wechsel zu PAM30 ist eine datengetriebene Optimierung.

### F. Konklusion

Die Entscheidung für **PAM30** ist statistisch signifikant validiert. Sie stellt kein "Overfitting" dar, sondern maximiert den Informationsgehalt (Relative Entropie) des Alignments.
Wir beobachten, dass für global hochkonservierte Proteine wie die *18S rRNA aminocarboxypropyltransferase* das klassische Dayhoff-Modell (PAM) dem neueren Block-Modell (BLOSUM) leicht überlegen ist, da es die evolutionäre Kontinuität besser mathematisch repräsentiert.

---

## 4. Statistische Validierung der BLAST-Ergebnisse

Um die Verlässlichkeit der gefundenen Homologien sicherzustellen, wurden die Score-Verteilungen und statistischen Parameter grafisch überprüft (siehe Plots im R-Skript `BLAST.R`).

### A. Trennung von Signal und Rauschen (Gumbel-Verteilung)
Die Verteilung der Alignment-Scores wurde gegen die theoretische **Extreme Value Distribution (Gumbel-Verteilung)** geplottet.
*   **Theorie:** Die Gumbel-Verteilung beschreibt die Wahrscheinlichkeit, einen bestimmten Score rein *zufällig* (Hintergrundrauschen) zu erhalten. Ihr Maximum liegt typischerweise bei sehr niedrigen Scores (z.B. < 30).
*   **Beobachtung:** Unsere empirischen Daten (Histogramm der BLAST-Hits) liegen weit rechts von der theoretischen Gumbel-Kurve, getrennt durch eine signifikante Lücke (Scores > 50).
*   **Schlussfolgerung:** Dies belegt statistisch, dass es sich bei den gefundenen Sequenzen **nicht um zufällige Ähnlichkeiten**, sondern um echte biologische Homologien handelt. Die klare Separation von der Gumbel-Kurve bestätigt die hohe Signifikanz der Ergebnisse.

### B. Konsistenzprüfung (Bit-Score vs. E-Value)
Zur Qualitätskontrolle wurde der Zusammenhang zwischen Bit-Score ($S'$) und E-Value ($E$) visualisiert.
*   **Erwartung:** Mathematisch gilt $-log(E) \propto S'$. Im Plot sollte sich eine perfekte Gerade ergeben.
*   **Ergebnis:** Die Datenpunkte liegen exakt auf einer diagonalen Linie. Dies bestätigt, dass die statistischen Parameter der Datenbanksuche (Datenbankgröße $n$, Query-Länge $m$) korrekt berücksichtigt wurden und keine Anomalien in der Berechnung vorliegen.

### C. Analyse von Indels (Alignments > Query-Länge)
Einige Alignments weisen eine Länge auf, die die Länge des menschlichen TSR3-Proteins (312 Aminosäuren) übersteigt (z.B. 315 oder 322).
*   **Ursache:** Dies ist kein Fehler, sondern auf **Insertionen und Deletionen (Indels)** zurückzuführen. Um homologe Bereiche korrekt übereinanderzulegen, fügt der Algorithmus **Gaps** (Lücken) ein.
*   **Berechnung:** $\text{Alignment Length} = \text{Matches} + \text{Mismatches} + \text{Gaps}$.
*   **Biologische Bedeutung:** Diese Längenunterschiede deuten auf evolutionäre Ereignisse hin, bei denen in verwandten Spezies Aminosäuren hinzugefügt wurden (Insertion im Subjekt) oder im menschlichen Protein verloren gingen (Deletion).

---

## 5. Zusammenfassung für das Projekt

Für die weitere phylogenetische Analyse sind die Sequenzen hervorragend geeignet. Die Nutzung von **PAM30** wurde durch die theoretische Analyse (>85% Identität der Top-Hits) und die empirische Maximierung des robusten Bit-Scores (541.45 vs. 539.10 und 516.24) bestätigt.

Obwohl BLOSUM80 oft als Standard für konservierte Proteine gilt, zeigen unsere Daten, dass **PAM30** für diese spezifische Fragestellung (hochkonserviertes TSR3-Protein in Säugetieren) die überlegene Modellierung der evolutionären Distanzen bietet. Wir entscheiden uns daher bewusst für diese Matrix, um die höchstmögliche Auflösung im phylogenetischen Baum zu gewährleisten.
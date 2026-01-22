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

Um die Wahl der Substitutionsmatrix über reine Heuristiken hinaus zu validieren, wurde eine quantitative Analyse durchgeführt. Diese kombiniert einen klassischen Bit-Score-Vergleich mit einer untersuchung der Informationsdichte (Relative Entropie), um das optimale evolutionäre Modell für das vorliegende Protein **TSR3** zu identifizieren.

### A. Theoretischer Hintergrund: Markov-Prozess vs. Block-Clustering

Die Konkurrenz zwischen PAM- und BLOSUM-Matrizen ist fundamental durch das zugrundeliegende evolutionäre Modell bestimmt:

1.  **PAM-Modell (Point Accepted Mutation):**
    PAM-Matrizen basieren auf einem **zeitkontinuierlichen Markov-Prozess**. Die Evolution wird als Folge unabhängiger Punktmutationen modelliert. Mathematisch ist die Übergangswahrscheinlichkeitsmatrix $M$ für eine Distanz $t$ definiert als:
    $$M_t = (M_1)^t$$
    Dies passt biologisch hervorragend zu essentiellen "Housekeeping"-Proteinen wie TSR3, die unter globalem Selektionsdruck stehen und sich langsam, aber stetig durch evolutionäre Drift verändern ("Global Alignment"-Charakteristik).

2.  **BLOSUM-Modell (Blocks Substitution Matrix):**
    BLOSUM basiert auf der Cluster-Analyse hochkonservierter lokaler Domänen ("Blocks") ohne explizites Zeitmodell. Da TSR3 jedoch über die gesamte Sequenzlänge hochkonserviert ist und nicht nur in Inseln, könnte der lokale Cluster-Ansatz von BLOSUM die globale Mutationsdynamik weniger präzise abbilden als der Markov-Ansatz.

### B. Empirische Ergebnisse I: Der Bit-Score Vergleich

Zunächst wurden die **Bit-Scores** ($S'$) der Top-30 Orthologe verglichen. Der Bit-Score ist matrix-unabhängig normalisiert und repräsentiert den reinen Informationsgehalt des Alignments in Bits:

**Tabelle 1: Vergleich der Top-30 Orthologe**

| Matrix       | Ø Identität | Ø Bit-Score | Delta (zu BLOSUM62) |
| :----------- | :---------- | :---------- | :------------------ |
| **BLOSUM62** | 85.55 %     | 516.24      | Referenz            |
| **PAM70**    | 85.59 %     | 530.66      | + 14.42             |
| **BLOSUM80** | 85.51 %     | 539.10      | + 22.86             |
| **BLOSUM90** | 85.53 %     | 541.24      | + 25.00             |
| **PAM30**    | 85.61 %     | 541.45      | + 25.21             |

*Befund:* Alle "strengeren" Matrizen liefern signifikant höhere Scores als der Standard BLOSUM62. Dies bestätigt, dass die Sequenzen hochkonserviert sind (>85% Identität). Rein nach Bit-Score scheinen PAM30 und BLOSUM90 die stärksten Kandidaten zu sein.

### C. Empirische Ergebnisse II: Validierung durch Relative Entropie

Ein rein score-basierter Vergleich birgt das Risiko des "Overfittings" (Wahl einer zu strengen Matrix). Um die **Modell-Passform** zu prüfen, haben wir die theoretische Relative Entropie ($H$) der Matrizen mit der tatsächlich erzielten Informationsdichte (Bits per Residue) verglichen.

Die theoretische Entropie $H$ einer Matrix gibt an, wie viel Information (in Bits) pro Position erwartet wird, wenn das Alignment der Ziel-Divergenz der Matrix entspricht.

**Analyse der Informationsdichte (siehe Abb. "Matrix-Passform"):**

1.  **BLOSUM62 (Under-fitting):**
    * *Theoretisches $H$:* ~0.7 Bits
    * *Beobachtet:* ~1.6 Bits
    * *Interpretation:* Die Matrix erwartet weit entfernte Sequenzen. Da unsere Daten jedoch hochkonserviert sind, wird die Erwartung massiv übertroffen. Die Matrix ist "zu weich" und trennt Signal nicht scharf genug vom Rauschen.

2.  **PAM30 (Over-fitting / Zu streng):**
    * *Theoretisches $H$:* ~2.57 Bits
    * *Beobachtet:* ~1.7 Bits
    * *Interpretation:* Hier liegt die theoretische Erwartung weit *über* dem beobachteten Wert. PAM30 erwartet fast identische Sequenzen. Da unsere Top-30 Hits zwar sehr ähnlich, aber nicht identisch sind, bestraft PAM30 die vorhandenen Variationen zu stark. Die Matrix ist für die vorliegende Datenstruktur zu restriktiv.

3.  **Der "Sweet Spot" (PAM70 & BLOSUM90):**
    * Bei **PAM70** und **BLOSUM90** decken sich das theoretische $H$ und die beobachteten Bits/Residue fast perfekt (~1.6 - 1.7 Bits).
    * Dies indiziert, dass diese Matrizen genau die evolutionäre Distanz modellieren, die unsere Daten tatsächlich aufweisen. Modell und Realität sind im Einklang.

### D. Konklusion und Matrix-Selektion

Basierend auf der kombinierten Analyse entscheiden wir uns für **PAM70**.

* **Gegen PAM30:** Obwohl PAM30 den höchsten nominalen Bit-Score liefert, zeigt die Entropie-Analyse, dass die Matrix theoretisch zu streng ist ($H_{theo} \gg H_{obs}$). Dies birgt die Gefahr, dass Alignments bei leicht divergenteren Sequenzen (außerhalb der Top 30) instabil werden oder künstlich verkürzt werden.
* **Gegen BLOSUM90:** BLOSUM90 zeigt eine ähnlich gute Passform der Entropie wie PAM70. Die Entscheidung fällt jedoch aufgrund der unter **(A)** beschriebenen theoretischen Überlegung zugunsten von PAM: Da es sich bei TSR3 um ein Protein mit kontinuierlicher evolutionärer Historie handelt, ist das Markov-Modell der PAM-Serie dem blockbasierten Ansatz vorzuziehen.

**Ergebnis:** Die **PAM70**-Matrix stellt den robustesten Kompromiss dar. Sie maximiert die Informationsausbeute, ohne (wie PAM30) zu streng zu strafen, und bildet die zugrundeliegende Mutationsdynamik (Markov-Prozess) biologisch adäquater ab als die BLOSUM-Serie.

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

Für die weitere phylogenetische Analyse sind die Sequenzen hervorragend geeignet. Die Nutzung von **PAM70** wurde durch die theoretische Analyse (>85% Identität der Top-Hits) und die empirische Maximierung des robusten Bit-Scores (diesen teil überarbeiten) bestätigt.

Obwohl BLOSUM80 oft als Standard für konservierte Proteine gilt, zeigen unsere Daten, dass **PAM70** für diese spezifische Fragestellung (hochkonserviertes TSR3-Protein in Säugetieren) die überlegene Modellierung der evolutionären Distanzen bietet. Wir entscheiden uns daher bewusst für diese Matrix, um die höchstmögliche Auflösung im phylogenetischen Baum zu gewährleisten.
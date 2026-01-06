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

## 3. Vorschlag für erweiterte stochastische Analysen

Um die Wahl der Matrix nicht nur auf "Gefühl" oder Standardregeln zu stützen, können folgende empirische Methoden angewandt werden:

### A. Information Content (Relative Entropie) Analyse
Nach Altschul (1991) ist eine Matrix dann optimal, wenn ihre relative Entropie $H$ der tatsächlichen Divergenz der gefundenen Sequenzen entspricht.
*   **Methode:** Man berechnet die durchschnittliche Identität der Top-Hits (hier ~85%).
*   **Vergleich:** Man wählt die Matrix, deren theoretische Ziel-Identität (Target Frequency) diesem Wert am nächsten kommt.
    *   BLOSUM62: ~62% Identität
    *   BLOSUM80: ~80% Identität
    *   PAM30:    ~90% Identität (kurze Distanzen)
*   **Fazit:** Da unsere Hits >>62% liegen, sind BLOSUM80 oder PAM30 mathematisch die korrekte Wahl.

### B. Bit-Score Distribution & E-Value Sensitivität
Der Bit-Score $S'$ ist matrix-unabhängig normalisiert und erlaubt den direkten Vergleich.
*   **Vorgehen:** Man vergleicht die Bit-Scores des Top-Hits (Homo sapiens) für alle drei Matrizen.
*   **Beobachtung:**
    *   **BLOSUM62:** 628
    *   **PAM30:** 647
    *   **BLOSUM80:** 656
*   **Interpretation:** Der Score steigt beim Wechsel von BLOSUM62 zu den strengeren Matrizen (PAM30, BLOSUM80) deutlich an. Dies bestätigt, dass das Standard-Modell (BLOSUM62) nicht optimal passt. Interessanterweise liefert **BLOSUM80** den höchsten Score, was darauf hindeutet, dass PAM30 möglicherweise bereits *zu* strikt für die leichten Variationen ist oder BLOSUM80 das "Sweet Spot" der Konservierung für dieses Protein am besten abbildet.

### C. Alignment-Längen und Gap-Statistik
Eine unpassende Matrix (z.B. eine für zu hohe Distanz wie BLOSUM45 auf identische Sequenzen) führt oft zu "über-extendierten" Alignments oder unnötigen Gaps ("Gap wander").
*   **Analyse:** Prüfen, ob sich die *Alignment Length* oder *Gaps* signifikant ändern.
*   **Ergebnis:** Die Gaps und Längen blieben in unserer Stichprobe stabil. Dies zeigt, dass die Homologie so stark ist, dass selbst eine "zu weiche" Matrix (BLOSUM62) das Alignment nicht zerstört, aber die Scores unter BLOSUM80/PAM30 präziser sind.

### D. Validierung der Matrix-Unabhängigkeit (Methodische Kritik)
Eine wissenschaftlich kritische Frage ist, ob die zur Entscheidung herangezogene *Percent Identity* selbst von der verwendeten Matrix (initial BLOSUM62) abhängt und somit das Ergebnis verfälscht.

*   **Theoretische Abhängigkeit:** Ja, die berechnete Identität $\frac{\text{Matches}}{\text{Länge}}$ hängt vom Alignment ab.
*   **Empirische Robustheit:** In diesem Projekt handelt es sich um einen **iterativen Validierungsprozess**.
    1.  **Schritt 1 (Heuristik):** Der initiale Scan mit BLOSUM62 zeigte bereits extrem hohe Identitäten.
    2.  **Stabilität:** Bei derart hoher Sequenzähnlichkeit ist das Alignment topologisch robust (Matches bleiben Matches).
    3.  **Schlussfolgerung:** Der "Fehler" durch die Matrixwahl bei der Identitätsberechnung ist vernachlässigbar. Der Wechsel zu BLOSUM80 ist daher eine **Optimierung der statistischen Power**, basierend auf der initialen Beobachtung.

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

Für die weitere phylogenetische Analyse sind die Sequenzen hervorragend geeignet. Die Nutzung von **BLOSUM80** wurde durch die theoretische Analyse (>80% Identität) und die empirische Maximierung des Bit-Scores (656 vs 628/647) bestätigt. Sie stellt das optimale Modell für die Analyse dieser hochkonservierten Proteinfamilie dar.
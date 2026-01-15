# Analyse-Leitfaden: Multiple Sequence Alignment (MSA)

Dieser Bericht dient als "Skizze" und Leitfaden für die wissenschaftliche Auswertung des MSA. Er verbindet biologische Interpretation mit algorithmischer Kritik.

## 1. Einleitung & Zielsetzung
Das MSA ist keine absolute Wahrheit, sondern eine **Hypothese der Homologie**.
*   **Ziel:** Anordnung der Sequenzen so, dass jede Spalte (Residue) evolutionär voneinander abstammt.
*   **Biologische Bedeutung:** Konservierte Regionen deuten oft auf funktionelle Wichtigkeit (aktive Zentren, Bindestellen) oder strukturelle Stabilität (Disulfidbrücken, Core-Packing) hin. Variable Regionen liegen oft an der Proteinoberfläche und sind toleranter gegenüber Mutationen.

## 2. Visuelle Inspektion & Qualitätssicherung
*Hier sollen die Beobachtungen eingetragen werden, die beim Betrachten der Alignment-Visualisierung (PDF/ggmsa) gemacht wurden.*

### A. Globale Konservierung
*   **Beobachtung:** Wie hoch ist die generelle Übereinstimmung über die gesamte Länge?
*   **Interpretation:** TSR3 ist ein essentielles Protein für die Ribosomenbiogenese.
    *   *Hypothese:* Wir erwarten eine hohe Konservierung, besonders in der katalytischen Domäne (Aminocarboxypropyltransferase).
    *   *Check:* Gibt es Bereiche, die über alle Spezies hinweg zu 100% identisch sind?

### B. Gap-Muster (Lücken)
*   **Beobachtung:** Wo treten Gaps (Insertionen/Deletionen) auf?
*   **Biologische Plausibilität:**
    *   Gaps sollten biologisch eher in "Loops" (Schlaufenverbindungen) auftreten, nicht mitten in Sekundärstrukturen (Alpha-Helices, Beta-Faltblätter).
    *   *Kritischer Blick:* Gibt es "treppenartige" Gaps am Anfang oder Ende? Dies sind oft Artefakte des Algorithmus.

## 3. Algorithmische Kritik: ClustalW
Das verwendete Verfahren (ClustalW) ist ein **progressiver Algorithmus**. Es ist wichtig, dessen mathematische Grenzen zu diskutieren.

### A. Der "Greedy"-Ansatz
*   **Funktionsweise:** ClustalW berechnet zuerst einen "Guide Tree" und aligniert dann schrittweise die ähnlichsten Sequenzen zuerst.
*   **Das Problem ("Once a gap, always a gap"):** Ein Fehler, der früh im Prozess (bei den ähnlichsten Sequenzen) gemacht wird, kann später nicht mehr korrigiert werden.
*   **Auswirkung auf deine Daten:** Wenn du weit entfernte Spezies hast, könnte ClustalW diese suboptimal alignieren, weil es durch die initialen Entscheidungen der nahen Verwandten "voreingenommen" ist.

### B. Vergleich mit iterativen Methoden (Diskussion)
*   Für eine fundierte wissenschaftliche Diskussion sollte erwähnt werden, dass modernere Algorithmen existieren:
    *   **MUSCLE / T-Coffee / MAFFT:** Diese sind **iterativ**. Sie bauen ein Alignment, bewerten es, und bauen es um, um eine mathematische Zielfunktion (Objective Function) zu maximieren.
    *   *Fazit:* Für die Top 30 Hits (die vermutlich recht ähnlich sind), ist ClustalW meist ausreichend ("good enough"). Bei sehr divergenten Sequenzen wäre es eine Fehlerquelle.

## 4. Identifikation funktioneller Domänen
*   **Aufgabe:** Suche nach spezifischen Motiven.
*   **Tipp:** Vergleiche dein Alignment mit Datenbanken wie *Pfam* oder *InterPro* für TSR3.
*   **Frage:** Sind die Regionen, die in der Literatur als wichtig für die 18S rRNA Modifikation beschrieben sind, in deinem Alignment konserviert?

## 5. Zusammenfassung der Ergebnisse
*   Ist das Alignment vertrauenswürdig genug als Basis für den phylogenetischen Baum?
*   Müssen bestimmte Sequenzen ausgeschlossen werden (z.B. Fragmente, die viel zu kurz sind und riesige Gaps erzeugen)?

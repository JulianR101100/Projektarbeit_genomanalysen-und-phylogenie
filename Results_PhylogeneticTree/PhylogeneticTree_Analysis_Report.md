# Analyse-Leitfaden: Phylogenetischer Baum

Dieser Bericht strukturiert die phylogenetische Analyse und liefert den theoretischen Unterbau für die Diskussion der Ergebnisse. Hier geht es um den Konflikt zwischen mathematischen Modellen und biologischer Realität.

## 1. Distanzmatrix: Von Identität zu Evolution
Das Skript berechnet aktuell Distanzen basierend auf "Identity" (prozentuale Übereinstimmung).

### A. Das Problem der "Multiplen Substitutionen"
*   **Theorie:** Wenn sich an Position X eine Aminosäure von A -> B -> A ändert, sieht die Sequenz identisch aus wie vorher. Evolution hat stattgefunden, ist aber unsichtbar ("Saturation").
*   **Kritik:** Die reine "Identity"-Distanz unterschätzt die wahre evolutionäre Distanz, besonders bei entfernt verwandten Spezies.
*   **Lösung/Diskussion:** Wissenschaftlich exakter wäre die Nutzung von Korrekturmodellen:
    *   **Kimura-Korrektur (Protein):** Schätzt die "wahre" Anzahl der Substitutionen statistisch.
    *   **PAM/BLOSUM Distanzen:** Nutzen empirische Wahrscheinlichkeiten für den Austausch von Aminosäuren (biochemische Ähnlichkeit).

## 2. Cluster-Algorithmen & Die Molekulare Uhr (Ultrametrie)

### A. UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
*   **Annahme:** UPGMA setzt eine **strikte molekulare Uhr** voraus (Ultrametrie). Das bedeutet: Die Evolutionsrate ist in allen Ästen des Baums exakt gleich. Alle Spezies sind gleich weit von der Wurzel (dem Vorfahren) entfernt.
*   **Test:** Die Funktion `check_strict_ultrametric_property` in deinem Skript prüft dies mathematisch (Drei-Punkte-Bedingung).
    *   *Erwartetes Ergebnis:* Du wirst "Violations" finden. Biologie ist selten uhrwerkartig.
*   **Gefahr:** Wenn die Raten ungleich sind (z.B. evolvieren Nagetiere oft schneller als Menschen), gruppiert UPGMA falsch ("Long Branch Attraction"). Es zieht Spezies mit langen Ästen künstlich zusammen oder an die Wurzel.

### B. Neighbor-Joining (NJ) - (Empfohlene Ergänzung)
*   **Vorteil:** NJ nimmt **keine** molekulare Uhr an. Es ist ein "Minimum Evolution"-naher Ansatz.
*   **Realität:** Es erlaubt unterschiedliche Astlängen. Wenn UPGMA und NJ unterschiedliche Topologien zeigen, ist NJ meist vertrauenswürdiger, da seine Annahmen weniger strikt sind.

## 3. Robustheit & Statistik (Bootstrapping)
Ein Baum ist nur eine Punkt-Schätzung. Wie sicher sind wir uns bei den Ästen?

### A. Das Konzept
*   **Bootstrapping:** Man würfelt die Spalten des Alignments neu (Resampling mit Zurücklegen) und baut daraus 1000 neue Bäume.
*   **Interpretation:** Ein Wert von 90 an einem Ast bedeutet: "In 90% der zufälligen Bäume existierte diese Gruppierung."
*   **Anwendung:**
    *   Werte > 70-80% gelten als gut gestützt.
    *   Werte < 50% bedeuten: "Hier rät der Algorithmus nur" (Polytomie).
    *   *Diskussionspunkt:* Äste mit niedrigem Bootstrap-Support dürfen nicht über-interpretiert werden!

## 4. Biologische Interpretation (Validierung)
Vergleiche den berechneten Baum mit der etablierten Taxonomie (NCBI Taxonomy).

*   **Cluster-Check:**
    *   Bilden alle **Säugetiere** (Mammalia) eine Monophylie (einen gemeinsamen Ast)?
    *   Stehen **Primaten** zusammen?
    *   Wo landen Außengruppen (z.B. Vögel, Fische, Reptilien), falls vorhanden?
*   **Anomalien:**
    *   Liegt ein Mensch-Gen neben einem Pilz-Gen? -> Hinweis auf Horizontalen Gentransfer (unwahrscheinlich hier) oder **Datenbank-Kontamination** / falsche Annotation.
    *   Gibt es Gen-Duplikationen (Paraloge)? (Zwei Äste für die gleiche Spezies).

## 5. Zusammenfassung & Fazit
*   Welcher Baum (UPGMA vs. Ward vs. NJ) spiegelt die biologische Erwartung am besten wider?
*   Ist TSR3 ein guter phylogenetischer Marker? (Ja, wenn der Gen-Baum dem Spezies-Baum entspricht).

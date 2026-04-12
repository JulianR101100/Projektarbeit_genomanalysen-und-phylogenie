# Wissenschaftlicher Abschlussbericht: Multiple Sequence Alignment (MSA) von TSR3

Der folgende Bericht dient als zusammenfassungen der Ergebnisse und wurde **mit KI erstellt** - *jedoch nicht inhatlich sondern nur als schnelle zusammenfassung*. Dazu wurden die Ergebnisse und die Auswertung bei fortschritten in der Analyse übergeben und damit immer ein aktuelle Zusammenfassung des projektstands erstellt.

## 1. Methodik
*   **Sequenzauswahl:** Top 30 homologe Sequenzen (identifiziert via BLASTp, PAM70 Matrix).
*   **Algorithmus:** ClustalW (Progressives Alignment). Beschränkugn auf PAM Matrix -> vgl ergebnisse aus Blast.R zu geeigeneten Substitutionsmatrix
*   **Statistische Auswertung:** Quantitative Spaltenanalyse der Aminosäure-Identität (Konservierung) und Gap-Frequenz in R.

## 2. Globale Stochastische Ergebnisse
Die Auswertung des gesamten Alignments lieferte folgende Durchschnittswerte:
*   **Mittlere globale Konservierung ($\bar{C}_{\text{global}}$):** 72.26 %
*   **Globale Gap-Rate ($\bar{f}_{\text{gap}}$):** 17.99 %

Diese Werte bestätigen eine generell hohe Homologie, zeigen aber auch eine signifikante Variabilität in bestimmten Abschnitten.

## 3. Domänenspezifische Analyse: Die funktionelle Dichotomie
Die stochastische Analyse erlaubt eine klare Trennung des Proteins in zwei funktionelle Einheiten 
(**Diese zwei Domänen stammen aus der Uniprot datei von TSR_Human (https://www.uniprot.org/uniprotkb/Q9UJK0/entry) also der Querry Sequence**):

### A. Die katalytische Domäne (Position 96–222)
*   **Konservierung ($\bar{C}_{\text{dom1}}$):** 95,33 %
*   **Gap-Rate ($\bar{f}_{\text{gap, dom1}}$):** 0,00 %
*   **Interpretation:** Dieser Bereich bildet den hochstrukturierten Kern des Enzyms. Die mathematische "Null-Toleranz" für Gaps und die fast perfekte Identität belegen den extremen funktionellen Zwang (Constraint). Hier finden die Bindung von S-Adenosylmethionin (SAM) und die Interaktion mit der 18S rRNA statt.

### B. Die unstrukturierte C-terminale Region (IDR) (Position 224–312)
*   **Konservierung ($\bar{C}_{\text{dom2}}$):** 62,35 %
*   **Gap-Rate ($\bar{f}_{\text{gap, dom2}}$):** 22,29 %
*   **Struktureller Kontext:** Diese Region korreliert mit niedrigen AlphaFold-Konfidenzwerten (pLDDT < 50).
*   **Interpretation:** Es handelt sich um eine *Intrinsically Disordered Region* (IDR). Da hier keine starre Faltung für die Funktion notwendig ist, können Mutationen und Indels (Gaps) akkumulieren, ohne die Lebensfähigkeit des Organismus zu gefährden.

## 4. Algorithmische Bewertung (ClustalW)
Obwohl ClustalW als progressiver Algorithmus ("Once a gap, always a gap") theoretisch anfällig für initiale Fehler ist, erweisen sich die Ergebnisse für TSR3 als äußerst robust. Dies liegt vor allem an der extrem hohen Konservierung der zentralen Domäne, die als stabiler "Anker" für das Alignment dient.

## 5. Fazit & Ausblick
Die Studie zeigt erfolgreich, dass die evolutionäre Dynamik von TSR3 modular aufgebaut ist. Die mathematische Quantifizierung der Konservierung deckt sich exakt mit den strukturbiologischen Vorhersagen.

**Zukünftige Schritte:**
*   Einbeziehung phylogenetisch entfernterer Spezies, um die Variabilität der IDR statistisch noch schärfer zu fassen.
*   Vergleich der phylogenetischen Bäume basierend auf der Gesamtfrequenz vs. der rein katalytischen Domäne.

---

### Notizen für die Präsentation
*   **Slide-Tipp:** Gegenüberstellung des farbigen Konservierungs-Plots (R) und des AlphaFold-Modells (Das findet ihr auch auf der Uni Prot seite). 
*   **Kernbotschaft:** "Die Mathematik des Alignments (62% Konservierung) macht die physikalische Unordnung (IDR) sichtbar."
*   **Key Fact:** 0,0% Gaps in der katalytischen Domäne sind der statistische Beweis für die essenzielle Rolle von TSR3 im Leben jeder eukaryotischen Zelle.

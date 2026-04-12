# **KI Bericht zur allgemeinen Vorrecherche: molekularen Charakterisierung, Pathologie und evolutionären Signifikanz des Gens TSR3 (18S rRNA Aminocarboxypropyltransferase)**


## **Executive Summary**

Dieser Forschungsbericht präsentiert eine umfassende Analyse des durch Ihre BLAST-Suche identifizierten Gens *TSR3* und seines korrespondierenden Proteinprodukts, der 18S rRNA Aminocarboxypropyltransferase. Das Protein, das lange Zeit als "Orphan"-Gen mit unbekannter Funktion galt, wurde in jüngster Zeit als der entscheidende enzymatische Akteur in der finalen Reifung der kleinen ribosomalen Untereinheit (40S) bei Eukaryoten identifiziert. Seine Hauptfunktion besteht in der Katalyse einer chemisch hochkomplexen Hypermodifikation – der Übertragung einer Aminocarboxypropyl-Gruppe auf ein prä-modifiziertes Nukleotid im ribosomalen Decoding-Center (Position U1248 im Menschen).

Die vorliegende Untersuchung geht weit über die bloße Identifikation hinaus. Sie beleuchtet die kritische Rolle von TSR3 als "molekularen Torwächter" der Translationsfidelity und analysiert die verheerenden Konsequenzen seiner Dysfunktion, die von spezifischen Anfälligkeiten in der Leukämogenese (NPM1-mutierte AML) bis hin zu komplexen neurologischen Entwicklungsstörungen im Rahmen des 16p13.3-Deletionssyndroms reichen. Besonderes Augenmerk liegt auf der evolutionären Perspektive: Der Bericht kontrastiert die eukaryotische Abhängigkeit von dieser komplexen Modifikation mit den eleganten alternativen Lösungsstrategien von Prokaryoten und extrem reduzierten Parasiten (Microsporidien), die ohne dieses Gen überleben. Diese vergleichende Analyse offenbart TSR3 nicht nur als essentielles Housekeeping-Gen, sondern als einen evolutionären Marker für die zunehmende Komplexität der Qualitätskontrolle in der eukaryotischen Proteinbiosynthese.

## **1\. Molekulare Identität und Enzymatische Mechanik**

### **1.1 Das Zielsubstrat: Die Hypermodifikation von Helix 31**

Die eindeutige Identifizierung von *TSR3* (in früheren Datenbankeinträgen oft als *C16orf42* oder *20S rRNA accumulation protein 3 homolog* geführt) verweist auf einen hochspezialisierten Prozess im Zytoplasma menschlicher Zellen.1 Das Ziel dieses Enzyms ist nicht irgendein Protein, sondern ein einzelnes, spezifisches Nukleotid innerhalb der 18S ribosomalen RNA (rRNA), dem strukturellen Kern der kleinen ribosomalen Untereinheit.

Im Menschen handelt es sich um das Uridin an Position 1248 (U1248), im Modellorganismus *Saccharomyces cerevisiae* um U1191.3 Dieses Nukleotid befindet sich an der Spitze von Helix 31 (h31), einer Struktur, die oft als der "Schnabel" (Beak) der 40S-Untereinheit bezeichnet wird. Diese Region ragt direkt in das Decoding-Center des Ribosoms hinein, jenen Ort, an dem die mRNA-Codons von den tRNA-Anticodons gelesen werden.

Die Modifikation an dieser Stelle ist chemisch so anspruchsvoll, dass sie nicht in einem einzigen Schritt erfolgen kann. Es handelt sich um eine sogenannte "Hypermodifikation", die in drei sequenziellen biochemischen Phasen abläuft, wobei TSR3 den finalen Schritt exekutiert:

1. **Pseudouridylierung:** Zunächst wird das Uridin durch den snoRNP-Komplex (in Hefe snR35) isomerisiert, wobei die glykosidische Bindung aufgebrochen und neu geknüpft wird, um Pseudouridin ($\\Psi$) zu bilden.4  
2. **Methylierung:** Das Enzym Nep1 (Emg1) methyliert dieses Pseudouridin an der N1-Position, wodurch 1-Methylpseudouridin ($m^1\\Psi$) entsteht. Mutationen in diesem Schritt führen beim Menschen bereits zum Bowen-Conradi-Syndrom, was die kritische Natur dieser Kaskade unterstreicht.6  
3. **Aminocarboxypropyl-Transfer:** Hier greift TSR3 ein. Es erkennt das vor-modifizierte $m^1\\Psi$ und überträgt eine 3-Amino-3-carboxypropyl-Gruppe (acp) auf die N3-Position.1

Das Endprodukt ist **1-Methyl-3-(3-amino-3-carboxypropyl)pseudouridin ($m^1acp^3\\Psi$)**. Diese Modifikation ist einzigartig, da sie dem Ribosom eine extrem sperrige, chemisch komplexe Seitenkette hinzufügt, die sowohl eine positive Ladung (Aminogruppe) als auch eine negative Ladung (Carboxylgruppe) trägt (zwitterionischer Charakter).

### **1.2 Strukturbiologie: Ein evolutionärer "Hack" der SPOUT-Superfamilie**

Die Kristallstrukturanalyse von TSR3-Homologen (z.B. aus dem Archaeon *Vulcanisaeta distributa*) hat spektakuläre Einblicke in die Evolution von Enzymen geliefert. TSR3 gehört strukturell zur **SPOUT-Superfamilie** (SpoU-TrmD) von Methyltransferasen.3 Diese Enzyme sind typischerweise dafür bekannt, eine Methylgruppe von S-Adenosylmethionin (SAM) auf RNA zu übertragen.

TSR3 stellt jedoch eine faszinierende Ausnahme dar. Obwohl es die charakteristische "Knoten"-Struktur der SPOUT-Klasse beibehält, bindet es den Kofaktor SAM in einer völlig neuartigen Orientierung. In klassischen Methyltransferasen wird SAM so gebunden, dass die Methylgruppe exponiert ist. Im katalytischen Zentrum von TSR3 hingegen ist SAM so positioniert, dass die Methylgruppe verborgen bleibt, während die lange Aminocarboxypropyl-Seitenkette des Methionin-Teils für den nukleophilen Angriff durch die rRNA präsentiert wird.3

**Wissenschaftliche Einordnung:** Dies ist ein Paradebeispiel für divergente Evolution, bei der ein bestehendes Proteingerüst (SPOUT-Fold) "repurposed" (zweckentfremdet) wurde, um eine völlig neue chemische Reaktion zu katalysieren – den Transfer einer C4-Kette statt einer C1-Gruppe. Dabei entsteht als Nebenprodukt nicht das übliche S-Adenosylhomocystein (SAH), sondern 5'-Methylthioadenosin (MTA).1

### **1.3 Lokalisation und Qualitätskontrolle**

Im Gegensatz zu den frühen rRNA-Modifikationen (wie 2'-O-Methylierung), die bereits im Nukleolus während der Transkription stattfinden, ist TSR3 im **Zytoplasma** lokalisiert.3 Diese räumliche Trennung ist von höchster funktioneller Relevanz. Sie impliziert, dass die $acp^3$-Modifikation einer der allerletzten Schritte in der Reifung der kleinen ribosomalen Untereinheit ist. TSR3 fungiert somit als "Qualitätskontroll-Posten". Erst wenn das Ribosom fast vollständig assembliert und in das Zytoplasma exportiert wurde, wird diese letzte, entscheidende chemische Gruppe hinzugefügt, um das Ribosom "scharf zu schalten" für die Translation.9

## **2\. Physiologische Relevanz und Funktion im Körper**

Die Relevanz dieses Proteins für den Körper erschließt sich aus der zentralen Rolle des Ribosoms. Ein Defekt in TSR3 legt nicht einfach nur ein Enzym lahm, sondern beeinträchtigt die fundamentale Fähigkeit der Zelle, den genetischen Code präzise in Proteine zu übersetzen.

### **2.1 Stabilisierung des Decoding-Centers**

Die Modifikation $m^1acp^3\\Psi$ befindet sich im Herzen des Ribosoms, direkt dort, wo die genetische Entscheidung über die Aminosäuresequenz fällt.

* **Strukturelle Integrität:** Die sperrige acp-Gruppe interagiert mit dem Rückgrat der rRNA und stabilisiert die Haarnadelstruktur von Helix 31\. Ohne diese Stabilisierung wäre diese Schleife zu flexibel ("wobbelig"), was die Präzision der tRNA-Bindung gefährden würde.10  
* **Magnesium-Koordination:** Die Carboxylgruppe der acp-Seitenkette ist in der Lage, Magnesiumionen ($Mg^{2+}$) zu koordinieren. Magnesium ist essenziell, um die negativen Ladungen der RNA-Phosphate abzuschirmen und eine dichte Packung der RNA-Struktur zu ermöglichen. Fehlt TSR3, fehlt diese ionische Brücke, was zu einer lokalen Destabilisierung führt.12

### **2.2 Der Rio2-Checkpoint: Ein molekularer Sicherheitsschalter**

Neuere Studien haben eine direkte mechanistische Verbindung zwischen TSR3 und dem Reifungsfaktor **Rio2** aufgedeckt.9 Rio2 ist eine ATPase/Kinase, die an die unreife 40S-Untereinheit bindet und verhindert, dass diese verfrüht mit der 60S-Untereinheit fusioniert oder die Translation startet.

Die Daten legen nahe, dass die Aktivität von TSR3 (die Installation der acp-Gruppe) das Signal für die **Ablösung von Rio2** ist.

* **Im Normalzustand:** TSR3 modifiziert Helix 31 \-\> Konformationsänderung \-\> Rio2 löst sich ab \-\> Ribosom ist reif für die Translation.  
* **Bei TSR3-Defekt:** Die Modifikation fehlt \-\> Rio2 bleibt gebunden oder bindet erneut \-\> Die 40S-Untereinheit bleibt in einem unreifen Komplex gefangen oder tritt als defekter Komplex in den Translationspool ein.9

Dies erklärt, warum TSR3-Defekte nicht sofort tödlich sind, aber zu einer langsamen Akkumulation von Translationsfehlern und einer ineffizienten Ribosomenbiogenese führen. Die Ribosomen sind zwar da, aber sie sind "verstopft" durch alte Reifungsfaktoren oder strukturell suboptimal.

## **3\. Pathologie: Krankheiten, Gendefekte und Klinische Phänotypen**

Obwohl *TSR3* in den klassischen Datenbanken (wie OMIM) noch nicht als alleiniger Verursacher eines benannten Syndroms mit eigenem Eintrag geführt wird, zeigt die tiefergehende Literaturrecherche eine signifikante Beteiligung an mehreren pathologischen Prozessen. Es ist entscheidend, hier zwischen somatischen Mutationen (Krebs) und konstitutionellen Defekten (Erbkrankheiten) zu unterscheiden.

### **3.1 Das 16p13.3-Deletionssyndrom (ATR-16 und assoziierte Phänotypen)**

Das Gen *TSR3* befindet sich auf dem kurzen Arm von Chromosom 16 (Locus 16p13.3).2 Diese Region ist klinisch hochrelevant, da sie häufig von chromosomalen Mikrodeletionen betroffen ist. Patienten mit einer Deletion in 16p13.3 zeigen ein komplexes Krankheitsbild, zu dem *TSR3* als haploinsuffizientes Gen beiträgt.

* **Der Phänotyp:** Patienten zeigen typischerweise eine intellektuelle Beeinträchtigung (Intellectual Disability, ID), Entwicklungsverzögerungen, Mikrozephalie (kleiner Kopfumfang) und charakteristische Gesichtszüge (Dysmorphien wie Hypertelorismus oder breite Nasenwurzel).13 Oft liegt auch eine Alpha-Thalassämie vor, da die benachbarten Hämoglobin-Gene (*HBA1/HBA2*) ebenfalls gelöscht sind.  
* **Die Rolle von TSR3:** Während Gene wie *CREBBP* (Rubinstein-Taybi-Syndrom) Teile des Phänotyps erklären, wird die neurologische Komponente zunehmend mit Defekten der Ribosomenbiogenese in Verbindung gebracht. Neuronen haben einen enorm hohen Bedarf an Proteinsynthese für das axonale Wachstum und die synaptische Plastizität. Eine Reduktion der TSR3-Dosis um 50% (Haploinsuffizienz) führt wahrscheinlich zu einem Engpass in der Ribosomenproduktion, der sich spezifisch in der neuronalen Entwicklung manifestiert.13 Dies passt zum Bild anderer Ribosomopathien, die häufig kraniofaziale und neurologische Defekte zeigen.

### **3.2 NPM1-mutierte Akute Myeloische Leukämie (AML): Eine therapeutische Achillesferse**

Eine der spannendsten Entdeckungen der letzten Jahre ist die Rolle von *TSR3* als essentielle Abhängigkeit (Dependency) in einer spezifischen Form von Blutkrebs.

* **Der Mechanismus:** Etwa 30% aller AML-Patienten tragen eine Mutation im Gen *NPM1* (Nucleophosmin). Diese Mutation (NPM1c) führt dazu, dass das Protein NPM1 aus dem Nukleolus in das Zytoplasma verlagert wird. NPM1 ist normalerweise ein Chaperon für die Ribosomenbiogenese. Sein Fehlen im Nukleolus setzt die Krebszellen unter enormen "nukleolären Stress".15  
* **Die TSR3-Abhängigkeit:** Um trotz dieses Stresses zu überleben und zu proliferieren, werden diese Leukämiezellen hyper-abhängig von den verbleibenden Reifungsfaktoren, insbesondere von TSR3. CRISPR-Screens haben gezeigt, dass ein Knockout von TSR3 in normalen Zellen toleriert wird, in NPM1-mutierten AML-Zellen jedoch den sofortigen Zelltod auslöst (synthetische Letalität).15  
* **Klinische Relevanz:** Dies macht TSR3 zu einem potenziellen Ziel für neue Krebstherapien. Inhibitoren, die die enzymatische Tasche von TSR3 blockieren, könnten spezifisch Leukämiezellen abtöten, während gesunde Körperzellen überleben. Studien zeigen bereits, dass die Hemmung von TSR3 die Zellen wieder sensitiv für Chemotherapeutika wie Venetoclax macht.16

### **3.3 Abgrenzung zu Skelettdysplasien (Wichtige Klärung)**

In der Literatur 28 taucht der Begriff "TSR3" auch im Kontext von Skeletterkrankungen wie der Weill-Marchesani-ähnlichen Dysplasie oder Geleophysischen Dysplasie auf (verursacht durch Mutationen in ADAMTSL2).  
Wichtige Unterscheidung: Hierbei bezieht sich "TSR3" auf die Thrombospondin Type 1 Repeat Domain 3, eine Proteinstrukturdomäne innerhalb des ADAMTSL2-Proteins, nicht auf das Gen TSR3 (Ribosome Maturation Factor). Diese terminologische Überschneidung ist eine häufige Falle in der Bioinformatik. Für Ihren Bericht ist es essentiell festzuhalten: Das Gen TSR3 verursacht keine Skelettdysplasie durch Matrix-Defekte, sondern wirkt primär zellulär über das Ribosom. Die in 28 beschriebenen Mutationen betreffen ein völlig anderes Protein.

### **3.4 Autoimmunität und GlycoRNA**

Eine neuartige Forschungslinie verbindet die Modifikation acp3U mit Autoimmunerkrankungen. Es wurde entdeckt, dass acp3U auf der Oberfläche von tRNAs (und potenziell auch rRNA-Fragmenten) als Ankerpunkt für Zucker-Ketten (Glykane) dient – sogenannte GlycoRNAs.17  
Diese GlycoRNAs befinden sich auf der Zelloberfläche. Wenn die acp3U-Modifikation dort "nackt" (ohne Zucker) exponiert wird oder durch Zellzerfall freigesetzt wird, kann sie vom angeborenen Immunsystem (über TLR7-Rezeptoren) als "fremd" erkannt werden. Dies könnte erklären, warum Defekte in RNA-modifizierenden Enzymen manchmal zu Autoimmunreaktionen führen (wie bei Lupus), da der Körper beginnt, seine eigene RNA anzugreifen.18

## **4\. Evolutionäre Konservierung und Strategien des "Fehlens"**

Ihre Frage nach den Strategien weit entfernter Lebewesen, die dieses Gen nicht besitzen, führt zu einer faszinierenden evolutionären Dichotomie.

### **4.1 Die Eukaryotische Lösung: Komplexität**

Alle Eukaryoten (vom Menschen bis zur Hefe) und viele Archaeen besitzen das TSR3-Gen und die entsprechende Modifikation $acp^3U$ in ihrer rRNA.3 Dies deutet darauf hin, dass diese komplexe Modifikation ein uraltes Merkmal ist, das wahrscheinlich schon im letzten gemeinsamen Vorfahren der Eukaryoten (LECA) vorhanden war. Sie dient der Feinabstimmung eines immer komplexer werdenden Ribosoms.

### **4.2 Die Bakterielle Lösung: Die "Methylierungs-Strategie"**

Hier wird es spannend: Bakterien (wie *E. coli*) besitzen **kein** Homolog von TSR3, das an der rRNA wirkt. Sie haben auch keine $acp^3U$-Modifikation in ihrer 16S rRNA (dem Gegenstück zur eukaryotischen 18S rRNA) an der entsprechenden Position.20

Wie kompensieren Bakterien das Fehlen?  
Anstatt einer einzigen, chemisch massiven und teuren Modifikation (acp3U) an der Spitze von Helix 31, nutzen Bakterien zwei kleinere Methylierungen an benachbarten Nukleotiden:

1. **m²G966:** Eine Methylgruppe am Guanosin (Position 966).  
2. **m⁵C967:** Eine Methylgruppe am Cytosin (Position 967).

Diese werden durch die Enzyme RsmD und RsmB installiert.21  
Die Strategie: Bakterien setzen auf "verteilte Stabilisierung". Zwei kleine hydrophobe Methylgruppen sorgen für eine ausreichende Versteifung der RNA-Schleife durch Basen-Stapelung (Base Stacking). Das ist energetisch günstiger und reicht für die etwas simplere Initiationsmaschinerie der Bakterien aus. Eukaryoten hingegen benötigen die "schwere Artillerie" der acp-Gruppe, vermutlich weil ihr Translationsstart (mit Dutzenden Initiationsfaktoren) mechanisch anspruchsvoller ist und eine rigidere Struktur im Decoding-Center verlangt.

### **4.3 Das TapT-Paradoxon**

Interessanterweise besitzen Bakterien durchaus die *chemische Fähigkeit*, acp3U herzustellen. Sie haben ein Gen namens **TapT** (oder YfiP), das entfernt mit TSR3 verwandt ist.20

* **Der Unterschied:** TapT modifiziert **tRNA** (an Position 47), aber niemals die rRNA.  
* **Evolutionäre Logik:** Bakterien nutzen die Stabilisierungskraft von acp3U selektiv für ihre "Adapter" (tRNAs), haben sich aber evolutionär dagegen entschieden, sie im "Maschinenraum" (Ribosom) zu verbauen. Dies zeigt, dass Evolution Module (wie die acp-Synthese) für völlig unterschiedliche Zwecke in verschiedenen Domänen des Lebens einsetzen kann.

### **4.4 Strategien reduzierter Genome: Der Fall der Microsporidien**

Microsporidien (z.B. *Encephalitozoon cuniculi*) sind parasitäre Pilze mit extrem reduzierten Genomen. Sie leben intrazellulär und haben alles "Unnötige" über Bord geworfen.

* **Das Problem:** Ihre Ribosomen sind radikal verkleinert, ihnen fehlen viele rRNA-Expansionssegmente, die bei anderen Eukaryoten für Stabilität sorgen.24  
* **Die Kompensations-Strategie:** Da sie sich energetisch teure Modifikationsenzyme wie TSR3 oft sparen oder deren Substrate (rRNA-Helices) degeneriert sind, nutzen sie **strukturelle Mimikry**.  
  * **Protein statt RNA:** Wo rRNA-Helices instabil sind, haben sich ribosomale Proteine (z.B. msL1) so entwickelt, dass sie den leeren Raum physisch ausfüllen und die Stabilität übernehmen, die normalerweise durch RNA-Struktur oder chemische Modifikation gewährleistet würde.26  
  * **Hibernation Factors:** Sie nutzen Faktoren wie Mdf1, um das Ribosom in Ruhephasen (Sporen) zu "verklammern" und so die Integrität ohne ständige biochemische Wartung zu sichern.24

## **5\. Zusammenfassung und Ausblick**

Das Gen *TSR3* ist ein faszinierendes Beispiel für die versteckte Komplexität unserer Zellen. Was als bloßes Ergebnis einer BLAST-Suche begann, entpuppt sich als zentraler Baustein der zellulären Qualitätskontrolle.

**Zentrale Erkenntnisse für Ihre Recherche:**

1. **Funktion:** TSR3 ist ein spezialisierter "Endfertiger" des Ribosoms. Es installiert eine chemisch einzigartige Antenne im Lesekopf der Proteinfabrik.  
2. **Krankheit:** Ein Ausfall führt nicht zum sofortigen Tod, sondern zu einem schleichenden Funktionsverlust (Haploinsuffizienz), der sich besonders in Geweben mit hoher Teilungsrate (Blut) oder komplexer Entwicklung (Gehirn) zeigt (16p13.3 Deletion, AML-Abhängigkeit).  
3. **Evolution:** Eukaryoten leisten sich den energetischen Luxus dieser Modifikation für maximale Präzision. Bakterien wählen die "Billig-Lösung" (einfache Methylierung). Parasiten mit minimalen Genomen ersetzen chemische Stabilität durch Protein-Stützräder.

Für die moderne Medizin ist TSR3 ein "Target of Interest": Seine Blockade könnte eine Achillesferse von Leukämien treffen, während seine Rolle bei Entwicklungsstörungen hilft, die komplexen Symptome von Patienten mit Chromosom-16-Defekten besser zu verstehen.

### **Tabellarische Übersicht: Strategien der Stabilisierung von Helix 31**

| Organismus-Domäne | Enzym | Substrat & Position | Modifikation | Strategie / Lösung |
| :---- | :---- | :---- | :---- | :---- |
| **Eukaryoten** (Mensch) | **TSR3** | 18S rRNA (U1248) | **m¹acp³Ψ** (Hypermodifikation) | Maximale strukturelle Stabilität und Ladungskontrolle für komplexe Regulation. |
| **Archaeen** (*Haloferax*) | **Tsr3** | 16S rRNA | **acp³U** | Nutzung derselben komplexen Chemie für Thermostabilität. |
| **Bakterien** (*E. coli*) | **Keines** (für rRNA) | 16S rRNA | **m²G \+ m⁵C** (Methylierung) | "Verteilte Stabilität": Zwei simple Methylgruppen ersetzen die komplexe acp-Gruppe. Kosteneffizient. |
| **Bakterien** (*E. coli*) | **TapT** | tRNA (Pos. 47\) | **acp³U** | Chemische Fähigkeit vorhanden, aber nur für tRNA-Stabilität genutzt. |
| **Microsporidien** | Reduziert/Variabel | Degenerierte rRNA | Protein-Substitution | Ersatz von RNA-Chemie durch stützende Proteine (msL1) oder Hibernationsfaktoren. |

*Quellennachweise:* 1

#### **Referenzen**

1. TSR3 \- 18S rRNA aminocarboxypropyltransferase \- Homo sapiens (Human) | UniProtKB, Zugriff am Januar 8, 2026, [https://www.uniprot.org/uniprotkb/Q9UJK0/entry](https://www.uniprot.org/uniprotkb/Q9UJK0/entry)  
2. TSR3 TSR3 ribosome maturation factor \[Homo sapiens (human ..., Zugriff am Januar 8, 2026, [https://www.ncbi.nlm.nih.gov/gene?Db=gene\&Cmd=DetailsSearch\&Term=115939](https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=115939)  
3. Ribosome biogenesis factor Tsr3 is the aminocarboxypropyl transferase responsible for 18S rRNA hypermodification in yeast and humans \- PubMed, Zugriff am Januar 8, 2026, [https://pubmed.ncbi.nlm.nih.gov/27084949/](https://pubmed.ncbi.nlm.nih.gov/27084949/)  
4. Ribosome biogenesis factor Tsr3 is the aminocarboxypropyl transferase responsible for 18S rRNA hypermodification in yeast and humans \- ResearchGate, Zugriff am Januar 8, 2026, [https://www.researchgate.net/publication/301327909\_Ribosome\_biogenesis\_factor\_Tsr3\_is\_the\_aminocarboxypropyl\_transferase\_responsible\_for\_18S\_rRNA\_hypermodification\_in\_yeast\_and\_humans](https://www.researchgate.net/publication/301327909_Ribosome_biogenesis_factor_Tsr3_is_the_aminocarboxypropyl_transferase_responsible_for_18S_rRNA_hypermodification_in_yeast_and_humans)  
5. Ribosome biogenesis factor Tsr3 is the aminocarboxypropyl transferase responsible for 18S rRNA hypermodification in yeast and humans \- PubMed Central, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC4872110/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4872110/)  
6. Effects of the Bowen-Conradi syndrome mutation in EMG1 on its nuclear import, stability and nucleolar recruitment | Human Molecular Genetics | Oxford Academic, Zugriff am Januar 8, 2026, [https://academic.oup.com/hmg/article/25/24/5353/2593613](https://academic.oup.com/hmg/article/25/24/5353/2593613)  
7. Mutation of a Gene Essential for Ribosome Biogenesis, EMG1, Causes Bowen-Conradi Syndrome | Request PDF \- ResearchGate, Zugriff am Januar 8, 2026, [https://www.researchgate.net/publication/26235618\_Mutation\_of\_a\_Gene\_Essential\_for\_Ribosome\_Biogenesis\_EMG1\_Causes\_Bowen-Conradi\_Syndrome](https://www.researchgate.net/publication/26235618_Mutation_of_a_Gene_Essential_for_Ribosome_Biogenesis_EMG1_Causes_Bowen-Conradi_Syndrome)  
8. S-Adenosylmethionine: more than just a methyl donor \- PMC \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC10491745/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10491745/)  
9. The modifying enzyme Tsr3 establishes the hierarchy of Rio kinase binding in 40S ribosome assembly \- PMC \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC8925970/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8925970/)  
10. A Comparison of the Crystal Structures of Eukaryotic and Bacterial SSU Ribosomal RNAs Reveals Common Structural Features in the Hypervariable Regions \- PMC, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC3364965/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3364965/)  
11. Selection of Peptides Targeting Helix 31 of Bacterial 16S Ribosomal RNA by Screening M13 Phage-Display Libraries \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC6259748/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6259748/)  
12. Biogenesis and functions of aminocarboxypropyluridine in tRNA \- PMC \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC6895100/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6895100/)  
13. Distal Partial Trisomy 15q26 and Partial Monosomy 16p13.3 in a 36-Year-Old Male with Clinical Features of Both Chromosomal Abnormalities \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC5463451/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5463451/)  
14. (PDF) A pure de novo 16p13.3 duplication and amplification in a patient with femoral hypoplasia, psychomotor retardation, heart defect, and facial dysmorphism-a case report and literature review of the partial 16p13.3 trisomy syndrome \- ResearchGate, Zugriff am Januar 8, 2026, [https://www.researchgate.net/publication/366714869\_A\_pure\_de\_novo\_16p133\_duplication\_and\_amplification\_in\_a\_patient\_with\_femoral\_hypoplasia\_psychomotor\_retardation\_heart\_defect\_and\_facial\_dysmorphism-a\_case\_report\_and\_literature\_review\_of\_the\_partial\_](https://www.researchgate.net/publication/366714869_A_pure_de_novo_16p133_duplication_and_amplification_in_a_patient_with_femoral_hypoplasia_psychomotor_retardation_heart_defect_and_facial_dysmorphism-a_case_report_and_literature_review_of_the_partial_)  
15. Posttranscriptional depletion of ribosome biogenesis factors engenders therapeutic vulnerabilities in NPM1-mutant AML | Blood | American Society of Hematology \- ASH Publications, Zugriff am Januar 8, 2026, [https://ashpublications.org/blood/article/146/10/1239/537926/Posttranscriptional-depletion-of-ribosome](https://ashpublications.org/blood/article/146/10/1239/537926/Posttranscriptional-depletion-of-ribosome)  
16. Breaking ribosomes to fight leukemia | Blood | American Society of Hematology, Zugriff am Januar 8, 2026, [https://ashpublications.org/blood/article/146/10/1155/547029/Breaking-ribosomes-to-fight-leukemia](https://ashpublications.org/blood/article/146/10/1155/547029/Breaking-ribosomes-to-fight-leukemia)  
17. The modified RNA base acp3U is an attachment site for N-glycans in glycoRNA \- Glycobiology Research & Training Center, Zugriff am Januar 8, 2026, [https://grtc.ucsd.edu/education/courses/The-modified-RNA-base-acp3U-is-an-attachment-site-for-N-glycans-in-glycoRNA.pdf](https://grtc.ucsd.edu/education/courses/The-modified-RNA-base-acp3U-is-an-attachment-site-for-N-glycans-in-glycoRNA.pdf)  
18. The modified RNA base acp3U is an attachment site for N-glycans in glycoRNA | Request PDF \- ResearchGate, Zugriff am Januar 8, 2026, [https://www.researchgate.net/publication/383323764\_The\_modified\_RNA\_base\_acp3U\_is\_an\_attachment\_site\_for\_N-glycans\_in\_glycoRNA](https://www.researchgate.net/publication/383323764_The_modified_RNA_base_acp3U_is_an_attachment_site_for_N-glycans_in_glycoRNA)  
19. Phenotypic characterization of yeast TSR3 deletion (Δtrs3) and human... \- ResearchGate, Zugriff am Januar 8, 2026, [https://www.researchgate.net/figure/Phenotypic-characterization-of-yeast-TSR3-deletion-Dtrs3-and-human-TSR3-depletion\_fig2\_301327909](https://www.researchgate.net/figure/Phenotypic-characterization-of-yeast-TSR3-deletion-Dtrs3-and-human-TSR3-depletion_fig2_301327909)  
20. Identification of the 3-amino-3-carboxypropyl (acp) transferase enzyme responsible for acp3U formation at position 47 in Escherichia coli tRNAs \- PubMed Central, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7026641/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026641/)  
21. Impact of methylations of m2G966/m5C967 in 16S rRNA on bacterial fitness and translation initiation \- PMC \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC3439901/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3439901/)  
22. Role of the Ribosomal P-Site Elements of m2G966, m5C967, and the S9 C-Terminal Tail in Maintenance of the Reading Frame during Translational Elongation in Escherichia coli \- PMC \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC3754560/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3754560/)  
23. Molecular basis of tRNA aminocarboxypropyl-transferase TapT for substrate recognition \- PMC \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12630136/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12630136/)  
24. Adaptation to genome decay in the structure of the smallest eukaryotic ribosome, Zugriff am Januar 8, 2026, [https://www.researchgate.net/publication/358269677\_Adaptation\_to\_genome\_decay\_in\_the\_structure\_of\_the\_smallest\_eukaryotic\_ribosome](https://www.researchgate.net/publication/358269677_Adaptation_to_genome_decay_in_the_structure_of_the_smallest_eukaryotic_ribosome)  
25. Muller's Ratchet and Ribosome Degeneration in the Obligate Intracellular Parasites Microsporidia \- Semantic Scholar, Zugriff am Januar 8, 2026, [https://pdfs.semanticscholar.org/1f31/1554aacd3b29435b7bf699be78308f34293c.pdf](https://pdfs.semanticscholar.org/1f31/1554aacd3b29435b7bf699be78308f34293c.pdf)  
26. A Conserved Ribosomal Protein Has Entirely Dissimilar Structures in Different Organisms, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC10764239/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10764239/)  
27. Ribosomopathies: Global process, tissue specific defects \- PMC \- NIH, Zugriff am Januar 8, 2026, [https://pmc.ncbi.nlm.nih.gov/articles/PMC4590025/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4590025/)  
28. DISEASES \- TSR3 \- JensenLab, Zugriff am Januar 8, 2026, [https://diseases.jensenlab.org/Entity?documents=10\&type1=9606\&id1=ENSP00000007390\&type2=-26\&id2=DOID:0111724](https://diseases.jensenlab.org/Entity?documents=10&type1=9606&id1=ENSP00000007390&type2=-26&id2=DOID:0111724)
# Library and Settings ----
rm(list = ls())
library(ggplot2)
# library(tidyr)
# library("tidyverse")
library(Biostrings)
library(pwalign) # Weil Biositrings manche sub-matritzen nicht mehr untersützt
library(seqinr)
library(ordinal)
library(dplyr)
library(stringr)
library(readr)
library(reshape2)


# Select Working Directory ----
wds <- data.frame(
  system = c("Julian Windows PC/Laptop", 
             "Linux Uni",
             "Johann",
             "Rphi"),
  path = c("C:/Users/julia/OneDrive/Integrated Life Siences/Genomanalysen und Phylogenie/Projektarbeit",
           "/home/stud/ha24vepa/Documents/Bash_Linux_Introduction_supplements/Genom_und_Phylogenie_kurs/...",
           "C:/Users/johan/Uni/Genomanalyse_Projekt/Projektarbeit",
           "C:/Users/darko/OneDrive/Dokumente/Genomanalysen_und_Phylogenie_Projekt"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(wds))) {
  if (dir.exists(wds$path[i])) {
    setwd(wds$path[i])
    message("Working Directory set for ", wds$system[i])
    break
  } else {
    message("Path for ", wds$system[i], " not found.")
  }
}
rm(wds, i)

# ==== Blast Ergebnisse Einlesen (verschiedene Substitutionsmatritzen) ====

# BLOSUM 62
blast_data_62 <- read.table("Data/query_seq_blast_BLOSUM62.blast")
colnames(blast_data_62) <- c("Accession_ID", "Alignment_length", "Percent_identity", 
                             "Gaps_number", "Score", "Bit_score", "E_value")

# BLOSUM80
blast_data_80 <- read.table("Data/query_seq_blast_BLOSUM80.blast")
colnames(blast_data_80) <- c("Accession_ID", "Alignment_length", "Percent_identity", 
                             "Gaps_number", "Score", "Bit_score", "E_value")

# BLOSUM90
blast_data_90 <- read.table("Data/query_seq_blast_BLOSUM90.blast")
colnames(blast_data_90) <- c("Accession_ID", "Alignment_length", "Percent_identity", 
                             "Gaps_number", "Score", "Bit_score", "E_value")

# PAM70
blast_data_70 <- read.table("Data/query_seq_blast_PAM70.blast")
colnames(blast_data_70) <- c("Accession_ID", "Alignment_length", "Percent_identity", 
                             "Gaps_number", "Score", "Bit_score", "E_value")

# PAM30
blast_data_30 <- read.table("Data/query_seq_blast_PAM30.blast")
colnames(blast_data_30) <- c("Accession_ID", "Alignment_length", "Percent_identity", 
                             "Gaps_number", "Score", "Bit_score", "E_value")


# ==== Daten Anotieren =====
annotate_blast_results <- function(blast_results, fasta_db_path = "Data/Rohdaten/database.fasta", 
                                   e_value_threshold, create_csv = FALSE) {
  
  # 1. Error fals es fehler mit datenbank etc gibt (is i guess besser als R fehler)
  if (!file.exists(fasta_db_path)) {
    stop("Fehler: Die Datenbank-Datei '", fasta_db_path, "' wurde nicht gefunden.")
  }
  
  # 2. Metadaten aus FASTA-Headern extrahieren
  # Wir lesen direkt nur die Header-Zeilen (starten mit ">")
  raw_headers <- grep("^>", readLines(fasta_db_path), value = TRUE)
  
  # Bereinigung: ">" entfernen
  headers_clean <- str_sub(raw_headers, 2)
  
  # Parsing: ID, Organismus und Beschreibung extrahieren
  # Erstellt direkt den Lookup-Dataframe
  annotation_df <- data.frame(
    Accession_ID = word(headers_clean, 1),
    Organism     = str_extract(headers_clean, '(?<="organism_name":")[^"]+'),
    Description  = str_extract(headers_clean, '(?<="description":")[^"]+'),
    stringsAsFactors = FALSE
  )
  
  # 3. Verknüpfung mit BLAST-Daten & Filterung
  annotated_data <- blast_results %>%
    left_join(annotation_df, by = "Accession_ID") %>%
    filter(E_value <= e_value_threshold) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, Accession_ID, Organism, Description, everything()) # Wichtige Spalten nach vorne ziehen
  
  # Info-Output für den User
  message(sprintf("Annotation abgeschlossen. %d Treffer verbleiben (E-Value <= %s).", 
                  nrow(annotated_data), as.character(e_value_threshold)))
  
  # 4. Optional: Speichern
  if (create_csv) {
    output_file <- "Results_BLAST/anotated_top_hits.csv"
    write_csv(annotated_data, output_file)
    message("Datei gespeichert unter: ", output_file)
  }
  
  return(annotated_data)
}

# Anotieren der Blast ergenisse mit Pfad zur FASTA-Datenbank und E-Value Schwellenwert
# Ihr könnte die Ergebnisse in einer csv/ exel Datei speichern wenn create_csv = TRUE
annotated_blast_results <- annotate_blast_results(
  blast_results = blast_data_62,
  e_value_threshold = 1e-160,
  create_csv = FALSE)

# Top 30 Ergenisse -> Achtet darauf, wenn e value zu nidrig ist, werden natürlich nicht ganz 30 gezeigt
top30_identity <- head(annotated_blast_results, 30)
top30_identity[, c("Organism", "Percent_identity", "Alignment_length", "E_value")]
# Es gibt 
blast_data_30[blast_data_30$Alignment_length > 312, ]
# Gen ist sehr gut konserviert.  nur zwei Primaten (9, 18) mit gaps (alignmentlen > 312)
top30_identity[top30_identity$Alignment_length > 312, c("Organism", "Percent_identity", "Alignment_length", "E_value")]

# ==== Blast Ergebnisse Visualisieren ====

#' Plot: Empirische Score-Verteilung vs. Theoretische Gumbel-Verteilung
#' @param lambda (Standard für BLOSUM62: 0.267)
#' @param K Karlin-Altschul Parameter (Standard für BLOSUM62: 0.041)
#' @param m Länge der Query-Sequenz (Standard TSR3: 312 AA)
#' @param n Größe der Datenbank in Residues (Standard Projekt-DB: ~156.222 AA)
plot_score_gumbel <- function(blast_data, lambda = 0.267, K = 0.041, m = 312, n = 156222) {
  
  # Bereinigung: Nur gültige Scores nutzen
  scores <- na.omit(blast_data$Score)
  
  # Theoretische Gumbel-Dichtefunktion definieren
  # Parameterberechnung nach Karlin-Altschul
  gumbel_density <- function(x) {
    # Location parameter u
    u <- (log(K * m * n)) / lambda
    # Scale parameter beta
    beta <- 1 / lambda
    
    # Gumbel PDF
    (1/beta) * exp(-(x - u)/beta) * exp(-exp(-(x - u)/beta))
  }
  
  # Plot mit ggplot2
  ggplot(data.frame(Score = scores), aes(x = Score)) +
    # 1. Histogramm auf Dichte-Skala normieren (damit es zur Kurve passt)
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "lightblue", color = "black", alpha = 0.6) +
    # 2. Theoretische Kurve
    stat_function(fun = gumbel_density, color = "red", linewidth = 1.2) +
    # Styling
    labs(
      title = "Verteilung der Raw-Scores vs. Gumbel-Theorie",
      subtitle = paste("Parameter: Lambda =", lambda, "| K =", K, "| DB-Größe =", n),
      x = "Raw Score (S)",
      y = "Dichte (Wahrscheinlichkeit)"
    ) +
    theme_minimal()
}

#' Plot: Bit-Score vs. E-Value (Signifikanz-Check)
#' @description Zeigt den logarithmischen Zusammenhang zwischen Score und E-Value.
#' Es sollte eine perfekte Gerade entstehen.
plot_bitscore_significance <- function(blast_data) {
  
  # Filtern: E-Values von 0 können nicht log-transformiert werden
  plot_data <- blast_data %>%
    filter(E_value > 0)
  
  ggplot(plot_data, aes(x = Bit_score, y = -log10(E_value))) +
    geom_point(alpha = 0.6, color = "darkblue") +
    labs(
      title = "Zusammenhang: Bit-Score vs. Statistische Signifikanz",
      subtitle = "Erwartung: Linearer Zusammenhang (-log10 E proportional zu Bit-Score)",
      x = "Bit Score (S')",
      y = "-log10(E-Value)"
    ) +
    theme_bw()
}

# 1. Gumbel-Verteilung Plotten (Test mit PAM30 für mehr "Rauschen")
# Wir nutzen blast_data_30, da hier kürzere/schwächere Alignments erwartet werden.
# "Problem" BLAST speichert nur signifikante Hits (datenbank sonst viel zu groß), liegen unsere Daten
# wahrscheinlich immer noch weit rechts vom theoretischen Gumbel-Maximum (Rauschen).

print(plot_score_gumbel(blast_data_30, lambda = 0.267, K = 0.041)) # Werte für PAM30 ggf. anpassen

# 2. Signifikanz
print(plot_bitscore_significance(annotated_blast_results))



# ==== Allgemiene Analyse für Wahl der Substitutionsmatrix ====
# 1. identitäten und allgemeiner Vergeleich
calc_mean_identity <- function(blast_data) {
  mean_identity <- mean(blast_data$Percent_identity)
  return(round(mean_identity, 2))
}
calc_target_identity <- function(blast_data, top_n = 30) {
  target_identity <- mean(head(blast_data$Percent_identity, top_n))
  return(round(target_identity, 2))
}

# mean_identity_all <- mean(blast_data_62$Percent_identity)
# target_identity <- mean(head(blast_data_62$Percent_identity, 30))

message("Durchschnittliche Identität (Alle Hits):",
        "\nBLOSUM62: ", calc_mean_identity(blast_data_62), "%",
        "\nBLOSUM80: ", calc_mean_identity(blast_data_80), "%",
        "\nBLOSUM90: ", calc_mean_identity(blast_data_90), "%",
        "\nPAM30:    ", calc_mean_identity(blast_data_30), "%",
        "\nPAM70:    ", calc_mean_identity(blast_data_70), "%")

message("Durchschnittliche Identität (Top 30 Hits):",
        "\nBLOSUM62: ", calc_target_identity(blast_data_62), "%",
        "\nBLOSUM80: ", calc_target_identity(blast_data_80), "%",
        "\nBLOSUM90: ", calc_target_identity(blast_data_90), "%",
        "\nPAM30:    ", calc_target_identity(blast_data_30), "%",
        "\nPAM70:    ", calc_target_identity(blast_data_70), "%")

get_mean_bitscore_top30 <- function(df) {
  # Wir nehmen Zeile 2 bis 30 (oder weniger, falls weniger Treffer existieren)
  n <- min(nrow(df), 30)
  return(round(mean(df$Bit_score[2:n]), 2))
}

message("Durchschnittlicher Bit-Score (Hits 2-30):",
        "\nBLOSUM62: ", get_mean_bitscore_top30(blast_data_62),
        "\nBLOSUM80: ", get_mean_bitscore_top30(blast_data_80),
        "\nBLOSUM90: ", get_mean_bitscore_top30(blast_data_90),
        "\nPAM30:    ", get_mean_bitscore_top30(blast_data_30),
        "\nPAM70:    ", get_mean_bitscore_top30(blast_data_70))

# ==== Entscheidung mittels (KL-Divergenz / Relative Entropie) ====
# 2. Definition der Theoretischen Entropie (H)
# Werte basierend auf Literatur (Altschul et al. / NCBI Dokumentation).
# H (in Bits) ist ein Maß für den Informationsgehalt der Matrix pro alignierter Position.
# Hohes H = Matrix für sehr ähnliche Sequenzen (kurze Evolution).
# Niedriges H = Matrix für entfernte Sequenzen (lange Evolution).
matrix_entropy <- data.frame(
  Matrix = c("PAM30", "PAM70", "BLOSUM90", "BLOSUM80", "BLOSUM62"),
  H_theoretical = c(2.57, 1.60, 1.18, 0.99, 0.70),
  stringsAsFactors = FALSE
)

# 2. Funktion zur Berechnung der Observed Bits per Residue
# Wir prüfen, wie viel Information wir tatsächlich pro Position extrahieren konnten.
calculate_observed_entropy <- function(blast_data, matrix_name) {
  # Wir betrachten die Top 30 Hits (Orthologe + evtl. Paraloge/Isoformen)
  # Um robust zu sein, nehmen wir die tatsächlichen Top 30 des jeweiligen Laufs.
  top30 <- head(blast_data, 30)
  
  # Berechnung: Bits per Residue = BitScore / Alignment Länge
  # Wenn die Matrix gut passt, sollte dieser Wert nahe an der theoretischen Entropie (H) liegen.
  # Ist er viel niedriger, bestraft die Matrix Mismatches/Gaps zu hart -> Matrix ist zu streng.
  top30$Bits_per_Residue <- top30$Bit_score / top30$Alignment_length
  
  # Durchschnitt über die Top 30
  mean_obs <- mean(top30$Bits_per_Residue, na.rm = TRUE)
  
  return(data.frame(
    Matrix = matrix_name,
    Observed_Bits_per_Residue = round(mean_obs, 3)
  ))
}

# Daten aggregieren
entropy_comparison <- rbind(
  calculate_observed_entropy(blast_data_30, "PAM30"),
  calculate_observed_entropy(blast_data_70, "PAM70"),
  calculate_observed_entropy(blast_data_90, "BLOSUM90"),
  calculate_observed_entropy(blast_data_80, "BLOSUM80"),
  calculate_observed_entropy(blast_data_62, "BLOSUM62")
)

# Merge mit theoretischen Werten
entropy_comparison <- merge(entropy_comparison, matrix_entropy, by = "Matrix")

# Differenz berechnen
entropy_comparison$Difference <- entropy_comparison$H_theoretical - entropy_comparison$Observed_Bits_per_Residue

# Vergleich: Observed vs. Theoretical Entropy (Bits per Residue)
entropy_comparison

# Interpretation / Entscheidungshilfe
# Wir suchen die Matrix, die:
# 1. Einen möglichst hohen Informationsgehalt (H) bietet (um feine Unterschiede aufzulösen).
# 2. Dabei aber keinen "Einbruch" bei den Observed Bits zeigt (Differenz sollte nicht riesig sein).
# PAM30 hat das höchste H (2.57). Wenn Observed ebenfalls hoch ist (~2.0+), 
# bestätigt das, dass die Sequenzen konserviert genug für diese strenge Matrix sind.

# Visualisierung
# Wir formen die Daten für ggplot um
plot_data <- pivot_longer(entropy_comparison, cols = c("H_theoretical", "Observed_Bits_per_Residue"), 
                          names_to = "Metric", values_to = "Bits")

# Faktor-Level ordnen für schöne Darstellung (nach theoretischer Härte)
plot_data$Matrix <- factor(plot_data$Matrix, levels = c("PAM30", "PAM70", "BLOSUM90", "BLOSUM80", "BLOSUM62"))

ggplot(plot_data, aes(x = Matrix, y = Bits, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Matrix-Passform: Theorie vs. Praxis",
    subtitle = "Vergleich der relativen Entropie (H) mit den tatsächlich erzielten Bits/Residue.\nIst Observed << Theoretical, ist die Matrix zu streng (viele Strafen).",
    y = "Bits per Residue",
    x = "Substitutionsmatrix",
    fill = "Metrik"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("H_theoretical" = "grey70", "Observed_Bits_per_Residue" = "dodgerblue3"),
                    labels = c("Theoretisches H (Max)", "Beobachtet (Top 30 Avg)"))

# FAZIT: 
# PAM30 bietet das höchste theoretische Potential (H=2.57).
# Da unser TSR3 Protein hoch konserviert ist, erreichen wir auch empirisch sehr hohe Werte.
# Wir nutzen also PAM70, um die maximale Information aus den Alignments zu ziehen.

# ==== Speichern der Top30 mit gewählter Substitutionsmatrix (PAM70) ====
# 1. Anotieren der PAM70 Ergebnisse
pam70_annotated <- annotate_blast_results(
  blast_results = blast_data_70,
  e_value_threshold = 1e-155,
  create_csv = FALSE
)

# 2. Top 30 selektieren
top30_pam70 <- head(pam70_annotated, 30)

# 3. Speichern der CSV
csv_output_path <- "Results_BLAST/top30_pam70_annotated.csv"
write_csv(top30_pam70, csv_output_path)
message("Top 30 PAM70 Ergebnisse gespeichert in: ", csv_output_path)

# 4. Sequenzen extrahieren und als FASTA speichern Extrahiere Sequenzen für MSA...

db_sequences <- readAAStringSet("Data/Rohdaten/database.fasta")

# Die Accession_IDs aus unseren Top 30 extrahieren
# Hinweis: Die Accession_ID im Dataframe (z.B. "9606_0:0048f7") muss exakt mit dem Namen im AAStringSet übereinstimmen.
# Die BLAST-Ergebnisse haben oft nur den ersten Teil vor dem Leerzeichen als ID, aber readAAStringSet liest den ganzen Header.
# Wir müssen sicherstellen, dass wir die richtigen Sequenzen finden.

# Wir bereinigen die Namen in db_sequences, damit sie matchen
# (Analog zur Logik in annotate_blast_results: ID ist das erste Wort)
names(db_sequences) <- word(names(db_sequences), 1)

# Filtern: Nur Sequenzen behalten, deren Name in unserer Top30 Liste ist
# Wir nutzen match, um auch die Reihenfolge beizubehalten (optional, aber nett)
target_ids <- top30_pam70$Accession_ID
matched_indices <- match(target_ids, names(db_sequences))

# Checken ob alle gefunden wurden
if(any(is.na(matched_indices))) {
  warning("Achtung: Nicht alle Top 30 IDs wurden in der Datenbank gefunden!")
}

top30_seqs_final <- db_sequences[na.omit(matched_indices)]

# Speichern
fasta_output_path <- "Results_MultibleSequenceAlign/top30_sequences.fasta"
writeXStringSet(top30_seqs_final, fasta_output_path)
message("Top 30 Sequenzen gespeichert in: ", fasta_output_path)


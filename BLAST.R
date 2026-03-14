# Library and Settings ----
rm(list = ls())
library(ggplot2)
library(tidyr)
library("tidyverse")
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
             "Raphi"),
  path = c("C:/Users/julia/OneDrive/Integrated Life Siences/Genomanalysen und Phylogenie/Projektarbeit",
           "/home/stud/ha24vepa/Documents/Bash_Linux_Introduction_supplements/Genom_und_Phylogenie_kurs/...",
           "C:/Users/johan/Uni/Genomanalyse_Projekt/Projektarbeit",
           "C:/Users/darko/OneDrive/Dokumente/Genomanalysen_und_Phylogenie_Projekt/Projektarbeit_genomanalysen-und-phylogenie"),
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

# ==== 0.1.Blast Ergebnisse Einlesen (verschiedene Substitutionsmatritzen) ====

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


# ==== 0.2.Daten Anotieren =====
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
    output_file <- "Results_BLAST/anotated_top30_Blosum90.csv"
    write_csv(annotated_data, output_file)
    message("Datei gespeichert unter: ", output_file)
  }
  
  return(annotated_data)
}

# Anotieren der Blast ergenisse mit Pfad zur FASTA-Datenbank und E-Value Schwellenwert
# Ihr könnte die Ergebnisse in einer csv/ exel Datei speichern wenn create_csv = TRUE
annotated_blast_results <- annotate_blast_results(
  blast_results = blast_data_70,
  e_value_threshold = 1e-160,
  create_csv = FALSE)

# Top 30 Ergenisse -> Achtet darauf, wenn e value zu nidrig ist, werden natürlich nicht ganz 30 gezeigt
top30_identity <- head(annotated_blast_results, 30)
top30_identity[, c("Organism", "Percent_identity", "Alignment_length", "E_value")]
# Es gibt 
blast_data_30[blast_data_30$Alignment_length > 312, ]
# Gen ist sehr gut konserviert.  nur zwei Primaten (9, 18) mit gaps (alignmentlen > 312)
top30_identity[top30_identity$Alignment_length > 312, c("Organism", "Percent_identity", "Alignment_length", "E_value")]

# ==== 1.  Blast Ergebnisse Visualisieren ====

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
#' @param blast_data
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

print(plot_score_gumbel(blast_data_30, lambda = 0.294, K = 0.110)) # Werte für PAM30 ggf. anpassen

# 2. Signifikanz
print(plot_bitscore_significance(annotated_blast_results))


# ==== 1.2.Estimate Alignment Score Distribution (Vergleich Empirische CDF mit Theoretischer CDF) ====
#anstatt histogramm vs. Gumbel

plot_ECDF_gumbel <- function(blast_data, lambda = 0.267, K = 0.041, m = 312, n = 156222) {

  # Nur gültige Scores
  scores_sorted <- sort(na.omit(blast_data$Score))
  
  # Empirische CDF
  N_hits <- length(scores_sorted)
  ecdf_values <- (1:N_hits) / N_hits
  
  plot(scores_sorted, ecdf_values,
       main = "Empirical CDF of BLAST Raw Scores vs. Theoretical Gumbel CDF",
       xlab = "Max Score",
       ylab = "CDF")
  
  # Gumbel-Parameter
  a_m <- log(K * m * n) / lambda
  b <- 1 / lambda
  
  # Theoretische Gumbel-CDF
  curve(pgumbel(x, loc = a_m, scale = b),
        add = TRUE, col = "red", lwd = 2)
}

# print(plot_ECDF_gumbel(blast_data_70, lambda = 0.294, K = 0.110))

#Die empirische CDF der beobachteten BLAST-Raw-Scores liegt deutlich rechts
#der theoretischen Gumbel-CDF des Nullmodells. 
#Dies ist erwartbar, da BLAST nicht alle zufälligen Paaralignments reportet, 
#sondern nur signifikante Treffer. Der Plot zeigt daher weniger einen 
#Goodness-of-fit des Nullmodells als vielmehr, dass die beobachteten TSR3-Hits 
#deutlich höhere Scores aufweisen als unter Zufall zu erwarten wäre.

# ==== 2.  Allgemiene Analyse für Wahl der Substitutionsmatrix ====
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

# ==== 3.  Entscheidung mittels (KL-Divergenz / Relative Entropie) ====
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
  # Wir betrachten die übergebenen Daten (z.B. Top 30 oder den gesamten Datensatz)
  # Berechnung: Bits per Residue = BitScore / Alignment Länge
  # Wenn die Matrix gut passt, sollte dieser Wert nahe an der theoretischen Entropie (H) liegen.
  blast_data$Bits_per_Residue <- blast_data$Bit_score / blast_data$Alignment_length
  
  # Durchschnitt über die Daten
  mean_obs <- mean(blast_data$Bits_per_Residue, na.rm = TRUE)
  
  return(data.frame(
    Matrix = matrix_name,
    Observed_Bits_per_Residue = round(mean_obs, 3)
  ))
}

# Daten aggregieren für Top 500
entropy_comparison <- rbind(
  calculate_observed_entropy(head(blast_data_30, 100), "PAM30"),
  calculate_observed_entropy(head(blast_data_70, 100), "PAM70"),
  calculate_observed_entropy(head(blast_data_90, 100), "BLOSUM90"),
  calculate_observed_entropy(head(blast_data_80, 100), "BLOSUM80"),
  calculate_observed_entropy(head(blast_data_62, 100), "BLOSUM62")
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
    title = "Matrix-Passform (Top ....): Theorie vs. Praxis",
    subtitle = "Vergleich der relativen Entropie (H) mit den tatsächlich erzielten Bits/Residue.\nIst Observed << Theoretical, ist die Matrix zu streng (viele Strafen).",
    y = "Bits per Residue",
    x = "Substitutionsmatrix",
    fill = "Metrik"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("H_theoretical" = "grey70", "Observed_Bits_per_Residue" = "dodgerblue3"),
                    labels = c("Theoretisches H (Max)", "Beobachtet (Top 30 Avg)"))


# ==== 4.  Validierung der Hintergrundfrequenzen (Compositional Check) ====

#' Berechnet die Hintergrundfrequenzen der Aminosäuren aus einem MSA oder einer Sequenz
calculate_background_frequencies <- function(x) {
  aa_list <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  if (inherits(x, "MsaAAMultipleAlignment") || inherits(x, "AAMultipleAlignment")) {
    raw_data <- as.matrix(unmasked(x))
  } else {
    raw_data <- as.character(x)
  }
  
  counts <- table(factor(raw_data[raw_data %in% aa_list], levels = aa_list))
  return(as.numeric(counts / sum(counts)))
}

# 1. Query-Sequenz einlesen
query_seq <- readAAStringSet("Data/Rohdaten/query_sequence.fasta")

# 2. Beobachtete Frequenzen (p_i) berechnen
as_freq <- calculate_background_frequencies(query_seq)

# 3. Standard-Frequenzen (Robinson-Robinson Modell, Grundlage vieler Matrizen)
std_freq <- c(A=0.0825, R=0.0553, N=0.0406, D=0.0545, C=0.0137, Q=0.0393, 
              E=0.0675, G=0.0707, H=0.0227, I=0.0596, L=0.0966, K=0.0584, 
              M=0.0242, F=0.0386, P=0.0470, S=0.0656, T=0.0534, W=0.0108, 
              Y=0.0292, V=0.0687)
# Korrektur kleiner Tippfehler in Namen (f -> F)
names(std_freq) <- toupper(names(std_freq))

# 4. Vergleich & Korrelation (Prüfung auf ungewöhnliche Komposition)
comp_check <- data.frame(AA = names(std_freq), Observed = as_freq, Standard = std_freq)
rho <- cor(comp_check$Observed, comp_check$Standard)

message("Die Korrelation der Aminosäure-Zusammensetzung beträgt: ", round(rho, 3))
comp_check_sorted <- comp_check %>% mutate(Diff = Observed - Standard) %>% arrange(desc(Diff))

# FAZIT: 
# PAM30 bietet das höchste theoretische Potential (H=2.57).
# Da unser TSR3 Protein hoch konserviert ist, erreichen wir auch empirisch sehr hohe Werte.
# Wir nutzen also PAM70, um die maximale Information aus den Alignments zu ziehen.

# ==== 5.  Speichern der Top30 mit gewählter Substitutionsmatrix (PAM70) ====
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


# ==== 6.  Validierung der Zielhäufigkeiten (q_ij Analyse) ====

#' Hilfsfunktion zum Einlesen von Matrix-Dateien (z.B. BLSOUM90.txt) mit KI erstellt
read_matrix_file <- function(file_path) {
  lines <- readLines(file_path)
  # Header-Zeile mit AA-Labels finden
  header_idx <- grep("^[ ]+A  B  C", lines)
  if (length(header_idx) == 0) header_idx <- grep("^[ ]+A  R  N", lines)
  
  labels <- strsplit(trimws(lines[header_idx]), "[ ]+")[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[grepl("^[A-Z*]", data_lines)]
  
  mat <- matrix(NA, nrow = length(data_lines), ncol = length(labels))
  for (i in seq_along(data_lines)) {
    nums <- strsplit(trimws(substr(data_lines[i], 2, nchar(data_lines[i]))), "[ ]+")[[1]]
    mat[i, ] <- as.numeric(nums)
  }
  rownames(mat) <- substr(data_lines, 1, 1)
  colnames(mat) <- labels
  return(mat)
}

#' Validiert die Wahl der Substitutionsmatrix durch Abgleich 
#' der theoretischen Zielhäufigkeiten mit den empirisch im MSA beobachteten Mutationen.
validate_substitution_matrix <- function(msa_obj, matrix_source = "PAM70") {
  
  message("Starte statistische Validierung für: ", 
          if(is.matrix(matrix_source)) "Manuelle Matrix" else matrix_source)
  
  msa_mat <- as.matrix(unmasked(msa_obj))
  aa_list <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  # Optimierte Zählung der Austausche (Frequenz-basiert statt paarweise Iteration)
  sub_counts <- matrix(0, nrow = 20, ncol = 20, dimnames = list(aa_list, aa_list))
  
  for (col in 1:ncol(msa_mat)) {
    residues <- msa_mat[, col]
    # Filterung auf valide Aminosäuren
    valid_residues <- residues[residues %in% aa_list]
    
    if (length(valid_residues) > 1) {
      # Frequenzen der AS in dieser Spalte
      counts <- table(factor(valid_residues, levels = aa_list))
      counts_vec <- as.numeric(counts)
      
      # Äußeres Produkt berechnet alle Kombinationen: counts[i] * counts[j]
      # Dies entspricht exakt der Anzahl der Paare (i, j)
      col_sub_matrix <- outer(counts_vec, counts_vec)
      
      # Korrektur der Diagonale: Ein Aminosäure-Zeichen kann nicht mit sich selbst 
      # an derselben Position gepaart werden (n * (n-1) statt n * n)
      diag(col_sub_matrix) <- counts_vec * (counts_vec - 1)
      
      sub_counts <- sub_counts + col_sub_matrix
    }
  }
  
  # Normalisierung zu Wahrscheinlichkeiten q_ij
  q_ij_obs <- sub_counts / sum(sub_counts)
  p_i_obs <- rowSums(q_ij_obs)
  
  # Berechnung der Log-Odds Verhältnisse (Observed Matrix)
  log_odds_obs <- matrix(NA, 20, 20, dimnames = list(aa_list, aa_list))
  for (i in aa_list) {
    for (j in aa_list) {
      if (q_ij_obs[i, j] > 0 && p_i_obs[i] > 0 && p_i_obs[j] > 0) {
        log_odds_obs[i, j] <- log2(q_ij_obs[i, j] / (p_i_obs[i] * p_i_obs[j]))
      }
    }
  }
  
  # Matrix-Quelle bestimmen
  if (is.matrix(matrix_source)) {
    sub_mat_theo <- matrix_source
  } else if (file.exists(matrix_source)) {
    sub_mat_theo <- read_matrix_file(matrix_source)
  } else {
    # Dynamisches Laden aus Biostrings/pwalign
    data(list = matrix_source, package = "Biostrings", envir = environment())
    sub_mat_theo <- get(matrix_source)
  }
  
  # Korrelation berechnen (nur untere Dreiecksmatrix inkl. Diagonale)
  common_aa <- intersect(aa_list, rownames(sub_mat_theo))
  obs_vector <- log_odds_obs[common_aa, common_aa][lower.tri(log_odds_obs[common_aa, common_aa], diag = TRUE)]
  theo_vector <- sub_mat_theo[common_aa, common_aa][lower.tri(sub_mat_theo[common_aa, common_aa], diag = TRUE)]
  
  valid_idx <- !is.na(obs_vector) & !is.na(theo_vector)
  df_comp <- data.frame(Observed = obs_vector[valid_idx], Theoretical = theo_vector[valid_idx])
  
  correlation <- cor(df_comp$Observed, df_comp$Theoretical)
  
  p <- ggplot(df_comp, aes(x = Theoretical, y = Observed)) +
    geom_point(alpha = 0.5, color = "dodgerblue4") +
    geom_smooth(method = "lm", color = "firebrick", linetype = "dashed", formula = y ~ x) +
    labs(title = paste("Validation:", if(is.matrix(matrix_source)) "Manual" else matrix_source),
         subtitle = paste("Pearson R =", round(correlation, 4)),
         x = "Theoretical Log-Odds (Matrix)",
         y = "Empirical Log-Odds (MSA)") +
    theme_minimal()
  
  return(list(correlation = correlation, plot = p, q_ij = q_ij_obs))
}

# ==== 7.  Kompakter Cross-Check Workflow (BLOSUM90 Beispiel) ====

library(msa)
# 1. MSA für BLOSUM90 on-the-fly erstellen
# Top 30 IDs aus BLOSUM90 BLAST Ergebnissen
top30_b90_ids <- head(blast_data_90$Accession_ID, 30)

# Sequenzen aus DB laden
db_seqs <- readAAStringSet("Data/Rohdaten/database.fasta")
names(db_seqs) <- word(names(db_seqs), 1)
top30_b90_seqs <- db_seqs[top30_b90_ids]

# Schnelles MSA (nur im Arbeitsspeicher)
msa_b90 <- msa(top30_b90_seqs, method = "ClustalW", substitutionMatrix = "blosum")
msa_pam70 <- readRDS("Results_MultibleSequenceAlign/msa_result.rds")

# 2. Kreuz-Validierung
# MSA(B90) gegen Matrix(B90)
val_b90_self <- validate_substitution_matrix(msa_b90, "Data/BLOSUM90.txt")

# MSA(B90) gegen Matrix(PAM70)
val_b90_vs_pam <- validate_substitution_matrix(msa_b90, "PAM70")

print(val_b90_self$plot)
print(val_b90_vs_pam$plot)

message("Ergebnis des Cross-Checks (MSA auf Basis von BLOSUM90):",
        "\nKorrelation mit BLOSUM90 Matrix: ", round(val_b90_self$correlation, 4),
        "\nKorrelation mit PAM70 Matrix:    ", round(val_b90_vs_pam$correlation, 4))


# ==== 8.  Globale Analyse (Gesamter Datensatz) ====

message("Starte globale Analyse für alle BLAST-Treffer (~500 Sequenzen)...")

# 1. Alle Treffer annotieren (E-Value <= 1e-10 für biologische Relevanz)
all_hits_annotated <- annotate_blast_results(
  blast_results = blast_data_70,
  e_value_threshold = 1e-10, 
  create_csv = FALSE
)

# 2. Globale Entropie-Analyse
global_entropy <- calculate_observed_entropy(all_hits_annotated, "PAM70 (Global)")
message("Globale beobachtete Entropie: ", global_entropy$Observed_Bits_per_Residue)

# 3. Sequenzen für den gesamten Datensatz extrahieren
all_ids <- all_hits_annotated$Accession_ID
db_all <- readAAStringSet("Data/Rohdaten/database.fasta")
names(db_all) <- word(names(db_all), 1)
all_seqs <- db_all[all_ids]

# 4. Globales MSA (Achtung: Rechenzeit ca. 1-3 Minuten)
msa_all <- msa(all_seqs, method = "ClustalW", substitutionMatrix = "blosum")

# 5. Globale Validierung der Zielhäufigkeiten
val_global_70 <- validate_substitution_matrix(msa_all, "PAM70")
val_global_90 <- validate_substitution_matrix(msa_all, "Data/BLOSUM90.txt")
print(val_global_70$plot)
print(val_global_90$plot)

# Weitere
val_global_80 <- validate_substitution_matrix(msa_all, "BLOSUM80")
val_global_62 <- validate_substitution_matrix(msa_all, "BLOSUM62")
val_global_30 <- validate_substitution_matrix(msa_all, "PAM30")

# 6. Globale Hintergrundfrequenzen (Validierung des Protein-Pools)
global_p_i <- calculate_background_frequencies(msa_all)
comp_check_global <- data.frame(AA = names(std_freq), Global_Observed = global_p_i, Standard = std_freq)
rho_global <- cor(comp_check_global$Global_Observed, comp_check_global$Standard)
message("Korrelation (Globaler Pool vs Std): ", round(rho_global, 3))


message(
  "Ergebnis des Cross-Checks (MSA auf Basis von BLOSUM90):",
  "\nKorrelation mit PAM70 Matrix: ", round(val_global_70$correlation, 4),
  "\nKorrelation mit PAM30 Matrix:    ", round(val_global_30$correlation, 4),
  "\nKorrelation mit BLOSUM90 Matrix: ", round(val_global_90$correlation, 4),
  "\nKorrelation mit BLOSUM80 Matrix:    ", round(val_global_80$correlation, 4),
  "\nKorrelation mit BLOSUM62 Matrix:    ", round(val_global_62$correlation, 4)
)


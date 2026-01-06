# Library and Settings ----
rm(list = ls())
library(ggplot2)
# library("tidyverse")
# library(Biostrings)
library(seqinr)
library(ordinal)
library(dplyr)
library(stringr)
library(readr)


# Select Working Directory ----
wds <- data.frame(
  system = c("Julian Windows PC/Laptop", 
             "Linux Uni",
             "Johann ...",
             "Rphi ..."),
  path = c("C:/Users/julia/OneDrive/Integrated Life Siences/Genomanalysen und Phylogenie/Projektarbeit",
           "/home/stud/ha24vepa/Documents/Bash_Linux_Introduction_supplements/Genom_und_Phylogenie_kurs/...",
           "Johann:Copy paste dein Working directory",
           "Rphi: Copy paste dein Working directory"),
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



# ==== Analyse zur Wahl der Substitutionsmatrix (BLOSUM62 vs BLOSUM80 vs PAM30) ====

# 1. Daten Import der anderen Blast Ergebnisse




# 2. Berechnung der durchschnittlichen Identität (basierend auf BLOSUM62 Standard)
# Wir betrachten die Top-Treffer (ähnlich zu unserer Auswahl für die Phylogenie),
# da diese die Signale sind, die wir optimal auflösen wollen.
mean_identity_all <- mean(blast_data_62$Percent_identity)
target_identity <- mean(head(blast_data_62$Percent_identity, 30))

message("Durchschnittliche Identität (Alle Hits BLOSUM62): ", round(mean_identity_all, 2), "%")
message("Durchschnittliche Identität (Top 30 Hits BLOSUM62): ", round(mean_identity_top30, 2), "%")

# 3. Entscheidungshilfe / Logik
# Theorie:
# BLOSUM62: ~62% Identität (Allrounder)
# BLOSUM80: ~80% Identität (für nähere Verwandte)
# PAM30:    <30 PAM (~70-90% Identität, für sehr kurze evolutionäre Distanzen)

target_identity


# 4. Vergleich der Bit-Scores für den Top-Hit (als Validierung)
# Wenn die Matrix besser zum Alignment passt, sollte der Bit-Score (Information content) steigen.
# Annahme: Der erste Hit ist der gleiche (z.B. Homo sapiens Query self-hit oder Top Ortholog)

score_62 <- head(blast_data_62$Bit_score, 1)
score_80 <- head(blast_data_80$Bit_score, 1)
score_30 <- head(blast_data_30$Bit_score, 1)
message("BLOSUM62: ", score_62)
message("BLOSUM80: ", score_80)
message("PAM30:    ", score_30)
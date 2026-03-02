# Library and Settings ----
rm(list = ls())
library(ggplot2)
library(Biostrings)
library(msa)
library(readr)
library(stringr)
library(ggmsa)
library(dplyr) 
library(stringr)


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
# ==== 1. Sequenzen & Metadaten Laden & Vorbereiten ====

seqs <- readAAStringSet("Results_MultibleSequenceAlign/top30_sequences.fasta")
meta <- read_csv("Results_BLAST/top30_pam70_annotated.csv", show_col_types = FALSE)

# Robust renaming: Ensure the order of metadata matches the order of sequences
# Create the new labels in the metadata
meta <- meta |>
  mutate(label = str_replace_all(Organism, " ", "_"),
         label = paste(label, Accession_ID, sep = "_"))

# Match the metadata to the sequences based on Accession_ID (which are the current names of seqs)
match_indices <- match(names(seqs), meta$Accession_ID)

# Assign new names only if a match is found (should be for all)
if(any(is.na(match_indices))) {
  warning("Some sequences could not be matched to metadata!")
}
names(seqs) <- meta$label[match_indices]

# ==== 2. Multiple Sequence Alignment (MSA) durchführen ====

message("Starte Multiple Sequence Alignment mit ClustalW...")

msa_res <- msa(seqs, method = "ClustalW")
# Kurze Vorschau in der Konsole
# msa_res

# ==== 3. Ergebnisse Exportieren & Speichern ====

writeXStringSet(unmasked(msa_res),
                "Results_MultibleSequenceAlign/msa_alignment.fasta")
saveRDS(msa_res,
        "Results_MultibleSequenceAlign/msa_result.rds")

# ==== 4. Variablen für Visualisierung ====

# Wir extrahieren die Sequenzen in ein neues Objekt
plot_seqs <- unmasked(msa_res)

# Wir gleichen die aktuellen, robusten Namen mit den Metadaten ab,
# um den reinen Organismus-Namen (ohne Accession ID) abzurufen
match_plot <- match(names(plot_seqs), meta$label)

# Reine Organismus-Namen laden und Leerzeichen durch Unterstriche ersetzen
clean_org_names <- str_replace_all(meta$Organism[match_plot], " ", "_")

# Zwingend für ggmsa: Namen eindeutig machen (hängt .1, .2 an Duplikate an)
names(plot_seqs) <- make.unique(clean_org_names)

# ==== 5. Visualisierung des Alignments ====
msa_plot <- ggmsa(plot_seqs,
                  font = NULL, 
                  color = "Chemistry_AA") +
  theme(axis.text.y = element_text(size = 8))
msa_plot

# Hochauflösender Export für die Präsentation
  # ggsave("Results_MultibleSequenceAlign/MSA_Visualization_HighRes.png", 
       # plot = msa_plot, 
       # width = 50,  
       # height = 8,  
       # dpi = 300, 
       # limitsize = FALSE)

message("Hochauflösendes Bild wurde gespeichert!")

# ==== 6a. Quantitative Konservierungs-Analyse & Statistik ====

message("Berechne quantitative Konservierungs-Metriken...")

# zur Berechnung der Spaltenstatistik
calculate_msa_stats <- function(msa_obj) {
  # In Matrix konvertieren für effizienten Zugriff
  msa_matrix <- as.matrix(unmasked(msa_obj))
  n_seq <- nrow(msa_matrix)
  aln_width <- ncol(msa_matrix)

  # Statistiken pro Spalte berechnen
  stats_list <- lapply(seq_len(aln_width), function(i) {
    col_data <- msa_matrix[, i]
    gaps <- sum(col_data == "-")

    # Nicht-Gap Zeichen für Konservierung (Mehrheitsprinzip)
    non_gaps <- col_data[col_data != "-"]

    # Konservierung: Wie viel Prozent der Sequenzen tragen die häufigste Aminosäure?
    conservation <- if (length(non_gaps) == 0) 0 else max(table(non_gaps)) / n_seq

    data.frame(
      position = i, 
      conservation = conservation, 
      gap_freq = gaps / n_seq
    )
  })

  return(do.call(rbind, stats_list))
}

msa_stats <- calculate_msa_stats(msa_res)

# Grundlegende deskriptive Statistik für die Interpretation
mean_cons <- mean(msa_stats$conservation)
mean_gaps <- mean(msa_stats$gap_freq)
message(sprintf("Analyse abgeschlossen: Mittlere Konservierung %.2f%%, Mittlere Gap-Rate %.2f%%", 
                mean_cons * 100, mean_gaps * 100))

# ==== 6b. Domänen-Mapping (UniProt -> MSA Koordinaten) ====

message("Mappe UniProt-Domänen der humanen Sequenz auf Alignment-Positionen...")

# 1. Humane Sequenz im Alignment finden (Sicherstellen, dass der Name exakt stimmt)
# Hier wird angenommen, dass in "names(plot_seqs)" der String "Homo_sapiens" vorkommt.
human_idx <- grep("Homo_sapiens", names(plot_seqs))
if(length(human_idx) == 0) {
  warning("Humane Sequenz nicht gefunden! Bitte Namensgebung prüfen.")
}

human_aln_seq <- as.character(plot_seqs[[human_idx[1]]])

# 2. Mapping-Funktion: Sequenz-Position (ungapped) -> Alignment-Position (gapped)
map_seq_to_aln <- function(aln_seq, seq_positions) {
  chars <- strsplit(aln_seq, "")[[1]]
  is_aa <- chars != "-"             # Logischer Vektor: TRUE für Aminosäuren, FALSE für Gaps
  seq_pos_mapping <- cumsum(is_aa)  # Kumulative Summe gibt die AS-Position an
  
  # Finde die Alignment-Spalte für die gesuchten AS-Positionen
  aln_positions <- sapply(seq_positions, function(p) match(p, seq_pos_mapping))
  return(aln_positions)
}

# UniProt Koordinaten für TSR3
uniprot_coords <- c(
  dom1_start = 96, dom1_end = 222,  # Ribosome biogenesis protein, C-terminal
  dom2_start = 224, dom2_end = 312  # Disordered Region
)

# Mappen
aln_coords <- map_seq_to_aln(human_aln_seq, uniprot_coords)

message(sprintf("Domäne 1 im Alignment: Spalte %d bis %d", aln_coords["dom1_start"], aln_coords["dom1_end"]))
message(sprintf("Domäne 2 im Alignment: Spalte %d bis %d", aln_coords["dom2_start"], aln_coords["dom2_end"]))

# ==== 6c. Bereichsspezifische Statistiken ====

# Statistik für Domäne 1 (Strukturiert)
dom1_stats <- msa_stats |>
  filter(position >= aln_coords["dom1_start"] & position <= aln_coords["dom1_end"]) |>
  summarise(mean_cons = mean(conservation), mean_gaps = mean(gap_freq))

# Statistik für Domäne 2 (Disordered)
dom2_stats <- msa_stats |>
  filter(position >= aln_coords["dom2_start"] & position <= aln_coords["dom2_end"]) |>
  summarise(mean_cons = mean(conservation), mean_gaps = mean(gap_freq))

message(sprintf("Auswertung Domäne 1 (Strukturiert): Konservierung %.2f%%, Gaps %.2f%%", 
                dom1_stats$mean_cons * 100, dom1_stats$mean_gaps * 100))
message(sprintf("Auswertung Domäne 2 (Disordered): Konservierung %.2f%%, Gaps %.2f%%", 
                dom2_stats$mean_cons * 100, dom2_stats$mean_gaps * 100))
# ==== 7. Visualisierung der Statistiken (mit Domänen-Markierung) ====

message("Erstelle erweitertes Statistik-Diagramm...")

stats_plot <- ggplot(msa_stats, aes(x = position)) +
  
  # Domänen-Hintergründe einzeichnen (vor den Balken, damit sie im Hintergrund liegen)
  annotate("rect", xmin = aln_coords["dom1_start"], xmax = aln_coords["dom1_end"], 
           ymin = -1, ymax = 1, alpha = 0.15, fill = "green") +
  annotate("rect", xmin = aln_coords["dom2_start"], xmax = aln_coords["dom2_end"], 
           ymin = -1, ymax = 1, alpha = 0.15, fill = "orange") +
  
  # Text-Labels für die Domänen hinzufügen
  annotate("text", x = mean(c(aln_coords["dom1_start"], aln_coords["dom1_end"])), 
           y = 1.05, label = "Ribosome Biogenesis Dom.", size = 3, fontface = "bold", color = "darkgreen") +
  annotate("text", x = mean(c(aln_coords["dom2_start"], aln_coords["dom2_end"])), 
           y = 1.05, label = "Disordered Region", size = 3, fontface = "bold", color = "darkorange") +
  
  # Konservierung als Balken (Blau)
  geom_bar(aes(y = conservation, fill = "Konservierung"), stat = "identity", alpha = 0.8) +
  # Gaps als "hängende" Balken (Rot) zur visuellen Trennung
  geom_bar(aes(y = -gap_freq, fill = "Gap-Frequenz"), stat = "identity", alpha = 0.8) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Konservierung" = "#2c7bb6", "Gap-Frequenz" = "#d7191c")) +
  scale_y_continuous(labels = abs, limits = c(-1, 1.1)) + # Y-Achse leicht erhöht für Labels
  labs(title = "MSA Konservierungs- und Gap-Profil von TSR3",
       subtitle = sprintf("Gesamt-Konservierung: %.1f%% | Domäne 1: %.1f%% | Domäne 2: %.1f%%", 
                          mean_cons*100, dom1_stats$mean_cons*100, dom2_stats$mean_cons*100),
       x = "Position im Alignment",
       y = "Frequenz (0.0 - 1.0)",
       fill = "Metrik") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

stats_plot
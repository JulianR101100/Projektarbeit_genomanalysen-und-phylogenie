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
# ==== 1. Sequenzen & Metadaten Laden & Vorbereiten ====

seqs <- readAAStringSet("Results_MultibleSequenceAlign/top30_sequences.fasta")
meta <- read_csv("Results_BLAST/top30_pam70_annotated.csv", show_col_types = FALSE)

# Probelm: IDs wichtig für fehlerfreies MSA aber auch später hinderlich für Lesbarkeit
# Umbenennung: Metadata -> Match order of Sequences
# Create new labels in the metadata
meta <- meta |>
  mutate(label = str_replace_all(Organism, " ", "_"),
         label = paste(label, Accession_ID, sep = "_"))

# Match the metadata to the sequences based on Accession_ID (which are the current names of seqs)
match_indices <- match(names(seqs), meta$Accession_ID)

# Assign new names only if a match is found (should be for all -> Assure warning)
if(any(is.na(match_indices))) {
  warning("Some sequences could not be matched to metadata!")
}
names(seqs) <- meta$label[match_indices]

# ==== 2. Multiple Sequence Alignment (MSA) durchführen ====

message("Starte Multiple Sequence Alignment mit ClustalW...")

msa_res <- msa(seqs, method = "ClustalW", substitutionMatrix = "pam")
# Kurze Vorschau in der Konsole
# msa_res

# ==== 3. Ergebnisse Exportieren & Speichern ====
# rds für R; fasta einfach so
writeXStringSet(unmasked(msa_res),
                "Results_MultibleSequenceAlign/msa_alignment.fasta")
saveRDS(msa_res,
        "Results_MultibleSequenceAlign/msa_result.rds")

# ==== 4. Variablen für Visualisierung ====

# extrahieren Sequenzen in ein neues Objekt
plot_seqs <- unmasked(msa_res)

# Wir gleichen die aktuellen, robusten Namen mit den Metadaten ab,
# um den reinen Organismus-Namen (ohne Accession ID) abzurufen
match_plot <- match(names(plot_seqs), meta$label)

# Reine Organismus-Namen laden und: Leerzeichen -> Unterstriche
clean_org_names <- str_replace_all(meta$Organism[match_plot], " ", "_")

# Namen eindeutig machen (hängt .1, .2 an Duplikate an), ggmsa funktioniert sonst nicht
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
  # -> Als Matrix effizienter
  msa_matrix <- as.matrix(unmasked(msa_obj))
  n_seq <- nrow(msa_matrix)
  aln_width <- ncol(msa_matrix)

  # Statistiken pro Spalte
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

# Mit KI unterstützung geschrieben

# 1. Humane Sequenz im Alignment finden
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

# UniProt Koordinaten für TSR3 (jetzt mit beiden Disordered Regions)
uniprot_coords <- c(
  ribo_start = 96,  ribo_end = 222,   # Ribosome biogenesis protein, C-terminal
  dis1_start = 1,   dis1_end = 21,    # Disordered Region (N-Terminus)
  dis2_start = 224, dis2_end = 312    # Disordered Region (C-Terminus)
)

# Mappen auf Alignment-Koordinaten
aln_coords <- map_seq_to_aln(human_aln_seq, uniprot_coords)

message(sprintf("Ribosom-Domäne im Alignment: Spalte %d bis %d", aln_coords["ribo_start"], aln_coords["ribo_end"]))
message(sprintf("Disordered Region 1: Spalte %d bis %d", aln_coords["dis1_start"], aln_coords["dis1_end"]))
message(sprintf("Disordered Region 2: Spalte %d bis %d", aln_coords["dis2_start"], aln_coords["dis2_end"]))

# ==== 6c. Bereichsspezifische Statistiken ====

# Statistik für die strukturierte Ribosom-Domäne
ribo_stats <- msa_stats |>
  filter(position >= aln_coords["ribo_start"] & position <= aln_coords["ribo_end"]) |>
  summarise(mean_cons = mean(conservation), mean_gaps = mean(gap_freq))

# Statistik für BEIDE Disordered Regions zusammen
disordered_stats <- msa_stats |>
  filter((position >= aln_coords["dis1_start"] & position <= aln_coords["dis1_end"]) |
         (position >= aln_coords["dis2_start"] & position <= aln_coords["dis2_end"])) |>
  summarise(mean_cons = mean(conservation), mean_gaps = mean(gap_freq))

message(sprintf("Auswertung Ribosom-Domäne: Konservierung %.2f%%, Gaps %.2f%%", 
                ribo_stats$mean_cons * 100, ribo_stats$mean_gaps * 100))
message(sprintf("Auswertung Disordered Regions (kombiniert): Konservierung %.2f%%, Gaps %.2f%%", 
                disordered_stats$mean_cons * 100, disordered_stats$mean_gaps * 100))

# ==== 7. Visualisierung der Statistiken (mit Domänen-Markierung) ====

message("Erstelle erweitertes Statistik-Diagramm...")

stats_plot <- ggplot(msa_stats, aes(x = position)) +

  # Hintergrund-Markierung: Ribosom-Domäne (Kräftiges Grün)
  annotate("rect", xmin = aln_coords["ribo_start"], xmax = aln_coords["ribo_end"], 
           ymin = -1, ymax = 1, alpha = 0.20, fill = "#4daf4a") +

  # Hintergrund-Markierung: Disordered Regions (Kräftiges Orange: #f37a00)
  annotate("rect", xmin = aln_coords["dis1_start"], xmax = aln_coords["dis1_end"], 
           ymin = -1, ymax = 1, alpha = 0.20, fill = "darkorange") +
  annotate("rect", xmin = aln_coords["dis2_start"], xmax = aln_coords["dis2_end"], 
           ymin = -1, ymax = 1, alpha = 0.20, fill = "darkorange") +

  # Text-Labels für die Bereiche
  annotate("text", x = mean(c(aln_coords["ribo_start"], aln_coords["ribo_end"])), 
           y = 1.05, label = "Ribosome Biogenesis Dom.", size = 3, fontface = "bold", color = "#2e7d32") +
  # Label mittig über der größeren Disordered Region
  annotate("text", x = mean(c(aln_coords["dis2_start"], aln_coords["dis2_end"])), 
           y = 1.05, label = "Disordered Region", size = 3, fontface = "bold", color = "darkorange") +

  # Konservierung als Balken (Blau)
  geom_bar(aes(y = conservation, fill = "Konservierung"), stat = "identity", alpha = 0.8) +
  # Gaps als "hängende" Balken (Rot) zur visuellen Trennung
  geom_bar(aes(y = -gap_freq, fill = "Gap-Frequenz"), stat = "identity", alpha = 0.8) +

  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Konservierung" = "#2c7bb6", "Gap-Frequenz" = "#d7191c")) +
  scale_y_continuous(labels = abs, limits = c(-1, 1.1)) + # Y-Achse leicht erhöht für Labels
  labs(title = "MSA Konservierungs- und Gap-Profil von TSR3",
       subtitle = sprintf("Gesamt-Konservierung: %.1f%% | Ribosom-Domäne: %.1f%% | Disordered (Gesamt): %.1f%%", 
                          mean_cons*100, ribo_stats$mean_cons*100, disordered_stats$mean_cons*100),
       x = "Position im Alignment",
       y = "Frequenz (0.0 - 1.0)",
       fill = "Metrik") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

stats_plot

ggsave(
  filename = "Alignment_PLOT+Regions.pdf", 
  plot = stats_plot,
  width = 7,
  height = 5,
  device = "pdf"
)

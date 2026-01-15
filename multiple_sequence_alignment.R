# Library and Settings ----
rm(list = ls())
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
meta <- read_csv("Results_BLAST/top30_pam30_annotated.csv", show_col_types = FALSE)

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
msa_res

# ==== 3. Ergebnisse Exportieren & Speichern ====

writeXStringSet(unmasked(msa_res),
                "Results_MultibleSequenceAlign/msa_alignment.fasta")
saveRDS(msa_res,
        "Results_MultibleSequenceAlign/msa_result.rds")

# ==== 4. Visualisierung des Alignments ====

# ggmsa ist nicht so easy man muss das iwi konvertieren, frag nicht wieso,
# Die ergebnisse daraus sind eh hesslich also würd ich schauen ob es ein extenres 
# tol gibt mit dem wir das schön machen

ggmsa(unmasked(msa_res),
      start = 1,
      end = 100,
      font = NULL, # Font NULL beschleunigt das Rendering bei vielen Sequenzen sieht aber besser aus ohne
      color = "Chemistry_AA") +
  geom_seqlogo() +
  theme(axis.text.y = element_text(size = 8)) # Beschriftung etwas kleiner


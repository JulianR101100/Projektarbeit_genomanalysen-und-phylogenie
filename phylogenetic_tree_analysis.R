# Library and Settings ----
rm(list = ls())
library(Biostrings) 
library(seqinr)
library(pwalign)
library(gplots) 
library(ape)
library(dendextend)
library(phytools)
library(msa)
library(readr)
library(stringr)
library(colorspace)

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
# ==== 1. Distanzberechnungen ====

# MSA laden und konvertieren
alignment_biostrings <- readAAMultipleAlignment("Results_MultibleSequenceAlign/msa_alignment.fasta", format = "fasta")
alignment_align <- msaConvert(alignment_biostrings, type = "seqinr::alignment")

# Distanzmatrix (Identität) berechnen
dist_msa <- dist.alignment(alignment_align, matrix = "identity")
dist_matrix <- as.matrix(dist_msa)
dist_matrix_sq <- dist_matrix^2

# Heatmap speichern
pdf("Results_PhylogeneticTree/Distance_Heatmap.pdf", width = 12, height = 12)
heatmap.2(dist_matrix_sq,
          trace = "none",
          col = rev(hcl.colors(100, palette = "heat")),
          margins = c(12, 12),
          cexRow = 0.5,
          cexCol = 0.5,
          dendrogram = "none",
          key = TRUE,
          denscol="black",
          densadj = 0.25,
          main = "Distanzmatrix Heatmap (Identität^2)")
dev.off()
message("Distanz-Heatmap gespeichert: Results_PhylogeneticTree/Distance_Heatmap.pdf")

# ==== 2. Ultrametrie und Hierarchisches Clustering ====

message("Prüfe Ultrametrie-Eigenschaft...")

# Funktion: Strikte Ultrametrie prüfen
check_strict_ultrametric_property <- function(dist_matrix, tol = 1e-12) {
  points <- rownames(dist_matrix)
  combinations <- combn(points, 3)
  
  tests <- 0
  violations <- 0
  violating_triplets <- character(0)
  
  for (i in 1:ncol(combinations)) {
    comb <- combinations[, i]
    A <- comb[1]; B <- comb[2]; C <- comb[3]
    
    d_AB <- dist_matrix[A, B]
    d_AC <- dist_matrix[A, C]
    d_BC <- dist_matrix[B, C]
    
    tests <- tests + 1
    d <- sort(c(d_AB, d_AC, d_BC))
    d1 <- d[1]; d2 <- d[2]; d3 <- d[3]
    
    # Kriterium: Die zwei größten Distanzen müssen gleich sein
    is_strict_ultra <- (abs(d3 - d2) <= tol) && (d2 > d1 + 1e-6)
    # d2 > d1 + tol liefert bei Wert 0,01 weniger Verletzungen als bei 0,1 
    #macht mathematisch keinen Sinn macht da tol einmal Strenger wird (d3-d3)<=tol
    #und einmnal entspannter bei d2>d1+tol (ist wiederspruch)
    #Toleranz sollte eher vergleich der distanzen Beeinflussen weniger ob d2>d1 ist 
    #(diese Bedingung muss egal wie für strickte Ultrametrie erfüllt sein)
    
    if (!is_strict_ultra) {
      violations <- violations + 1L
      violating_triplets <- c(violating_triplets, paste(A, B, C, sep = "-"))
    }
  }
  return(list(
    tests = tests,
    violations = violations,
    strict_ultrametric = (violations == 0L),
    violating_triplets = violating_triplets))
}



# ==== 2.2. Ultrametrie und Hierarchisches Clustering ====

check_ultrametric_property <- function(dist_matrix, tol = 1e-12) {
  points <- rownames(dist_matrix)
  combinations <- combn(points, 3)
  
  tests <- 0L
  violations <- 0L
  violating_triplets <- character(0)
  
  for (i in 1:ncol(combinations)) {
    comb <- combinations[, i]
    A <- comb[1]; B <- comb[2]; C <- comb[3]
    
    # Extract distances
    d_AB <- dist_matrix[A, B]
    d_AC <- dist_matrix[A, C]
    d_BC <- dist_matrix[B, C]
    
    tests <- tests + 1L
    
    # Ultrametric inequalities (all must hold)
    ok <- (d_AB <= max(d_AC, d_BC) + tol) &&
      (d_AC <= max(d_AB, d_BC) + tol) &&
      (d_BC <= max(d_AB, d_AC) + tol)
    
    if (!ok) {
      violations <- violations + 1L
      violating_triplets <- c(violating_triplets, paste(A, B, C, sep = "-"))
    }
  }
  
  return(list(
    tests = tests,
    violations = violations,
    ultrametric = (violations == 0L),
    violating_triplets = violating_triplets))
}


# ======== 2.3. Ultrametric checks for multiple tolerances ========

dist_msa_sq <- as.dist(dist_matrix_sq)

tolerances <- c(0.1, 0.03, 0.02, 0.01, 0.005, 0.001)
out_file <- "Results_PhylogeneticTree/Ultrametric_Check.txt"
n_examples <- 5   # number of violating triplets to display

# Datei neu anlegen / überschreiben
cat("Ultrametric property checks for multiple tolerances\n",
    "===================================================\n\n",
    file = out_file)

for (tol in tolerances) {
  
  strict_res <- check_strict_ultrametric_property(dist_matrix_sq, tol = tol)
  nonstrict_res <- check_ultrametric_property(dist_matrix_sq, tol = tol)
  
  # Ergebnisse in Datei schreiben
  cat("Tolerance:", tol, "\n",
      "------------------------------\n",
      file = out_file, append = TRUE)
  
  cat("Strict ultrametric test:\n",
      "Tests: ", strict_res$tests, "\n",
      "Violations: ", strict_res$violations, "\n",
      "Strict ultrametric: ", strict_res$strict_ultrametric, "\n",
      file = out_file, append = TRUE, sep = "")
  
  # if (length(strict_res$violating_triplets) > 0) {
  #   n_show <- min(n_examples, length(strict_res$violating_triplets))
  #   cat("First violating triplets (strict):\n",
  #       paste(strict_res$violating_triplets[1:n_show], collapse = "\n"),
  #       "\n",
  #       file = out_file, append = TRUE)
  # } else {
  #   cat("First violating triplets (strict): none\n",
  #       file = out_file, append = TRUE)
  # }
  
  cat("\nNon-strict ultrametric test:\n",
      "Tests: ", nonstrict_res$tests, "\n",
      "Violations: ", nonstrict_res$violations, "\n",
      "Non-strict ultrametric: ", nonstrict_res$ultrametric, "\n",
      file = out_file, append = TRUE, sep = "")
  
  # if (length(nonstrict_res$violating_triplets) > 0) {
  #   n_show <- min(n_examples, length(nonstrict_res$violating_triplets))
  #   cat("First violating triplets (non-strict):\n",
  #       paste(nonstrict_res$violating_triplets[1:n_show], collapse = "\n"),
  #       "\n",
  #       file = out_file, append = TRUE)
  # } else {
  #   cat("First violating triplets (non-strict): none\n",
  #       file = out_file, append = TRUE)
  # }
  
  cat("\n\n", file = out_file, append = TRUE)
  
  message("Tolerance ", tol,
          ": strict violations = ", strict_res$violations,
          ", non-strict violations = ", nonstrict_res$violations)
}

#Auskommentierte Bereiche für Angabe der Violating Triplets
message("Ultrametric checks saved to: ", out_file)

# ======== 2.4. Hierarchisches Clustering (UPGMA, WPGMA, Ward)========
hc_upgma <- hclust(dist_msa_sq, method = "average")
hc_wpgma <- hclust(dist_msa_sq, method = "mcquitty")
hc_ward <- hclust(dist_msa_sq, method = "ward.D2")

# ==== 3. Phylogenetische Bäume aus hierarchischem Clustering ====

phylo_tree_upgma <- as.phylo(hc_upgma)
phylo_tree_wpgma <- as.phylo(hc_wpgma)
phylo_tree_ward <- as.phylo(hc_ward)

# Bäume plotten und speichern
pdf("Results_PhylogeneticTree/Phylogenetic_Trees.pdf", width = 10, height = 8)

plot(phylo_tree_upgma, cex = 0.6, main = "UPGMA phylogenetischer Baum")
plot(phylo_tree_wpgma, cex = 0.6, main = "WPGMA phylogenetischer Baum")
plot(phylo_tree_ward, cex = 0.6, main = "Ward.D2 phylogenetischer Baum")

dev.off()
message("Phylogenetische Bäume gespeichert: Results_PhylogeneticTree/Phylogenetic_Trees.pdf")

# Tanglegrams (Baumvergleiche) speichern
pdf("Results_PhylogeneticTree/Tree_Comparisons.pdf", width = 12, height = 8)

# UPGMA vs WPGMA
phylo_list_upgma_wpgma <- dendlist(as.dendrogram(hc_upgma), as.dendrogram(hc_wpgma))
phylo_list_upgma_wpgma %>%
  dendextend::untangle() %>%
  tanglegram(main_left = "UPGMA", main_right = "WPGMA",
             lwd = 2, edge.lwd = 2, common_subtrees_color_branches = TRUE,
             columns_width = c(5, 1, 5), margin_inner = 12, lab.cex = 0.6)

# UPGMA vs Ward
phylo_list_upgma_ward <- dendlist(as.dendrogram(phylo_tree_upgma), as.dendrogram(phylo_tree_ward))
phylo_list_upgma_ward %>%
  dendextend::untangle() %>%
  tanglegram(main_left = "UPGMA", main_right = "Ward",
             lwd = 2, edge.lwd = 2, common_subtrees_color_branches = TRUE,
             columns_width = c(5, 1, 5), margin_inner = 12, lab.cex = 0.6)

dev.off()
message("Tanglegrams gespeichert: Results_PhylogeneticTree/Tree_Comparisons.pdf")


# ==== 4. Robustheit ====

# Funktion: Robustheit berechnen
calculate_robustness <- function(tree, dist_matrix, rho) {
  cophenetic_distances <- cophenetic(tree)
  max_diff <- max(abs(cophenetic_distances - dist_matrix))
  
  num_tips <- length(tree$tip.label)
  num_nodes <- tree$Nnode
  internal_nodes <- (num_tips + 1):(num_tips + num_nodes)
  internal_edges <- which(tree$edge[, 2] %in% internal_nodes)
  internal_edge_lengths <- tree$edge.length[internal_edges]
  
  if(length(internal_edge_lengths) > 0) {
      min_edge_length <- min(internal_edge_lengths)
      robustness_truth <- max_diff <= rho * min_edge_length
      est_rho <- max_diff/min_edge_length
      return(list(truth=robustness_truth, est_rho=est_rho, max_dist=max_diff, w0=min_edge_length))
  } else {
      return(list(truth=NA, est_rho=NA, max_dist=max_diff, w0=NA))
  }
}

# Robustheit berechnen (UPGMA)
robustness_upgma <- calculate_robustness(phylo_tree_upgma, dist_matrix_sq, 1)

# Ergebnisse anhängen
cat("\nRobustheit des UPGMA Baums:\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Wahrheit (Truth):", robustness_upgma$truth, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Geschätztes Rho (Est Rho):", robustness_upgma$est_rho, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)

message("Robustheitsprüfung UPGMA abgeschlossen: Results_PhylogeneticTree/Ultrametric_Check.txt")

#Robustgeit berechnen (WPGMA)
robustness_wpgma <- calculate_robustness(phylo_tree_wpgma, dist_matrix_sq, 1)

#Ergebnisse anhängen
cat("\nRobustheit des WPGMA Baums:\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Wahrheit (Truth):", robustness_wpgma$truth, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Geschätztes Rho (Est Rho):", robustness_wpgma$est_rho, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)

message("Robustheitsprüfung WPGMA abgeschlossen: Results_PhylogeneticTree/Ultrametric_Check.txt")

#Robustheit berechnen (WARD)
robustness_ward <- calculate_robustness(phylo_tree_ward, dist_matrix_sq, 1)

#Ergebnisse anhängen
cat("\nRobustheit des WARD Baums:\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Wahrheit (Truth):", robustness_ward$truth, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Geschätztes Rho (Est Rho):", robustness_ward$est_rho, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)

message("Robustheitsprüfung WARD abgeschlossen: Results_PhylogeneticTree/Ultrametric_Check.txt")


#========  5. Four-Point-Condition =======
#Funktion Four-Point-Condition 

check_four_point_condition <- function(dist_matrix, tol=1e-12) {
  points <- rownames(dist_matrix)
  combinations <- combn(points, 4)
  
  tests <- 0
  violations <- 0
  violating_quadruplets <- character(0)
  
  for (i in 1:ncol(combinations)) {
    comb <- combinations[, i]
    A <- comb[1]; B <- comb[2]; C <- comb[3]; D <- comb[4]
    
    d_AB <- dist_matrix[A, B]
    d_AC <- dist_matrix[A, C]
    d_AD <- dist_matrix[A, D]
    d_BC <- dist_matrix[B, C]
    d_BD <- dist_matrix[B, D]
    d_CD <- dist_matrix[C, D]
    
    tests <- tests + 1
    ok <- (d_AB+d_CD <= max(d_AC+d_BD, d_AD+d_BC) + tol) &&
      (d_AC+d_BD <= max(d_AB+d_CD, d_AD+d_BC) + tol) &&
      (d_AD+d_BC <= max(d_AB+d_CD, d_AC+d_BD) + tol)
    
    if (!ok) {
      violations <- violations + 1
    }
  }
  return(list(
    tests = tests,
    violations = violations,
    additive = (violations == 0L),
    violating_quadruplets = violating_quadruplets))
}

#Four-Point-Condition prüfen für dist_matrix_sq

Four_Point_Condition_check <- check_four_point_condition(dist_matrix_sq, tol=0.1)

#Ergebnisse anhängen
cat("\nErgebnis Four-Point-Condition:\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Tests:", Four_Point_Condition_check$tests, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Additivität:", Four_Point_Condition_check$additive, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)
cat("Verletzungen:", Four_Point_Condition_check$violations, "\n", file = "Results_PhylogeneticTree/Ultrametric_Check.txt", append = TRUE)

message("Four-Point-Condition getestet: Results_PhylogeneticTree/Ultrametric_Check.txt")

#========= 6. NJ Tree =========
#Erstellen und speichern des Baums

pdf("Results_PhylogeneticTree/NJ_tree.pdf", width = 12, height = 12)
nj_tree <- nj(dist_matrix_sq)
plot(nj_tree, cex = 0.5,
     type = "phylogram",
     direction = "rightwards",
     main = "Neighbor-Joining phylogenetic tree")
edgelabels(round(nj_tree$edge.length,2),
           bg="white", col="black", frame = "none",
           cex=0.5, adj = c(0.5,0))
dev.off()
message("NJ_tree gespeichert: Results_PhylogeneticTree/NJ_tree.pdf")

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

k_all=6

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

#großteil der "Verschönerungsfunktionen" mit Hilfe von Ki erstellt

# =========================
# 1) MSA laden + Distanz
# =========================
alignment_biostrings <- readAAMultipleAlignment(
  "Results_MultibleSequenceAlign/msa_alignment.fasta",
  format = "fasta"
)
alignment_align <- msaConvert(alignment_biostrings, type = "seqinr::alignment")

dist_msa <- dist.alignment(alignment_align, matrix = "identity")
dist_matrix <- as.matrix(dist_msa)
dist_matrix_sq <- dist_matrix^2
dist_msa_sq <- as.dist(dist_matrix_sq)

# =========================
# 2) hclust Objekte
# =========================
hc_upgma <- hclust(dist_msa_sq, method = "average")
hc_wpgma <- hclust(dist_msa_sq, method = "mcquitty")
hc_ward  <- hclust(dist_msa_sq, method = "ward.D2")

# =========================
# 3) Hilfsfunktionen
# =========================
shorten_label <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  sapply(parts, function(p) paste(p[1:min(2, length(p))], collapse = "_"))
}

# Edge-Farben nach Cluster-Zugehörigkeit
color_tree_by_clusters <- function(phy, cl_membership, palette) {
  n_tips  <- length(phy$tip.label)
  n_edges <- nrow(phy$edge)
  
  descendants <- function(node) {
    kids <- phy$edge[phy$edge[, 1] == node, 2]
    out <- integer(0)
    for (k in kids) {
      if (k <= n_tips) out <- c(out, k) else out <- c(out, descendants(k))
    }
    out
  }
  
  edge_cols <- rep("black", n_edges)
  
  for (e in 1:n_edges) {
    child <- phy$edge[e, 2]
    tips <- if (child <= n_tips) child else descendants(child)
    tip_names <- phy$tip.label[tips]
    
    cl_ids <- unique(cl_membership[tip_names])
    cl_ids <- cl_ids[!is.na(cl_ids)]
    
    if (length(cl_ids) == 1) edge_cols[e] <- palette[cl_ids] else edge_cols[e] <- "black"
  }
  
  edge_cols
}

# Cluster-Namen (du kannst die Texte frei anpassen!)
make_cluster_legend_names <- function(k) {
  # Default-Namen: Cluster 1..k
  paste0("Cluster ", 1:k)
}

# =========================
# 4) Plot-Funktion für UPGMA/WPGMA/WARD (mit Legende)
# =========================
plot_cluster_colored_tree <- function(hc, method_name, out_pdf,
                                      k = 6,
                                      palette = c("deepskyblue3", "forestgreen", "red", "orange", "magenta", "cyan"),
                                      tip_cex = 0.85,
                                      legend_names = NULL) {
  
  phy <- ladderize(as.phylo(hc), right = FALSE)
  phy$tip.label <- shorten_label(phy$tip.label)
  
  # Cluster aus hc
  hc_labels_short <- shorten_label(hc$labels)
  cl <- cutree(hc, k = k)
  names(cl) <- hc_labels_short
  
  # Edge-Farben
  edge_cols <- color_tree_by_clusters(phy, cl, palette)
  
  # Legendennamen
  if (is.null(legend_names)) legend_names <- make_cluster_legend_names(k)
  
  # PDF
  pdf(out_pdf, width = 15, height = 15)
  par(mar = c(7, 14, 3, 2))  # unten extra Platz für Legende
  
  plot(phy,
       type = "phylogram",
       direction = "rightwards",
       use.edge.length = TRUE,
       edge.color = edge_cols,
       cex = tip_cex,
       main = method_name)
  
  # Scale-Bar (sichtbar innerhalb)
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  scale_x <- lp$x.lim[1] + 0.02 * diff(lp$x.lim)
  scale_y <- lp$y.lim[1] + -0.02 * diff(lp$y.lim)  # <-- sicher sichtbar
  add.scale.bar(x = scale_x, y = scale_y, cex = 0.9)
  
  # Legende direkt darunter (bottomleft)

  legend("topleft",
         legend = legend_names[1:k_all],
         col = palette[1:k],
         lwd = 3,
         bty = "n",
         cex = 0.85)
  
  dev.off()
  message("Saved: ", out_pdf)
}

# =========================
# 5) Ausführen: UPGMA/WPGMA/WARD
# =========================
out_dir <- "Results_PhylogeneticTree"

plot_cluster_colored_tree(
  hc_upgma, "UPGMA",
  file.path(out_dir,paste0( "UPGMA_cluster_colored_k_",k_all,".pdf")),
  k = k_all,
  tip_cex = 1,
  legend_names = c(
    "Feuchtnasenprimaten ",   # deepskyblue3
    "Primaten ",                               # forestgreen
    "Nagetiere ",             # red
    "Biber (Castor) – separater Nagetier-Cluster",         # orange
    "Paarhufer + Wal ", # magenta
    "Raubtiere + Flossenfüßer "      # cyan
  )
)

plot_cluster_colored_tree(
  hc_wpgma, "WPGMA",
  file.path(out_dir,paste0( "WPGMA_cluster_colored_k_",k_all,".pdf")),
  k = k_all,
  tip_cex = 1,
  legend_names = c(
    "Feuchtnasenprimaten ",   # deepskyblue3
    "Primaten ",                               # forestgreen
    "Nagetiere ",             # red
    "Biber (Castor) – separater Nagetier-Cluster",         # orange
    "Paarhufer + Wal ", # magenta
    "Raubtiere + Flossenfüßer "      # cyan
    )
)
edgelabels(
  round(phy$edge.length, 3),
  frame = "none",
  cex = 0.6,
  adj = c(0.5, -0.2)
)
plot_cluster_colored_tree(
  hc_ward, "Ward",
  file.path(out_dir, paste0("WARD_cluster_colored_k_",k_all,".pdf")),
  k = k_all,
  tip_cex = 1,
  legend_names = c(
    "Feuchtnasenprimaten ",   # deepskyblue3
    "Primaten ",                               # forestgreen
    "Nagetiere ",             # red
    "Biber (Castor) – separater Nagetier-Cluster",         # orange
    "Paarhufer + Wal ", # magenta
    "Raubtiere + Flossenfüßer "      # cyan
  )
)

# =========================
# 6) NJ Baum (cluster-colored) + Legende
# =========================
nj_tree <- nj(dist_msa_sq)
nj_tree <- ladderize(nj_tree, right = FALSE)
nj_tree$tip.label <- shorten_label(nj_tree$tip.label)

# Cluster-Zuordnung aus Distanzmatrix (damit du "vergleichbare" Cluster bekommst)
k_nj <- k_all
palette_nj <- c("deepskyblue3", "forestgreen", "red","orange", "magenta", "cyan")
legend_names_nj <- c(
  "Feuchtnasenprimaten ",   # deepskyblue3
  "Primaten ",                               # forestgreen
  "Nagetiere ",             # red
  # orange "Biber (Castor) – separater Nagetier-Cluster",
  "Paarhufer + Wal ", # magenta
  "Raubtiere + Flossenfüßer "      # cyan
)

hc_groups <- hclust(dist_msa_sq, method = "mcquitty")
cl_nj <- cutree(hc_groups, k = k_nj)
names(cl_nj) <- shorten_label(hc_groups$labels)
cl_nj <- cl_nj[nj_tree$tip.label]

# Castor demselben Cluster wie die anderen Nager zuordnen
rodent_cluster <- cl_nj["Sciurus_vulgaris"]
cl_nj["Castor_canadensis"] <- rodent_cluster

edge_cols_nj <- color_tree_by_clusters(nj_tree, cl_nj, palette_nj)

pdf(file.path(out_dir,paste0( "NJ_tree_cluster_colored_k_",k_all,".pdf")), width = 14, height = 10)
par(mar = c(7, 16, 3, 2))

plot(nj_tree,
     type = "phylogram",
     direction = "rightwards",
     use.edge.length = TRUE,
     edge.color = edge_cols_nj,
     cex = 1,
     main = "Neighbor-Joining phylogenetic tree (cluster-colored)")
edgelabels(
       round(nj_tree$edge.length, 3),
       frame = "none",
       cex = 0.6,
       adj = c(0.5, -0.2)
     )

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
scale_x <- lp$x.lim[1] + 0.02 * diff(lp$x.lim)
scale_y <- lp$y.lim[1] + -0.02 * diff(lp$y.lim)
add.scale.bar(x = scale_x, y = scale_y, cex = 0.9)


legend("topleft",
       legend = legend_names_nj[1:5],
       col = palette_nj[c(1, 2, 3, 5, 6)],
       lwd = 3,
       bty = "n",
       cex = 0.85)

dev.off()
message("NJ_tree gespeichert: Results_PhylogeneticTree/NJ_tree_cluster_colored_k=",k_all,".pdf")


#---- schöne Heatmap ----

# MSA laden und konvertieren
alignment_biostrings <- readAAMultipleAlignment("Results_MultibleSequenceAlign/msa_alignment.fasta", format = "fasta")
alignment_align <- msaConvert(alignment_biostrings, type = "seqinr::alignment")

# Distanzmatrix (Identität) berechnen
dist_msa <- dist.alignment(alignment_align, matrix = "identity")
dist_matrix <- as.matrix(dist_msa)
dist_matrix_sq <- dist_matrix^2

short_names <- shorten_label(rownames(dist_matrix_sq))
rownames(dist_matrix_sq) <- short_names
colnames(dist_matrix_sq) <- short_names

# Heatmap speichern
pdf("Results_PhylogeneticTree/Distance_Heatmap.pdf", width = 12, height = 12)
heatmap.2(dist_matrix_sq,
          trace = "none",
          col = rev(hcl.colors(100, palette = "heat")),
          margins = c(12, 12),
          cexRow = 0.9,
          cexCol = 0.9,
          dendrogram = "none",
          key = TRUE,
          denscol="black",
          densadj = 0.25,
          main = "Distanzmatrix Heatmap (Identität^2)")
dev.off()
message("Distanz-Heatmap gespeichert: Results_PhylogeneticTree/Distance_Heatmap.pdf")

# =========================
# 7) Baum-Vergleiche
#    - Tanglegrams für UPGMA/WPGMA/WARD
#    - NJ vs UPGMA/WPGMA/WARD
# =========================

# Lokale Kopien mit kurzen Labels
hc_upgma_cmp <- hc_upgma
hc_wpgma_cmp <- hc_wpgma
hc_ward_cmp  <- hc_ward

hc_upgma_cmp$labels <- shorten_label(hc_upgma_cmp$labels)
hc_wpgma_cmp$labels <- shorten_label(hc_wpgma_cmp$labels)
hc_ward_cmp$labels  <- shorten_label(hc_ward_cmp$labels)

# hclust -> phylo für Vergleich mit NJ
tree_upgma_cmp <- ladderize(as.phylo(hc_upgma_cmp), right = FALSE)
tree_wpgma_cmp <- ladderize(as.phylo(hc_wpgma_cmp), right = FALSE)
tree_ward_cmp  <- ladderize(as.phylo(hc_ward_cmp),  right = FALSE)

# NJ-Baum als Vergleichsobjekt
tree_nj_cmp <- nj_tree
tree_nj_cmp$tip.label <- shorten_label(tree_nj_cmp$tip.label)

# Ausgabe-PDF
cmp_pdf <- file.path(out_dir, "Tree_Comparisons_and_NJ.pdf")
pdf(cmp_pdf, width = 14, height = 13)

# -----------------------------
# Tanglegram: UPGMA vs WPGMA
# -----------------------------
dendextend::dendlist(
  as.dendrogram(hc_upgma_cmp),
  as.dendrogram(hc_wpgma_cmp)
) %>%
  dendextend::untangle(method = "step2side") %>%
  dendextend::tanglegram(
    main_left  = "UPGMA",
    main_right = "WPGMA",
    lwd = 2,
    edge.lwd = 2,
    common_subtrees_color_branches = TRUE,
    columns_width = c(5, 1, 5),
    margin_inner = 10,
    lab.cex = 0.75
  )

# -----------------------------
# Tanglegram: UPGMA vs WARD
# -----------------------------
dendextend::dendlist(
  as.dendrogram(hc_upgma_cmp),
  as.dendrogram(hc_ward_cmp)
) %>%
  dendextend::untangle(method = "step2side") %>%
  dendextend::tanglegram(
    main_left  = "UPGMA",
    main_right = "WARD",
    lwd = 2,
    edge.lwd = 2,
    common_subtrees_color_branches = TRUE,
    columns_width = c(5, 1, 5),
    margin_inner = 10,
    lab.cex = 0.75
  )

# -----------------------------
# Tanglegram: WPGMA vs WARD
# -----------------------------
dendextend::dendlist(
  as.dendrogram(hc_wpgma_cmp),
  as.dendrogram(hc_ward_cmp)
) %>%
  dendextend::untangle(method = "step2side") %>%
  dendextend::tanglegram(
    main_left  = "WPGMA",
    main_right = "WARD",
    lwd = 2,
    edge.lwd = 2,
    common_subtrees_color_branches = TRUE,
    columns_width = c(5, 1, 5),
    margin_inner = 10,
    lab.cex = 0.75
  )

# -----------------------------
# NJ vs UPGMA
# -----------------------------
obj_nj_upgma <- phytools::cophylo(tree_nj_cmp, tree_upgma_cmp, rotate = TRUE)
plot(obj_nj_upgma,
     link.type = "curved",
     fsize = 0.75,
     ftype = "i",
     lwd = 2)
title("NJ vs UPGMA")

# -----------------------------
# NJ vs WPGMA
# -----------------------------
obj_nj_wpgma <- phytools::cophylo(tree_nj_cmp, tree_wpgma_cmp, rotate = TRUE)
plot(obj_nj_wpgma,
     link.type = "curved",
     fsize = 0.75,
     ftype = "i",
     lwd = 2)
title("NJ vs WPGMA")

# -----------------------------
# NJ vs WARD
# -----------------------------
obj_nj_ward <- phytools::cophylo(tree_nj_cmp, tree_ward_cmp, rotate = TRUE)
plot(obj_nj_ward,
     link.type = "curved",
     fsize = 0.75,
     ftype = "i",
     lwd = 2)
title("NJ vs WARD")

dev.off()
message("Baumvergleiche gespeichert: ", cmp_pdf)





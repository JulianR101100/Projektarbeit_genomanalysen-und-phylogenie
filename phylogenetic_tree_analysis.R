# phylogenetic_tree_analysis.R
# constructs distance matrices and phylogenetic trees from MSA results.

library(seqinr)
library(gplots)
library(ape)
library(msa)

# 1. Load MSA result
# We can load the RDS saved from the previous script
# Or read the aligned FASTA if preferred.
# Using the msa object provides better conversion options.
if(file.exists("msa_result.rds")) {
  my_msa <- readRDS("msa_result.rds")
} else {
  stop("msa_result.rds not found. Please run multiple_sequence_alignment.R first.")
}

# 2. Convert to seqinr format for distance matrix calculation
msa_seqinr <- msaConvert(my_msa, type = "seqinr::alignment")

# 3. Calculate Distance Matrix
# Using dist.alignment from seqinr
dist_matrix <- dist.alignment(msa_seqinr, matrix = "identity")
# Convert to standard matrix for heatmap
mat <- as.matrix(dist_matrix)

# 4. Visualize with Heatmap
# Using heatmap.2 from gplots
# Adjust margins to fit labels
heatmap.2(mat, 
          main = "Distance Matrix Heatmap",
          trace = "none", 
          margins = c(10, 10),
          col = bluered(75))

# 5. Check Properties (Four-point & Ultrametric)
# Instructions refer to functions defined in lecture.
# Placeholders:
# check_four_point(dist_matrix)
# check_ultrametric(dist_matrix)

# 6. Tree Construction
# Hierarchical Clustering (UPGMA, NJ, etc.)
# 'hclust' uses the distance matrix

# UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
tree_upgma <- hclust(dist_matrix, method = "average")
plot(tree_upgma, main = "Phylogenetic Tree (UPGMA)")

# Neighbor Joining (requires 'ape' package usually, or 'nj' function)
# dist_matrix is a 'dist' object, suitable for ape::nj
tree_nj <- nj(dist_matrix)
plot(tree_nj, main = "Phylogenetic Tree (Neighbor Joining)")

print("Analysis complete. Check plots.")

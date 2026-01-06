# multiple_sequence_alignment.R
# Performs Multiple Sequence Alignment (MSA) using ClustalW.

library(msa)
library(Biostrings)
# library(ggmsa) # Optional for visualization

# 1. Load the sequences extracted from the BLAST step
sequences <- readAAStringSet("top30_sequences.fasta")

# 2. Perform MSA
# Method: ClustalW as per instructions
print("Running ClustalW alignment...")
my_msa <- msa(sequences, method = "ClustalW")

print(my_msa)

# 3. Save the alignment
# Converting to a standard format (e.g., FASTA) allows other tools to use it
# We can also save the msa object itself for R
saveRDS(my_msa, "msa_result.rds")

# Export to FASTA for external viewers or subsequent steps if needed
# msaPrettyPrint output can be complex, but for simple export:
unmasked_alignment <- unmasked(my_msa)
writeXStringSet(unmasked_alignment, "msa_alignment.fasta")

print("MSA completed and saved to 'msa_alignment.fasta' and 'msa_result.rds'")

# 4. Visualization (Optional)
# if(require(ggmsa)) {
#   ggmsa(unmasked_alignment, start = 1, end = 50, color = "Chemistry_AA")
# }

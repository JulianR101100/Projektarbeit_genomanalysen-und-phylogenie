# Projektarbeit Genomik und Phylogenie

## Rohdaten Ordner

In diesem Ordner könnt ihr die zwei datein von StudOn reinkopieren (die Rohdaten), einmal die Datenbank gegen die wir blasten und einmal die querry Sequence.
Die bash skripte erzeugen dann die ergebnisse die aber auch schon gespeichert sind in zb `query_seq_blast_BLOSUM62.blast`. Die werden dann auch in den R skripten verwendet also ich glaube theoretisch müsstet ihr die bash skripte nicht ausführen ig. Hier noch mal zur übersicht der "Arbeitsaufrtag"

---

## Arbeitsauftrag

### BLAST

You find a file called database.fasta on StudOn. This is the database fasta file you will use for the BLAST search. You will also find a file called query_sequence.fasta. This is the query sequence you will use for the BLAST search. Your first job is to turn the database.fasta into a BLAST database using the makeblastdb command in the terminal. Afterwards, run BLAST in the Linux terminal using the bash code we have also used in course:

```bash
blastp -query <QUERY_FILENAME> \
-db <DATABASE_FILENAME> \
-outfmt "7 delim= sacc length pident gaps score bitscore evalue" > <OUTPUT_FILENAME>
```

Interpret the results. What can you conclude from the output file? Which species does the sequence belong to? Which gene is behind the protein sequence? Based on the biological function, would you expect the protein sequence to be conserved or not? Are all the hits in your search orthologs or would you suspect paralogs?

Extract the top 30 hits from the BLAST search. Based on the sequence identity percentage of these top hits and in all hits, would you re-run BLAST with a different substitution matrix? If yes, which one would you choose and why? Remember, there is no clear-cut "one is right everything else is wrong" answer to this question. You should justify your choice based on the information you have and your opinion about the fitting of the model to your research question. Also, if you deem another substitution matrix more appropriate, you should re-run BLAST with this matrix and extract the top 30 hits again!

### Multiple sequence alignment

After you extracted the top 30 hits, extract the according sequences from the database fasta file. Before conducting MSA, think about whether it makes sense to rename your sequences prior to subsequent analyses.

Perform multiple sequence alignment. Use the msa function from the msa package with the ClustalW method. You can also use the ggmsa package to visualize the alignment. Alternatively, you can use a different visualization tool like the alignment viewer.

### Phylogenetic tree analysis

Given the MSA results, construct a distance matrix and visualize it using a heatmap. You can use the dist.alignment function from the seqinr package to calculate the distance matrix. Then, use the heatmap.2 function from the gplots package to visualize the distance matrix. What do you see? Do you already observe any patterns regarding the conservation of the protein sequence? What can you say about the evolutionary relationship between the sequences? Don't forget to check out the manual page for the function to see how to toggle the parameters to get the best visualization for your purpose.

After constructing the distance matrix, check the four-point condition and the ultrametric property. Use the functions we defined in the lecture. Also, construct trees with the clustering algorithms we discussed in the lecture and used in the practical course. Compare the trees. Browse the literature to find out which clustering algorithm is most suitable for your research question. What do you observe? What can you conclude from the trees? Do your results overlap with published results?

Finally, discuss the limitations of your approach and make suggestions for future research. What would you do differently if you had more time, resources, or additional data? Also check you trees for robustness as seen in the lecture and discuss the results.
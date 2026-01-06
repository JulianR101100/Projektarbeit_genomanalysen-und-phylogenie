#!/bin/bash

# 1. Weitere Datein sollten geordnet in verschiedenen Ordnern liegen (-> database_blast_BLOSUM62, logs_BLOSUM62)
mkdir -p database_blast_BLOSUM62 logs_BLOSUM62

# 2. Create BLAST database from database.fasta
# Wie im Skript: "turn the database.fasta into a BLAST database using the makeblastdb command"
# Database datein in Ordner database_blast speichern
makeblastdb -in Rohdaten/database.fasta -dbtype prot \
-out database_blast_BLOSUM62/database \
> logs_BLOSUM62/makeblastdb.log 2>&1

# 3. Run BLAST search with all target sequences
# sacc gibt nur ID der Treffer aus (alternativ stitle ->Sequenzname/Beschreibung; sseqid -> volle Sequenz-ID)
blastp -query Rohdaten/query_sequence.fasta \
-db database_blast_BLOSUM62/database \
-outfmt "7 delim=	 sacc length pident gaps score bitscore evalue" \
-matrix BLOSUM62 \
> query_seq_blast_BLOSUM62.blast

echo "BLAST search mit BlOSUM62 allen möglichen Zielsequenzen durchgeführt und Ergebnisse in 'query_seq_blast_BLOSUM62.blast' gespeichert."
# head -n 100 query_seq_blast_BLOSUM62.blast
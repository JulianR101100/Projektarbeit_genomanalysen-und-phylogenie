#!/bin/bash

# 1. Weitere Datein sollten geordnet in verschiedenen Ordnern liegen (-> database_blast, logs)
mkdir -p datebase_blast_top1000 logs_top1000

# 2. Create BLAST database from database.fasta
# Wie im Skript: "turn the database.fasta into a BLAST database using the makeblastdb command"
# Database datein in Ordner datebase_blast_top1000 speichern
makeblastdb -in Rohdaten/database.fasta -dbtype prot \
-out datebase_blast_top1000/database \
> logs_top1000/makeblastdb.log 2>&1

# 3. Run BLAST search with max_target_seqs set to 1000
# sacc gibt nur ID der Treffer aus (alternativ stitle ->Sequenzname/Beschreibung; sseqid -> volle Sequenz-ID)
blastp -query Rohdaten/query_sequence.fasta \
-db datebase_blast_top1000/database \
-max_target_seqs 1000 \
-outfmt "7 delim=	 sacc length pident gaps score bitscore evalue" \
> query_top1000_seq_blast.blast

echo "BLAST sarch mit maximal 1000 Zielsequenzen durchgeführt und Ergebnisse in 'query_top1000_seq_blast.blast' gespeichert."
# head -n 30 query_top1000_seq_blast.blast
#!/bin/bash

# 1. Weitere Datein sollten geordnet in verschiedenen Ordnern liegen (-> database_blast_PAM30, logs_PAM30)
mkdir -p database_blast_PAM30 logs_PAM30

# 2. Create BLAST database from database.fasta
# Wie im Skript: "turn the database.fasta into a BLAST database using the makeblastdb command"
# Database datein in Ordner database_blast speichern
makeblastdb -in Rohdaten/database.fasta -dbtype prot \
-out database_blast_PAM30/database \
> logs_PAM30/makeblastdb.log 2>&1

# 3. Run BLAST search with all target sequences
# sacc gibt nur ID der Treffer aus (alternativ stitle ->Sequenzname/Beschreibung; sseqid -> volle Sequenz-ID)
blastp -query Rohdaten/query_sequence.fasta \
-db database_blast_PAM30/database \
-outfmt "7 delim=	 sacc length pident gaps score bitscore evalue" \
-matrix PAM30 \
> query_seq_blast_PAM30.blast

echo "BLAST search mit PAM30 allen möglichen Zielsequenzen durchgeführt und Ergebnisse in 'query_seq_blast_PAM30.blast' gespeichert."

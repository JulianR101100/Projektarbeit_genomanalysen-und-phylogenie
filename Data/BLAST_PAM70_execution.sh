#!/bin/bash

# 1. Ordner anlegen
mkdir -p database_blast_PAM70 logs_PAM70

# 2. BLAST-Datenbank erstellen (falls noch nicht geschehen)
makeblastdb -in Rohdaten/database.fasta -dbtype prot \
-out database_blast_PAM70/database \
> logs_PAM70/makeblastdb.log 2>&1

# 3. BLAST mit PAM70
# Hinweis: Standard-Gapkosten für PAM70 sind meist -gapopen 10 -gapextend 1
blastp -query Rohdaten/query_sequence.fasta \
-db database_blast_PAM70/database \
-outfmt "7 delim=    sacc length pident gaps score bitscore evalue" \
-matrix PAM70 \
-gapopen 10 -gapextend 1 \
> query_seq_blast_PAM70.blast

echo "BLAST search mit PAM70 abgeschlossen."

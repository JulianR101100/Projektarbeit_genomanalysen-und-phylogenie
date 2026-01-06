#!/bin/bash

# Die Datei top_hits_summary.csv wurde erfolgreich erstellt und enthält nun die gewünschten Informationen (Organismus und Beschreibung).
  # Beispielausgabe (erste Zeilen):
  # 1 rank,subject_id,organism,description,length,pident,bitscore,evalue
  # 2 1,9606_0:0048f7,Homo sapiens,18S rRNA aminocarboxypropyltransferase,312,100.000,628,0.0
  # 3 2,9544_0:004ef9,Macaca mulatta,ribosome biogenesis protein TSR3 homolog,312,95.192,597,0.0
  # 4 3,9509_0:0037a9,Ateles geoffroyi,,312,91.667,572,0.0
  #(Hinweis: Manche Einträge haben keine Description im Original-Header, wie bei Platz 3 zu sehen).

set -euo pipefail

### PARAMETER ###
BLAST_FILE="query_1000_seq_blast.blast"
FASTA="Rohdaten/database.fasta"
TOP_N=50
OUT_CSV="top_hits_summary.csv"

TMP_HITS=$(mktemp)

### 1. Top-N Treffer extrahieren (Header entfernen)
grep -m "$TOP_N" -v '^#' "$BLAST_FILE" \
  | nl -w1 -s $'\t' \
  > "$TMP_HITS"

### 2. CSV-Header
echo "rank,subject_id,organism,description,length,pident,bitscore,evalue" \
  > "$OUT_CSV"

### 3. Zusammenführen und Annotieren via FASTA-Header
while IFS=$'\t' read -r rank sacc length pident gaps score bitscore evalue; do

    # Headerzeile aus der FASTA-Datei suchen (JSON-Metadaten)
    # Wir suchen nach ">SACC" (mit oder ohne folgendes Leerzeichen)
    header=$(grep -F ">$sacc" "$FASTA" | head -n 1 || true)

    if [ -n "$header" ]; then
        # Extrahiere Organism Name und Description aus dem JSON-String
        organism=$(echo "$header" | sed -n 's/.*"organism_name":"\([^"]*\)".*/\1/p')
        desc=$(echo "$header" | sed -n 's/.*"description":"\([^"]*\)".*/\1/p' | sed 's/,/;/g')
    else
        organism="Unknown"
        desc="Description not found"
    fi

    echo "$rank,$sacc,$organism,$desc,$length,$pident,$bitscore,$evalue" \
      >> "$OUT_CSV"

done < "$TMP_HITS"

### 4. Cleanup
rm -f "$TMP_HITS"

echo "CSV erfolgreich erzeugt: $OUT_CSV"
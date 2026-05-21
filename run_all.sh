#!/bin/bash

FASTQ_DIR="data/fastq"
REFERENCES=("data/J29B_expected.fasta" "data/J30B_expected.fasta")
GAPOPEN=10
GAPEXTEND=0.5
FRAGMENT_LENGTH=296

TMPCONFIG=$(mktemp /tmp/config_XXXXXX.txt)

for fastq in "$FASTQ_DIR"/*.fastq; do
    for ref in "${REFERENCES[@]}"; do
        cat > "$TMPCONFIG" <<EOF
reference_file = $ref
fragments_file = $fastq
gapopen = $GAPOPEN
gapextend = $GAPEXTEND
fragment_length = $FRAGMENT_LENGTH
EOF
        echo "Running: $(basename "$fastq") vs $(basename "$ref")"
        python main.py "$TMPCONFIG"
    done
done

rm "$TMPCONFIG"
echo "All done. Results saved in results/"

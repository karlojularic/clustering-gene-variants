#!/bin/bash

FASTQ_DIR="data/fastq"
REFERENCES=("data/J29B_expected.fasta" "data/J30B_expected.fasta")
GAPOPEN=10
GAPEXTEND=0.5
FRAGMENT_LENGTH=296

TMPCONFIG=$(mktemp /tmp/config_XXXXXX.txt)

MIN_READS=100

for fastq in "$FASTQ_DIR"/J*.fastq; do
    count=$(.venv/bin/python3 -c "
from Bio import SeqIO
print(sum(1 for r in SeqIO.parse('$fastq', 'fastq') if len(r.seq) == $FRAGMENT_LENGTH))
")
    if [ "$count" -lt "$MIN_READS" ]; then
        echo "Skipping $(basename "$fastq") — only $count reads of length $FRAGMENT_LENGTH"
        continue
    fi

    for ref in "${REFERENCES[@]}"; do
        cat > "$TMPCONFIG" <<EOF
reference_file = $ref
fragments_file = $fastq
gapopen = $GAPOPEN
gapextend = $GAPEXTEND
fragment_length = $FRAGMENT_LENGTH
EOF
        echo "Running: $(basename "$fastq") vs $(basename "$ref")"
        .venv/bin/python main.py "$TMPCONFIG"
    done
done

rm "$TMPCONFIG"
echo "All done. Results saved in results/"

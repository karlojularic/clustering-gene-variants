#!/usr/bin/env python3
import argparse
from Bio import SeqIO

from distance import distance_matrix

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('fragments', type=str)

    args = parser.parse_args()

    reference = SeqIO.parse(args.reference, 'fasta')
    fragments = SeqIO.parse(args.fragments, 'fastq')

    stored_fragments = {}    # Of length 296
    total = 0

    for fragment in fragments:
        total += 1
        if len(fragment.seq) == 296:
            stored_fragments[fragment.id] = str(fragment.seq)

    print(f"Number of fragments stored: {len(stored_fragments)}")
    print(f"Total number of fragments: {total}")

    matrix = distance_matrix(stored_fragments)
    print(f"Distance matrix:")
    for i, row in enumerate(matrix):
        print(f"Fragment {i}: {row}")

if __name__ == '__main__':
    main()
#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import spoa

from distance import distance_matrix
from clustering import greedy_clustering

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

    conseunsus, msa = spoa.poa(stored_fragments.values())

    matrix = distance_matrix(msa)
    print(f"Distance matrix:")
    for i, row in enumerate(matrix):
        print(f"Fragment {i}: {row}")

    clusters = greedy_clustering(msa, matrix)
    new_clusters = {}
    for cluster_id, cluster in clusters.items():
        if len(cluster) >= 10:              # cluster has to have at least 10 reads
            new_clusters[cluster_id] = cluster

    for cluster_id, cluster in new_clusters.items():
        print(f"Cluster {cluster_id}: {cluster}")

    
if __name__ == '__main__':
    main()
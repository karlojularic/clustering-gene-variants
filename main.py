#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import spoa

from distance import distance_matrix
from clustering import greedy_clustering, cluster_representative

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('fragments', type=str)

    args = parser.parse_args()

    reference = SeqIO.parse(args.reference, 'fasta')
    fragments = SeqIO.parse(args.fragments, 'fastq')

    stored_fragments = {}
    total = 0

    for fragment in fragments:
        total += 1
        if len(fragment.seq) == 296:
            stored_fragments[fragment.id] = str(fragment.seq)

    print(f"Number of fragments stored: {len(stored_fragments)}")
    print(f"Total number of fragments: {total}")

    if len(stored_fragments) == 0:
        print("No fragments of length 296 found.")
        return

    consensus, msa = spoa.poa(stored_fragments.values())

    matrix = distance_matrix(msa)
    print(f"Distance matrix:")
    for i, row in enumerate(matrix):
        print(f"Fragment {i}: {row}")

    clusters = greedy_clustering(msa, matrix)
    new_clusters = {}
    for cluster_id, cluster in clusters.items():
        if len(cluster) >= 10:              # cluster has to have at least 10 reads
            new_clusters[cluster_id] = cluster

    if len(new_clusters) == 0:
        print("No clusters with at least 10 fragments found.")
        return

    fragment_ids = list(stored_fragments.keys())
    fragment_sequences = list(stored_fragments.values())

    with open("representatives.fasta", "w", encoding="utf-8") as output_file:
        for cluster_id, cluster in new_clusters.items():
            representative_index = cluster_representative(cluster, matrix)
            representative_id = fragment_ids[representative_index]
            representative_sequence = fragment_sequences[representative_index]

            print(f"Cluster {cluster_id}: {cluster}")
            print(f"Representative: {representative_id}")

            output_file.write(f">cluster_{cluster_id}_{representative_id}\n")
            output_file.write(f"{representative_sequence}\n")

    print("Representatives written to representatives.fasta")

if __name__ == '__main__':
    main()
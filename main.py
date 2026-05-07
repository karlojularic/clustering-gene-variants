#!/usr/bin/env python3
import argparse
import subprocess
from Bio import SeqIO
import spoa

from distance import distance_matrix
from clustering import greedy_clustering, cluster_representative

def parsing_config(config_file):

    config = {}

    with open(config_file, 'r') as file:
        for line in file:
            if '=' in line:
                key, value = line.split('=', 1)
                config[key.strip()] = value.strip()
    return config


def emboss_testing(representatives_file, reference_file, gapopen, gapextend):
    result = subprocess.run(
        [
            "needle",
            representatives_file,
            reference_file,
            "-gapopen", str(gapopen),
            "-gapextend", str(gapextend),
            "-outfile", "representative_vs_expected.txt",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"needle failed:\n{result.stderr}")
    else:
        print("Alignment written to representative_vs_expected.txt")
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()

    config = parsing_config(args.config_file)
    reference_file = config['reference_file']
    fragments_file = config['fragments_file']
    gapopen = config['gapopen']
    gapextend = config['gapextend']
    fragment_length = config['fragment_length']

    reference = SeqIO.parse(reference_file, 'fasta')
    fragments = SeqIO.parse(fragments_file, 'fastq')

    stored_fragments = {}
    total = 0

    for fragment in fragments:
        total += 1
        if len(fragment.seq) == int(fragment_length):
            stored_fragments[fragment.id] = str(fragment.seq)

    print(f"Number of fragments stored: {len(stored_fragments)}")
    print(f"Total number of fragments: {total}")

    if len(stored_fragments) == 0:
        print(f"No fragments of length {fragment_length} found.")
        return

    consensus, msa = spoa.poa(stored_fragments.values())

    matrix = distance_matrix(msa)

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
    emboss_testing("representatives.fasta", reference_file, gapopen, gapextend)

if __name__ == '__main__':
    main()
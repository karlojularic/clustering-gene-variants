#!/usr/bin/env python3
import argparse
import os
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


def emboss_testing(representatives_file, reference_file, gapopen, gapextend, output_file):
    result = subprocess.run(
        [
            "needle",
            representatives_file,
            reference_file,
            "-gapopen", str(gapopen),
            "-gapextend", str(gapextend),
            "-outfile", output_file,
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"needle failed:\n{result.stderr}")
    else:
        print(f"Alignment written to {output_file}")
    


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

    sample_name = os.path.splitext(os.path.basename(fragments_file))[0]
    ref_name = os.path.splitext(os.path.basename(reference_file))[0]
    output_prefix = f"results/{sample_name}_vs_{ref_name}"
    os.makedirs("results", exist_ok=True)

    reference = SeqIO.parse(reference_file, 'fasta')
    fragments = SeqIO.parse(fragments_file, 'fastq')

    all_fragments = {f.id: str(f.seq) for f in fragments}
    total = len(all_fragments)

    target_length = int(fragment_length)
    stored_fragments = {}

    for fid, seq in all_fragments.items():
        if abs(len(seq) - target_length) <= 5:          # length is not capped at the set number
            stored_fragments[fid] = seq

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

    representatives_file = f"{output_prefix}_representatives.fasta"
    with open(representatives_file, "w", encoding="utf-8") as output_file:
        for cluster_id, cluster in new_clusters.items():
            representative_index = cluster_representative(cluster, matrix)
            representative_id = fragment_ids[representative_index]
            representative_sequence = fragment_sequences[representative_index]

            print(f"Cluster {cluster_id}: {cluster}")
            print(f"Representative: {representative_id}")

            output_file.write(f">cluster_{cluster_id}_{representative_id}\n")
            output_file.write(f"{representative_sequence}\n")

    print(f"Representatives written to {representatives_file}")

    for cluster_id, cluster in new_clusters.items():
        representative_index = cluster_representative(cluster, matrix)
        representative_sequence = fragment_sequences[representative_index]
        representative_id = fragment_ids[representative_index]

        cluster_fasta = f"{output_prefix}_cluster{cluster_id}_representative.fasta"
        with open(cluster_fasta, "w", encoding="utf-8") as f:
            f.write(f">cluster_{cluster_id}_{representative_id}\n")
            f.write(f"{representative_sequence}\n")

        alignment_file = f"{output_prefix}_cluster{cluster_id}_alignment.txt"
        emboss_testing(cluster_fasta, reference_file, gapopen, gapextend, alignment_file)

if __name__ == '__main__':
    main()
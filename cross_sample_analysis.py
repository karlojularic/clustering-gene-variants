#!/usr/bin/env python3

import os
import re
import csv
import subprocess
import tempfile
import spoa
from Bio import SeqIO
from distance import distance_matrix
from clustering import greedy_clustering, cluster_representative


def collect_representatives(results_dir):
    reps = []
    pattern = re.compile(r"(.+)_vs_J29B_expected_(cluster\d+)_representative\.fasta")

    for filename in sorted(os.listdir(results_dir)):
        m = pattern.match(filename)
        if not m:
            continue
        sample, cluster = m.group(1), m.group(2)
        filepath = os.path.join(results_dir, filename)
        for record in SeqIO.parse(filepath, "fasta"):
            reps.append((sample, cluster, str(record.seq)))

    return reps


def run_needle(query_fasta, reference_fasta, output_file, gapopen, gapextend):
    result = subprocess.run(
        ["needle", query_fasta, reference_fasta,
         "-gapopen", str(gapopen), "-gapextend", str(gapextend),
         "-outfile", output_file],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"needle error: {result.stderr}")


def parse_needle_identity(filepath):
    results = {}
    current_ref = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            m = re.match(r"^# 2: (.+)", line)
            if m:
                current_ref = m.group(1)
            m = re.match(r"^# Identity:\s+\d+/\d+ \(([\d.]+)%\)", line)
            if m and current_ref:
                results[current_ref] = float(m.group(1))
                current_ref = None
    return results


def main():
    results_dir = "results"
    ref_j29 = "data/J29B_expected.fasta"
    ref_j30 = "data/J30B_expected.fasta"
    gapopen = 10
    gapextend = 0.5
    output_csv = "cross_sample_results.csv"
    variant_fasta = "variant_representatives.fasta"

    print("Collecting cluster representatives...")
    reps = collect_representatives(results_dir)
    if not reps:
        print(f"No files found in {results_dir}. Run run_all.sh first.")
        return

    print(f"Found {len(reps)} representatives from {len(set(r[0] for r in reps))} samples.")

    sequences = [seq for _, _, seq in reps]
    labels = [(sample, cluster) for sample, cluster, _ in reps]

    print("Running SPOA MSA on all representatives...")
    _, msa = spoa.poa(sequences)

    print("Computing distance matrix...")
    matrix = distance_matrix(msa)

    print("Greedy clustering of representatives...")
    raw_clusters = greedy_clustering(msa, matrix)

    variants = []
    for variant_id, members in raw_clusters.items():
        unique_samples = sorted(set(labels[i][0] for i in members))
        rep_index = cluster_representative(members, matrix)

        variants.append({
            "variant_id": variant_id,
            "rep_sequence": sequences[rep_index],
            "sample_count": len(unique_samples),
            "samples": unique_samples,
        })

    print(f"\nFound {len(variants)} unique variants across {len(reps)} clusters.\n")

    with open(variant_fasta, "w") as f:
        for v in variants:
            f.write(f">variant_{v['variant_id']}\n{v['rep_sequence']}\n")
    print(f"Variant representatives written to {variant_fasta}")

    print("Aligning variants against reference alleles...")

    for v in variants:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
            tmp.write(f">variant_{v['variant_id']}\n{v['rep_sequence']}\n")
            tmp_path = tmp.name

        for ref_name, ref_file in [("J29B", ref_j29), ("J30B", ref_j30)]:
            out_file = os.path.join(results_dir, f"variant_{v['variant_id']}_vs_{ref_name}_alignment.txt")
            run_needle(tmp_path, ref_file, out_file, gapopen, gapextend)
            identities = parse_needle_identity(out_file)
            if identities:
                best_allele = max(identities, key=identities.get)
                v[f"best_{ref_name}_allele"] = best_allele
                v[f"best_{ref_name}_identity"] = identities[best_allele]
            else:
                v[f"best_{ref_name}_allele"] = ""
                v[f"best_{ref_name}_identity"] = 0.0

        os.unlink(tmp_path)

    fieldnames = [
        "variant_id", "sample_count", "samples",
        "best_J29B_allele", "best_J29B_identity",
        "best_J30B_allele", "best_J30B_identity",
    ]
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for v in sorted(variants, key=lambda x: -x["sample_count"]):
            writer.writerow({
                "variant_id": v["variant_id"],
                "sample_count": v["sample_count"],
                "samples": ";".join(v["samples"]),
                "best_J29B_allele": v.get("best_J29B_allele", ""),
                "best_J29B_identity": v.get("best_J29B_identity", ""),
                "best_J30B_allele": v.get("best_J30B_allele", ""),
                "best_J30B_identity": v.get("best_J30B_identity", ""),
            })

    print(f"\nResults written to {output_csv}")


if __name__ == "__main__":
    main()

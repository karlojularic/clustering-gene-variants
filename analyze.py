#!/usr/bin/env python3
import os
import re
import csv

def parse_alignment_file(filepath):
    alignments = []
    current = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            m = re.match(r"^# 1: (.+)", line)
            if m:
                current["query"] = m.group(1)
            m = re.match(r"^# 2: (.+)", line)
            if m:
                current["reference"] = m.group(1)
            m = re.match(r"^# Identity:\s+(\d+)/(\d+) \(([\d.]+)%\)", line)
            if m:
                current["identity"] = float(m.group(3))
            if line.startswith("#===") and "query" in current and "identity" in current:
                alignments.append(current)
                current = {}
    if "query" in current and "identity" in current:
        alignments.append(current)
    return alignments


def main():
    rows = []

    for filename in sorted(os.listdir("results")):
        if not filename.endswith("_alignment.txt"):
            continue
        if "_cluster" not in filename:
            continue

        m = re.match(r"(.+)_vs_(.+)_(cluster\d+)_alignment\.txt", filename)
        if not m:
            continue
        sample, ref, cluster = m.group(1), m.group(2), m.group(3)

        filepath = os.path.join("results", filename)
        alignments = parse_alignment_file(filepath)

        if not alignments:
            continue

        best = max(alignments, key=lambda a: a["identity"])

        row = {
            "sample": sample,
            "reference_db": ref,
            "cluster": cluster,
            "best_allele": best.get("reference", ""),
            "best_identity": best["identity"],
        }
        for a in alignments:
            row[a.get("reference", "unknown")] = a["identity"]

        rows.append(row)

    if not rows:
        print("No results found.")
        return

    allele_cols = []
    for row in rows:
        for k in row:
            if k not in ("sample", "reference_db", "cluster", "best_allele", "best_identity") and k not in allele_cols:
                allele_cols.append(k)

    fieldnames = ["sample", "reference_db", "cluster", "best_allele", "best_identity"] + allele_cols

    output_csv = "analysis_results.csv"
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    current_sample = None
    for row in rows:
        if row["sample"] != current_sample:
            current_sample = row["sample"]
            print(f"\n{current_sample} (ref: {row['reference_db']})")
        print(f"  {row['cluster']}: {row['best_allele']} ({row['best_identity']}%)")
    
    analyze_global_variants()


def sequence_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        return 1.0

    differences = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            differences += 1

    return differences / len(seq1)


def analyze_global_variants():
    threshold = 0.05
    cluster_rows = []

    for filename in sorted(os.listdir("results")):
        if not filename.endswith("_clusters.csv"):
            continue

        filepath = os.path.join("results", filename)

        with open(filepath, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row["cluster_size"] = int(row["cluster_size"])
                cluster_rows.append(row)

    if not cluster_rows:
        print("No cluster summary files found for global variant analysis.")
        return

    variants = []

    for row in cluster_rows:
        sequence = row["representative_sequence"]
        added = False

        for variant in variants:
            distance = sequence_distance(sequence, variant["representative_sequence"])

            if distance < threshold:
                variant["clusters"].append(row)
                variant["total_reads"] += row["cluster_size"]
                added = True
                break

        if not added:
            variants.append({
                "variant_id": f"variant_{len(variants) + 1}",
                "representative_sequence": sequence,
                "clusters": [row],
                "total_reads": row["cluster_size"]
            })

    total_reads = sum(v["total_reads"] for v in variants)

    with open("global_variant_frequencies.csv", "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "variant_id",
            "total_reads",
            "frequency_percent",
            "number_of_clusters",
            "samples",
            "representative_sequence"
        ]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for variant in variants:
            samples = sorted(set(row["sample"] for row in variant["clusters"]))
            frequency = (variant["total_reads"] / total_reads) * 100

            writer.writerow({
                "variant_id": variant["variant_id"],
                "total_reads": variant["total_reads"],
                "frequency_percent": round(frequency, 2),
                "number_of_clusters": len(variant["clusters"]),
                "samples": ";".join(samples),
                "representative_sequence": variant["representative_sequence"]
            })

    print(f"\nGlobal variant analysis:")
    print(f"Found {len(variants)} global variants.")
    print("Written: global_variant_frequencies.csv")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import os
import re
import csv

RESULTS_DIR = "results"

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

    for filename in sorted(os.listdir(RESULTS_DIR)):
        if not filename.endswith("_alignment.txt"):
            continue
        if "_cluster" not in filename:
            continue

        # parse filename: SAMPLE_vs_REF_clusterN_alignment.txt
        m = re.match(r"(.+)_vs_(.+)_(cluster\d+)_alignment\.txt", filename)
        if not m:
            continue
        sample, ref, cluster = m.group(1), m.group(2), m.group(3)

        filepath = os.path.join(RESULTS_DIR, filename)
        alignments = parse_alignment_file(filepath)

        if not alignments:
            continue

        # find best matching allele
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

    # collect all allele column names
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

    print(f"Saved {len(rows)} rows to {output_csv}")

    # print summary per sample
    print("\n--- Summary (best allele per cluster) ---")
    current_sample = None
    for row in rows:
        if row["sample"] != current_sample:
            current_sample = row["sample"]
            print(f"\n{current_sample} (ref: {row['reference_db']})")
        print(f"  {row['cluster']}: {row['best_allele']} ({row['best_identity']}%)")


if __name__ == "__main__":
    main()

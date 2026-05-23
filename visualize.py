#!/usr/bin/env python3
import csv
import collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

CSV_FILE = "analysis_results.csv"
OUTPUT_PDF = "analysis_report.pdf"

def load_data():
    rows = []
    with open(CSV_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["best_identity"] = float(row["best_identity"])
            rows.append(row)
    return rows


def plot_identity_distribution(rows, pdf):
    # only J30B reference (cleaner allele names)
    data = [r["best_identity"] for r in rows if "J30B" in r["reference_db"]]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(data, bins=20, color="steelblue", edgecolor="white")
    ax.set_xlabel("Identičnost s najbližim alelom (%)")
    ax.set_ylabel("Broj klastera")
    ax.set_title("Distribucija identičnosti klastera s referentnim alelima (J30B referenca)")
    ax.axvline(80, color="red", linestyle="--", label="80% prag")
    ax.legend()
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_allele_frequency(rows, pdf):
    for ref_label, ref_key in [("J29B referenca", "J29B"), ("J30B referenca", "J30B")]:
        subset = [r for r in rows if ref_key in r["reference_db"]]
        counter = collections.Counter(r["best_allele"] for r in subset)

        alleles = list(counter.keys())
        counts = [counter[a] for a in alleles]

        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(alleles, counts, color="steelblue", edgecolor="white")
        ax.set_xlabel("Alel")
        ax.set_ylabel("Broj klastera koji se vežu na taj alel")
        ax.set_title(f"Učestalost alela kroz sve uzorke — {ref_label}")
        ax.bar_label(bars)
        plt.xticks(rotation=20, ha="right")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def plot_clusters_per_sample(rows, pdf):
    for ref_label, ref_key in [("J29B referenca", "J29B"), ("J30B referenca", "J30B")]:
        subset = [r for r in rows if ref_key in r["reference_db"]]
        counter = collections.Counter(r["sample"] for r in subset)

        samples = sorted(counter.keys(), key=lambda s: s.split("_")[0] + s)
        counts = [counter[s] for s in samples]
        short_names = [s.replace("_CE_IonXpress_", "_") for s in samples]

        fig, ax = plt.subplots(figsize=(12, 5))
        bars = ax.bar(short_names, counts, color="steelblue", edgecolor="white")
        ax.set_xlabel("Uzorak")
        ax.set_ylabel("Broj pronađenih klastera")
        ax.set_title(f"Broj klastera po uzorku — {ref_label}")
        ax.bar_label(bars)
        plt.xticks(rotation=45, ha="right", fontsize=7)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def plot_identity_per_sample(rows, pdf):
    for ref_label, ref_key in [("J29B referenca", "J29B"), ("J30B referenca", "J30B")]:
        subset = [r for r in rows if ref_key in r["reference_db"]]

        # group by sample, take max identity per sample
        sample_identities = collections.defaultdict(list)
        for r in subset:
            sample_identities[r["sample"]].append(r["best_identity"])

        samples = sorted(sample_identities.keys())
        max_identities = [max(sample_identities[s]) for s in samples]
        short_names = [s.replace("_CE_IonXpress_", "_") for s in samples]

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.bar(short_names, max_identities, color="steelblue", edgecolor="white")
        ax.set_xlabel("Uzorak")
        ax.set_ylabel("Najviša identičnost (%)")
        ax.set_title(f"Najviša identičnost po uzorku — {ref_label}")
        ax.axhline(80, color="red", linestyle="--", label="80% prag")
        ax.set_ylim(60, 100)
        ax.legend()
        plt.xticks(rotation=45, ha="right", fontsize=7)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def main():
    rows = load_data()
    print(f"Učitano {len(rows)} redova iz {CSV_FILE}")

    with PdfPages(OUTPUT_PDF) as pdf:
        plot_identity_distribution(rows, pdf)
        plot_allele_frequency(rows, pdf)
        plot_clusters_per_sample(rows, pdf)
        plot_identity_per_sample(rows, pdf)

    print(f"PDF spremljen: {OUTPUT_PDF}")


if __name__ == "__main__":
    main()

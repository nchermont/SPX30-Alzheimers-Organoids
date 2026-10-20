"""Freeze the gene annotation used by the original paper figures.

The original analysis annotated genes three times, with whatever versions
were online at the time (August 2025):
  * biomaRt (Ensembl `external_gene_name`) for the interaction and
    normalised-count tables;
  * org.Hs.eg.db (SYMBOL) for the per-group ISS-vs-GROUND tables;
  * a second biomaRt query, from a different Ensembl release, for the two
    ALL/ tables (about 1% of their names differ from the first query).
None of these versions can be queried again, so this script extracts the three
mappings from the published tables into data/gene_annotation.tsv, which
workflow/scripts/paper_figures.R then uses offline.

Columns of data/gene_annotation.tsv:
  ensembl         Ensembl gene id
  symbol_biomart  biomaRt external_gene_name; "NA" = id not returned by biomaRt,
                  "" = returned without a name
  symbol_orgdb    org.Hs.eg.db SYMBOL; "NA" = no symbol
  symbol_all      name used in the ALL/ tables (Ensembl id when unnamed);
                  "NA" = gene not in those tables

Ensembl ids written in place of missing org.Hs.eg.db symbols are treated as
missing.

Usage:
  python scripts/freeze_gene_annotation.py --published path/to/paper_figures \
      [--out data/gene_annotation.tsv]
"""

import argparse
import csv
from pathlib import Path

# the 70-sample interaction tables cover every gene of the 46-sample ones
BIOMART_TABLES = [
    "interaction_APP_vs_WTAD_ground_vs_iss.csv",
    "interaction_APP_vs_WTAD_ground_vs_iss (1).csv",
]
ALL_TABLES = [
    "ALL/DESeq2_ISS_vs_GROUND_results_ordered.csv",
    "ALL/DESeq2_AD_vs_WTAD_results_ordered.csv",
]
ORGDB_TABLES = [
    "AD-WTAD/ad_iss_vs_ground_all.csv",
    "AD-WTAD/wtad_iss_vs_ground_all.csv",
    "APP/app_iss_vs_ground_all.csv",
    "PSEN1/psen1_iss_vs_ground_all.csv",
    "PSEN2/psen2_iss_vs_ground_all.csv",
]


def rows(published, name):
    with open(published / name, newline="") as fh:
        yield from csv.DictReader(fh)


def main(published, out_file):
    biomart = {}
    for name in BIOMART_TABLES:
        for r in rows(published, name):
            biomart.setdefault(r["ensembl_gene_id"], r["external_gene_name"])

    orgdb = {}
    for name in ORGDB_TABLES:
        for r in rows(published, name):
            symbol = r["symbol"]
            if symbol.startswith("ENSG"):
                continue
            previous = orgdb.setdefault(r["ensembl"], symbol)
            assert previous == symbol, f"{r['ensembl']}: {previous} != {symbol}"

    # ALL tables repeat a gene once per biomaRt name; keep the first
    all_names = {}
    for name in ALL_TABLES:
        for r in rows(published, name):
            all_names.setdefault(r["ensembl"], r["gene_symbol"])

    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w", newline="") as fh:
        out = csv.writer(fh, delimiter="\t", lineterminator="\n")
        out.writerow(["ensembl", "symbol_biomart", "symbol_orgdb", "symbol_all"])
        for gene in sorted(set(biomart) | set(orgdb)):
            out.writerow([gene, biomart.get(gene, "NA"), orgdb.get(gene, "NA"),
                          all_names.get(gene, "NA")])
    print(f"{len(set(biomart) | set(orgdb))} genes "
          f"({len(biomart)} biomaRt, {len(orgdb)} org.Hs.eg.db)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--published", type=Path, required=True,
                        help="directory with the published tables")
    parser.add_argument("--out", type=Path, default=Path("data/gene_annotation.tsv"))
    args = parser.parse_args()
    main(args.published, args.out)

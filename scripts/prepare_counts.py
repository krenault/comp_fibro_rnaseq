#!/usr/bin/env python3
"""Combine species prepDE matrices, collapse to one-to-one ENSG orthologs,
restrict to TPM samples.

Usage:
  python3 prepare_counts.py --counts-dir DIR --tpm FILE --out FILE
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd

EXPERIMENTS = ("glucose", "hypoxia", "temperature")
ENSG_RE = re.compile(r"(ENSG\d+)")

def discover_species(base: Path) -> list[Path]:
    dirs = []
    for p in sorted(base.iterdir()):
        if p.is_dir() and any((p / e / "gene_count_matrix.csv").exists() for e in EXPERIMENTS):
            dirs.append(p)
    return dirs

def load_species(species_dir: Path) -> pd.DataFrame | None:
    parts = []
    for exp in EXPERIMENTS:
        f = species_dir / exp / "gene_count_matrix.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, index_col=0)
        parts.append(df)
    if not parts:
        return None
    merged = pd.concat(parts, axis=1)
    return merged

def combine(base: Path) -> pd.DataFrame:
    stacked = []
    for d in discover_species(base):
        mat = load_species(d)
        if mat is None:
            print(f"  skip {d.name}")
            continue
        print(f"  {d.name}: {mat.shape[0]:,} genes × {mat.shape[1]} samples")
        stacked.append(mat)
    if not stacked:
        raise SystemExit(f"No gene_count_matrix.csv files under {base}")
    out = pd.concat(stacked, axis=0)
    out.index.name = "gene_id"
    return out.reset_index()


def collapse_orthologs(df: pd.DataFrame, ref_genes: set[str]) -> pd.DataFrame:
    gene_col = df.columns[0]
    df = df.copy()
    df["base_ensg"] = df[gene_col].astype(str).map(
        lambda x: m.group(1) if (m := ENSG_RE.search(x)) else x
    )
    count_cols = [c for c in df.columns if c not in {gene_col, "base_ensg"}]
    keep = df[df["base_ensg"].isin(ref_genes)]
    agg = keep.groupby("base_ensg", sort=False)[count_cols].sum()
    agg.index.name = gene_col
    print(f"  collapsed {len(keep):,} rows → {len(agg):,} orthologs × {len(count_cols)} samples")
    missing = ref_genes - set(agg.index.astype(str))
    if missing:
        print(f"  {len(missing):,} reference genes absent from counts")
    return agg

def mask_na(counts: pd.DataFrame, tpm: pd.DataFrame) -> pd.DataFrame:
    shared_samples = [c for c in counts.columns if c in tpm.columns]
    dropped = set(counts.columns) - set(shared_samples)
    if dropped:
        print(f"  dropping {len(dropped)} count samples not in TPM")
    counts = counts[shared_samples].astype(float)
    tpm = tpm.reindex(index=counts.index, columns=shared_samples)
    mask = tpm.isna() & (counts == 0)
    n_nonzero_kept = int((tpm.isna() & (counts > 0)).sum().sum())
    out = counts.where(~mask, np.nan)
    print(f"  masked {int(mask.sum().sum()):,} zeros→NA; kept {n_nonzero_kept:,} nonzero where TPM is NA")
    print(f"  final: {out.shape[0]:,} genes × {out.shape[1]} samples")
    return out

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts-dir", type=Path, required=True,
                   help="Directory of <species>/<experiment>/gene_count_matrix.csv")
    p.add_argument("--tpm", type=Path, required=True,
                   help="One-to-one ortholog TPM matrix (first column = ENSG IDs)")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    print("1. combine species matrices")
    combined = combine(args.counts_dir)

    print("2. collapse to one-to-one ENSG orthologs")
    ref = set(pd.read_csv(args.tpm, usecols=[0]).iloc[:, 0].astype(str))
    print(f"  {len(ref):,} orthologs in TPM")
    one2one = collapse_orthologs(combined, ref)

    print("3. filter samples to TPM and mask orthology NAs")
    tpm = pd.read_csv(args.tpm, index_col=0)
    final = mask_na(one2one, tpm)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

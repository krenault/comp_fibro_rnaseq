#!/usr/bin/env python3
"""Sanity checks for the ortholog count matrix. Exit 1 on FAIL."""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ENSG_RE = re.compile(r"^(ENSG\d+)")


def extract_ensg(gene_id: str) -> str | None:
    m = ENSG_RE.match(str(gene_id))
    return m.group(1) if m else None


def check_pre_merge(path: Path) -> list[str]:
    print(f"\n=== pre-merge: ≤1 nonzero row per ENSG per sample ({path.name}) ===")
    if not path.exists():
        return [f"SKIP: {path}"]
    df = pd.read_csv(path, index_col=0)
    ensg = pd.Series([extract_ensg(x) for x in df.index], index=df.index)
    fails = []
    for col in df.columns:
        active = df[col].astype(float)
        active = active[active > 0]
        if active.empty:
            continue
        dup = ensg.loc[active.index]
        dup = dup[dup.duplicated(keep=False)]
        if len(dup):
            fails.append(f"  FAIL {col}: {len(dup)} nonzero rows share an ENSG")
    if fails:
        print("\n".join(fails[:10]))
        return fails
    print(f"  PASS: {df.shape[1]} samples")
    return []


def check_sum_eq_max(path: Path, tpm: Path | None) -> list[str]:
    print(f"\n=== sum vs max aggregation ({path.name}) ===")
    if not path.exists():
        return [f"SKIP: {path}"]
    df = pd.read_csv(path, index_col=0)
    ensg = pd.Series([extract_ensg(x) for x in df.index], index=df.index)
    ok = ensg.notna()
    df, ensg = df.loc[ok], ensg.loc[ok]
    if tpm and tpm.exists():
        ref = set(pd.read_csv(tpm, usecols=[0]).iloc[:, 0].astype(str))
        keep = ensg.isin(ref)
        df, ensg = df.loc[keep], ensg.loc[keep]
    summed, maxed = df.groupby(ensg).sum(), df.groupby(ensg).max()
    n_diff = ((summed - maxed).abs() > 0).any(axis=1).sum()
    if n_diff:
        msg = [f"  FAIL: {n_diff} ENSGs where sum != max"]
        print(msg[0])
        return msg
    print(f"  PASS: {len(summed):,} ENSGs")
    return []


def check_final(path: Path) -> list[str]:
    print(f"\n=== final matrix ({path.name}) ===")
    if not path.exists():
        return [f"SKIP: {path}"]
    df = pd.read_csv(path, index_col=0)
    print(f"  {df.shape[0]:,} genes × {df.shape[1]} samples")
    bad = [x for x in df.index[:500] if not str(x).startswith("ENSG")]
    if bad:
        return [f"  FAIL: non-ENSG index (e.g. {bad[0]})"]
    n_na = int(df.isna().sum().sum())
    print(f"  NA cells: {n_na:,}")
    if n_na == 0:
        print("  WARN: no NAs — ortholog masking may not have been applied")
    print("  PASS")
    return []


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--combined", type=Path, help="pre-collapse stacked counts")
    p.add_argument("--final", type=Path, required=True, help="NA-masked one2one matrix")
    p.add_argument("--tpm", type=Path)
    args = p.parse_args()
    fails: list[str] = []
    if args.combined:
        fails += check_pre_merge(args.combined)
        fails += check_sum_eq_max(args.combined, args.tpm)
    fails += check_final(args.final)
    real = [f for f in fails if "FAIL" in f]
    if real:
        print(f"\nFAILED: {len(real)} check(s)")
        sys.exit(1)
    print("\nAll checks passed (or skipped).")


if __name__ == "__main__":
    main()

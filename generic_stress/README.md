# Public generic cell stress vs Drive fibroblasts

Scripts only. Four-class public cellular stress signature (ASTRA + extra cell
GEO). Drive fibroblasts are held out of discovery, then scored with
class-holdout GSEA and per-species HMP.

**Methods:** [METHODS.md](METHODS.md). Drive counting/DESeq2 is the parent
pipeline in [krenault/comp_fibro_rnaseq](https://github.com/krenault/comp_fibro_rnaseq).
Tables, plots, and DGE CSVs are written locally and are not in git.

## Run

```bash
# Put pairwise Drive DESeq2 CSVs in drive_dge/ (or set DRIVE_DGE)
# optional: cp config.example.env .env
bash scripts/generic_stress/00_run.sh          # skip ingest
bash scripts/generic_stress/00_run.sh --ingest # re-download ASTRA + extra GEO
```

| Script | What |
|--------|------|
| `00_run.sh` | runner |
| `01_ingest.R` | ASTRA + extra cell GEO (`01_astra`, `01_geo`, `01_metadata`) |
| `02_public_hmp.R` | mixed-effect class models + HMP (1,344 genes) |
| `03_drive_holdout.R` | class-holdout GSEA on Drive |
| `04_plots.R` | figures in `plots/mammalian_generic_stress/` |
| `05_species_hmp.R` | per-species Drive HMP, conservation, ORA |
| `06_reml_celltype.R` | optional REML / cell-type Jaccard |

R packages: `R_PACKAGES.txt`. Pull Drive DGE: `bash scripts/pull_drive_dge.sh`.

## Layout

```
scripts/generic_stress/                 analysis path (00–06)
scripts/shared/                         HMP helper, plot theme, ROOT, GEO helpers
scripts/pull_drive_dge*.sh              sync drive_dge/ from Google Drive
```

Local (gitignored): `drive_dge/`, `results/`, `plots/`, `raw/`.

Paths are relative to this folder (`FLEX_ROOT`). Override Drive CSVs with
`DRIVE_DGE`.

## Rules

- Drive is held out of discovery. Do not rewrite `drive_dge/`.
- Glucose is never a Drive comparison.
- Do not add more same-class GEO or rebuild the primary list as fibroblast-only unless that is a new analysis.

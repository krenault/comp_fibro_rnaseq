# Public generic cell stress vs Drive fibroblasts

Scripts only for git. Four-class public cellular stress signature (ASTRA + extra
cell GEO). Drive fibroblasts are held out of discovery, then scored with
class-holdout GSEA and per-species HMP.

**Methods:** [METHODS.md](METHODS.md). Drive counting/DESeq2 is the parent
pipeline in [krenault/comp_fibro_rnaseq](https://github.com/krenault/comp_fibro_rnaseq)
(this folder lives at `generic_stress/` there). Tables, plots, and DGE CSVs
are local and are not committed.

## Run

```bash
bash scripts/generic_stress/00_run.sh          # skip ingest
bash scripts/generic_stress/00_run.sh --ingest # re-download ASTRA + extra GEO
```

| Script | What |
|--------|------|
| `00_run.sh` | runner |
| `01_ingest.R` | ASTRA + extra cell GEO (`01_astra`, `01_geo`, `01_metadata`) |
| `02_public_hmp.R` | mixed-effect class models + HMP (pan-cell and fibroblast-only) |
| `03_drive_holdout.R` | class-holdout GSEA on Drive (both lists) |
| `04_plots.R` | figures in `plots/mammalian_generic_stress/` |
| `05_species_hmp.R` | per-species Drive HMP, conservation, ORA |
| `06_reml_celltype.R` | REML / cell-type Jaccard (skip with `SKIP_REML=1`) |
| `07_sensitivity_compare.R` | drop GSE328685 vs primary (`00_run.sh --sens-328685`) |
| `share/hmp_combine.R` | standalone HMP; copy this file to share the combination rule |

R packages: `R_PACKAGES.txt`. Pull Drive DGE: `bash scripts/pull_drive_dge.sh`.

## Layout

```
scripts/generic_stress/                 analysis path (00–07)
scripts/shared/                         plot theme, ROOT, GEO helpers (HMP wrapper)
share/hmp_combine.R                     standalone HMP, safe to send alone
scripts/pull_drive_dge*.sh              sync drive_dge/ from Google Drive
drive_dge/  results/  plots/  raw/      local only (gitignored)
```

## Rules

- Drive is held out of discovery. Do not rewrite `drive_dge/`.
- Glucose is never a Drive comparison.
- Do not add more same-class GEO or rebuild the primary list as fibroblast-only unless that is a new analysis.

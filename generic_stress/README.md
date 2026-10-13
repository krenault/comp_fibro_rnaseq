# Public generic cell stress

Four-class cellular stress signature from public data (ASTRA and additional
cell GEO: heat, hypoxia, oxidative stress, DNA damage). Fibroblast-panel
contrasts are scored after that signature is fixed, by class-holdout GSEA
and a per-species harmonic mean.

**Methods:** [METHODS.md](METHODS.md). Counting and DESeq2 for the fibroblast
panel are the parent pipeline in
[krenault/comp_fibro_rnaseq](https://github.com/krenault/comp_fibro_rnaseq).

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
| `03_drive_holdout.R` | class-holdout GSEA on the FH experiments (both lists) |
| `04_plots.R` | figures in `plots/mammalian_generic_stress/` |
| `05_species_hmp.R` | per-species FH experiment HMP, conservation, ORA |
| `06_reml_celltype.R` | REML / cell-type Jaccard (skip with `SKIP_REML=1`) |
| `07_sensitivity_compare.R` | drop GSE328685 vs primary (`00_run.sh --sens-328685`) |
| `share/hmp_combine.R` | harmonic-mean combination on a table of tests |

R packages: `R_PACKAGES.txt`.

## Layout

```
scripts/generic_stress/                 analysis path (00–07)
scripts/shared/                         plot theme, ROOT, GEO helpers (HMP wrapper)
share/hmp_combine.R                     harmonic-mean combination
```

## Design

- The fibroblast panel is not used to choose the public signature.
- Glucose is not a comparison in that scoring.
- The primary signature is the pan-cell four-class list. A fibroblast-only list is a separate run.

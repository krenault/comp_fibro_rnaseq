# Generic stress (`00_` … `06_`)

Drive is held out of discovery. Glucose is never a Drive comparison.
`drive_dge/` (or `$DRIVE_DGE`) is the held-out query and is not modified here.

```
bash scripts/generic_stress/00_run.sh              # skip ingest
bash scripts/generic_stress/00_run.sh --ingest     # re-download ASTRA + extra GEO

Rscript scripts/generic_stress/01_ingest.R
Rscript scripts/generic_stress/02_public_hmp.R
Rscript scripts/generic_stress/03_drive_holdout.R
Rscript scripts/generic_stress/04_plots.R
Rscript scripts/generic_stress/05_species_hmp.R
Rscript scripts/generic_stress/06_reml_celltype.R  # optional
```

Outputs: `results/mammalian_generic_stress/` and `plots/mammalian_generic_stress/`.

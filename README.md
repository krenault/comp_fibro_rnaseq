# drive_rnaseq

Private pipeline: StringTie GTFs → official `prepDE.py` → one-to-one ENSG counts → DESeq2.

Production DGE is **pairwise only**: one DESeq2 fit per contrast (e.g. 2.5 mM vs 8 mM on those samples alone; a separate fit for 30 mM vs 8 mM). Design is `~ individual + treatment`, gene filter **≥5 counts in ≥3 samples within the contrast**, NAs filled with **0**.

`--mode multi` is an optional comparison (all three levels of an axis in one fit, which is what Dylan used). It is **not** the Drive / paper analysis.

## Layout

```
scripts/prepDE.py            official StringTie helper (not ours)
scripts/process_stringtie.sh copy GTFs from Unity, run prepDE
scripts/prepare_counts.py    stack species, collapse orthologs, TPM sample filter, NA mask
scripts/run_deseq2.R         DESeq2 (default `--mode pairwise`)
scripts/validate.py          ortholog-merge sanity checks
```

No count matrices are stored in git. Point `DATA_DIR` at local data.

## Setup

```bash
cp config.example.env .env   # then edit; or export the same variables
# needs: python3 + pandas/numpy, R + DESeq2/tidyverse, ssh access to Unity
```

`config.example.env`:

- `DATA_DIR` — local directory for per-species `gene_count_matrix.csv` and the combined matrix
- `TPM` — one-to-one ortholog TPM CSV (first column = ENSG IDs); used to choose samples and mask orthology gaps
- `PREPDE_LENGTH` — fragment length passed to prepDE (`75` matches the Drive production run, which omitted `-l` and used the prepDE default)

## Run

```bash
set -a && source .env && set +a

# 1. StringTie GTFs → gene counts (read-only on Unity)
bash scripts/process_stringtie.sh

# 2. Combine, collapse ENSG orthologs, filter to TPM samples, mask zeros→NA
python3 scripts/prepare_counts.py \
  --counts-dir "$DATA_DIR" \
  --tpm "$TPM" \
  --out "$DATA_DIR/all_species_all_treatments_one2one_filtered_with_NA.csv"

python3 scripts/validate.py \
  --final "$DATA_DIR/all_species_all_treatments_one2one_filtered_with_NA.csv" \
  --tpm "$TPM"

# 3. DESeq2 (production settings)
Rscript scripts/run_deseq2.R \
  --counts "$DATA_DIR/all_species_all_treatments_one2one_filtered_with_NA.csv" \
  --outdir "$DATA_DIR/DESeq_results_singleside" \
  --mode pairwise \
  --na fill0
```

### NA treatments

`--na` is applied **within the samples of each fit**, then the count filter runs:

| mode | behavior |
|---|---|
| `fill0` | NA → 0 (production; DESeq2 cannot take NA) |
| `drop_mixed` | drop genes that are NA in some but not all samples; remaining NA → 0 |
| `drop_any` | drop any gene with ≥1 NA in that fit |

On the production glucose contrasts, **all three modes give identical DEG lists** (Jaccard = 1). Every NA is species-wide (the gene is NA in all samples of that contrast), so those genes already fail the ≥5-in-≥3 filter after `fill0`. There were **zero mixed-NA genes** (NA in some samples but not others). See `results/na_treatment_deseq.csv`.

### Optional: three-level (multi) model

Not production. Fits 2.5 / 8 / 30 mM together (or the three hypoxia / temperature levels), then extracts the two contrasts. Dylan’s DESeq used this form.

```bash
Rscript scripts/run_deseq2.R --counts ... --outdir ... --mode multi --na fill0
```

## Sample names

`{species}_{individual}_{temp}_{glucose}_{hypoxia}`  
example: `squirrel_OK14_37C_2.5mM_0`

Axes hold the other treatments at baseline: glucose at 37C / 0 hypoxia; hypoxia at 37C / 8mM; temperature at 8mM / 0 hypoxia.

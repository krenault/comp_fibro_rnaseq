# comp_fibro_rnaseq

StringTie GTFs → official `prepDE.py` → one-to-one ENSG counts → DESeq2.

Production DGE is **pairwise only**: one DESeq2 fit per contrast (e.g. 2.5 mM vs 8 mM on those samples alone; a separate fit for 30 mM vs 8 mM). Design is `~ individual + treatment`, gene filter **≥5 counts in ≥3 samples within the contrast.**

## Layout

```
scripts/prepDE.py            official StringTie helper 
scripts/process_stringtie.sh copy GTFs from Unity, run prepDE
scripts/prepare_counts.py    stack species, collapse orthologs, TPM sample filter, NA mask
scripts/run_deseq2.R         DESeq2
```

## Run

```bash
set -a && source .env && set +a

# 1. 
bash scripts/process_stringtie.sh

# 2. 
python3 scripts/prepare_counts.py \
  --counts-dir "$DATA_DIR" \
  --tpm "$TPM" \
  --out "$DATA_DIR/all_species_all_treatments_one2one_filtered_with_NA.csv"

python3 scripts/validate.py \
  --final "$DATA_DIR/all_species_all_treatments_one2one_filtered_with_NA.csv" \
  --tpm "$TPM"

# 3. 
Rscript scripts/run_deseq2.R \
  --counts "$DATA_DIR/all_species_all_treatments_one2one_filtered_with_NA.csv" \
  --outdir "$DATA_DIR/DESeq_results_singleside" \
  --mode pairwise \
  --na fill0
```
## Sample names

`{species}_{individual}_{temp}_{glucose}_{hypoxia}`  
example: `squirrel_OK14_37C_2.5mM_0`

Axes hold the other treatments at baseline: glucose at 37C / 0 hypoxia; hypoxia at 37C / 8mM; temperature at 8mM / 0 hypoxia.

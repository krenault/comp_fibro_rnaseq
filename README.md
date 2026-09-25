# comp_fibro_rnaseq

Reproducible RNA-seq count prep and DESeq2 for a multi-species fibroblast culture panel (glucose / temperature / hypoxia).

**Pipeline:** StringTie GTFs → official `prepDE.py` → one-to-one human ortholog (ENSG) counts → pairwise DESeq2.

Production DGE is **pairwise**: one DESeq2 fit per contrast (e.g. 2.5 mM vs 8 mM on those samples alone). Design is `~ individual + treatment`. Gene filter: **≥5 counts in ≥3 samples within the contrast.**

## Layout

```
scripts/prepDE.py            StringTie helper (Pertea lab)
scripts/process_stringtie.sh copy GTFs from a remote host, run prepDE
scripts/prepare_counts.py    stack species, collapse orthologs, TPM sample filter
scripts/run_deseq2.R         pairwise DESeq2
config.example.env           copy to .env (gitignored)
```

## Setup

```bash
cp config.example.env .env
# edit DATA_DIR, TPM, HOST, BASE
set -a && source .env && set +a
```

Needs: Python 3 + pandas; R with `DESeq2`, `tidyverse`; SSH access to the GTF host for step 1.

## Run

```bash
# 1. Remote GTFs → local gene_count_matrix.csv per species × experiment
bash scripts/process_stringtie.sh

# 2. Combine species, collapse to one-to-one ENSG, intersect TPM samples
python3 scripts/prepare_counts.py \
  --counts-dir "$DATA_DIR" \
  --tpm "$TPM" \
  --out "$DATA_DIR/all_species_all_treatments_one2one_filtered.csv"

# 3. Pairwise DESeq2
Rscript scripts/run_deseq2.R \
  --counts "$DATA_DIR/all_species_all_treatments_one2one_filtered.csv" \
  --outdir "$DATA_DIR/DESeq_results_pairwise" \
  --mode pairwise
```

## Sample names

`{species}_{individual}_{temp}_{glucose}_{hypoxia}`  
example: `squirrel_OK14_37C_2.5mM_0`

Axes hold other treatments at baseline: glucose at 37C / no hypoxia; hypoxia at 37C / 8 mM; temperature at 8 mM / no hypoxia.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use this pipeline in a paper, cite the associated study and this repository URL.

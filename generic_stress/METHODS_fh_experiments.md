# FH experiments — RNA-seq — counting, ortholog merge, pairwise DESeq2

This is **upstream of** the generic-stress analysis. Production scripts live
in the parent repository:
[krenault/comp_fibro_rnaseq](https://github.com/krenault/comp_fibro_rnaseq)
(`scripts/process_stringtie.sh`, `prepare_counts.py`, `run_deseq2.R`).

`drive_dge/` holds the pairwise DESeq2 tables for the FH experiments (or set
`DRIVE_DGE`). Scripts in `scripts/generic_stress/` do not overwrite them.

---

## Design used by the generic-stress query

- Pairwise DESeq2, one fit per contrast; design `~ individual + treatment`
- Gene filter: ≥5 counts in ≥3 samples **within the contrast**
- Human one-to-one orthologs (ENSG, version stripped)
- Axes: glucose at 37°C / no hypoxia; hypoxia at 37°C / 8 mM; temperature at 8 mM / no hypoxia
- **Generic-stress scoring drops glucose** (2.5 mM and 30 mM vs 8 mM)

Contrasts kept: heat 41°C vs 37°C, cold 32°C vs 37°C, hypoxia 6 h and 24 h vs 0 h.

This analysis builds on [dbioinfo/fibro](https://github.com/dbioinfo/fibro).
Production DGE here is prepDE, not featureCounts.

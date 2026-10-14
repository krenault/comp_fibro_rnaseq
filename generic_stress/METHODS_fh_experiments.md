# FH experiments — RNA-seq — counting, ortholog merge, pairwise DESeq2

This is **upstream of** the generic-stress analysis. Production scripts live
in the parent repository:
[krenault/comp_fibro_rnaseq](https://github.com/krenault/comp_fibro_rnaseq)
(`fh_dge/scripts/process_stringtie.sh`, `prepare_counts.py`, `run_deseq2.R`).

The scripts write DESeq2 tables locally under `fh_dge/results/`. Scripts in
`scripts/generic_stress/` read that folder, or `$DRIVE_DGE`, and do not overwrite it.

- `fh_dge/results/pairwise/` — one fit per contrast. This is the production set. The
  generic-stress scripts read it (or `$DRIVE_DGE`).
- `fh_dge/results/all_levels/` — one fit per axis, with every level of that axis in
  the same model.
- `fh_dge/results/hypoxia_6H_vs_24H_pairwise/` — hypoxia 6 h versus 24 h as its own
  fit, using only those two time points. Not the all-levels model, and not
  used in the generic-stress score.
- `fh_dge/results/hypoxia_6H_vs_24H_pairwise_fgsea/` — FGSEA of that contrast.
  Ranking is `sign(log2FC) * -log10(p)`, same gene sets as the other FH
  experiment enrichments (Hallmark, Reactome, KEGG, GO biological process,
  PANTHER). A positive NES is enriched among genes higher at 6 h.

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

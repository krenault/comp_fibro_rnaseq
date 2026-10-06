# Drive RNA-seq — counting, ortholog merge, pairwise DESeq2

This is **upstream of** the generic-stress analysis. Production scripts live
in the parent repository:
[krenault/comp_fibro_rnaseq](https://github.com/krenault/comp_fibro_rnaseq)
(`scripts/process_stringtie.sh`, `prepare_counts.py`, `run_deseq2.R`).

`drive_dge/` here is a copy of those pairwise DESeq2 tables (or set
`DRIVE_DGE`). Scripts in `scripts/generic_stress/` never overwrite it.

---

## Design used by the generic-stress query

- Pairwise DESeq2, one fit per contrast; design `~ individual + treatment`
- Gene filter: ≥5 counts in ≥3 samples **within the contrast**
- Human one-to-one orthologs (ENSG, version stripped)
- Axes: glucose at 37°C / no hypoxia; hypoxia at 37°C / 8 mM; temperature at 8 mM / no hypoxia
- **Generic-stress scoring drops glucose** (2.5 mM and 30 mM vs 8 mM) and drops
  rousette (no temperature or hypoxia files)

Contrasts kept: heat 41°C vs 37°C, cold 32°C vs 37°C, hypoxia 6 h and 24 h vs 0 h.

Pull a shared copy: `bash scripts/pull_drive_dge.sh`
([Google Drive folder](https://drive.google.com/drive/folders/1-cCvd4b5SyubpE2uWA3Kl3Eg4sMb3huo)).

Dylan’s public pipeline is [dbioinfo/fibro](https://github.com/dbioinfo/fibro).
Production DGE here is prepDE, not featureCounts.

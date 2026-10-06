# Methods — public generic cell-stress signature vs Drive fibroblasts

Scripts: `scripts/generic_stress/` (`00_`–`06_`).
Figures: `plots/mammalian_generic_stress/`.
Tables: `results/mammalian_generic_stress/`.

Drive DGE production (StringTie / prepDE / DESeq2) is the parent pipeline
([krenault/comp_fibro_rnaseq](https://github.com/krenault/comp_fibro_rnaseq);
[METHODS_drive_dge.md](METHODS_drive_dge.md)). Query tables: `drive_dge/` or
`$DRIVE_DGE`.

The design is: build a four-class cellular stress signature from public data,
hold Drive out of discovery, then ask which Drive responses recover that
signature and which do not.

---

## 1. Design

Public discovery and Drive are independent experiments. Drive libraries are
never used to define the generic list. Glucose is never a Drive comparison.
Rousette is dropped (no temperature or hypoxia files). In-vivo GEO is not in
discovery.

The public panel has heat and hypoxia as *classes*. Drive heat or hypoxia
recovering the *full* public list is therefore a matched-class positive
control, not evidence of a class-independent CSR. Honest tests drop the
matched class from the gene set before scoring Drive (section 6).

Join key throughout: human Ensembl gene ID (ENSG, version stripped). Mouse
GEO is mapped to human one-to-one orthologs before meta-analysis.

---

## 2. Drive query (held out)

Cultured dermal fibroblasts from 13 mammals; 12 species after dropping
rousette (human, gelada, rat, 13-lined ground squirrel, little brown bat,
honey badger, seal, dolphin, whale, rhino, dromedary camel, bactrian camel).
DESeq2 pairwise, design `~ individual + treatment`, human one-to-one
orthologs. Contrasts used here:

- heat: 41°C vs 37°C
- cold: 32°C vs 37°C
- hypoxia: 6 h and 24 h vs 0 h

Glucose (2.5 mM and 30 mM vs 8 mM) is excluded from every Drive axis in this
analysis. Files: `drive_dge/` (48 non-glucose, non-rousette contrasts).

---

## 3. Public discovery panel

Four stressor classes: heat, hypoxia, oxidative (H₂O₂), DNA damage (UV / IR).

**ASTRA.** Human cell-stress atlas, Zenodo 15885686 (Quintanilla et al.).
Tables `datasets.tsv` and `differential_expression.tsv`. Restricted to the
four classes above. Each ASTRA `record_id` is one contrast (dose / time / cell
type kept uncollapsed at ingest, then collapsed by GSE in the mixed model).

**Extra GEO.** Targeted fill for mouse cells and a few human fibroblast
series with usable counts or abundance tables. This was a catalog plus
manual GEO search, not an API sweep. Series that were Drive itself, in vivo,
cold, OGD/R, or otherwise off-class were recorded as rejected in
`results/mammalian_generic_stress/other_species_candidates.tsv`.

Locked discovery size: **108 contrasts, 59 GSE, 86 ASTRA + 22 GEO; 90 human
/ 18 mouse**. Class split: oxidative 38, heat 27, hypoxia 22, DNA damage 21.
About one fifth of units are fibroblast-like; the primary signature is not
restricted to fibroblasts (hypoxia would drop to three GSE).

Author-supplied DE-only tables and n = 1 abundance contrasts have their
standard errors inflated 1.5× so they cannot dominate ASTRA.

---

## 4. Effect sizes and mixed-effects meta-analysis

Adapted from Tyshkovskiy, Gladyshev et al., *Cell* 2023 (mammalian aging
meta-analysis), applied to cellular stress.

1. Per contrast, effect size (EC) is log2 fold-change. Standard error is
   taken from the DE table when present, otherwise backed out from the
   two-sided *p*-value via the normal quantile; if neither exists, a
   conservative MAD-based SE is used.
2. Pairwise Spearman correlations on the union of each contrast’s top 250
   genes (denoised similarity; diagnostic, not the signature rule).
3. Multiple Deming regression rescales ECs across contrasts (L-BFGS-B, 10
   random starts).
4. Per gene, two-stage random-effects: inverse-variance collapse of records
   within GSE, then `metafor::rma` (REML, DerSimonian–Laird or equal-effects
   fallback) across GSE. A gene needs ≥3 contrasts in a class model.
5. Class signatures (heat, hypoxia, oxidative, DNA damage), plus human-only,
   mouse-only, and global models, are genes with BH-adjusted mixed-model
   *p* < 0.05.

The *Cell* paper’s full `rma.mv` with crossed tissue/species intercepts is
too slow for ~15k genes × 15 signatures. Univariate REML on GSE-level
estimates still uses SE. REML vs DL concordance and fibroblast/cancer
Jaccard overlap are optional sensitivity (`scripts/generic_stress/06_reml_celltype.R`;
figure `reml_celltype_overlap`).

---

## 5. Generic signature (harmonic mean *p*)

HMP is the harmonic mean *p*-value (Wilson, *PNAS* 2019), used as in the
*Cell* 2023 STAR Methods. It combines the four *class* mixed-model
*p*-values for each gene, not the raw contrast *p*-values.

Direction-aware two-sided combination: one-sided *p* for up and for down are
formed from each class *p* and the sign of that class EC; HMP is computed on
each side; the reported two-sided *p* is twice the smaller of the two
(capped at 1). Implementation: `harmonicmeanp::p.hmp` when installed,
otherwise a Landau-tail approximation (`scripts/shared/hmp_twosided.R`).

A gene is **generic** if:

- BH-adjusted HMP < 0.05
- measured in ≥ 3 of 4 classes
- ≥ 80% of those classes agree on sign

Locked list: **1,344 genes (737 up, 607 down)**. A stricter both-species
list (independent human HMP and mouse HMP, same direction) is 104 genes.

This is not “genes that move in every stressor.” Class-restricted landmarks
(HSPA1A/B, HIF1A, HMOX1, VEGFA, FASN) fail the HMP / sign-agreement rule
and are not labelled generic.

---

## 6. Pathway analysis of the public generic list

Genes are ranked by −log10(HMP) × sign(mean class EC). GSEA
(`fgseaMultilevel`) uses MSigDB Hallmark, Reactome, and GO biological
process (min size 15, max 400). Table:
`results/mammalian_generic_stress/cell_meta_gsea_hallmark_reactome_gobp.csv`.
Locked figure `public_pathways` plots the three collections in one ranking (shape =
collection); redundant cell-cycle GO terms are collapsed for display.

---

## 7. Honest Drive scoring (class hold-out)

Leave-one-class-out HMP, using the already-fit class mixed models (no
refitting):

| Gene set | Classes combined | Min classes | Scored on |
|----------|------------------|-------------|-----------|
| full | heat + hypoxia + oxidative + DNA damage | 3 | Drive cold |
| no_heat | hypoxia + oxidative + DNA damage | 2 | Drive heat |
| no_hypoxia | heat + oxidative + DNA damage | 2 | Drive hypoxia 6 h and 24 h |
| unmatched | oxidative + DNA damage only | 2 | all Drive exposures (shared core) |

Same BH HMP < 0.05 and ≥ 80% sign as the primary list.

Each Drive species × exposure is ranked by −log10(*p*) × sign(LFC) and
tested by GSEA against the matching hold-out up/down sets. Call: NES > 0 and
FDR < 0.05 = induced the public up-set; NES < 0 and FDR < 0.05 = opposite.

Cold is the clean test (no public cold in discovery). Heat after dropping
public heat, and hypoxia after dropping public hypoxia, are the honest
matched-class tests. Unmatched (oxidative + DNA damage only) is a sensitivity
that holds the gene set fixed across every Drive contrast; it is not the
primary signature.

---

## 8. Conserved Drive genes vs the public generic

Within each Drive exposure, a gene is conserved if |median LFC| ≥ 0.25
across species, called in ≥ 5 species, and ≥ 80% of those species share
direction. Overlap with the public 1,344 is the generic slice; the remainder
is Drive-residual. Figure `conserved_overlap`.

---

## 9. Per-species Drive cross-perturbation HMP

Same combination rule as section 5, applied *within one Drive species*
across heat, cold, and hypoxia (24 h preferred; 6 h used only if 24 h is
missing for that gene). Glucose dropped.

A gene is a species cross-perturbation hit if BH HMP < 0.05, ≥ 3 arms
measured, and ≥ 80% sign agreement.

Conservation of those hits across the 12 Drive species: count how many
species call the gene, and require ≥ 80% of those calls to share direction.
Headline conserved core: **≥ 5 species**. Over-representation
(`fgsea::fora`) of that core, split up/down, uses Hallmark + Reactome +
GO BP with universe = genes tested in the species HMP.

6 h and 24 h hypoxia are not independent classes, so the species HMP uses
one hypoxia arm (24 h preferred).

---

## 10. Figures

`plots/mammalian_generic_stress/`

| File | Question |
|------|----------|
| discovery_classes | What is in public discovery? |
| public_pathways | Pathways of the 1,344-gene generic list |
| gene_universe | Gene-universe overlap with Drive (not an ID problem) |
| conserved_overlap | Conserved Drive genes: generic vs residual |
| holdout_scorecard | Drive vs class-held-out public up-genes. Strip (colour) = which public classes enter that HMP set: heat vs hypoxia+H2O2+UV; hypoxia vs heat+H2O2+UV; cold vs all four (includes heat) |
| unmatched_ox_dna_scorecard | Same H2O2+UV gene set scored on every Drive exposure |
| cold_forest | GSEA NES of Drive cold vs the 737 public generic-up genes |
| human_vs_others_holdout_nes | GSEA NES of each Drive contrast vs the held-out public generic-up set; human vs the other 11 mammals |
| reml_celltype_overlap | Jaccard of alternative generic lists vs the 1,344 |
| volcano_<species> | Drive volcano per species. Labels = DESeq2 padj < 0.05 AND public generic |
| species_hmp_vs_public | Numerator = genes in the public 1,344; denominator = that species' heat+cold+hypoxia HMP list |
| species_hmp_conserved_ora | ORA of the ≥5-species core |

House style: `theme_pub`, `treatment_colors`, blue–white–red.

---

## 11. What this analysis does not do

- Drive is not used to discover the generic list.
- Glucose is not a Drive comparison.
- In-vivo GEO is not in discovery.
- The primary signature is not fibroblast-only.
- Vote / median-z is not the paper method.
- Adding more same-class GEO, or pulling in ER stress, would rewrite the
  1,344-gene list; that is a new analysis, not a sensitivity of this one.

---

## 12. How to run

```bash
# Skip ingest if ASTRA + extra GEO are already on disk
bash scripts/generic_stress/00_run.sh

# Re-download ASTRA / extra GEO
bash scripts/generic_stress/00_run.sh --ingest
```

`01` ingest, `02` signature, `03` Drive query, `04` figures, `05`
per-species Drive HMP. Optional: `scripts/generic_stress/06_reml_celltype.R`.

Paths are relative to this analysis folder. Override with `FLEX_ROOT` and
`DRIVE_DGE` (`config.example.env`). Software: see `R_PACKAGES.txt`
(tidyverse, DESeq2, metafor, fgsea, msigdbr, harmonicmeanp, AnnotationDbi /
org.Hs.eg.db, babelgene).

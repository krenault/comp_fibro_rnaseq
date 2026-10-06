#!/usr/bin/env python3
"""Publication metadata for mammalian generic-stress discovery.

One row per analysis unit used in 02_public_hmp.R (ASTRA experiment_id, or GEO driver).
Drive is held out. Other-species candidates are written separately.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ["FLEX_ROOT"]) if os.environ.get("FLEX_ROOT") else Path(__file__).resolve().parents[2]
OUT = ROOT / "results/mammalian_generic_stress"
ASTRA = ROOT / "results/astra/comparison_catalog.csv"


def lineage(cell_type: str) -> str:
    t = (cell_type or "").lower()
    rules = [
        ("fibroblast-derived neuron", "neuron (iPSC/fibroblast-derived)"),
        ("cortical brain organoids", "organoid (neural)"),
        ("ipsc-3d neural", "organoid (neural)"),
        ("tissue-engineered cartilage", "organoid (cartilage)"),
        ("astrocyte", "astrocyte"),
        ("podocyte", "podocyte"),
        ("pancreatic_endocrine", "endocrine"),
        ("htert-hpne", "epithelial (pancreas)"),
        ("c2c12", "myoblast"),
        ("mesc", "pluripotent stem"),
        ("esc", "pluripotent stem"),
        ("nsc", "neural stem"),
        ("msc", "mesenchymal stem"),
        ("cd34", "hematopoietic stem"),
        ("wi-38", "fibroblast"),
        ("mrc5", "fibroblast"),
        ("imr90", "fibroblast"),
        ("hdf", "fibroblast"),
        ("hlf", "fibroblast"),
        ("hpf", "fibroblast"),
        ("nih3t3", "fibroblast (immortalized)"),
        ("mef", "fibroblast (primary/MEF)"),
        ("fibroblast", "fibroblast"),
        ("huvec", "endothelial"),
        ("haec", "endothelial"),
        ("ibmec", "endothelial"),
        ("pcs-100-013", "endothelial"),
        ("rpe", "epithelial (RPE)"),
        ("arpe", "epithelial (RPE)"),
        ("hela", "epithelial (cancer)"),
        ("mcf7", "epithelial (cancer)"),
        ("hct116", "epithelial (cancer)"),
        ("mda-mb-231", "epithelial (cancer)"),
        ("lncap", "epithelial (cancer)"),
        ("bewo", "epithelial (cancer)"),
        ("pterygium", "epithelial (cancer)"),
        ("hek293", "epithelial (embryonal kidney / cancer)"),
        ("rcc4", "epithelial (cancer)"),
        ("u2os", "osteosarcoma"),
        ("tk6", "lymphoblastoid"),
        ("lcl", "lymphoblastoid"),
        ("ramos", "lymphoid (cancer)"),
        ("tib-152", "lymphoid (cancer)"),
        ("atcc-bxs0114", "mesenchymal"),
    ]
    for key, lab in rules:
        if key in t:
            return lab
    return "other"


def cell_state_label(group: str) -> str:
    g = (group or "").lower()
    if "cancer" in g:
        return "cancer"
    if "organoid" in g:
        return "organoid"
    if "non-diseased" in g or "reference" in g:
        return "non-diseased"
    return group or ""


def uniq_join(vals) -> str:
    seen = []
    for v in vals:
        v = (v or "").strip()
        if v and v not in seen:
            seen.append(v)
    return "; ".join(seen)


FIELDS = [
    "experiment_id",
    "in_discovery",
    "source",
    "gse_id",
    "pmid",
    "species",
    "common_name",
    "stress_class",
    "perturbation",
    "cell_type",
    "cell_lineage",
    "tissue_of_origin",
    "cell_state",
    "genotype",
    "is_nondiased_wt",
    "is_fibroblast_like",
    "n_records_collapsed",
    "n_controls",
    "n_treated",
    "n_samples",
    "de_method",
    "notes",
]


def row(**kw):
    out = {k: kw.get(k, "") for k in FIELDS}
    for k in ("n_records_collapsed", "n_controls", "n_treated", "n_samples"):
        if out[k] is None:
            out[k] = ""
        else:
            out[k] = str(out[k])
    ct = str(out["cell_type"])
    lin = out["cell_lineage"] or lineage(ct)
    out["cell_lineage"] = lin
    out["is_fibroblast_like"] = "TRUE" if lin.startswith("fibroblast") else "FALSE"
    out["in_discovery"] = str(out["in_discovery"]).upper() if out["in_discovery"] != "" else "TRUE"
    return out


def astra_rows():
    with ASTRA.open() as f:
        recs = list(csv.DictReader(f))
    by = defaultdict(list)
    for r in recs:
        by[r["experiment_id"]].append(r)
    out = []
    for eid, rs in sorted(by.items()):
        n_ctrl = sum(int(r["No_of_controls"] or 0) for r in rs)
        n_trt = sum(int(r["No_of_treated"] or 0) for r in rs)
        n_samp = sum(int(r["No_of_samples"] or 0) for r in rs)
        wt = all(r["is_nondiased_wt"] == "TRUE" for r in rs)
        out.append(
            row(
                experiment_id=eid,
                in_discovery=True,
                source="ASTRA",
                gse_id=rs[0]["GEO_id"],
                pmid="" if rs[0]["PMID"] in ("None", "NA", "") else rs[0]["PMID"],
                species="Homo sapiens",
                common_name="human",
                stress_class=rs[0]["stress_class"],
                perturbation=uniq_join(r["treatment_type"] for r in rs),
                cell_type=rs[0]["Cell_type"],
                tissue_of_origin=rs[0]["Tissue"],
                cell_state=cell_state_label(rs[0]["Cell_group"]),
                genotype=rs[0]["Genotype"],
                is_nondiased_wt="TRUE" if wt else "FALSE",
                n_records_collapsed=len(rs),
                n_controls=n_ctrl,
                n_treated=n_trt,
                n_samples=n_samp,
                de_method="ASTRA edgeR (Zenodo DE table)",
                notes="Time/dose records collapsed to median z within experiment_id"
                if len(rs) > 1
                else "",
            )
        )
    return out


GEO_DISC = [
    dict(
        experiment_id="human_GSE281533_hyp",
        gse_id="GSE281533",
        species="Homo sapiens",
        common_name="human",
        stress_class="hypoxia",
        perturbation="hypoxia vs normoxia (WT, ignore DPP4-KD arms)",
        cell_type="HLF lung fibroblasts",
        cell_lineage="fibroblast",
        tissue_of_origin="Lung",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
        notes="Not in ASTRA; fibroblast hypoxia match for Drive",
    ),
    dict(
        experiment_id="human_GSE139963_hyp",
        gse_id="GSE139963",
        species="Homo sapiens",
        common_name="human",
        stress_class="hypoxia",
        perturbation="hypoxia vs PBS (no TGFβ arms)",
        cell_type="HPF pulmonary fibroblasts",
        cell_lineage="fibroblast",
        tissue_of_origin="Lung",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="log2FPKM_LFC",
        notes="Weaker than GSE281533; FPKM LFC only",
    ),
    dict(
        experiment_id="human_GSE225095_H2O2",
        gse_id="GSE225095",
        species="Homo sapiens",
        common_name="human",
        stress_class="oxidative",
        perturbation="H2O2 vs proliferating (RNA-seq arms)",
        cell_type="IMR90",
        cell_lineage="fibroblast",
        tissue_of_origin="Lung",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
        notes="Not in ASTRA; RNA arms only (not Ribo-seq)",
    ),
    dict(
        experiment_id="human_GSE262772_20Gy",
        gse_id="GSE262772",
        species="Homo sapiens",
        common_name="human",
        stress_class="dna_damage",
        perturbation="20 Gy gamma vs control, 24h",
        cell_type="HDF primary dermal fibroblasts",
        cell_lineage="fibroblast",
        tissue_of_origin="Skin",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
        notes="Same series as held-out bat radiation; human/mouse used in discovery",
    ),
    dict(
        experiment_id="mouse_GSE98906_heat",
        gse_id="GSE98906",
        species="Mus musculus",
        common_name="mouse",
        stress_class="heat",
        perturbation="heat shock vs untreated",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="log2RPKM_LFC",
    ),
    dict(
        experiment_id="mouse_GSE98906_H2O2",
        gse_id="GSE98906",
        species="Mus musculus",
        common_name="mouse",
        stress_class="oxidative",
        perturbation="H2O2 vs untreated",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="log2RPKM_LFC",
    ),
    dict(
        experiment_id="mouse_GSE65636_HS45",
        gse_id="GSE65636",
        species="Mus musculus",
        common_name="mouse",
        stress_class="heat",
        perturbation="HS 45 min vs control (no tunicamycin)",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
    ),
    dict(
        experiment_id="mouse_GSE65636_HS135",
        gse_id="GSE65636",
        species="Mus musculus",
        common_name="mouse",
        stress_class="heat",
        perturbation="HS 135 min vs control (no tunicamycin)",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
    ),
    dict(
        experiment_id="mouse_GSE262772_20Gy",
        gse_id="GSE262772",
        species="Mus musculus",
        common_name="mouse",
        stress_class="dna_damage",
        perturbation="20 Gy gamma vs control, 24h",
        cell_type="primary fibroblasts",
        cell_lineage="fibroblast",
        tissue_of_origin="Skin",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
    ),
    dict(
        experiment_id="mouse_GSE215293_heat",
        gse_id="GSE215293",
        species="Mus musculus",
        common_name="mouse",
        stress_class="heat",
        perturbation="43C 1h vs 37C (WT MEF)",
        cell_type="MEF",
        cell_lineage="fibroblast (primary/MEF)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="log2FPKM_LFC",
        notes="WT arms only",
    ),
    dict(
        experiment_id="mouse_GSE243284_heat",
        gse_id="GSE243284",
        species="Mus musculus",
        common_name="mouse",
        stress_class="heat",
        perturbation="45C 10 min vs no HS (siCtrl)",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="siCtrl",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
        n_controls=2,
        n_treated=2,
        n_samples=4,
        notes="siCtrl only; knockdown arms unused",
    ),
    dict(
        experiment_id="mouse_GSE197536_HS_2S",
        gse_id="GSE197536",
        species="Mus musculus",
        common_name="mouse",
        stress_class="heat",
        perturbation="HS_2S vs HS control",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2 (RSEM expected counts)",
    ),
    dict(
        experiment_id="mouse_GSE197536_HS_8M",
        gse_id="GSE197536",
        species="Mus musculus",
        common_name="mouse",
        stress_class="heat",
        perturbation="HS_8M vs HS control",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2 (RSEM expected counts)",
    ),
    dict(
        experiment_id="mouse_GSE197536_H2O2_2h",
        gse_id="GSE197536",
        species="Mus musculus",
        common_name="mouse",
        stress_class="oxidative",
        perturbation="H2O2 2h vs control",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2 (RSEM expected counts)",
    ),
    dict(
        experiment_id="mouse_GSE197536_H2O2_7h",
        gse_id="GSE197536",
        species="Mus musculus",
        common_name="mouse",
        stress_class="oxidative",
        perturbation="H2O2 7h vs control",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2 (RSEM expected counts)",
    ),
    dict(
        experiment_id="mouse_GSE281390_UVA",
        gse_id="GSE281390",
        species="Mus musculus",
        common_name="mouse",
        stress_class="dna_damage",
        perturbation="UVA 10 J/cm2 vs mock",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="log2FPKM_LFC",
    ),
    dict(
        experiment_id="mouse_GSE66286_UVC",
        gse_id="GSE66286",
        species="Mus musculus",
        common_name="mouse",
        stress_class="dna_damage",
        perturbation="UVC 80 J/m2, 6h vs 0h",
        cell_type="NIH 3T3",
        cell_lineage="fibroblast (immortalized)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="author_DESeq (significant genes only)",
        notes="Author table is not genome-wide; ~1800 genes",
    ),
    dict(
        experiment_id="mouse_GSE114086_GM_hyp",
        gse_id="GSE114086",
        species="Mus musculus",
        common_name="mouse",
        stress_class="hypoxia",
        perturbation="2% O2 vs 21% in growth medium",
        cell_type="C2C12",
        cell_lineage="myoblast",
        tissue_of_origin="Skeletal muscle",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="log2Cuff_LFC",
        notes="Cufflinks fractional counts; not DESeq2",
    ),
    dict(
        experiment_id="mouse_GSE114086_DM_hyp",
        gse_id="GSE114086",
        species="Mus musculus",
        common_name="mouse",
        stress_class="hypoxia",
        perturbation="2% O2 vs 21% in differentiation medium",
        cell_type="C2C12",
        cell_lineage="myoblast",
        tissue_of_origin="Skeletal muscle",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="log2Cuff_LFC",
    ),
    dict(
        experiment_id="mouse_GSE171871_mESC_hyp",
        gse_id="GSE171871",
        species="Mus musculus",
        common_name="mouse",
        stress_class="hypoxia",
        perturbation="hypoxia H4 vs normoxia N4",
        cell_type="mESC",
        cell_lineage="pluripotent stem",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        de_method="DESeq2",
    ),
    dict(
        experiment_id="mouse_GSE205381_hyp",
        gse_id="GSE205381",
        species="Mus musculus",
        common_name="mouse",
        stress_class="hypoxia",
        perturbation="hypoxia 10h vs control",
        cell_type="MEF",
        cell_lineage="fibroblast (primary/MEF)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        n_controls=1,
        n_treated=1,
        n_samples=2,
        de_method="log2abund_n1",
        notes="Unreplicated; keep with a sensitivity flag",
    ),
    dict(
        experiment_id="mouse_GSE205381_ox",
        gse_id="GSE205381",
        species="Mus musculus",
        common_name="mouse",
        stress_class="oxidative",
        perturbation="oxidative 3h vs control",
        cell_type="MEF",
        cell_lineage="fibroblast (primary/MEF)",
        tissue_of_origin="Embryo",
        cell_state="non-diseased",
        genotype="WT",
        is_nondiased_wt="TRUE",
        n_controls=1,
        n_treated=1,
        n_samples=2,
        de_method="log2abund_n1",
        notes="Unreplicated; keep with a sensitivity flag",
    ),
]


CAND_FIELDS = [
    "gse_id",
    "species",
    "common_name",
    "cell_type",
    "stress_class",
    "system",
    "recommendation",
    "already_processed",
    "reason",
]


CANDIDATES = [
    dict(
        gse_id="GSE262772",
        species="Myotis lucifugus; Eptesicus fuscus; Eonycteris spelaea; Artibeus jamaicensis",
        common_name="4 bats",
        cell_type="primary fibroblasts",
        stress_class="dna_damage",
        system="cell culture",
        recommendation="hold_out_loso",
        already_processed="TRUE",
        reason="Best extra-species cell RNA-seq we already have. Human+mouse of this series are in discovery; bats stay held out as a leave-one-clade universality check (same design as Drive fibroblasts).",
    ),
    dict(
        gse_id="Drive panel",
        species="Rattus norvegicus + 12 other mammals",
        common_name="rat and Drive species",
        cell_type="primary fibroblasts",
        stress_class="heat; cold; hypoxia; glucose",
        system="cell culture",
        recommendation="hold_out_query",
        already_processed="TRUE",
        reason="Drive is the query, not discovery. Dumping Drive rat into the meta would circularize the reviewer test.",
    ),
    dict(
        gse_id="GSE274106",
        species="Canis lupus familiaris",
        common_name="dog",
        cell_type="DH82 macrophage-like",
        stress_class="heat",
        system="cell culture",
        recommendation="candidate_if_counts",
        already_processed="FALSE",
        reason="Only clean non-human/mouse cultured-cell heat RNA-seq that turned up. Macrophage-like, not fibroblast. Needs GEO counts or SRA re-quant and dog→human orthologs. One contrast does not make a third discovery species.",
    ),
    dict(
        gse_id="GSE189108",
        species="Rattus norvegicus",
        common_name="rat",
        cell_type="engineered NRVM cardiac tissue",
        stress_class="hypoxia",
        system="engineered tissue",
        recommendation="do_not_add_now",
        already_processed="FALSE",
        reason="Counts exist, but it is engineered cardiac tissue in an oxygen-gradient device, not a monolayer CSR assay. Would mix system with ASTRA/Drive cells.",
    ),
    dict(
        gse_id="GSE166073",
        species="Rattus norvegicus",
        common_name="rat",
        cell_type="brain microvascular endothelial cells",
        stress_class="dna_damage",
        system="cell culture",
        recommendation="candidate_needs_sra",
        already_processed="FALSE",
        reason="X-ray 5/20 Gy in cultured BMEC. Cell, not fibroblast. Check whether a gene-count table is deposited; otherwise SRA. Still one class in one cell type.",
    ),
    dict(
        gse_id="GSE209531",
        species="Rattus norvegicus",
        common_name="rat",
        cell_type="PC12",
        stress_class="mixed (OGD/R)",
        system="cell culture",
        recommendation="exclude_mixed",
        already_processed="FALSE",
        reason="Oxygen-glucose deprivation/reperfusion is hypoxia plus nutrient plus reperfusion, not a single ASTRA class.",
    ),
    dict(
        gse_id="GSE282023",
        species="Mus musculus",
        common_name="mouse",
        cell_type="immortalized fibroblasts",
        stress_class="hypoxia",
        system="cell culture",
        recommendation="candidate_add_mouse",
        already_processed="FALSE",
        reason="Would fill the mouse fibroblast-hypoxia gap (current mouse hypoxia is C2C12, mESC, unreplicated MEF). Not a new species.",
    ),
    dict(
        gse_id="GSE296910",
        species="Mesocricetus auratus",
        common_name="Syrian hamster",
        cell_type="skeletal muscle satellite cells",
        stress_class="cold",
        system="cell culture",
        recommendation="out_of_class",
        already_processed="TRUE",
        reason="Already processed and usable vs Drive cold, but cold is not one of the four ASTRA discovery classes (heat/hypoxia/oxidative/dna_damage).",
    ),
    dict(
        gse_id="GSE93935",
        species="Ictidomys tridecemlineatus",
        common_name="13-lined ground squirrel",
        cell_type="iPSC-derived neurons",
        stress_class="cold",
        system="cell culture",
        recommendation="demote",
        already_processed="TRUE",
        reason="Squirrel CIRBP/RBM3 biology was judged suspect; also cold, not an ASTRA class.",
    ),
    dict(
        gse_id="GSE334758",
        species="Equus caballus",
        common_name="horse",
        cell_type="skeletal muscle myoblasts",
        stress_class="hypoxia",
        system="cell culture",
        recommendation="exclude",
        already_processed="TRUE",
        reason="Chronic 72h muscle hypoxia anti-correlates genome-wide with Drive fibroblast hypoxia.",
    ),
    dict(
        gse_id="GSE261409; GSE261939; GSE198577; GSE240251; GSE49485",
        species="sheep; bat; camel; pig; mole rat/rat",
        common_name="various",
        cell_type="tissue (in vivo)",
        stress_class="hypoxia; heat; dehydration",
        system="in vivo",
        recommendation="exclude_in_vivo",
        already_processed="TRUE",
        reason="In vivo (or bedGraph-only). Kept out of a cellular CSR meta by design.",
    ),
    dict(
        gse_id="(none found)",
        species="Macaca mulatta / fascicularis",
        common_name="macaque",
        cell_type="cultured cells",
        stress_class="heat; hypoxia; H2O2; UV",
        system="cell culture",
        recommendation="insufficient",
        already_processed="FALSE",
        reason="No macaque cultured-cell RNA-seq with a clean ASTRA-class contrast and deposited counts turned up.",
    ),
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    disc = astra_rows()
    for g in GEO_DISC:
        disc.append(
            row(
                in_discovery=True,
                source="GEO",
                n_records_collapsed=1,
                **g,
            )
        )
    # stable sort: species, class, source, gse
    disc.sort(key=lambda r: (r["common_name"], r["stress_class"], r["source"], r["gse_id"], r["experiment_id"]))
    path = OUT / "discovery_metadata.tsv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(disc)

    cand_path = OUT / "other_species_candidates.tsv"
    with cand_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CAND_FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(CANDIDATES)

    n = len(disc)
    n_fib = sum(r["is_fibroblast_like"] == "TRUE" for r in disc)
    n_hum = sum(r["common_name"] == "human" for r in disc)
    n_mou = sum(r["common_name"] == "mouse" for r in disc)
    lineages = defaultdict(int)
    states = defaultdict(int)
    for r in disc:
        lineages[r["cell_lineage"]] += 1
        states[f"{r['common_name']}|{r['cell_state']}"] += 1
    print(f"Wrote {path}  n={n}  human={n_hum} mouse={n_mou}  fibroblast-like={n_fib}")
    print("cell_state:", dict(states))
    print("lineage (top):")
    for k, v in sorted(lineages.items(), key=lambda kv: -kv[1]):
        print(f"  {v:3d}  {k}")
    print(f"Wrote {cand_path}")


if __name__ == "__main__":
    main()

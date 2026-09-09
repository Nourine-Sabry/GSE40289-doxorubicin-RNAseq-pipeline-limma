# Full-Dataset Analysis

## Methodology:

### Design:

This analysis extends the 16hr-only analysis to all 24 samples in GSE40289, spanning genotype (WT / Top2b-knockout) × treatment (dox / PBS) × time (16hr / 72hr) — 8 groups of 3 replicates each. Genotype, treatment, and time were all parsed from the `description.1` field, which encodes all three consistently (e.g. *"16 hr Top2b Del Ctr1"*, *"72 hr WT Dox2"*). Preprocessing (control-probe removal, log-transform, quantile normalization, low-expression filtering) is identical to the 16hr-only analysis and leaves 28,922 probes across 24 arrays.

A cell-means design (`~0 + group`, 8 columns) was fit, and nine contrasts were tested:

| Contrast | Definition |
|---|---|
| `dox_wt_16` / `dox_top2b_16` | dox vs pbs, at 16hr, within each genotype |
| `dox_wt_72` / `dox_top2b_72` | dox vs pbs, at 72hr, within each genotype |
| `interaction_16` / `interaction_72` | genotype × treatment interaction, at each timepoint |
| `time_effect_wt` / `time_effect_top2b` | does the dox response change between 16hr and 72hr, within each genotype |
| `three_way` | genotype × treatment × time (does the genotype-dependent dox response itself shift over time) |

An omnibus moderated F-test (`topTable(coef = c("dox_wt_16", "dox_top2b_16", "dox_wt_72", "dox_top2b_72"))`) was also run to identify genes responding to doxorubicin in *any* genotype/timepoint combination, jointly.

### Comparison to Zhang et al. 2012

Two things are worth mentioning before comparing results from the paper vs from my analysis:

1. **Different statistical pipeline.** The paper processed microarray data through Agilent Feature Extraction Software and Rosetta Resolver, applying a flat >1.5-fold threshold across 59,306 probes (unfiltered, uncollapsed). This analysis uses `limma`'s moderated t-statistics on 28,922 filtered probes with quantile normalization. Some numeric mismatch against the paper's reported fold-changes is expected from pipeline differences alone, independent of whether the biology reproduces.
2. **Different enrichment method.** The paper used Ingenuity Pathway Analysis (IPA), a commercial curated pathway database, on genes passing a >1.5-fold change threshold. This analysis uses KEGG pathways via `msigdbr`/`fgsea`, a rank-based method that doesn't use a hard fold-change cutoff.

Given that, the fairest comparison is **direction and genotype-specificity**, i.e., does WT show the reported effect and does the knockout not?

## Results

### Differential expression across all contrasts

| Contrast | Down | Up | Total DE | % of 28,922 |
|---|---:|---:|---:|---:|
| `dox_wt_16` | 3,689 | 3,251 | 6,940 | 24.0% |
| `dox_top2b_16` | 0 | 0 | 0 | 0% |
| `dox_wt_72` | 2,226 | 1,377 | 3,603 | 12.5% |
| `dox_top2b_72` | 0 | 0 | 0 | 0% |
| `interaction_16` | 1,088 | 1,644 | 2,732 | 9.4% |
| `interaction_72` | 0 | 1 | 1 | ~0% |
| `time_effect_wt` | 849 | 642 | 1,491 | 5.2% |
| `time_effect_top2b` | 0 | 0 | 0 | 0% |
| `three_way` | 36 | 37 | 73 | 0.25% |

The pattern holds at both timepoints separately: doxorubicin produces a large transcriptional response in WT hearts and *zero* FDR-significant genes in Top2b-knockout hearts, at both 16hr and 72hr independently. The `time_effect_wt` result (1,491 genes) shows the WT dox response itself changes shape between 16hr and 72hr — consistent with the paper's description of the 16hr and 72hr IPA profiles being *"dramatically shifted"* relative to each other — while `time_effect_top2b` is null, meaning there's no measurable dox response at either timepoint in the knockout to shift in the first place. The `three_way` interaction (73 genes) is small but non-zero, consistent with a genuinely genotype-dependent time-course rather than a flat absence of any
higher-order structure.

### Omnibus F-test: genes responding anywhere

The top 10 genes by the joint F-test across all four dox-vs-control contrasts include *Bax* (BCL2-associated X protein) at rank 10, which is notable because Bax is one of the genes the paper explicitly calls out as upregulated in WT following doxorubicin. Pulling its per-contrast logFC directly from that same table:

| Contrast | Bax logFC |
|---|---:|
| `dox_wt_16` | **+1.08** |
| `dox_top2b_16` | −0.14 |
| `dox_wt_72` | +0.26 |
| `dox_top2b_72` | +0.08 |

Much like the original paper, there is a clear, timepoint-appropriate increase in WT at 16hr (+1.08 log2FC, ~2.1-fold) that is essentially absent in the knockout (−0.14) and has faded by 72hr in WT as well, consistent with Bax's role in the acute apoptotic response rather than the later mitochondrial-biogenesis phase of the phenotype.

---

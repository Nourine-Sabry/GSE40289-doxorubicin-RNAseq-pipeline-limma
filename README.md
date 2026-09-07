# Mechanism of doxorubicin-induced cardiotoxicity (GSE40289)

![Banner](images/BannerImg.jpg)

## About the project:

This project was completed as part of the ['Differential Expression Analysis with limma in R'](https://app.datacamp.com/learn/courses/differential-expression-analysis-with-limma-in-r) course on DataCamp. The project covers preprocessing, differential expression analysis, and pathway enrichment analysis using the limma package.

## About the dataset:

Doxorubicin is a drug used to treat certain types of cancer. It has been found to cause cardiotoxicity as a side effect.  [Zhang et al. (2012)](https://www.nature.com/articles/nm.2919) investigated the hypothesis that doxorubicin damages heart cells by binding to the protein topoisomerase-II beta (Top2b). The experiment followed a 2x2 factorial experimental design where two types of mice were used: genetically normal wild type mice, and Top2b null mice which had cardiomyocyte-specific Top2b deletion. The mice were treated with doxorubicin or with PBS as a control solution. If doxorubicin requires Top2b to exert cardiotoxicity, the Top2b null mice should not be affected by doxorubicin treatment. We will test this hypothesis using the dataset generated from this study [GSE40289](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40289), which contains expression levels of 29,532 genes and 12 mice, with 3 replicates for each combination of the two factors. 

## Research question:

## Methodology:

### Obtaining the data:

- Data were obtained from GEO accession [GSE40289](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40289), an Agilent two-color mouse microarray study (platform GPL13912) comparing cardiac gene expression in wild-type mice and mice with a cardiomyocyte-specific knockout of *Top2b* (DNA topoisomerase II beta), each injected with doxorubicin (25 mg/kg, i.p.) or PBS. This analysis uses the 16-hour timepoint samples only (n = 3 biological replicates per genotype × treatment group, 12 samples total); the series also includes a 72-hour timepoint, which was excluded to keep a single, clean time window.
- Raw data and platform annotation were retrieved directly from GEO with `GEOquery::getGEO()`. Genotype and treatment were assigned from each sample's free-text `title` field (e.g. *"WT Cardiomyocytes 16hr after PBS rep1"*) rather than from the positional `characteristics_ch1.N` columns, since GEOquery numbers those columns by the order metadata tags happen to appear per sample, and that order was not consistent for every sample in this series.

### Preprocessing:
 
Agilent control probes were removed. Expression values were log-transformed and quantile-normalized across arrays (`limma::normalizeBetweenArrays`). Probes with a mean expression at or below zero after normalization were filtered out, leaving 28,965 probes across 12 arrays.

### Differential expression analysis:
 
- Samples were grouped by genotype × treatment (`wt.pbs`, `wt.dox`, `top2b.pbs`, `top2b.dox`). A linear model was fit without an intercept (`limma::lmFit`) and three contrasts were tested with moderated t-statistics (`contrasts.fit` → `eBayes`):
 
| Contrast | Definition | Question |
|---|---|---|
| `dox_wt` | wt.dox − wt.pbs | Doxorubicin response in wild-type hearts |
| `dox_top2b` | top2b.dox − top2b.pbs | Doxorubicin response in Top2b-knockout hearts |
| `interaction` | (top2b.dox − top2b.pbs) − (wt.dox − wt.pbs) | Does the doxorubicin response depend on genotype? |
 
- Genes were called differentially expressed at BH-adjusted p < 0.05 (`decideTests`, default settings).


### Gene set enrichment analysis
 
- Enrichment was assessed by GSEA. GSEA instead uses the complete ranked gene list and asks whether a pathway's genes are skewed toward one end of it. Genes were ranked by the moderated t-statistic from `eBayes`, which reflects both effect size and estimation confidence. GSEA was run with `fgsea` against the KEGG legacy pathway collection from `msigdbr` (mouse, `collection = "C2"`, `subcollection = "CP:KEGG_LEGACY"`, minSize = 10, maxSize = 500).
 
- Two visualization approaches were produced in parallel for each figure type: base `limma`/base R plotting functions, and equivalent `ggplot2` versions. Note that the `ggplot2` volcano plots use a nominal (uncorrected) p-value cutoff (p < 0.05, |log2FC| > 1) purely for point coloring, whereas the `decideTests` calls used for the DEG counts below use BH-adjusted p-values; the two are not counting the same "significant" set and should not be directly compared gene-for-gene.

---
 
## Results:
 
### Sample structure and normalization:
 
- Density plots of raw log-intensities showed a modest right-shift in wild-type arrays relative to Top2b-knockout arrays prior to normalization; this resolved into near-identical distributions across all 12 arrays after quantile normalization, as expected, and this shape was preserved after low-expression filtering — indicating the normalization behaved as intended and no arrays required exclusion on QC grounds.
 
- MDS of the filtered, normalized data separated samples cleanly along the first dimension (53% of variance) by genotype, with wild-type and Top2b-knockout samples forming two distinct clusters (Fig. MDS-genotype). This confirms genotype as the dominant source of transcriptional variance in this dataset, consistent with a global regulatory role for Top2b in cardiomyocytes independent of doxorubicin exposure. Direct inspection of *Top2b* probe intensity confirmed clear separation between genotypes in the expected direction, validating that the knockout samples are correctly labeled and the knockdown is transcriptionally detectable.
 
### Differential expression analysis:
 
| Contrast | DE genes (FDR < 0.05) | % of 28,965 tested |
|---|---:|---:|
| `dox_wt` | 8,108 | 28.0% |
| `dox_top2b` | 0 | 0% |
| `interaction` | 3,278 | 11.3% |

- Wild-type hearts showed a large transcriptional difference to doxorubicin 16 hours post-injection, with over a quarter of tested genes significantly differentially expressed. In contrast, Top2b-knockout hearts showed no genes reaching FDR-significance for the same doxorubicin treatment. The p-value distribution for `dox_wt` was strongly right-skewed (excess of small p-values, as expected under true signal), while the `dox_top2b` distribution was close to uniform, the shape expected under the null and independent confirmation that the zero-DEG result is not an artifact of the significance cutoff. The `interaction` contrast was significant for 3,278 genes, of which the large majority (3,207) overlapped with the `dox_wt` gene set, indicating that most of the genotype-dependent difference in doxorubicin response is explained by genes that respond in wild-type hearts but fail to respond in the knockout, rather than by a distinct interaction-specific gene program.
 
- This pattern reproduces the central finding of the original GSE40289 study: cardiomyocyte-specific deletion of Top2b blunts the transcriptional response to doxorubicin.
 
### Gene set enrichment:
 
- GSEA against KEGG pathways in `dox_wt` showed upregulation of p53 signaling and DNA replication pathways, alongside downregulation of cardiac muscle contraction, calcium signaling, and cardiomyopathy-related pathways (hypertrophic, dilated, arrhythmogenic right ventricular). This is consistent with the established mechanism of doxorubicin cardiotoxicity: Top2b-mediated DNA double-strand breaks trigger a p53-driven damage response, which in turn suppresses genes required for normal cardiomyocyte contractile function.
 
- GSEA on `dox_top2b` identified coordinated shifts in several pathways: upregulation of peroxisome, fatty acid metabolism, and branched-chain amino acid degradation pathways, and downregulation of complement/coagulation, leukocyte transendothelial migration, and cell adhesion pathways. Biologically, this suggests the doxorubicin challenge is not entirely inert in Top2b-knockout hearts; a milder, non-p53-driven metabolic and inflammatory response persists, but the acute DNA-damage transcriptional program seen in wild-type hearts is specifically absent, consistent with Top2b's proposed role as the trigger for that program.

## Discussion:

## Acknowledgement:

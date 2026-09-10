# Analysis 2: Reproducing the Original Paper

## About this analysis

This analysis uses all 24 samples in GSE40289, spanning genotype, treatment, and both timepoints (16hr and 72hr), in an attempt to reproduce what Zhang et al. (2012) actually reported, using a completely different software pipeline than the one they used?

The paper used Agilent Feature Extraction Software, Rosetta Resolver, and Ingenuity Pathway Analysis (IPA), a commercial curated pathway database. This analysis uses `limma`, quantile normalization, and KEGG pathways through `msigdbr` and `fgsea`. The goal is for the two pipelines to land on the same biological conclusion despite being built from different tools and different statistical assumptions.

## Methodology:

### Design:

A cell means design (`~0 + group`, 8 columns, one per genotype by treatment by time combination) was fit, and nine contrasts were tested:

| Contrast | Definition |
|---|---|
| `dox_wt_16` / `dox_top2b_16` | dox vs pbs, at 16hr, within each genotype |
| `dox_wt_72` / `dox_top2b_72` | dox vs pbs, at 72hr, within each genotype |
| `interaction_16` / `interaction_72` | genotype by treatment interaction, at each timepoint |
| `time_effect_wt` / `time_effect_top2b` | does the dox response change between 16hr and 72hr, within each genotype |
| `three_way` | genotype by treatment by time (does the genotype dependent dox response itself shift over time) |

An omnibus moderated F-test (`topTable(coef = c("dox_wt_16", "dox_top2b_16", "dox_wt_72", "dox_top2b_72"))`) was also run, to find genes responding to doxorubicin in any genotype and timepoint combination jointly, rather than in one specific contrast.

This model is fit once, across all 24 samples together, and the timepoint specific contrasts are pulled out of that single fit.

### Comparison to Zhang et al. 2012:

The main difference between the original paper and this reanalysis:

1. First, the statistics are different. The paper applied a flat 1.5 fold change threshold across 59,306 uncollapsed probes, processed through Rosetta Resolver. This analysis uses `limma`'s moderated t-statistics on 28,922 filtered, quantile normalized probes. Some numeric mismatch against the paper's exact fold changes is expected from the pipeline difference alone, independent of whether the biology holds up.

2. Second, the enrichment method is different. The paper used IPA, a commercial pathway database with its own curated, proprietary pathway definitions, on genes passing that 1.5 fold cutoff. This analysis uses KEGG pathways through `fgsea`, a rank based method that doesn't use a hard cutoff at all.

Given all that, the fairest comparison is direction and genotype specificity: does wild type show the reported effect, and does the knockout not, rather than expecting exact fold change values.

---

## Results

### Differential expression across all contrasts:

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

The pattern from analysis 1 holds at both timepoints when checked independently: doxorubicin drives a large transcriptional response in wild type hearts and zero FDR significant genes in the knockout, at 16hr and at 72hr separately. `time_effect_wt` (1,491 genes) shows the wild type dox response itself changes character between 16hr and 72hr, which matches the paper's own description of the two timepoints' pathway profiles being dramatically different from each other. `time_effect_top2b` is null, meaning there's no measurable dox response at either timepoint in the knockout for that response to shift in the first place. The `three_way` interaction consisted of 73 genes, which fits a genotype dependent time course.

### MDS and clustering:

![](/images/mds_all.jpg)

Samples separate along the first MDS dimension almost entirely by genotype, not by time or treatment. This is the same pattern seen in the 16 hour only analysis: genotype is the largest source of transcriptional variance in this dataset.

![](/images/heatmap_all.jpg)

The top 50 most variable genes cluster samples into two genotype clusters, with time and treatment structure visible as a secondary pattern within each block.

### Volcano plots:

![](/images/volc_wt_16.jpg)

![](/images/volc_wt_72.jpg)

Both wild type contrasts show a broad, roughly symmetric spread of up and down regulated genes, with the 16hr response noticeably larger than the 72hr one, matching the DEG counts above.

### Pathway enrichment:

![](/images/gsea_dox_wt_16.jpg)

At 16hr, GSEA in wild type shows the same pattern as the 16 hour only analysis: p53 signaling and DNA replication up, cardiac contractile and cardiomyopathy related pathways down. The top hit, KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS, is again the strongest signal by a wide margin (padj = 3.2 x 10⁻¹⁷). 

![](/images/gsea_dox_wt_72.jpg)

At 72hr, wild type shows a top hit of KEGG_DRUG_METABOLISM_CYTOCHROME_P450 (padj = 8.9 x 10⁻¹²), along with oxidative phosphorylation and several cardiomyopathy related pathways among the most significantly downregulated. This is a match to the paper's 72 hour finding, defective mitochondrial function and oxidative phosphorylation, specifically in wild type.

## Discussion:

- Integration of the complete time course reveals that the doxorubicin response in wild-type hearts comprises two temporally distinct phases with divergent biological themes. At 16 hours, the transcriptional landscape is dominated by DNA damage and apoptotic signaling, characterized by upregulation of p53 pathway genes and specific induction of Bax. Bax encodes a pro-apoptotic protein that promotes mitochondrial outer membrane permeabilization, and its induction aligns with progression of the DNA damage response to the point of apoptotic commitment in a subset of cells. By 72 hours, this acute genotoxic signature had substantially attenuated (Bax log-fold change decreased from +1.08 at 16 hours to +0.26 at 72 hours in wild-type), and a distinct transcriptional program emerged: significant downregulation of genes involved in oxidative phosphorylation (OXPHOS) and mitochondrial function, concomitant with upregulation of drug metabolism and xenobiotic clearance pathways.

- The early p53-mediated apoptotic response represents the immediate cellular reaction to acute genotoxicity. Subsequent suppression of mitochondrial and OXPHOS genes likely reflects a downstream consequence, wherein surviving cardiomyocytes exhibit diminished capacity for ATP production via oxidative phosphorylation. Given the exceptional dependence of cardiomyocytes on mitochondrial ATP generation relative to most other cell types, such metabolic reprogramming constitutes a plausible molecular basis for the progressive contractile weakness (reduced ejection fraction) observed clinically following repeated doxorubicin dosing, despite the present dataset encompassing only a single acute exposure.

- Genotype dependence was evident across both temporal phases, providing further evidence to the central role of Top2b in this cascade. Bax induction at 16 hours was essentially restricted to wild-type hearts, and OXPHOS suppression at 72 hours was markedly more pronounced and statistically robust in wild-type compared to knockout. Collectively, these findings support a two-step causal model: Top2b is required for initiation of the DNA damage response, and this primary genotoxic event subsequently drives mitochondrial dysfunction. Disruption of the initial link effectively prevents propagation of the downstream mitochondrial phenotype, at least to a comparable magnitude

- Notably, a modest but detectable OXPHOS suppression signal was observed in Top2b knockout hearts at 72 hours, albeit substantially weaker than in wild-type. Several non-mutually exclusive biological interpretations merit consideration. First, Top2b-independent mechanisms of doxorubicin action, most notably reactive oxygen species generation via redox cycling, may exert mild selective pressure on mitochondrial gene expression even in the absence of Top2b-mediated DNA damage. Second, this signal could represent a downstream echo of the limited transcriptional changes still detectable in knockout hearts at 16 hours (interaction and time-effect contrasts were not entirely null at this early time point), rather than a truly independent mitochondrial effect. Third, given the limited statistical power inherent to three biological replicates per group, a portion of this signal may reflect biological variability coincidentally aligning with a large, coherent gene set. The present data cannot definitively distinguish among these possibilities. Nevertheless, the findings support an interpretation of the original claim: Top2b is the dominant driver of doxorubicin-induced mitochondrial dysfunction, though it does not account for the entirety of the observed phenotype.

## References:

1. Gabani, M., Castañeda, D., Nguyen, Q. M., Choi, S. K., Chen, C., Mapara, A., Kassan, A., Gonzalez, A. A., Khataei, T., Ait-Aissa, K., & Kassan, M. (2021). Association of Cardiotoxicity With Doxorubicin and Trastuzumab: A Double-Edged Sword in Chemotherapy. Cureus, 13(9), e18194. https://doi.org/10.7759/cureus.18194

2. Zhang et al. Transcriptomic profiling reveals p53 as a key regulator of doxorubicin-induced cardiotoxicity. Cell Death Discovery (2019).

3. Kelly, C., Kiltschewskij, D.J., Leong, A.J. et al. Identifying common pathways for doxorubicin and carfilzomib-induced cardiotoxicities: transcriptomic and epigenetic profiling. Sci Rep 15, 4395 (2025). https://doi.org/10.1038/s41598-025-87442-5

4. Guo, Y., Tang, Y., Lu, G., & Gu, J. (2023). p53 at the Crossroads between Doxorubicin-Induced Cardiotoxicity and Resistance: A Nutritional Balancing Act. Nutrients, 15(10), 2259. https://doi.org/10.3390/nu15102259

## Acknowledgement

This analysis builds on the pipeline developed in analysis 1, which itself began as an exercise from DataCamp's "Differential Expression Analysis with Limma in R" course. This extension to the full dataset, the contrast design across genotype, treatment, and time, and the direct comparison against Zhang et al. (2012) were carried out independently as a public reanalysis project.

# Mechanism of doxorubicin-induced cardiotoxicity: Reanalysis of GSE40289

![Banner](images/BannerImg.jpg)

## Background:

Doxorubicin is one of the most effective chemotherapy drugs in use, but its clinical value is limited by dose dependent cardiotoxicity. Zhang et al. (2012) proposed that this toxicity is not simply a side effect of redox cycling, but is instead driven by doxorubicin binding topoisomerase II beta (Top2b) in cardiomyocytes. They tested this by comparing wild type mice to mice with a cardiomyocyte specific Top2b deletion, treated with either doxorubicin or PBS, and profiled cardiac gene expression by microarray at 16 and 72 hours after injection. The hypothesis was that if Top2b is required for the toxicity, knocking it out in heart muscle should blunt or remove the transcriptional damage response that doxorubicin normally triggers.

The paper's data is public as GEO accession [GSE40289](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40289). This repo is a from scratch reanalysis of that data in R, done in two stages with two different goals.

## Structure

**[Analysis 1: DataCamp reproduction](./analysis1_results.md)**
Uses the 16 hour timepoint only (12 samples, the 2x2 genotype by treatment design), matching the scope of a limma teaching exercise from DataCamp's ["Differential Expression Analysis with limma in R"](https://app.datacamp.com/learn/courses/differential-expression-analysis-with-limma-in-r) course.

**[Analysis 2: reproducing the original paper](./analysis2_results.md)**
Extends the pipeline to the full dataset, all 24 samples across genotype, treatment, and both timepoints (16hr and 72hr). The aim here was to see how closely I can reproduce Zhang et al.'s results.

## Data and tools:

- **Data**: GEO GSE40289, Agilent two color mouse microarray (platform GPL13912)
- **Core pipeline**: R, `GEOquery`, `limma`, `Biobase`
- **Visualization**: `ggplot2`, `pheatmap`
- **Pathway enrichment**: `fgsea`, `msigdbr` (KEGG legacy gene sets)
- **Reference paper**: Zhang, S. et al. "Identification of the molecular basis of doxorubicin-induced cardiotoxicity." *Nature Medicine* 18, 1639-1642 (2012). [doi:10.1038/nm.2919](https://doi.org/10.1038/nm.2919)

## A note on how these two analyses relate:

They share the same underlying data and the same core preprocessing logic (control probe removal, log transform, quantile normalization, low expression filtering), but they are not the same experiment scaled up. Analysis 1 fits a model on the 16 hour samples alone. Analysis 2 fits a single model across all 24 samples and pulls timepoint specific contrasts out of it. Because `limma`'s empirical Bayes moderation borrows information across every gene in the fitted model, the exact significance counts for a "16 hour dox in WT" contrast are not identical between the two analyses even though the contrast formula is the same in both. This is expected and is discussed in more detail in the analysis 2 write up, but it's worth knowing going in: these are two related but distinct models, not one analysis with an extra data source bolted on.

## Acknowledgement:

This project began as an exercise from DataCamp's "Differential Expression Analysis with Limma in R" course. Analysis 1 stays close to that course's structure and scope. Analysis 2 goes beyond it, extending the pipeline to the full dataset and comparing it against the original published findings. All code, interpretation, and the comparison against the paper were done independently as a public reanalysis project.

#Investigating doxorubicin-induced cardiotoxicity

library(GEOquery)
library(limma)
library(Biobase)

gset <- getGEO("GSE40289", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL13912", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

eset <- gset

# Clean phenotype column names
colnames(pData(eset)) <- make.names(colnames(pData(eset)))

# Check new names
colnames(pData(eset))

# Subset 16 hr samples
table(pData(eset)[,"description.1"])
idx_16h <- grepl("16 hr", pData(eset)[,"description.1"])
eset <- eset[, idx_16h]
ncol(eset)
table(pData(eset)[,"description.1"])

titles <- pData(eset)$title

pData(eset)$gen_short   <- ifelse(grepl("Top2b", titles, ignore.case = TRUE),
                                  "top2b", "wt")
pData(eset)$treat_short <- ifelse(grepl("doxorubicin", titles, ignore.case = TRUE),
                                  "dox", "pbs")
pData(eset)$group <- paste(pData(eset)$gen_short, pData(eset)$treat_short, sep = ".")

table(pData(eset)$gen_short, pData(eset)[,"characteristics_ch1.3"])
table(pData(eset)$treat_short, pData(eset)[,"characteristics_ch1.4"])

table(fData(eset)$CONTROL_TYPE)
eset <- eset[fData(eset)$CONTROL_TYPE == "FALSE", ]
dim(exprs(eset))

# Log transform
exprs(eset) <- log(exprs(eset))
plotDensities(eset,  group = pData(eset)[,"gen_short"], legend = "topright")
# Quantile normalize
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset,  group = pData(eset)[,"gen_short"], legend = "topright")
# Determine the genes with mean expression level greater than 0
keep <- rowMeans(exprs(eset)) > 0
sum(keep)
# Filter the genes
eset <- eset[keep]
plotDensities(eset, group = pData(eset)[,"gen_short"], legend = "topright")

# Find the row which contains Top2b expression data
top2b <- which(fData(eset)[,"GENE_SYMBOL"] == "Top2b")

# Plot Top2b expression versus genotype
boxplot(exprs(eset)[top2b, ] ~ pData(eset)[, "gen_short"],
        main = fData(eset)[top2b, ])

# Plot principal components labeled by genotype
plotMDS(eset, labels = pData(eset)[,"gen_short"], gene.selection = "common")

# Plot principal components labeled by treatment
plotMDS(eset, labels = pData(eset)[,"treat_short"], gene.selection = "common")

group <- factor(pData(eset)$group)

# Create design matrix with no intercept
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Count the number of samples modeled by each coefficient
colSums(design)

# Create a contrasts matrix
cm <- makeContrasts(dox_wt = wt.dox - wt.pbs,
                    dox_top2b = top2b.dox - top2b.pbs,
                    interaction = (top2b.dox - top2b.pbs) - (wt.dox - wt.pbs),
                    levels = design)

# View the contrasts matrix
cm

# Fit the model
fit <- lmFit(eset, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2)
summary(results)

# Create a Venn diagram
vennDiagram(results)

# Obtain the summary statistics for the contrast dox_wt
stats_dox_wt <- topTable(fit2, coef = "dox_wt", number = nrow(fit2),
                         sort.by = "none")
# Obtain the summary statistics for the contrast dox_top2b
stats_dox_top2b <- topTable(fit2, coef = "dox_top2b", number = nrow(fit2),
                            sort.by = "none")
# Obtain the summary statistics for the contrast interaction
stats_interaction <- topTable(fit2, coef = "interaction", number = nrow(fit2),
                              sort.by = "none")

# Create histograms of the p-values for each contrast
hist(stats_dox_wt[,"P.Value"])
hist(stats_dox_top2b[,"P.Value"])
hist(stats_interaction[,"P.Value"])

# Extract the gene symbols
gene_symbols <- fit2$genes[,"GENE_SYMBOL"]

# Create a volcano plot for the contrast dox_wt
volcanoplot(fit2, coef = "dox_wt", highlight = 5, names = gene_symbols)

# Create a volcano plot for the contrast dox_top2b
volcanoplot(fit2, coef = "dox_top2b", highlight = 5, names = gene_symbols)

# Create a volcano plot for the contrast interaction
volcanoplot(fit2, coef = "interaction", highlight = 5, names = gene_symbols)

# Extract the entrez gene IDs
entrez <- fit2$genes[,"GENE_ID"]

# Test for enriched KEGG Pathways for contrast dox_wt
enrich_dox_wt <- kegga(fit2, coef = "dox_wt", geneid = entrez, species = "Mm")

# View the top 5 enriched KEGG pathways
topKEGG(enrich_dox_wt, number = 5)

# Test for enriched KEGG Pathways for contrast interaction
enrich_interaction <- kegga(fit2, coef = "interaction", geneid = entrez, species = "Mm")

# View the top 5 enriched KEGG pathways
topKEGG(enrich_interaction, number = 5)

#########################################
# New plots + GSEA
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(msigdbr)
library(fgsea)

theme_set(theme_minimal(base_size = 13))

# Alternative plots:

## Density plot (replaces plotDensities):
plot_densities_gg <- function(eset, group_col, title = "") {
  df <- as.data.frame(exprs(eset))
  df$probe <- rownames(df)
  df_long <- melt(df, id.vars = "probe", variable.name = "sample", value.name = "expr")
  grp <- pData(eset)[as.character(df_long$sample), group_col]
  df_long$group <- grp
  
  ggplot(df_long, aes(x = expr, color = group, group = sample)) +
    geom_density(linewidth = 0.5, alpha = 0.7) +
    labs(title = title, x = "Expression", y = "Density", color = "Group") +
    scale_color_brewer(palette = "Set2")
}

plot_densities_gg(eset, "gen_short", "Density: log-transformed")

## MDS plot (replaces plotMDS):
plot_mds_gg <- function(eset, label_col, title = "") {
  mds <- plotMDS(eset, plot = FALSE, gene.selection = "common")
  df <- data.frame(
    x = mds$x, y = mds$y,
    label = pData(eset)[, label_col]
  )
  ggplot(df, aes(x = x, y = y, color = label, label = label)) +
    geom_point(size = 3) +
    geom_text_repel(show.legend = FALSE, size = 3.2) +
    labs(title = title, x = "Leading logFC dim 1", y = "Leading logFC dim 2", color = label_col) +
    scale_color_brewer(palette = "Set1")
}

plot_mds_gg(eset, "gen_short", "MDS by genotype")
plot_mds_gg(eset, "treat_short", "MDS by treatment")

## Top2b boxplot:
top2b_df <- data.frame(
  expr = exprs(eset)[top2b, ],
  genotype = pData(eset)$gen_short
)
ggplot(top2b_df, aes(x = genotype, y = expr, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, size = 2) +
  labs(title = "Top2b expression by genotype", x = NULL, y = "log expression") +
  scale_fill_brewer(palette = "Pastel1") +
  theme(legend.position = "none")

## Volcano plot (replaces limma::volcanoplot)
plot_volcano_gg <- function(stats_table, title, fc_col = "logFC",
                            p_col = "P.Value", label_col = "GENE_SYMBOL",
                            p_thresh = 0.05, fc_thresh = 1, top_n_label = 8) {
  df <- stats_table
  df$sig <- "NS"
  df$sig[df[[p_col]] < p_thresh & df[[fc_col]] >  fc_thresh] <- "Up"
  df$sig[df[[p_col]] < p_thresh & df[[fc_col]] < -fc_thresh] <- "Down"
  
  top_genes <- df %>%
    filter(sig != "NS") %>%
    arrange(.data[[p_col]]) %>%
    head(top_n_label)
  
  ggplot(df, aes(x = .data[[fc_col]], y = -log10(.data[[p_col]]), color = sig)) +
    geom_point(alpha = 0.5, size = 1.3) +
    geom_text_repel(data = top_genes, aes(label = .data[[label_col]]),
                    size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    labs(title = title, x = "log2 fold change", y = "-log10(p-value)", color = NULL)
}

plot_volcano_gg(stats_dox_wt, "dox_wt")
plot_volcano_gg(stats_dox_top2b, "dox_top2b")
plot_volcano_gg(stats_interaction, "interaction")

#################################################
# GSEA

options(timeout = 300)

msig_kegg_mm <- tryCatch(
  msigdbr(species = "Mus musculus", db_species = "MM",
          collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  error = function(e) {
    message("Mouse-native fetch failed (", conditionMessage(e),
            "), falling back to human-ortholog-mapped gene sets.")
    msigdbr(species = "Mus musculus",
            collection = "C2", subcollection = "CP:KEGG_LEGACY")
  }
)

stopifnot(nrow(msig_kegg_mm) > 0)

saveRDS(msig_kegg_mm, "msig_kegg_mm.rds")

id_col <- intersect(c("entrez_gene", "ncbi_gene"), colnames(msig_kegg_mm))
stopifnot(length(id_col) >= 1)
id_col <- id_col[1]
message("Using column '", id_col, "' as the gene ID for pathway_list.")

pathway_list <- split(as.character(msig_kegg_mm[[id_col]]), msig_kegg_mm$gs_name)
stopifnot(length(pathway_list) > 0)

run_gsea <- function(stats_table, entrez_col_vec, pathways, title) {
  ranks <- stats_table$t
  names(ranks) <- as.character(entrez_col_vec)
  ranks <- ranks[!is.na(names(ranks)) & names(ranks) != "" & !duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  
  fg <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
  fg <- fg[order(fg$padj), ]
  
  top20 <- head(fg, 20)
  p <- ggplot(top20, aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70")) +
    labs(title = title, x = NULL, y = "Normalized Enrichment Score (NES)", fill = "padj < 0.05")
  
  list(results = fg, plot = p, ranks = ranks)
}

gsea_dox_wt      <- run_gsea(stats_dox_wt, entrez, pathway_list, "GSEA: dox_wt (KEGG)")
gsea_dox_top2b   <- run_gsea(stats_dox_top2b, entrez, pathway_list, "GSEA: dox_top2b (KEGG)")
gsea_interaction <- run_gsea(stats_interaction, entrez, pathway_list, "GSEA: interaction (KEGG)")

gsea_dox_wt$plot
gsea_dox_top2b$plot
gsea_interaction$plot

tidy_gsea <- function(fg) {
  fg$leadingEdge <- sapply(fg$leadingEdge, paste, collapse = ",")
  as.data.frame(fg)
}

head(tidy_gsea(gsea_dox_wt$results), 10)
head(tidy_gsea(gsea_dox_top2b$results), 10)
head(tidy_gsea(gsea_interaction$results), 10)

top_pathway <- gsea_dox_wt$results$pathway[1]
plotEnrichment(pathway_list[[top_pathway]], gsea_dox_wt$ranks) +
  labs(title = top_pathway)
# Doxorubicin-induced cardiotoxicity: full dataset, all groups
# GSE40289: genotype (wt/top2b) x treatment (dox/pbs) x time (16hr/72hr)
# 8 groups x 3 replicates = 24 samples

library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(msigdbr)
library(fgsea)

theme_set(theme_minimal(base_size = 13))

# Obtaining the data

gset <- getGEO("GSE40289", GSEMatrix = TRUE, getGPL = TRUE)
idx  <- if (length(gset) > 1) grep("GPL13912", attr(gset, "names")) else 1
eset <- gset[[idx]]

colnames(pData(eset)) <- make.names(colnames(pData(eset)))

desc <- pData(eset)$description.1
pData(eset)$time_short   <- ifelse(grepl("16 hr", desc), "16hr", "72hr")
pData(eset)$gen_short    <- ifelse(grepl("Top2b", desc, ignore.case = TRUE), "top2b", "wt")
pData(eset)$treat_short  <- ifelse(grepl("Dox", desc, ignore.case = TRUE), "dox", "pbs")
pData(eset)$group        <- paste(pData(eset)$gen_short, pData(eset)$treat_short,
                                   pData(eset)$time_short, sep = ".")

table(pData(eset)$group)

eset <- eset[fData(eset)$CONTROL_TYPE == "FALSE", ]

# Preprocessing

exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
eset <- eset[rowMeans(exprs(eset)) > 0]

# -----------------------------------------------------------------------------
# 3. Design and contrasts -- 8-group cell-means model
# -----------------------------------------------------------------------------

group  <- factor(pData(eset)$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
colSums(design)

cm <- makeContrasts(
  dox_wt_16      = wt.dox.16hr    - wt.pbs.16hr,
  dox_top2b_16   = top2b.dox.16hr - top2b.pbs.16hr,
  dox_wt_72      = wt.dox.72hr    - wt.pbs.72hr,
  dox_top2b_72   = top2b.dox.72hr - top2b.pbs.72hr,

  interaction_16 = (top2b.dox.16hr - top2b.pbs.16hr) - (wt.dox.16hr - wt.pbs.16hr),
  interaction_72 = (top2b.dox.72hr - top2b.pbs.72hr) - (wt.dox.72hr - wt.pbs.72hr),

  time_effect_wt    = (wt.dox.72hr    - wt.pbs.72hr)    - (wt.dox.16hr    - wt.pbs.16hr),
  time_effect_top2b = (top2b.dox.72hr - top2b.pbs.72hr) - (top2b.dox.16hr - top2b.pbs.16hr),

  three_way = ((top2b.dox.72hr - top2b.pbs.72hr) - (wt.dox.72hr - wt.pbs.72hr)) -
              ((top2b.dox.16hr - top2b.pbs.16hr) - (wt.dox.16hr - wt.pbs.16hr)),

  levels = design
)

fit  <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

any_dox_response <- topTable(fit2, coef = c("dox_wt_16", "dox_top2b_16",
                                             "dox_wt_72", "dox_top2b_72"),
                              number = nrow(fit2), sort.by = "F")
head(any_dox_response, 10)

stats_dox_wt_16     <- topTable(fit2, coef = "dox_wt_16",    number = nrow(fit2), sort.by = "none")
stats_dox_top2b_16  <- topTable(fit2, coef = "dox_top2b_16", number = nrow(fit2), sort.by = "none")
stats_dox_wt_72      <- topTable(fit2, coef = "dox_wt_72",    number = nrow(fit2), sort.by = "none")
stats_dox_top2b_72   <- topTable(fit2, coef = "dox_top2b_72", number = nrow(fit2), sort.by = "none")

gene_symbols <- fit2$genes[, "GENE_SYMBOL"]
entrez       <- fit2$genes[, "GENE_ID"]

# MDS and heatmap across all 8 groups

mds <- plotMDS(eset, plot = FALSE, gene.selection = "common")
mds_df <- data.frame(x = mds$x, y = mds$y, pData(eset)[, c("gen_short", "treat_short", "time_short")])

ggplot(mds_df, aes(x = x, y = y, color = gen_short, shape = time_short)) +
  geom_point(size = 3) +
  labs(title = "MDS: all 24 samples", x = "Dim 1", y = "Dim 2",
       color = "Genotype", shape = "Time") +
  theme_minimal(base_size = 13)

# Heatmap
top_var_genes <- order(apply(exprs(eset), 1, var), decreasing = TRUE)[1:50]
heatmap_mat <- exprs(eset)[top_var_genes, ]
row_labels <- fData(eset)[top_var_genes, "GENE_SYMBOL"]
row_labels[row_labels == "" | is.na(row_labels)] <- rownames(heatmap_mat)[row_labels == "" | is.na(row_labels)]
rownames(heatmap_mat) <- make.unique(row_labels)

annotation_col <- pData(eset)[, c("gen_short", "treat_short", "time_short")]
pheatmap(heatmap_mat, scale = "row",
         annotation_col = annotation_col, show_rownames = TRUE, fontsize_row = 6,
         main = "Top 50 most variable genes across all 24 samples")

# Volcano plots for the four dox-vs-control contrasts

plot_volcano_gg <- function(stats_table, title, fc_col = "logFC", p_col = "P.Value",
                             label_col = "GENE_SYMBOL", p_thresh = 0.05,
                             fc_thresh = 1, top_n_label = 8) {
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

plot_volcano_gg(stats_dox_wt_16, "dox_wt_16")
plot_volcano_gg(stats_dox_wt_72, "dox_wt_72")

# Pathway enrichment (fgsea, KEGG via msigdbr)

options(timeout = 300)

msig_kegg_mm <- tryCatch(
  msigdbr(species = "Mus musculus", db_species = "MM",
          collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  error = function(e) {
    message("Mouse-native fetch failed, falling back to human-ortholog-mapped gene sets.")
    msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  }
)

id_col <- intersect(c("entrez_gene", "ncbi_gene"), colnames(msig_kegg_mm))[1]
pathway_list <- split(as.character(msig_kegg_mm[[id_col]]), msig_kegg_mm$gs_name)

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

gsea_dox_wt_16    <- run_gsea(stats_dox_wt_16,    entrez, pathway_list, "GSEA: dox_wt_16 (KEGG)")
gsea_dox_top2b_16 <- run_gsea(stats_dox_top2b_16, entrez, pathway_list, "GSEA: dox_top2b_16 (KEGG)")
gsea_dox_wt_72     <- run_gsea(stats_dox_wt_72,    entrez, pathway_list, "GSEA: dox_wt_72 (KEGG)")
gsea_dox_top2b_72  <- run_gsea(stats_dox_top2b_72, entrez, pathway_list, "GSEA: dox_top2b_72 (KEGG)")

gsea_dox_wt_16$plot
gsea_dox_wt_72$plot

tidy_gsea <- function(fg) {
  fg$leadingEdge <- sapply(fg$leadingEdge, paste, collapse = ",")
  as.data.frame(fg)
}

head(tidy_gsea(gsea_dox_wt_16$results), 10)
head(tidy_gsea(gsea_dox_top2b_16$results), 10)
head(tidy_gsea(gsea_dox_wt_72$results), 10)
head(tidy_gsea(gsea_dox_top2b_72$results), 10)

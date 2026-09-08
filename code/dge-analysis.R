# Title : RAD52 is a therapeutic vulnerability for acquired PARP inhibitor resistance in ovarian cancer
# Code Authors : Prasanth Kumar Thuthika, Yukihide Ota

# libraries
library(edgeR)
library(limma)
library(readxl)
library(writexl)
library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tidyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)
library(forcats)

# Package versions
# limma 3.68.4
# edgeR 4.10.1
# msigdbr 26.1.0
# clusterProfiler 4.20.0
# org.Mm.eg.db 3.23.0


# differential gene expression analysis
# all gene count data
gene_counts <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/input/all.gene_counts.xlsx")
# sample columns
selected_samples <- c(
  "sample.Scr-UT_1", "sample.Scr-UT_2", "sample.Scr-UT_3",
  "sample.KO1-UT_1", "sample.KO1-UT_2", "sample.KO1-UT_3",
  "sample.Scr_ola1", "sample.Scr_ola2", "sample.Scr_ola3",
  "sample.KO1-ola1", "sample.KO1-ola2", "sample.KO1-ola3"
)
# count matrix
count_matrix <- as.matrix(gene_counts[, selected_samples])
rownames(count_matrix) <- gene_counts$ensembl_gene_id
# sample groups
group <- factor(c(
  rep("Scr_UT", 3),
  rep("KO1_UT", 3),
  rep("Scr_Ola", 3),
  rep("KO1_Ola", 3)
), levels = c("Scr_UT", "KO1_UT", "Scr_Ola", "KO1_Ola"))

# dgelist object
dge <- DGEList(counts = count_matrix, group = group)

# Filter low-expression genes
keep <- rowSums(cpm(dge) >= 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalization
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# voomWithQualityWeights
v <- voomWithQualityWeights(dge, design, plot = TRUE)

# Fit the linear model
fit <- lmFit(v, design)

# contrasts for comparisons
contrast_matrix <- makeContrasts(
  KO1_UT_vs_Scr_UT = KO1_UT - Scr_UT,
  KO1_Ola_vs_Scr_Ola = KO1_Ola - Scr_Ola,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

res_KO1_UT_vs_Scr_UT <- topTable(fit2, coef = "KO1_UT_vs_Scr_UT", number = Inf, sort.by = "none")
res_KO1_Ola_vs_Scr_Ola <- topTable(fit2, coef = "KO1_Ola_vs_Scr_Ola", number = Inf, sort.by = "none")

res_KO1_UT_vs_Scr_UT$ensembl_gene_id <- rownames(res_KO1_UT_vs_Scr_UT)
res_KO1_Ola_vs_Scr_Ola$ensembl_gene_id <- rownames(res_KO1_Ola_vs_Scr_Ola)

res1 <- res_KO1_UT_vs_Scr_UT[, c("ensembl_gene_id", "logFC", "P.Value", "adj.P.Val")]
res2 <- res_KO1_Ola_vs_Scr_Ola[, c("ensembl_gene_id", "logFC", "P.Value", "adj.P.Val")]

colnames(res1) <- c("ensembl_gene_id", "KO1_UT_vs_Scr_UT_log2FC", "KO1_UT_vs_Scr_UT_pval", "KO1_UT_vs_Scr_UT_padj")
colnames(res2) <- c("ensembl_gene_id", "KO1_Ola_vs_Scr_Ola_log2FC", "KO1_Ola_vs_Scr_Ola_pval", "KO1_Ola_vs_Scr_Ola_padj")

merged_results <- merge(res1, res2, by = "ensembl_gene_id")
final_results <- merge(gene_counts, merged_results, by = "ensembl_gene_id")

# Save results
write_xlsx(final_results, "rna_seq_results.xlsx")

# Calculate Z-Score Expression values
gene_counts <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/rna_seq_results.xlsx")

# Sample columns
sample_cols <- c(
  "sample.Scr-UT_1", "sample.Scr-UT_2", "sample.Scr-UT_3",
  "sample.KO1-UT_1", "sample.KO1-UT_2", "sample.KO1-UT_3",
  "sample.Scr_ola1", "sample.Scr_ola2", "sample.Scr_ola3",
  "sample.KO1-ola1", "sample.KO1-ola2", "sample.KO1-ola3"
)

count_matrix <- as.matrix(gene_counts[, sample_cols])
rownames(count_matrix) <- gene_counts$ensembl_gene_id

group <- factor(c(rep("Scr_UT", 3), rep("KO1_UT", 3), rep("Scr_Ola", 3), rep("KO1_Ola", 3)))
dge <- DGEList(counts = count_matrix, group = group)
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

v <- voomWithQualityWeights(dge, design, plot = TRUE)

# Z-score calculation
z_scores <- t(scale(t(v$E)))
z_scores_df <- as.data.frame(z_scores)

z_scores_df$ensembl_gene_id <- rownames(z_scores_df)
z_scores_df <- z_scores_df[, c(ncol(z_scores_df), 1:(ncol(z_scores_df) - 1))]
colnames(z_scores_df)[-1] <- paste0("z.", colnames(z_scores_df)[-1])

merged_df <- gene_counts %>%
  left_join(z_scores_df, by = "ensembl_gene_id")

# Save
write.xlsx(merged_df, file = "Z-Scores_expression.xlsx", overwrite = TRUE)
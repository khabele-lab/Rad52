# Title : RAD52 is a therapeutic vulnerability for acquired PARP inhibitor resistance in ovarian cancer
# Code Authors : Prasanth Kumar Thuthika, Yukihide Ota

# Hallmark - Gene Set Enrichment Analysis
# KO1_UT_vs_Scr_UT
results <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/rna_seq_results.xlsx")
results <- results %>%
  mutate(rank_metric = KO1_UT_vs_Scr_UT_log2FC) # logFC ranking

ranked_list <- results %>%
  dplyr::select(external_gene_name, rank_metric) %>%
  filter(!is.na(rank_metric)) %>%
  arrange(desc(rank_metric))
ranked_list <- ranked_list %>%
  filter(!is.na(external_gene_name), external_gene_name != "") %>%
  arrange(desc(rank_metric))
ranked_list_unique <- ranked_list %>%
  distinct(external_gene_name, .keep_all = TRUE)
ranked_vector <- ranked_list_unique$rank_metric
names(ranked_vector) <- ranked_list_unique$external_gene_name

# mouse hallmark gene sets
msig_hallmark_mouse <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

gsea_result <- GSEA(
  geneList = ranked_vector,
  TERM2GENE = msig_hallmark_mouse[, c("gs_name", "gene_symbol")],
  pvalueCutoff = 0.05,
  eps = 0,
  nPerm = 10000,
  verbose = TRUE
)
gsea_df <- as.data.frame(gsea_result)
View(gsea_df)  
write_xlsx(gsea_df, path = "GSEA_Hallmark_KO1_UT_vs_Scr_UT.xlsx")

# KO1_Ola_vs_Scr_Ola
results <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/rna_seq_results.xlsx")
results <- results %>%
  mutate(rank_metric = KO1_Ola_vs_Scr_Ola_log2FC) # logFC ranking

ranked_list <- results %>%
  dplyr::select(external_gene_name, rank_metric) %>%
  filter(!is.na(rank_metric)) %>%
  arrange(desc(rank_metric))
ranked_list <- ranked_list %>%
  filter(!is.na(external_gene_name), external_gene_name != "") %>%
  arrange(desc(rank_metric))
ranked_list_unique <- ranked_list %>%
  distinct(external_gene_name, .keep_all = TRUE)
ranked_vector <- ranked_list_unique$rank_metric
names(ranked_vector) <- ranked_list_unique$external_gene_name

# mouse hallmark gene sets
msig_hallmark_mouse <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

gsea_result <- GSEA(
  geneList = ranked_vector,
  TERM2GENE = msig_hallmark_mouse[, c("gs_name", "gene_symbol")],
  pvalueCutoff = 0.05,
  eps = 0,
  nPerm = 10000,
  verbose = TRUE
)
gsea_df <- as.data.frame(gsea_result)
View(gsea_df)  
write_xlsx(gsea_df, path = "GSEA_Hallmark_KO1_Ola_vs_Scr_Ola.xlsx")

# Supplementary Fig: S9-g
gseaplot2(
  gsea_result,
  geneSetID    = "HALLMARK_G2M_CHECKPOINT",
  title        = "GSEA: G2M Checkpoint",
  color        = "green3",
  base_size    = 10,
  rel_heights  = c(1.5, 0.3, 0.5),
  subplots     = 1:3,
  pvalue_table = FALSE
)

gseaplot2(
  gsea_result,
  geneSetID    = "HALLMARK_E2F_TARGETS",
  title        = "GSEA: E2F Targets",
  color        = "green3",
  base_size    = 10,
  rel_heights  = c(1.5, 0.3, 0.5),
  subplots     = 1:3,
  pvalue_table = FALSE
)

# Title : RAD52 is a therapeutic vulnerability for acquired PARP inhibitor resistance in ovarian cancer
# Code Authors : Prasanth Kumar Thuthika, Yukihide Ota
# Note: Adjustments were done after generating the figures using ggplot2

# Figure - 5a
results <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/rna_seq_results.xlsx")
colnames(results)
head(results)

logfc_cutoff <- 0.5
fdr_cutoff   <- 0.05

results <- results %>%
  mutate(
    neglog10p = -log10(KO1_Ola_vs_Scr_Ola_padj),
    status = case_when(
      KO1_Ola_vs_Scr_Ola_padj < fdr_cutoff & KO1_Ola_vs_Scr_Ola_log2FC >  logfc_cutoff ~ "Up",
      KO1_Ola_vs_Scr_Ola_padj < fdr_cutoff & KO1_Ola_vs_Scr_Ola_log2FC < -logfc_cutoff ~ "Down",
      TRUE ~ "NotSig"
    )
  )

target_genes <- c("Xrcc3", "Polq", "Rad51c", "Pold4", "Brcc3", "Fen1")

results$label <- ifelse(results$external_gene_name %in% target_genes, 
                        results$external_gene_name, NA)

labeled_genes <- subset(results,
                        !is.na(label) &
                          !is.na(KO1_Ola_vs_Scr_Ola_padj) &
                          !is.na(KO1_Ola_vs_Scr_Ola_log2FC) &
                          KO1_Ola_vs_Scr_Ola_padj < fdr_cutoff &
                          abs(KO1_Ola_vs_Scr_Ola_log2FC) > logfc_cutoff)

max_y <- max(-log10(results$KO1_Ola_vs_Scr_Ola_padj), na.rm = TRUE)
y_max <- ceiling(max_y + 1.5)

ggplot(results, aes(x = KO1_Ola_vs_Scr_Ola_log2FC,
                    y = -log10(KO1_Ola_vs_Scr_Ola_padj),
                    color = status)) +
  geom_point(alpha = 1.5, size = 2.5) +
  
  scale_color_manual(values = c(
    "Up"     = "firebrick",
    "Down"   = "blue",
    "NotSig" = "gray"
  )) +
  
  geom_point(data = labeled_genes,
             shape = 21, color = "black", fill = NA,
             size = 2.5, stroke = 1.0) + 
  
  geom_label_repel(
    data = labeled_genes,
    aes(label = label),
    fontface = "bold",
    color = "black",
    fill = "white",
    label.size = 1.5,
    label.padding = 0.5,
    box.padding = 1.0,
    point.padding = 0.1,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    force = 3,
    force_pull = 1.0,
    min.segment.length = 0,
    size = 3.5
  ) +
  
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 100, by = 4)) +
  
  theme_minimal(base_size = 10) +
  labs(
    title = "KO1_Ola_vs_Scr_Ola",
    x = "log2FC",
    y = expression(bold(-log[10]~"(adj-p)"))
  ) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
             linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(fdr_cutoff),
             linetype = "dashed", color = "gray50") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.text = element_text( size = 10, face = "bold"),
    legend.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.1),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.1)
  )

# Figure - 5b
data <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/Z-Scores_expression.xlsx")

polq_row <- data %>% 
  dplyr::filter(external_gene_name == "Polq")

polq_z_scores <- polq_row %>%
  dplyr::select(starts_with("z.")) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Zscore") %>%
  mutate(
    Group = case_when(
      grepl("Scr-UT", Sample) ~ "Scr_Veh",
      grepl("KO1-UT", Sample) ~ "KO#1_Veh",
      grepl("Scr_ola", Sample) ~ "Scr_Ola",
      grepl("KO1-ola", Sample) ~ "KO#1_Ola",
      TRUE ~ "Other"
    )
  )


polq_z_scores$Group <- factor(polq_z_scores$Group,
                              levels = c("Scr_Veh", "Scr_Ola", "KO#1_Veh", "KO#1_Ola"))

shape_map <- c("Scr_Veh" = 16, "Scr_Ola" = 15, "KO#1_Veh" = 17, "KO#1_Ola" = 18)
color_map <- c("Scr_Veh" = "black", "Scr_Ola" = "red", "KO#1_Veh" = "blue", "KO#1_Ola" = "darkgreen")

means_df <- polq_z_scores %>%
  group_by(Group) %>%
  summarise(mean_z = mean(Zscore), .groups = "drop") %>%
  mutate(x = as.numeric(Group))

# Plot
ggplot(polq_z_scores, aes(x = Group, y = Zscore, color = Group, shape = Group)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 1) +
  geom_point(size = 3, position = position_jitter(width = 0.15)) +
  geom_segment(data = means_df,
               aes(x = x - 0.2, xend = x + 0.2, y = mean_z, yend = mean_z),
               inherit.aes = FALSE,
               color = "black",
               linewidth = 0.8) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(
    x = NULL,
    y = "Expression Z-score"
  ) +
  ggtitle(expression(italic("Polq"))) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(angle = 90, hjust = 1, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    legend.position = "none",
    plot.title = element_text(face = "italic", hjust = 0.5)
  )

# Supplementary Figure S9 - b

results <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/Z-Scores_expression.xlsx")

z_score_cols <- results %>%
  dplyr::select(starts_with("z."))

z_matrix <- t(z_score_cols)
colnames(z_matrix) <- results$ensembl_gene_id 
rownames(z_matrix) <- gsub("z\\.", "", rownames(z_matrix))

# run PCA
pca <- prcomp(z_matrix, center = TRUE, scale. = TRUE)

# plot PCA
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)

pca_df$Group <- case_when(
  str_detect(pca_df$Sample, "Scr-UT") ~ "Scr_Veh",
  str_detect(pca_df$Sample, "KO1-UT") ~ "KO#1_Veh",
  str_detect(pca_df$Sample, "Scr_ola") ~ "Scr_Ola",
  str_detect(pca_df$Sample, "KO1-ola") ~ "KO#1_Ola",
  TRUE ~ "Other"
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c(
    "Scr_Veh" = "black",
    "KO#1_Veh" = "blue",
    "Scr_Ola" = "red",
    "KO#1_Ola" = "green"
  )) +
  scale_fill_manual(values = c(
    "Scr_Veh" = "black",
    "KO#1_Veh" = "blue",
    "Scr_Ola" = "red",
    "KO#1_Ola" = "green"
  )) +
  labs(
    title = "PCA",
    x = paste0("PC1 – ", round(summary(pca)$importance[2, 1] * 100, 2), " % variance"),
    y = paste0("PC2 – ", round(summary(pca)$importance[2, 2] * 100, 2), " % variance")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(face = "bold", size = 12),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    panel.grid = element_blank(),        
    panel.border = element_blank()
  )

# Supplementary Figure  S9 - e
gsea_df <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/GSEA_Hallmark_KO1_UT_vs_Scr_UT.xlsx")
colnames(gsea_df)
head(gsea_df)

top_gsea <- gsea_df %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  mutate(
    Description = fct_reorder(Description, NES),
    direction   = ifelse(NES > 0, "Positive", "Negative")
  )

ggplot(top_gsea, aes(x = NES, y = reorder(Description, NES), fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c("Positive" = "red", "Negative" = "blue"),
    name   = ""
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "KO1_UT_vs_Scr_UT",
    x     = "Normalized Enrichment Score",
    y     = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title.x       = element_text(face = "bold", size = 8),
    axis.title.y       = element_text(face = "bold", size = 8),
    axis.text.x        = element_text(face = "bold", size = 8),
    axis.text.y        = element_text(face = "bold", size = 8),
    legend.title       = element_text(face = "bold", size = 8),
    legend.text        = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    axis.line.x        = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x       = element_line(color = "black", linewidth = 0.5),
    legend.position = "none"
  )

# Supplementary Figure S9 - f
gsea_df <- read_excel("~/Desktop/Khabele_Lab/Yuki-Rad52/output/GSEA_Hallmark_KO1_Ola_vs_Scr_Ola.xlsx")
colnames(gsea_df)
head(gsea_df)

top_gsea <- gsea_df %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  mutate(
    Description = fct_reorder(Description, NES),
    direction   = ifelse(NES > 0, "Positive", "Negative")
  )

ggplot(top_gsea, aes(x = NES, y = reorder(Description, NES), fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c("Positive" = "red", "Negative" = "blue"),
    name   = ""
  ) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "KO1_Ola_vs_Scr_Ola",
    x     = "Normalized Enrichment Score",
    y     = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title.x       = element_text(face = "bold", size = 8),
    axis.title.y       = element_text(face = "bold", size = 8),
    axis.text.x        = element_text(face = "bold", size = 8),
    axis.text.y        = element_text(face = "bold", size = 8),
    legend.title       = element_text(face = "bold", size = 8),
    legend.text        = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    axis.line.x        = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x       = element_line(color = "black", linewidth = 0.5),
    legend.position = "none"
  )

# Supplementary figure S9 - i
# The code below is the same code used for Polq (Fig. 5B).
# Use the same code for the following genes by changing the gene name:
# Ccnd1, Ccna2, Ccne1, Pole, Smc4, and Cdc6.

data <- read_excel("Z-Scores_expression.xlsx")

polq_row <- data %>% 
  dplyr::filter(external_gene_name == "Cdc6")

polq_z_scores <- polq_row %>%
  dplyr::select(starts_with("z.")) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Zscore") %>%
  mutate(
    Group = case_when(
      grepl("Scr-UT", Sample) ~ "Scr_Veh",
      grepl("KO1-UT", Sample) ~ "KO#1_Veh",
      grepl("Scr_ola", Sample) ~ "Scr_Ola",
      grepl("KO1-ola", Sample) ~ "KO#1_Ola",
      TRUE ~ "Other"
    )
  )


polq_z_scores$Group <- factor(polq_z_scores$Group,
                              levels = c("Scr_Veh", "Scr_Ola", "KO#1_Veh", "KO#1_Ola"))

shape_map <- c("Scr_Veh" = 16, "Scr_Ola" = 15, "KO#1_Veh" = 17, "KO#1_Ola" = 18)
color_map <- c("Scr_Veh" = "black", "Scr_Ola" = "red", "KO#1_Veh" = "blue", "KO#1_Ola" = "darkgreen")

means_df <- polq_z_scores %>%
  group_by(Group) %>%
  summarise(mean_z = mean(Zscore), .groups = "drop") %>%
  mutate(x = as.numeric(Group))

# Plot
ggplot(polq_z_scores, aes(x = Group, y = Zscore, color = Group, shape = Group)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 1) +
  geom_point(size = 3, position = position_jitter(width = 0.15)) +
  geom_segment(data = means_df,
               aes(x = x - 0.2, xend = x + 0.2, y = mean_z, yend = mean_z),
               inherit.aes = FALSE,
               color = "black",
               linewidth = 0.8) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(
    x = NULL,
    y = "Expression Z-score"
  ) +
  ggtitle(expression(italic("Cdc6"))) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(angle = 90, hjust = 1, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    legend.position = "none",
    plot.title = element_text(face = "italic", hjust = 0.5)
  )

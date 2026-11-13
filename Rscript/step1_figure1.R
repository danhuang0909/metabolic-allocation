# Figure 1 and Supplementary Figures 1-2.
# Run this script from the project root with Rscript.

project_dir <- normalizePath(
  Sys.getenv("METABOLISM_PROJECT_DIR", unset = getwd()),
  mustWork = TRUE
)
source(file.path(project_dir, "Rscript", "step0_config.R"))
setwd(project_dir)


# 1. Libraries and environment

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(grid)
library(tiff)
library(RColorBrewer)
library(ggsignif)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggtext)

# getFCresult, getcorSpearmanForall, plot_pathway_boxplot

# 2. Shared palettes

class_colors    <- c("Normal" = "#00BFFF", "Tumor" = "#EE7621")
dataset_colors  <- c("CRC-SG1" = "#7FDBB3", "KUL3" = "#F5AD75", "KUL5" = "#85A4E6",
                     "Pelka" = "#F593CE", "SMC" = "#CAF77B")
method_colors   <- c("Absolute pathway activity" = "#74A9CF", "Pathway allocation" = "#EF8A62")
method_colors_lb <- setNames(method_colors, c("Absolute pathway\nactivity", "Pathway\nallocation"))
celltype_colors <- c("Epithelial" = "#BEBADA", "B cell" = "#8DD3C7", "T cell" = "#B3DE69",
                     "Myeloid" = "#80B1D3", "Fibroblast" = "#FB8072",
                     "Endothelial" = "#EEE685", "Plasma" = "#FDB462")
celltype_pal_set3 <- brewer.pal(8, "Set3")[c(1:4, 8, 5:7)]
epi_color <- "#C41E3A"

known_celltypes <- c("B cell", "T cell", "Myeloid", "Epithelial",
                     "Fibroblast", "Endothelial", "Plasma")

# 3. Helper functions

# Add UMAP axis arrows and labels to a plot
add_umap_axes <- function(p, umap, arrow_len = 5, label_size = 4.7) {
  x0 <- min(umap$UMAP_1)
  y0 <- min(umap$UMAP_2)
  p +
    annotate("segment", x = x0, y = y0, xend = x0 + arrow_len, yend = y0,
             colour = "black", linewidth = 1, arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + arrow_len,
             colour = "black", linewidth = 1, arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text", x = x0 + 4, y = y0 - 1.5, label = "UMAP 1", size = label_size) +
    annotate("text", x = x0 - 1.5, y = y0 + 4, label = "UMAP 2", size = label_size, angle = 90)
}

# Match cell type names inside sample ids
get_celltype <- function(ids) {
  sapply(ids, function(x) {
    for (ct in known_celltypes) if (grepl(ct, x, fixed = TRUE)) return(ct)
    NA
  })
}

# Parse ids of the form "<dataset>_<condition>_<cell type>_<patient>"
extract_sample_meta <- function(ids) {
  do.call(rbind, lapply(ids, function(x) {
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    data.frame(dataset = parts[1], condition = parts[2],
               celltype = paste(parts[3:(length(parts) - 1)], collapse = " "),
               patient = parts[length(parts)])
  }))
}

# Parse ids of the form "<dataset>_<cell type>_<patient>" (FC matrices, no condition)
extract_fc_meta <- function(ids) {
  do.call(rbind, lapply(ids, function(x) {
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    data.frame(dataset = parts[1],
               celltype = paste(parts[2:(length(parts) - 1)], collapse = " "),
               patient = parts[length(parts)])
  }))
}

# Per-cell-type pairwise Spearman correlations between pathway log2FCs
calc_pairwise_cor <- function(fc_df, celltype_col = "celltype") {
  pathway_cols <- setdiff(names(fc_df), celltype_col)
  do.call(rbind, lapply(unique(fc_df[[celltype_col]]), function(ct) {
    cor_mat <- cor(fc_df[fc_df[[celltype_col]] == ct, pathway_cols],
                   method = "spearman", use = "pairwise.complete.obs")
    data.frame(celltype = ct,
               pair_id = seq_len(sum(upper.tri(cor_mat))),
               rho = cor_mat[upper.tri(cor_mat)])
  }))
}

# Long-format pathway allocation plus per-sample Normal/Tumor rank differences
prepare_rank_data <- function(pathway_mean_df) {
  data <- cbind(as.data.frame(pathway_mean_df),
                extract_sample_meta(rownames(pathway_mean_df)))
  pathway_cols <- setdiff(names(data), c("dataset", "condition", "celltype", "patient"))

  all_long <- data %>%
    pivot_longer(all_of(pathway_cols), names_to = "pathway", values_to = "Gp") %>%
    mutate(condition = factor(condition, levels = c("Normal", "Tumor")))

  rank_change <- all_long %>%
    group_by(celltype, patient, condition) %>%
    mutate(rank = rank(-Gp)) %>%
    ungroup() %>%
    dplyr::select(dataset, celltype, patient, pathway, condition, rank) %>%
    pivot_wider(names_from = condition, values_from = rank) %>%
    mutate(rank_diff = Tumor - Normal, abs_rank_diff = abs(rank_diff)) %>%
    filter(!is.na(rank_diff))  # drop samples with only one condition

  list(data = data, pathway_cols = pathway_cols,
       all_long = all_long, rank_change = rank_change)
}

# X-axis labels with Epithelial highlighted in bold red (needs element_markdown)
epi_x_labels <- function(x) {
  ifelse(x == "Epithelial", "<b style='color:#C41E3A'>Epithelial</b>", x)
}

# Violin plots of a QC metric per patient, one panel per dataset
plot_vln_by_dataset <- function(meta, metric, ylab) {
  patient_counts <- meta %>%
    distinct(dataset, patient) %>%
    count(dataset, name = "patient_count")
  n_ds <- nrow(patient_counts)

  plots <- lapply(seq_len(n_ds), function(i) {
    ds <- patient_counts$dataset[i]
    p <- ggplot(meta %>% filter(dataset == ds),
                aes(x = patient, y = .data[[metric]], fill = class_new)) +
      geom_violin(width = 0.8, position = position_dodge(width = 0.8), scale = "width") +
      theme_hd_minimal_plain() +
      theme(axis.text.x = element_blank(),
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0.5, "cm"),
            plot.margin = margin(t = 3, l = if (i == 1) 3 else 0,
                                 b = 3, r = if (i == n_ds) 3 else 0)) +

      labs(title = ds, x = if (i == 4) "Patient" else "", y = if (i == 1) ylab else "") +
      scale_fill_manual(name = "Class", values = class_colors)
    if (i != 1) {
      p <- p + theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.y = element_blank())
    }
    p
  })
  ggarrange(plotlist = plots, nrow = 1, common.legend = TRUE,
            legend = "right", widths = patient_counts$patient_count)
}

# Paired-line plot of pathway allocation (Normal vs Tumor), highlighting rank shifts
plot_pathway_shift <- function(df, facet_by, title, free_y = FALSE) {
  ggplot(df, aes(x = condition, y = Gp, group = pathway)) +
    geom_line(data = df %>% filter(!highlight),
              color = "#E8E8E8", alpha = 0.6, linewidth = 0.5) +
    geom_point(data = df %>% filter(!highlight),
               aes(color = condition), size = 1.5, alpha = 0.4) +
    geom_line(data = df %>% filter(highlight),
              color = "#E41A1C", alpha = 0.8, linewidth = 1.3) +
    geom_point(data = df %>% filter(highlight),
               aes(color = condition), size = 2.5, alpha = 0.8) +
    geom_text_repel(data = df %>% filter(condition == "Normal", !highlight),
                    aes(label = pathway), direction = "y", nudge_x = -0.35, size = 2.5,
                    color = "#666666", segment.linewidth = 0.2, segment.color = "#CCCCCC",
                    max.overlaps = 15, min.segment.length = 0) +
    scale_color_manual(values = c("Normal" = "#74A9CF", "Tumor" = "#EF8A62"), guide = "none") +
    facet_wrap(as.formula(paste("~", facet_by)), nrow = 1,
               scales = if (free_y) "free_y" else "fixed") +
    scale_y_log10(labels = label_number()) +
    labs(x = NULL, y = "Global pathway allocation", title = title) +
    theme_hd_minimal_plain() +
    theme(strip.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, hjust = 0.5),
          panel.spacing = unit(1.5, "lines"))
}

# Format a Spearman correlation annotation from two vectors
cor_label <- function(x, y) {
  ct <- cor.test(x, y, method = "spearman")
  p_str <- if (ct$p.value < 2.2e-16) "p < 2.2e-16" else paste("p =", format(ct$p.value, digits = 2))
  paste0("R = ", round(unname(ct$estimate), 2), ", ", p_str)
}

# 4. Data loading

meta_data_all <- readRDS("./data/all_dataset_SCT_filter_meata_data.rds")
combined.sct_filter <- readRDS("./data/all_meta_6cohorts_NGpaper_integrated_result_meta_genes_filtered.rds")
all_gene_pathway <- readRDS("./data/all_genes_pathway_infor_select.rds")
allcells_pathwayMeanSCTMainClass <- readRDS("./data/allcells_pathwayMeanSCTMainClass.Rds")
allcells_pathwayMeanSCTMainClassNormalized <- readRDS("./data/allcells_pathwayMeanSCTMainClassNormalizedByMetacount.Rds")
pathway_mean_patientClassNormalized<-getAggregateDate2(allcells_pathwayMeanSCTMainClassNormalized,meta_data_all[meta_data_all$Cell_type!="Mast",])
pathway_mean_ClassNormalized<-getAggregateDate(allcells_pathwayMeanSCTMainClassNormalized,meta_data_all[meta_data_all$Cell_type!="Mast",])


# Add the proper readRDS()/load() calls here, otherwise a fresh knit will fail:
all_dataset_SCT_filter <- readRDS("./data/all_dataset_SCT_filter.rds")

# 5. Fig 1a - flowchart

img <- rasterGrob(readTIFF("./figures/figure1a_02.tif", native = TRUE), interpolate = TRUE)

p1_flowchart <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
  geom_blank() +
  annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

# 6. Fig 1b - UMAP by cell type

umap <- combined.sct_filter@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type = combined.sct_filter@meta.data$Cell_type,
        class_new = combined.sct_filter@meta.data$class_new)

cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p2_umap <- DimPlot(combined.sct_filter, label = FALSE, group.by = "Cell_type") +
  ggtitle("Metabolic genes") +
  scale_color_manual(name = "Cell type", values = celltype_pal_set3) +
  theme_void(base_size = 14, base_family = "Helvetica")
p2_umap <- add_umap_axes(p2_umap, umap) +
  geom_text(aes(label = cell_type, x = UMAP_1, y = UMAP_2, size = I(4.7), fontface = "plain"),
            data = cell_type_med) +
  theme(plot.title = element_text(lineheight = 0.8, size = 14, hjust = 0.5),
        plot.margin = unit(c(0, 0, 0, 1), "cm"),
        legend.position = "none")

# 7. Correlation heatmap data (used in Supp Fig 1)

metagenes <- unique(all_gene_pathway$gene)
all_genes <- rownames(all_dataset_SCT_filter@assays$SCT@counts)
raw_data <- as.matrix(all_dataset_SCT_filter@assays$SCT@counts[all_genes[all_genes %in% metagenes], ])
merge_metaSCT_raw <- all_dataset_SCT_filter@meta.data[colnames(raw_data), ]

# Raw mean expression per cell type x dataset x class
mean_genes_raw <- aggregate(t(raw_data),
                            list(s = as.character(merge_metaSCT_raw$Cell_type),
                                 p = merge_metaSCT_raw$dataset,
                                 type = merge_metaSCT_raw$class_new), mean)
mean_genes_df_raw <- mean_genes_raw[, -c(1:3)]
rownames(mean_genes_df_raw) <- paste(mean_genes_raw[, 1], mean_genes_raw[, 2], mean_genes_raw[, 3], sep = "_")
cor_raw <- plotCorHeatmap3(mean_genes_df_raw, "Correlation")

# Per-gene metabolic fraction (expression / total metabolic counts per cell)

mean_genes_df <- aggregate(
  t(sweep(raw_data, 2, colSums(raw_data), `/`)),
  list(s = as.character(merge_metaSCT_raw$Cell_type),
       p = merge_metaSCT_raw$dataset,
       type = merge_metaSCT_raw$class_new), mean)
rownames(mean_genes_df) <- do.call(paste, c(mean_genes_df[, 1:3], sep = "_"))
cor_norm <- plotCorHeatmap2(mean_genes_df[, -c(1:3)], "Per-gene metabolic fraction")
cor_norm


saveRDS(list(mean_genes_df_raw = mean_genes_df_raw, cor_raw = cor_raw,
             mean_genes_df = mean_genes_df, cor_meta = cor_norm),
        file = "./data/heatmap.rds")
cor_result <- readRDS("./data/heatmap.rds")

# 8. Fig 1c - total metabolic load per cell type

meta_aggre2 <- aggregate(meta_data_all$metagene_count,
                         list(patient = meta_data_all$patient,
                              type = meta_data_all$class_new,
                              cell = meta_data_all$Cell_type,
                              dataset = meta_data_all$dataset), mean)

p_mat <- ggplot(meta_aggre2[meta_aggre2$cell != "Mast", ], aes(x = type, y = x, fill = type)) +
  geom_jitter(aes(color = dataset), size = 2, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE) +
  geom_signif(comparisons = list(c("Normal", "Tumor")), map_signif_level = TRUE, color = "black") +
  theme_hd_minimal_plain() +
  theme(strip.text.y = element_text(angle = 0),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14),
        legend.position = "bottom") +
  scale_fill_manual(name = "Sample Type", values = class_colors) +
  scale_color_manual(name = "Dataset", values = dataset_colors) +
  guides(fill = "none", color = guide_legend(title = "Dataset")) +
  ylab("Total metabolic load") +
  facet_wrap(~cell, ncol = 1, strip.position = "right") +
  coord_flip()

# 9. Tumor-vs-Normal log2FC of total metabolic load

metabolic_activity <- aggregate(meta_data_all$metagene_count,
                                list(dataset = meta_data_all$dataset,
                                     cell = meta_data_all$Cell_type,
                                     patient = meta_data_all$patient,
                                     type = meta_data_all$class_new), mean)
metabolic_activity_normal <- metabolic_activity[metabolic_activity$type == "Normal" & metabolic_activity$cell != "Mast", ]
metabolic_activity_tumor  <- metabolic_activity[metabolic_activity$type == "Tumor"  & metabolic_activity$cell != "Mast", ]
rownames(metabolic_activity_normal) <- paste(metabolic_activity_normal$dataset, metabolic_activity_normal$cell,
                                             metabolic_activity_normal$patient, sep = "_")
rownames(metabolic_activity_tumor)  <- paste(metabolic_activity_tumor$dataset, metabolic_activity_tumor$cell,
                                             metabolic_activity_tumor$patient, sep = "_")
# Pairs missing a matched Normal sample become NA and are dropped
metabolic_activity_tumor$fc <- log2(metabolic_activity_tumor$x /
                                      metabolic_activity_normal[rownames(metabolic_activity_tumor), "x"])
metabolic_activity_tumor <- metabolic_activity_tumor[!is.na(metabolic_activity_tumor$fc), ]

# Pathway-level FC results (absolute and normalized), Mast cells excluded
meta_data <- meta_data_all[meta_data_all$Cell_type != "Mast", ]
fc_mainclass <- getFCresult(allcells_pathwayMeanSCTMainClass, meta_data)
fc_mainclassNormalized <- getFCresult(allcells_pathwayMeanSCTMainClassNormalized, meta_data)

xx_reSCTClass <- getcorSpearmanForall(metabolic_activity_tumor[, "fc", drop = FALSE], fc_mainclass$FC)
xx_reSCTClassNormalized <- getcorSpearmanForall(metabolic_activity_tumor[, "fc", drop = FALSE], fc_mainclassNormalized$FC)

# 10. Fig 1e - total load FC vs pathway FC (Epithelial)

total_fc_df <- metabolic_activity_tumor
total_fc_df$sample_id <- rownames(total_fc_df)

pathway_focus <- "Canonical amino acid metabolism"
build_plot_df <- function(fc_result) {
  df <- data.frame(Canonical_fc = fc_result$FC[, pathway_focus],
                   sample_id = rownames(fc_result$FC), stringsAsFactors = FALSE)
  df <- merge(total_fc_df[, c("sample_id", "cell", "dataset", "fc")], df, by = "sample_id")
  df[is.finite(df$fc) & is.finite(df$Canonical_fc), ]
}
plot_epi_abs <- build_plot_df(fc_mainclass)           %>% filter(cell == "Epithelial") %>% mutate(version = "Absolute pathway activity")
plot_epi_sp  <- build_plot_df(fc_mainclassNormalized) %>% filter(cell == "Epithelial") %>% mutate(version = "Pathway allocation")
plot_combined <- rbind(plot_epi_abs, plot_epi_sp)

label_df <- data.frame(
  version = c("Absolute pathway activity", "Pathway allocation"),
  label = c(cor_label(plot_epi_abs$fc, plot_epi_abs$Canonical_fc),
            cor_label(plot_epi_sp$fc, plot_epi_sp$Canonical_fc)),
  x = -0.2, y = 0.8,
  color = c("darkred", "darkgreen")
)

p <- ggplot(plot_combined, aes(x = fc, y = Canonical_fc)) +
  geom_point(aes(color = dataset), size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.8, linetype = "dashed") +
  geom_text(data = label_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 4, color = label_df$color) +
  facet_wrap(~version, ncol = 2, scales = "fixed") +
  scale_color_manual(values = dataset_colors) +
  labs(x = "Total metabolic load log2FC (T vs N) in Epithelial cells",
       y = "Canonical amino acid metabolism\nlog2FC (T vs N)",
       color = "Dataset") +
  theme_hd_minimal_plain() +
  theme(legend.position = "bottom")
p

# 11. Fig 1d - slope chart of pathway-load correlations

cor_by_method <- rbind(
  transform(xx_reSCTClass, method = "Absolute pathway\nactivity"),
  transform(xx_reSCTClassNormalized, method = "Pathway\nallocation")
)
mean_data <- cor_by_method %>%
  group_by(cell, method) %>%
  summarise(mean_rho = mean(rho.rho), .groups = "drop")
cell_order <- mean_data %>%
  filter(method == "Absolute pathway\nactivity") %>%
  arrange(desc(mean_rho)) %>% pull(cell)
mean_data$cell <- factor(mean_data$cell, levels = cell_order)

p_slope <- ggplot(mean_data, aes(x = method, y = mean_rho, group = cell)) +
  geom_line(color = "#BDBDBD", alpha = 0.5, linewidth = 0.8) +
  geom_point(aes(color = method), size = 3.5) +
  scale_color_manual(values = method_colors_lb, guide = "none") +
  geom_text_repel(data = mean_data %>% filter(method == "Absolute pathway\nactivity"),
                  aes(label = cell), nudge_x = -0.3, size = 3.2) +
  labs(x = NULL, y = "Mean correlation coefficient\nof log2FC(T vs N)") +
  ggtitle("Global pathway correlation\nwith total metabolic load") +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_slope

# 12. Inter-pathway correlation per cell type (bar + boxplot versions)

ap_fc <- as.data.frame(fc_mainclass$FC)           %>% mutate(celltype = get_celltype(rownames(.))) %>% filter(!is.na(celltype))
sp_fc <- as.data.frame(fc_mainclassNormalized$FC) %>% mutate(celltype = get_celltype(rownames(.))) %>% filter(!is.na(celltype))

ap_pairwise <- calc_pairwise_cor(ap_fc) %>% mutate(method = "Absolute pathway activity")
sp_pairwise <- calc_pairwise_cor(sp_fc) %>% mutate(method = "Pathway allocation")
all_pairwise <- rbind(ap_pairwise, sp_pairwise) %>%
  mutate(method = factor(method, levels = names(method_colors)))

cell_order_pw <- ap_pairwise %>%
  group_by(celltype) %>% summarise(median_rho = median(rho), .groups = "drop") %>%
  arrange(desc(median_rho)) %>% pull(celltype)
all_pairwise$celltype <- factor(all_pairwise$celltype, levels = cell_order_pw)

# Mean-correlation bar chart (exploratory; the boxplot below is the one used)
cor_summary <- all_pairwise %>%
  group_by(celltype, method) %>%
  summarise(mean_cor = mean(rho, na.rm = TRUE), .groups = "drop")
cell_order_bar <- cor_summary %>%
  filter(method == "Absolute pathway activity") %>%
  arrange(desc(mean_cor)) %>% pull(celltype)
cor_summary$celltype <- factor(cor_summary$celltype, levels = cell_order_bar)

p_bar <- ggplot(cor_summary, aes(fill = method, y = mean_cor, x = celltype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = method_colors) +
  labs(x = NULL, y = "Mean correlation coefficient\nof log2FC(T vs N)", fill = "Method") +
  ggtitle("Inter-pathway correlation (global)") +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(nrow = 2))
p_bar

# Boxplot of all pairwise correlations (used in Supp Fig 1)
p_supp <- ggplot(all_pairwise, aes(x = celltype, y = rho, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.8), size = 1.2, alpha = 0.4, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = method_colors) +
  labs(x = NULL, y = "Mean correlation coefficient\nof log2FC(T vs N)", fill = "Method") +
  ggtitle("Inter-pathway correlation (global)") +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(fill = guide_legend(nrow = 2))
p_supp

# 13. Fig 1f - mean absolute pathway change per cell type

ap_fc_full <- cbind(as.data.frame(fc_mainclass$FC), extract_fc_meta(rownames(fc_mainclass$FC)))
sp_fc_full <- cbind(as.data.frame(fc_mainclassNormalized$FC), extract_fc_meta(rownames(fc_mainclassNormalized$FC)))
pathway_cols <- setdiff(names(ap_fc_full), c("dataset", "celltype", "patient"))

combined <- inner_join(
  ap_fc_full %>% pivot_longer(all_of(pathway_cols), names_to = "pathway", values_to = "ap_fc"),
  sp_fc_full %>% pivot_longer(all_of(pathway_cols), names_to = "pathway", values_to = "sp_fc"),
  by = c("dataset", "celltype", "patient", "pathway")
) %>%
  group_by(celltype, pathway) %>%
  summarise(ap_fc = mean(ap_fc, na.rm = TRUE), sp_fc = mean(sp_fc, na.rm = TRUE), .groups = "drop")

mean_fc_by_celltype <- combined %>%
  group_by(celltype) %>%
  summarise(mean_ap_fc = mean(abs(ap_fc), na.rm = TRUE),
            mean_sp_fc = mean(abs(sp_fc), na.rm = TRUE), .groups = "drop")

mean_fc_long <- mean_fc_by_celltype %>%
  pivot_longer(c(mean_ap_fc, mean_sp_fc), names_to = "method", values_to = "mean_fc") %>%
  mutate(method = ifelse(method == "mean_ap_fc", "Absolute pathway\nactivity", "Pathway\nallocation"),
         celltype = factor(celltype, levels = arrange(mean_fc_by_celltype, desc(mean_sp_fc))$celltype))

x_colors <- ifelse(levels(mean_fc_long$celltype) == "Epithelial", epi_color, "black")
x_faces  <- ifelse(levels(mean_fc_long$celltype) == "Epithelial", "bold", "plain")

p2a_bar <- ggplot(mean_fc_long, aes(x = celltype, y = mean_fc, fill = celltype)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  facet_wrap(~method, ncol = 2, scales = "free_y") +
  labs(y = expression("Mean absolute change of global pathways")) +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = x_colors, face = x_faces),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
p2a_bar

# 14. Fig 1g - allocation heatmap (Normal) + log2FC heatmap

normal_mean <- pathway_mean_ClassNormalized %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(grepl("_Normal_", sample)) %>%
  mutate(celltype = sub(".*_Normal_", "", sample)) %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("celltype") %>%
  t()  # pathways x cell types
normal_mean <- normal_mean[rownames(normal_mean) != "sample", ]

log2fc_t <- as.matrix(t(fc_mainclassNormalized$FC_sum))  # pathways x cell types
log2fc_t <- log2fc_t[rownames(normal_mean), colnames(normal_mean)]

# Order columns by mean |log2FC| (desc), rows by Epithelial |log2FC| (desc)
mean_abs_fc <- colMeans(abs(log2fc_t))
celltypes_reordered <- names(sort(mean_abs_fc, decreasing = TRUE))
pathway_order <- names(sort(abs(log2fc_t[, "Epithelial"]), decreasing = TRUE))
normal_mean <- normal_mean[pathway_order, celltypes_reordered]
log2fc_t <- log2fc_t[pathway_order, celltypes_reordered]
mean_abs_fc <- mean_abs_fc[celltypes_reordered]

normal_log <- log2(normal_mean * 10^4)
col_normal <- colorRamp2(
  quantile(normal_log, c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE),
  c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#F16913", "#D94801"))
fc_max <- max(abs(as.vector(log2fc_t)), na.rm = TRUE)
col_fc <- colorRamp2(
  c(-fc_max, -fc_max * 0.6, -fc_max * 0.3, 0, fc_max * 0.3, fc_max * 0.6, fc_max),
  c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B"))

top_anno <- HeatmapAnnotation(
  `Mean\n|Log2FC|` = anno_simple(mean_abs_fc,
                                col = colorRamp2(c(0, max(mean_abs_fc)), c("#F7F7F7", "#B2182B")),
                                height = unit(0.5, "cm")),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 14),
  annotation_name_side = "right",
  height = unit(0.7, "cm"))

column_name_colors <- ifelse(celltypes_reordered == "Epithelial", epi_color, "black")
column_name_faces <- ifelse(celltypes_reordered == "Epithelial", "bold", "plain")
column_gp <- gpar(fontsize = 14, col = column_name_colors, fontface = column_name_faces)
rownames(normal_log)<-gsub("Non-canonical amino acids Metabolism","Non-canonical amino acids metabolism",rownames(normal_log))
ht1 <- Heatmap(normal_log, name = "Normal", col = col_normal,
               cluster_columns = FALSE, cluster_rows = FALSE,
               show_row_names = TRUE, show_column_names = TRUE,
               row_names_side = "left", row_names_gp = gpar(fontsize = 14),
               column_names_gp = column_gp,
               column_title = "Normal(Allocation)", column_title_gp = gpar(fontsize = 14),
               width = unit(4.2, "cm"), height = unit(8, "cm"))

ht2 <- Heatmap(log2fc_t, name = "Log2FC", col = col_fc, top_annotation = top_anno,
               cluster_columns = FALSE, cluster_rows = FALSE,
               show_row_names = FALSE, show_column_names = TRUE,
               column_names_gp = column_gp,
               column_title = "Tumor vs Normal", column_title_gp = gpar(fontsize = 14),
               width = unit(5.5, "cm"), height = unit(8, "cm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 val <- log2fc_t[i, j]
                 if (abs(val) > 0.12) {
                   grid.text(sprintf("%.2f", val), x, y,
                             gp = gpar(fontsize = 7.5,
                                       col = ifelse(abs(val) > fc_max * 0.5, "white", "black")))
                 }
               })

ht_grob <- grid.grabExpr(draw(ht1 + ht2))

# 15. Pathway rank-shift analysis (patient level)

rank_data <- prepare_rank_data(pathway_mean_patientClassNormalized)
data <- rank_data$data
pathway_cols <- rank_data$pathway_cols
all_long <- rank_data$all_long
rank_change <- rank_data$rank_change

# Per cell type x pathway: proportion of patients with |rank shift| > 2
proportion_matrix <- rank_change %>%
  group_by(celltype, pathway) %>%
  summarise(n_total = n(),
            n_sig = sum(abs_rank_diff > 2),
            proportion = n_sig / n_total,
            .groups = "drop")

# Pathway order: ascending Epithelial proportion; cell type order: descending sum of proportions
pathway_order <- proportion_matrix %>%
  filter(celltype == "Epithelial") %>%
  arrange(proportion) %>% pull(pathway)
celltype_mean <- proportion_matrix %>%
  group_by(celltype) %>%
  summarise(mean_prop = sum(proportion), .groups = "drop") %>%
  arrange(desc(mean_prop))
celltype_order <- celltype_mean$celltype

proportion_matrix <- proportion_matrix %>%
  mutate(pathway = factor(pathway, levels = pathway_order),
         celltype = factor(celltype, levels = celltype_order))

# 16. Fig 1h - diverging bar chart (Epithelial vs other cell types)

diverging_data <- proportion_matrix %>%
  group_by(pathway) %>%
  summarise(epi_value = proportion[celltype == "Epithelial"],
            max_other_value = max(proportion[celltype != "Epithelial"]),
            max_other_celltype = celltype[celltype != "Epithelial"][which.max(proportion[celltype != "Epithelial"])],
            .groups = "drop") %>%
  arrange(epi_value) %>%
  mutate(pathway = factor(pathway, levels = pathway),
         epi_pct = epi_value * 100,
         other_pct = -max_other_value * 100)
n_pathways <- length(levels(diverging_data$pathway))

right_df <- proportion_matrix %>%
  group_by(celltype) %>%
  summarise(sum_prop = sum(proportion, na.rm = TRUE), .groups = "drop") %>%
  arrange(sum_prop) %>%
  mutate(y_scaled = row_number(),
         bg_color = {
           v <- sum_prop
           if (length(v) == 1 || max(v) == min(v)) rep("#FEE2E2", length(v))
           else rgb(colorRamp(c("#FEE2E2", epi_color))((v - min(v)) / (max(v) - min(v))) / 255)
         },
         name_color = ifelse(celltype == "Epithelial", epi_color, "black"),
         name_face = ifelse(celltype == "Epithelial", "bold", "plain"))
n_ct <- nrow(right_df)

p_main <- ggplot(diverging_data) +
  geom_col(aes(x = other_pct, y = pathway),
           fill = "#D9D9D9", width = 0.6, color = "black", linewidth = 0.3) +
  geom_col(aes(x = epi_pct, y = pathway, fill = epi_pct),
           width = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = "#E8E4F3", high = "#5E5493", guide = "none") +
  geom_text(aes(x = other_pct - 2, y = pathway, label = max_other_celltype),
            hjust = 1, size = 2.8, fontface = "bold", color = "#333333") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  scale_y_discrete(expand = expansion(add = c(0, 2))) +
  scale_x_continuous(labels = function(x) paste0(abs(x), "%"),
                     limits = c(-105, 75), breaks = seq(-100, 100, 25), expand = c(0, 0)) +
  labs(x = "% samples with priority disruption", y = "Global pathways") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, b = 10, l = 10, r = 10)) +
  annotate("text", x = 23, y = n_pathways + 1,
           label = "Epithelial", size = 4.5, fontface = "bold", color = "#5E5493") +
  annotate("text", x = -50, y = n_pathways + 1,
           label = "Max in other cell types", size = 4.5, fontface = "bold", color = "#666666")

p_right <- ggplot(right_df, aes(x = 1, y = y_scaled)) +
  geom_tile(aes(fill = I(bg_color)), width = 1.3, color = "white", linewidth = 3) +
  geom_text(aes(label = round(sum_prop, 2)), size = 4) +
  geom_text(aes(x = -0.05, y = y_scaled, label = celltype,
                color = I(name_color), fontface = name_face),
            hjust = 1, size = 4) +
  scale_x_continuous(limits = c(-3.5, 1.8), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, n_ct + 2.5), expand = c(0, 0)) +
  annotate("text", x = 0, y = n_ct + 1, label = "Sum proportion", size = 4) +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(l = 5, r = 5, b = 20))

p_percent <- p_main + p_right + plot_layout(widths = c(2.3, 1))
p_percent

# 17. Fig 1i - example patient C157 (Epithelial)

c157_long <- data %>%
  filter(patient == "C157", celltype == "Epithelial") %>%
  pivot_longer(all_of(pathway_cols), names_to = "pathway", values_to = "Gp") %>%
  mutate(condition = factor(condition, levels = c("Normal", "Tumor")))

c157_highlight <- c157_long %>%
  group_by(celltype, condition) %>%
  mutate(rank = rank(-Gp)) %>%
  ungroup() %>%
  select(celltype, pathway, condition, rank) %>%
  pivot_wider(names_from = condition, values_from = rank) %>%
  mutate(rank_diff = abs(Tumor - Normal), highlight = rank_diff > 2) %>%
  select(celltype, pathway, rank_diff, highlight)

summary_data <- c157_long %>%
  left_join(c157_highlight, by = c("celltype", "pathway")) %>%
  group_by(celltype, pathway, condition, highlight) %>%
  summarise(mean_Gp = mean(Gp, na.rm = TRUE), .groups = "drop")

# Order pathways by Normal allocation (desc); shorten one long pathway name
pathway_order_c157 <- summary_data %>%
  filter(condition == "Normal") %>%
  arrange(desc(mean_Gp)) %>% pull(pathway)
summary_data$pathway <- factor(summary_data$pathway, levels = pathway_order_c157)
summary_data$pathway <- recode(summary_data$pathway,
                               "Metabolism of cofactors and vitamins" = "Cofactors and vitamins")

p_c157 <- ggplot(summary_data, aes(x = condition, y = mean_Gp, group = pathway)) +
  geom_line(data = summary_data %>% filter(!highlight),
            color = "#D9D9D9", alpha = 0.5, linewidth = 0.6) +
  geom_point(data = summary_data %>% filter(!highlight),
             aes(fill = condition), size = 2.5, shape = 21, color = "#D9D9D9", stroke = 0.3) +
  geom_line(data = summary_data %>% filter(highlight),
            color = "#E41A1C", alpha = 0.8, linewidth = 1.2) +
  geom_point(data = summary_data %>% filter(highlight),
             aes(fill = condition), size = 3.5, shape = 21, color = "#E41A1C", stroke = 0.5) +
  scale_fill_manual(values = c("Normal" = "#74A9CF", "Tumor" = "#EF8A62"), guide = "none") +
  geom_text_repel(data = summary_data %>% filter(condition == "Normal", highlight),
                  aes(label = pathway), direction = "y", nudge_x = -0.3, size = 4,
                  color = "#666666", fontface = "bold",
                  segment.linewidth = 0.2, segment.color = "#BDBDBD") +
  scale_y_log10(labels = label_number()) +
  labs(x = NULL, y = "Global pathway allocation",
       title = "Example: Pathway priority shift\nin epithelial cells (Patient C157)") +
  theme_hd_minimal_plain()
p_c157

# 18. Assemble Figure 1

pall1_0 <- ggarrange(p2_umap + guides(color = "none"), p_mat,
                     labels = letters[2:3], nrow = 2, ncol = 1,
                     font.label = list(size = 18, color = "black", face = "bold"),
                     heights = c(1, 1.8))
pall1_1 <- ggarrange(p1_flowchart, pall1_0,
                     labels = c(letters[1], ""), nrow = 1, ncol = 2,
                     font.label = list(size = 18, color = "black", face = "bold"),
                     widths = c(1.6, 1))
pall1_2 <- ggarrange(p_slope, p, ht_grob,
                     labels = letters[c(4:5, 7)], nrow = 1, ncol = 3,
                     font.label = list(size = 18, color = "black", face = "bold"),
                     widths = c(1, 1.7, 2.3))
pall1_3 <- ggarrange(p2a_bar, p_percent, p_c157,
                     labels = letters[c(6, 8, 9)], nrow = 1, ncol = 3,
                     font.label = list(size = 18, color = "black", face = "bold"),
                     widths = c(1.2, 2.3, 1))
pall1_4 <- ggarrange(pall1_1, pall1_2, pall1_3, nrow = 3, ncol = 1, heights = c(2, 1, 1))

save_tiff(pall1_4, "figure1.tiff", width = 19, height = 21)

# 19. Supp Fig 1e - marker gene dotplot

# Marker list is written directly in display order (no index-based reordering)
marker_genes <- list(
  "B cell"      = c("MS4A1", "CD19"),
  "Plasma"      = c("CD79A", "MZB1", "JCHAIN"),
  "T cell"      = c("CD3E", "CD3D", "CD3G"),
  "Epithelial"  = c("EPCAM", "MUC1", "CDH1"),
  "Mast"        = c("HDC", "CPA3", "TPSAB1"),
  "Myeloid"     = c("FCGR3A", "CD14", "CD163"),
  "Endothelial" = c("PECAM1", "VWF", "CLDN5"),
  "Fibroblast"  = c("FBLN1", "PDGFRA", "FGF7")
)

all_dataset_SCT_filter$cell_clustered <- factor(all_dataset_SCT_filter$cell_clustered,
                                                levels = names(marker_genes))


suppressMessages(
  trace(ggplot2::facet_grid, quote({
    if (!is.null(facets)) { rows <- facets; facets <- NULL }
  }), print = FALSE)
)


p2_dot <- DotPlot(object = all_dataset_SCT_filter, features = marker_genes,
                  group.by = "cell_clustered", cols = c("lightgrey", "#EE7621")) +
  theme_hd_minimal_plain() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5,size = 10),
        strip.text = element_text(angle = 60),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title = element_blank(),
        legend.position = "right") +
  ggtitle("") +
  guides(size = guide_legend(title = "Percent\nexpressed"))
p2_dot$guides$colour$title <- "Average\nexpression"


# 20. Supp Fig 1a-c - QC violin plots

p_vln1 <- plot_vln_by_dataset(all_dataset_SCT_filter@meta.data, "nFeature_RNA", "Number of detected genes")
p_vln2 <- plot_vln_by_dataset(all_dataset_SCT_filter@meta.data, "nCount_RNA", "UMI counts")
p_vln3 <- plot_vln_by_dataset(all_dataset_SCT_filter@meta.data, "percent.mt", "Percentage of\nmitochondrial genes")

# 21. Supp Fig 1d - UMAP colored by patient

umap_patient <- combined.sct_filter@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(dataset = combined.sct_filter@meta.data$dataset,
        patient = combined.sct_filter@meta.data$patient)


patient_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffbb78", "#ff9896",
                    "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#f5b0d3",
                    "#c7c7c7", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffbb78",
                    "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
                    "#f5b0d3", "#c7c7c7", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                    "#ffbb78", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d",
                    "#9edae5", "#f5b0d3")

p_patient_umap <- ggplot(umap_patient, aes(x = UMAP_1, y = UMAP_2, color = patient)) +
  geom_point(size = 0.001) +
  theme_hd_plain() +
  facet_wrap(~dataset, nrow = 1) +
  scale_color_manual(values = patient_colors) +
  guides(color = guide_legend(nrow = 4, override.aes = list(size = 4))) +
  xlab("UMAP 1") + ylab("UMAP 2")

# 22. Supp Fig 1g-h - correlation boxplots

cor_all <- rbind(
  transform(xx_reSCTClass, method = "Absolute pathway activity"),
  transform(xx_reSCTClassNormalized, method = "Pathway allocation")
) %>% mutate(method = factor(method, levels = names(method_colors)))

cell_order_cor <- cor_all %>%
  filter(method == "Absolute pathway activity") %>%
  group_by(cell) %>% summarise(median_rho = median(rho.rho), .groups = "drop") %>%
  arrange(desc(median_rho)) %>% pull(cell)
cor_all$cell <- factor(cor_all$cell, levels = cell_order_cor)

p_cor_all <- ggplot(cor_all, aes(x = cell, y = rho.rho, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.8), size = 1.5, alpha = 0.5, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = method_colors) +
  labs(x = NULL, y = "Mean correlation coefficient\nof log2FC(T vs N)", fill = "Method") +
  ggtitle("Global pathway correlation\nwith total metabolic load") +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  guides(fill = guide_legend(nrow = 2))
p_cor_all

# 23. Assemble Supplementary Figure 1

p_s13 <- ggarrange(p_cor_all, p_supp, nrow = 2, ncol = 1,
                   labels = letters[7:8],
                   font.label = list(size = 18, color = "black", face = "bold"),
                   common.legend = TRUE, legend = "bottom")
ps1_1 <- ggarrange(p_vln1, p_vln2, p_vln3, nrow = 3, ncol = 1,
                   labels = letters[1:3],
                   font.label = list(size = 18, color = "black", face = "bold"))
ps1_2 <- ggarrange(p2_dot, cor_result$cor_meta$pp, p_s13, nrow = 1, ncol = 3,
                   labels = c(letters[5:6], ""),
                   font.label = list(size = 18, color = "black", face = "bold"),
                   widths = c(1.4, 1.8, 1))
ps1_4 <- ggarrange(ps1_1, p_patient_umap, ps1_2, nrow = 3, ncol = 1,
                   labels = c("", letters[4], ""),
                   font.label = list(size = 18, color = "black", face = "bold"),
                   heights = c(2.7, 1.4, 2))

save_tiff(ps1_4, "supplementaryfigure1.tiff", width = 19, height = 22)

# 24. Supp Fig 2a - proportion heatmap with "Sum proportion" row

annotation_row <- celltype_mean %>%
  mutate(pathway = "Sum proportion", proportion = mean_prop) %>%
  dplyr::select(celltype, pathway, proportion)

combined_data <- bind_rows(proportion_matrix %>% dplyr::select(celltype, pathway, proportion),
                           annotation_row) %>%
  mutate(pathway = factor(pathway, levels = c(pathway_order, "Sum proportion")),
         celltype = factor(celltype, levels = celltype_order))

# Manual color scales (identity fill): white -> red for pathways, red -> purple for the sum row
pathway_color <- function(x) {
  intensity <- pmin(x / 0.8, 1)
  rgb(1 - intensity * (1 - 196 / 255),
      1 - intensity * (1 - 30 / 255),
      1 - intensity * (1 - 58 / 255))
}
sum_color <- function(x) {
  if (x <= 0.5) {  # same scaling as pathway_color for continuity
    return(pathway_color(x))
  }
  if (x <= 1.0) {  # 0.5-1.0: deepen red
    t <- (x - 0.4) / 0.6
    return(rgb(196 / 255 - t * 30 / 255, 30 / 255 - t * 10 / 255, 58 / 255 - t * 10 / 255))
  }
  if (x <= 1.5) {  # 1.0-1.5: red -> red-purple
    t <- (x - 1.0) / 0.5
    return(rgb(166 / 255 - t * 30 / 255, 20 / 255 + t * 30 / 255, 48 / 255 + t * 60 / 255))
  }
  if (x <= 2.0) {  # 1.5-2.0: red-purple -> light purple
    t <- (x - 1.5) / 0.5
    return(rgb(136 / 255 - t * 30 / 255, 50 / 255 + t * 50 / 255, 108 / 255 + t * 60 / 255))
  }
  t <- pmin((x - 2.0) / 0.6, 1)  # 2.0-2.6: light purple, capped
  rgb(106 / 255 - t * 15 / 255, 100 / 255 + t * 20 / 255, 168 / 255 + t * 20 / 255)
}

combined_data <- combined_data %>%
  mutate(tile_color = ifelse(pathway == "Sum proportion",
                             sapply(proportion, sum_color),
                             sapply(proportion, pathway_color)))

y_face <- c(rep("plain", n_pathways), "bold")
y_color <- c(rep("black", n_pathways), "#8B0000")
celltype_levels <- levels(combined_data$celltype)
x_colors <- ifelse(celltype_levels == "Epithelial", epi_color, "black")
x_faces <- ifelse(celltype_levels == "Epithelial", "bold", "plain")

p_combined <- ggplot(combined_data, aes(x = celltype, y = pathway)) +
  geom_tile(aes(fill = tile_color), color = "white", linewidth = 0.5) +
  scale_fill_identity() +
  geom_text(aes(label = round(proportion, 2),
                fontface = ifelse(pathway == "Sum proportion", "bold", "plain")),
            size = 3.2, color = "black") +
  geom_hline(yintercept = n_pathways + 0.5, color = "black", linewidth = 1.2) +
  labs(x = "Cell Type", title = "Sample proportion with pathway rank change >2") +
  theme_hd_minimal_plain() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = x_colors, face = x_faces),
        axis.text.y = element_text(size = 11, face = y_face, color = y_color,
                                   margin = margin(r = 5), hjust = 1),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 20))
p_combined


# 26. Cosine distance of pathway allocation (T vs N)

cosine_dist <- all_long %>%
  group_by(dataset, celltype, patient, condition) %>%
  mutate(Gp_norm = Gp / sum(Gp)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(dataset, celltype, patient, pathway),
              names_from = condition, values_from = Gp_norm) %>%
  filter(!is.na(Normal), !is.na(Tumor)) %>%
  group_by(dataset, celltype, patient) %>%
  summarise(cosine_sim = sum(Normal * Tumor) / (sqrt(sum(Normal^2)) * sqrt(sum(Tumor^2))),
            .groups = "drop") %>%
  mutate(cosine_dist = 1 - cosine_sim)

celltype_order_median <- cosine_dist %>%
  group_by(celltype) %>%
  summarise(median_dist = median(cosine_dist), .groups = "drop") %>%
  arrange(desc(median_dist)) %>% pull(celltype)
cosine_dist$celltype <- factor(cosine_dist$celltype, levels = celltype_order_median)

p_cosine_final <- ggplot(cosine_dist, aes(x = celltype, y = cosine_dist,
                                          fill = celltype == "Epithelial")) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.4) +
  geom_jitter(aes(color = dataset), width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("TRUE" = epi_color, "FALSE" = "#BDBDBD"), guide = "none") +
  scale_color_manual(values = dataset_colors, guide = "none") +
  scale_x_discrete(labels = epi_x_labels) +
  scale_y_log10() +
  labs(x = NULL, y = "Cosine distance (T vs N)",
       title = "Divergence of global pathway allocation") +
  theme_hd_plain() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        legend.position = "right")
p_cosine_final

# 27. Key pathway boxplots (Epithelial)

p_box <- plot_pathway_boxplot(
  data = as.data.frame(pathway_mean_patientClassNormalized),
  key_pathways = c("Energy metabolism", "Detoxification pathways", "Carbohydrate metabolism",
                   "Nucleotide metabolism", "Transport pathways", "Xenobiotics metabolism"),
  celltype = "Epithelial",
  facet_ncol = 6,
  y_label = "Global pathway allocation",
  title = "Epithelial cell",
  log_scale = TRUE
)
p_box

# 28. Proportion of significantly changed pathways per method

count_sig <- function(sig_df) {
  sig_df %>%
    mutate(sig = wcx_p < 0.05) %>%
    group_by(cell) %>%
    summarise(n_sig = sum(sig), prop_sig = mean(sig), .groups = "drop")
}
gp_sig <- count_sig(fc_mainclassNormalized$all_path_sigClass)
ap_sig <- count_sig(fc_mainclass$all_path_sigClass)

compare_sig <- inner_join(
  gp_sig %>% dplyr::select(cell, gp_n_sig = n_sig, gp_prop = prop_sig),
  ap_sig %>% dplyr::select(cell, ap_n_sig = n_sig, ap_prop = prop_sig),
  by = "cell")
compare_sig

custom_order <- c("Epithelial", "Endothelial", "Fibroblast", "Myeloid", "Plasma", "T cell", "B cell")
compare_prop_long <- compare_sig %>%
  pivot_longer(c(gp_prop, ap_prop), names_to = "method", values_to = "prop_sig") %>%
  mutate(method = ifelse(method == "gp_prop", "Pathway allocation", "Absolute pathway activity"),
         cell = factor(cell, levels = custom_order))


# cell type order; compute both style vectors from the order used here
x_colors <- ifelse(levels(compare_prop_long$cell) == "Epithelial", epi_color, "black")
x_faces <- ifelse(levels(compare_prop_long$cell) == "Epithelial", "bold", "plain")

p_prop <- ggplot(compare_prop_long, aes(x = cell, y = prop_sig, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Absolute pathway activity" = "#00BFFF",
                               "Pathway allocation" = "#EE7621")) +
  scale_y_continuous(limits = c(0, 1.15), labels = percent, expand = c(0, 0)) +
  labs(y = "Proportion of significantly changed\nglobal pathways (T vs N)",
       title = "Allocation reveals reprogramming\nmasked by total activity change") +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = x_colors, face = x_faces),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2))
p_prop

# 29. Supp Fig 2c - non-epithelial cells of patient C157


c157_all <- data %>%
  filter(patient == "C157") %>%
  pivot_longer(all_of(pathway_cols), names_to = "pathway", values_to = "Gp") %>%
  mutate(condition = factor(condition, levels = c("Normal", "Tumor")))

c157_all_highlight <- c157_all %>%
  group_by(celltype, condition) %>%
  mutate(rank = rank(-Gp)) %>%
  ungroup() %>%
  dplyr::select(celltype, pathway, condition, rank) %>%
  pivot_wider(names_from = condition, values_from = rank) %>%
  mutate(highlight = abs(Tumor - Normal) > 2) %>%
  dplyr::select(celltype, pathway, highlight)

c157_all <- c157_all %>% left_join(c157_all_highlight, by = c("celltype", "pathway"))

c157_other <- c157_all %>%
  filter(celltype != "Epithelial") %>%
  mutate(celltype = factor(celltype, levels = setdiff(celltype_levels, "Epithelial")))

p_s2_labeled <- plot_pathway_shift(
  c157_other, facet_by = "celltype", free_y = TRUE,
  title = "Non-epithelial cells show minimal priority disruption (Patient C157)")
p_s2_labeled

# 30. Supp Fig 2b - priority shift in epithelial cells of other patients

extra_patients <- c("C165", "C114", "SC044", "C110")

s3_data <- data %>%
  filter(celltype == "Epithelial", patient %in% extra_patients) %>%
  pivot_longer(all_of(pathway_cols), names_to = "pathway", values_to = "Gp") %>%
  mutate(condition = factor(condition, levels = c("Normal", "Tumor")))

s3_highlight <- s3_data %>%
  group_by(patient, condition) %>%
  mutate(rank = rank(-Gp)) %>%
  ungroup() %>%
  dplyr::select(patient, pathway, condition, rank) %>%
  pivot_wider(names_from = condition, values_from = rank) %>%
  mutate(highlight = abs(Tumor - Normal) > 2) %>%
  dplyr::select(patient, pathway, highlight)

s3_data <- s3_data %>% left_join(s3_highlight, by = c("patient", "pathway"))

p_s3_labeled <- plot_pathway_shift(
  s3_data, facet_by = "patient",
  title = "Example: Pathway priority shift\nin epithelial cells of different patients")
p_s3_labeled

# 31. Assemble Supplementary Figure 2


spall2 <- ggarrange(p_prop, p_cosine_final,p_box,
                      labels = letters[1:3], nrow = 1, ncol = 3,
                      font.label = list(size = 18, color = "black", face = "bold"),widths = c(1,1,2))
spall2_0 <- ggarrange(p_combined, p_s3_labeled,
                      labels = letters[4:5], nrow = 1, ncol = 2,
                      font.label = list(size = 18, color = "black", face = "bold"))
spall2_1 <- ggarrange(spall2 ,spall2_0, p_s2_labeled,
                      labels = c("","", letters[6]), nrow = 3, ncol = 1,
                      font.label = list(size = 18, color = "black", face = "bold"))

save_tiff(spall2_1, "supplementaryfigure2.tiff", width = 18, height = 20)

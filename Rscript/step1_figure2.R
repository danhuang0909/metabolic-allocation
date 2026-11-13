# Figure 2 and Supplementary Figures 3-8.
# Run this script from the project root with Rscript.

project_dir <- normalizePath(
  Sys.getenv("METABOLISM_PROJECT_DIR", unset = getwd()),
  mustWork = TRUE
)
source(file.path(project_dir, "Rscript", "step0_config.R"))
setwd(project_dir)


# ==============================================================================
# Pathway reallocation analysis: epithelial focus + per-cell-type networks
# ==============================================================================

# ---- 1. Setup ----------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggrepel)

panel_label_style <- list(size = 18, color = "black", face = "bold")

nt_colors <- c("Purine synthesis" = "#D6604D",
               "Pyrimidine synthesis" = "#B2182B",
               "Other nucleotide metabolism" = "#999999",
               "Purine catabolism" = "#2166AC",
               "Pyrimidine catabolism" = "#4393C3")

line_pal <- c("Downregulated module" = "#2980b9",
              "Energy & lipid module" = "#F72585",
              "Glycan & transport detox module" = "#e67e22")

# Wrapper around draw_pathway_network_gg for one cell type of a stats result
draw_celltype_network <- function(celltype, stats_result, r_threshold = 0.5,
                                  seed = 115, filter_kegg_class = NULL,
                                  show_kegg_legend = TRUE, show_label = TRUE,
                                  label_size = 3.2, p_cutoff = 0.05,
                                  remove_isolated = TRUE) {
  draw_pathway_network_gg(
    expr_matrix = stats_result$cell_fc_list[[celltype]],
    anno_df = all_gene_pathway,
    fc_df = stats_result$sig_data_all %>% filter(cell == celltype),
    cell_type = celltype,
    r_threshold = r_threshold,
    seed = seed,
    filter_kegg_class = filter_kegg_class,
    show_label = show_label,
    label_size = label_size,
    show_kegg_legend = show_kegg_legend,
    p_cutoff = p_cutoff,
    remove_isolated = remove_isolated
  )
}

# ---- 2. Data loading and FC computation --------------------------------------

allcells_pathwayMeanSCTMainClassNormalized <- readRDS("./data/allcells_pathwayMeanSCTMainClassNormalizedByMetacount.Rds")
allcells_pathwayMeanSCTNormalized <- readRDS("./data/allcells_pathwayMeanSCTNormalizedByMetacount.Rds")
allcells_pathwayMeanSCT <- readRDS("./data/allcells_pathwayMeanSCT.Rds")
all_gene_pathway <- readRDS("./data/all_genes_pathway_infor_select.rds")

pathwayRatioNormalized <- readRDS(file.path(data_dir, "pathwayRatioNormalized.rds"))

meta_data <- readRDS("./data/all_dataset_SCT_filter_meata_data.rds")
meta_data <- meta_data[meta_data$Cell_type != "Mast", ]

fc_SCT <- getFCresult(allcells_pathwayMeanSCT, meta_data)
fc_SCTNormalized <- getFCresult(allcells_pathwayMeanSCTNormalized, meta_data)
fc_pathwayRatioNormalized <- getFCresult(pathwayRatioNormalized, meta_data)


pathway_localclassNormalized <- getAggregateDate2(allcells_pathwayMeanSCTNormalized, meta_data)
pathway_localclassNormalized_normal <- t(pathway_localclassNormalized[grepl("Normal", rownames(pathway_localclassNormalized)), ])
colnames(pathway_localclassNormalized_normal) <- gsub("_Normal", "", colnames(pathway_localclassNormalized_normal))

pathway_mean_patientRatioNormalized <- getAggregateDate2(pathwayRatioNormalized, meta_data)
pathway_mean_patientRatioNormalized_normal <- t(pathway_mean_patientRatioNormalized[grepl("Normal", rownames(pathway_mean_patientRatioNormalized)), ])
colnames(pathway_mean_patientRatioNormalized_normal) <- gsub("_Normal", "", colnames(pathway_mean_patientRatioNormalized_normal))

sp_result <- prepare_pathway_stats(fc_data = fc_SCTNormalized, classification = all_gene_pathway)
ap_result <- prepare_pathway_stats(fc_data = fc_SCT, classification = all_gene_pathway)


epi_fc_sp <- fc_SCTNormalized$FC[grep("_Epithelial_", rownames(fc_SCTNormalized$FC)), , drop = FALSE]
sig_data_sp <- sp_result$sig_data_all
epi_data <- sig_data_sp %>% filter(cell == "Epithelial")

# ---- 3. Fig 2a: KEGG class transcriptional changes (Epithelial) --------------

res_epi <- plot_kegg_class_boxplot_unified(
  sig_data = sp_result$sig_data_all,
  cell_type = "Epithelial",
  single_cell_mode = TRUE,
  top_n_label = 10
)
p_change <- res_epi$plot

# ---- 4. Fig 2c: intra-class correlation (Epithelial) --------------------------

p_epi_cor <- plot_intra_class_correlation(
  stats_result = sp_result,
  cell_type = "Epithelial",
  facet_all_cells = FALSE
)
p_cor <- p_epi_cor$boxplot


# ---- 6. Fig 2b: carbohydrate pathway boxplots (Epithelial) --------------------

p_path <- plot_pathway_boxplot(
  data = pathway_localclassNormalized,
  key_pathways = c("Glycolysis/gluconeogenesis", "Pyruvate metabolism",
                   "Butanoate metabolism", "Citric acid cycle"),
  celltype = "Epithelial",
  facet_ncol = 4,
  y_label = "Local per-gene pathway allocation",
  title = "Carbohydrate metabolism in epithelial cells",
  log_scale = TRUE,
  name_nrow = FALSE
)
print(p_path)

# ---- 7. Fig 2d: carbohydrate subnetwork (Epithelial) --------------------------

p_net <- draw_celltype_network("Epithelial", sp_result,
                               r_threshold = 0.4, seed = 42,
                               filter_kegg_class = "Carbohydrate metabolism",
                               show_kegg_legend = FALSE)
print(p_net)

# ---- 8. Fig 2e/g/h: pathway FC scatter plots ----------------------------------

p_scatter_carbon <- plot_pathway_scatter(
  epi_fc = epi_fc_sp,
  path_x = "Citric acid cycle",
  path_y = "Pyruvate metabolism",
  title = "Within-class Correlation in epithelial cells",
  label_x = expression("Citric Acid Cycle" ~ Log[2] * "FC"),
  label_y = expression("Pyruvate metabolism" ~ Log[2] * "FC"),
  label_pos = c(-1.2, 1)
)

p1 <- plot_pathway_scatter(
  epi_fc = epi_fc_sp,
  path_x = "Citric acid cycle",
  path_y = "Transport, lysosomal",
  title = "Inter-class Correlation in epithelial cells",
  label_x = expression("Citric Acid Cycle" ~ Log[2] * "FC"),
  label_y = expression("Lysosomal Transport" ~ Log[2] * "FC"),
  label_pos = c(-1.5, 2)
)

p2 <- plot_pathway_scatter(
  epi_fc = epi_fc_sp,
  path_x = "Glycolysis/gluconeogenesis",
  path_y = "Pyrimidine synthesis",
  title = "Inter-class Correlation in epithelial cells",
  label_x = expression("Glycolysis/Gluconeogenesis" ~ Log[2] * "FC"),
  label_y = expression("Pyrimidine Synthesis" ~ Log[2] * "FC"),
  label_pos = c(-0.5, 2.5)
)

# ---- 9. Fig 2f: full epithelial network with module outlines ------------------

p_all <- draw_celltype_network("Epithelial", sp_result)
print(p_all$p)

# The module outline coordinates in ./result/network_position.rds were created
# interactively with locator() on this plot (downregulated / energy-lipid /


module_outline <- readRDS("./result/network_position.rds")

# Overlay module outlines
p_final <- p_all$p +
  geom_polygon(data = module_outline,
               aes(x = x, y = y, group = group, color = group),
               fill = NA, linetype = "dashed", linewidth = 0.82, alpha = 0.42,
               inherit.aes = FALSE) +
  scale_color_manual(values = line_pal, guide = "none")

# Add module annotations with arrows
annotation_df <- data.frame(
  label = c("Tissue-specialized metabolism (Down)",
            "Central carbon allocation (Up)",
            "Glycan remodeling (Mixed)"),
  group = c("Downregulated module", "Energy & lipid module", "Glycan & transport detox module"),
  x = c(2.8, 3, -10),
  y = c(5, -7, 5.5),
  arrow_x = c(2.7, 2, -9.3),
  arrow_y = c(4.1, -6, 4.7)
)

p_final <- p_final +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label, color = group),
            size = 3.8, fontface = "bold", inherit.aes = FALSE) +
  geom_segment(data = annotation_df,
               aes(x = x, y = y, xend = arrow_x, yend = arrow_y, color = group),
               arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.6,
               inherit.aes = FALSE) +
  scale_color_manual(values = line_pal, guide = "none")

print(p_final)

# ---- 10. Assemble Figure 2 -----------------------------------------------------

pall3_0 <- ggarrange(p_change, p_path,
                     labels = letters[1:2], nrow = 1, ncol = 2,
                     font.label = panel_label_style)
pall3_1 <- ggarrange(p_cor + theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 10)), p_net$p,
                     labels = letters[3:4], nrow = 1, ncol = 2,
                     font.label = panel_label_style, widths = c(1.5, 1.5))
pall3_2 <- ggarrange(p_scatter_carbon,  p2,p1,
                     labels = letters[c(5, 7:8)], nrow = 3, ncol = 1,
                     font.label = panel_label_style,
                     common.legend = TRUE, legend = "bottom")
pall3_3 <- ggarrange(pall3_1, p_final,
                     labels = c("", letters[6]), nrow = 2, ncol = 1,
                     font.label = panel_label_style, heights = c(1, 2))
pall3_4 <- ggarrange(pall3_3, pall3_2, nrow = 1, ncol = 2, widths = c(3, 1))
pall3_5 <- ggarrange(pall3_0, pall3_4, nrow = 2, ncol = 1, heights = c(1.2, 3))

save_tiff(pall3_5, "figure2.tiff", width = 18, height = 23)

# ==============================================================================
# Supplementary figures
# ==============================================================================

# ---- 11. Supp Fig 2: non-epithelial KEGG class changes + correlations ---------

cell_types <- unique(sig_data_sp$cell)
non_epi_types <- cell_types[cell_types != "Epithelial"]

res_non_epi <- plot_kegg_class_boxplot_unified(
  sig_data = sig_data_sp,
  single_cell_mode = FALSE,
  cell_types_facet = non_epi_types,
  top_n_label = 8,
  title = "Non-epithelial cell populations: \nPathway transcriptional changes across classes"
)
p_change_non_epi <- res_non_epi$plot + coord_flip()

p_multi_cor_non_epi <- plot_intra_class_correlation(
  stats_result = sp_result,
  facet_all_cells = TRUE,
  cell_list_facet = non_epi_types
)
p_cor_non_epi <- p_multi_cor_non_epi$boxplot

ps3_1 <- ggarrange(p_change_non_epi, p_cor_non_epi,
                   labels = letters[1:2], nrow = 1, ncol = 2,
                   font.label = panel_label_style)
save_tiff(ps3_1, "supplementaryfigure3.tiff", width = 17, height = 22)

# ---- 12. Supp Fig 5: intra-class correlation, absolute activity (all cells) ---

p_multi_cor_ap <- plot_intra_class_correlation(
  stats_result = ap_result,
  facet_all_cells = TRUE,
  cell_list_facet = cell_types
)
p_cor_ap <- p_multi_cor_ap$boxplot + facet_wrap(~cell, nrow = 2)
save_tiff(p_cor_ap, "supplementaryfigure5.tiff", width = 17, height = 17)

# ---- 13. Supp Figs 4 & 6: pathway correlation heatmaps (Epithelial) -----------

kegg_order <- epi_data %>%
  group_by(kegg_class) %>%
  summarise(median_fc = median(abs(mean_log2fc)), .groups = "drop") %>%
  arrange(desc(median_fc)) %>%
  pull(kegg_class)

local_to_kegg <- all_gene_pathway %>%
  select(final_name, kegg_class) %>%
  distinct() %>%
  deframe()

kegg_colors <- c("Energy metabolism" = "#c42e2b",
                 "Carbohydrate metabolism" = "#e67e22",
                 "Nucleotide metabolism" = "#2980b9",
                 "Protein metabolism" = "#27ae60",
                 "Canonical amino acid metabolism" = "#8e44ad",
                 "Lipid metabolism" = "#d81b87",
                 "Detoxification pathways" = "#EEA2AD",
                 "Xenobiotics metabolism" = "#5c3a1e",
                 "Metabolism of cofactors and vitamins" = "#BFEFFF",
                 "Transport pathways" = "#EEEE00",
                 "Glycan metabolism" = "#d0d0d0",
                 "Non-canonical amino acids Metabolism" = "#40E0D0",
                 "Miscellaneous" = "#707070",
                 "Other" = "#303030")

pp_sp <- plot_pathway_correlation_heatmap(
  fc_matrix = sp_result$cell_fc_list$Epithelial,
  sig_data = sp_result$sig_data_all,
  cell_type = "Epithelial",
  title_x = "allocation value",
  p_col = "wcx_p_pair",
  p_threshold = 0.05,
  classification = local_to_kegg,
  class_colors = kegg_colors,
  class_order = kegg_order,
  cluster_within = TRUE
)

pp_ap <- plot_pathway_correlation_heatmap(
  fc_matrix = ap_result$cell_fc_list$Epithelial,
  sig_data = ap_result$sig_data_all,
  cell_type = "Epithelial",
  title_x = "absolute activity",
  p_col = "wcx_p_pair",
  p_threshold = 0.05,
  classification = local_to_kegg,
  class_colors = kegg_colors,
  class_order = kegg_order,
  cluster_within = TRUE
)

save_tiff(pp_sp$ht, "supplementaryfigure4.tiff", width = 17, height = 14)
save_tiff(pp_ap$ht, "supplementaryfigure6.tiff", width = 17, height = 14)

# ---- 14. Supp Figs 7 & 8: per-cell-type networks -------------------------------

# SP networks for all non-epithelial cell types (identical parameters)
sp_net <- lapply(
  setNames(nm = c("Endothelial", "Fibroblast", "Myeloid", "Plasma", "B cell", "T cell")),
  draw_celltype_network,
  stats_result = sp_result
)

# AP networks for Epithelial at two correlation thresholds.


ap_net_epi <- lapply(
  setNames(nm = c(0.5, 0.6)),
  function(r) draw_celltype_network("Epithelial", ap_result, r_threshold = as.numeric(r))
)
ap_net_epi <- lapply(ap_net_epi, function(res) {
  res$p <- res$p + ggtitle("Epithelial inter-class pathway co-reallocation network\n(Absolute pathway activity)")
  res
})

ps7_1 <- ggarrange(ap_net_epi$`0.5`$p, ap_net_epi$`0.6`$p,
                   sp_net$Endothelial$p, sp_net$Fibroblast$p,
                   labels = letters[1:4], nrow = 2, ncol = 2,
                   font.label = panel_label_style)
save_tiff(ps7_1, "supplementaryfigure7.tiff", width = 20, height = 20)

ps8_1 <- ggarrange(sp_net$Myeloid$p, sp_net$Plasma$p, sp_net$`B cell`$p, sp_net$`T cell`$p,
                   labels = letters[1:4], nrow = 2, ncol = 2,
                   font.label = panel_label_style)
save_tiff(ps8_1, "supplementaryfigure8.tiff", width = 20, height = 20)

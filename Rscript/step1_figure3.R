# Figure 3 and Supplementary Figure 9.
# Run this script from the project root with Rscript.

project_dir <- normalizePath(
  Sys.getenv("METABOLISM_PROJECT_DIR", unset = getwd()),
  mustWork = TRUE
)
source(file.path(project_dir, "Rscript", "step0_config.R"))
setwd(data_dir)


# 1. Libraries and environment

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(pheatmap)
library(igraph)
library(FactoMineR)
library(grid)
library(tiff)

# 2. Shared palettes and helpers

dataset_colors <- c("CRC-SG1" = "#7FDBB3", "KUL3" = "#F5AD75", "KUL5" = "#85A4E6",
                    "Pelka" = "#F593CE", "SMC" = "#CAF77B")
masm_colors_2 <- c("1" = "#73b873", "2" = "#f8b360", "3" = "#4a7acc", "4" = "#b98fd8",
                 "5" = "#e86868", "6" = "#d8b99c", "7" = "#69b8b8")
module_colors <- c("Glycan remodeling(Mixed)" = "#e6ab02",
                   "Tissue-specialized catabolism(Down)" = "#2166ac",
                   "Central carbon allocation(Up)" = "#b2182b",
                   "Other & Heme & transport" = "gray")
cor_palette <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
panel_label_style <- list(size = 18, color = "black", face = "bold")


run_epi_pca <- function(mat, title, n_x=4,standardize = FALSE, ...) {
if (standardize) mat <- standardize_rowname(mat)
plot_epi_pca(epi_all_mat = mat, patient_infor = patient_infor,
ncp = n_x, top_var_n = 12, plot_title = title, ...)
}

standardize_rowname <- function(mat) {
    rn_tbl <- tibble(old = rownames(mat)) %>%
        mutate(group = str_extract(old, "Normal|Tumor"),
               id_clean = str_remove(old, "_Normal|_Tumor"),
               new = paste0(id_clean, "_", group))
    `rownames<-`(mat, rn_tbl$new)
}


# 3. Data loading

meta_data <- readRDS("all_meta_data.rds")
patient_infor <- readRDS("single_cell_patient_infor_all_samples.rds")
# Per-patient pseudobulk CPM of all genes, and KEGG pathway gene sets
cpm_all <- readRDS("sum_all_genes_count_rawRNA_cpm.rds")
kegg_path_genes <- readRDS("all_kegg_pathway_genes.rds")

all_gene_pathway <- readRDS("all_genes_pathway_infor_select.rds")
allcells_pathwayMeanSCT <- readRDS("allcells_pathwayMeanSCTNormalizedByMetacount.Rds")
result_df <- readRDS("sum_all_genes_count_rawRNA_test.rds")
# Precomputed: per-gene Spearman correlation with pathway reallocation
# (originally: getcorSpearmanForall(cpm_all_select_signal_fcAN, fc_pathwayMeanSCT$FC))
temp_corSCTnormalized <- readRDS("./signalgenes_correlation_pathwaySCTnormalized.Rds")

# Signaling-pathway genes (KEGG "signal" entries), excluding metabolic genes
kegg_signal <- kegg_path_genes[grepl("signal", kegg_path_genes$description), ]
cpm_all_select_signal <- cpm_all[, colnames(cpm_all) %in% kegg_signal$symbol &
                                   !colnames(cpm_all) %in% unique(all_gene_pathway$gene)]
cpm_all_select_signal_fcAN<-getFCdataSinglelogAllnormaleach(cpm_all_select_signal,str_split_fixed(rownames(cpm_all_select_signal),"_",4)[,3])

fc_pathwayMeanSCT <- getFCresult(allcells_pathwayMeanSCT, meta_data)
sp_result <- prepare_pathway_stats(fc_data = fc_pathwayMeanSCT,
                                   classification = all_gene_pathway)


result_df <- result_df[result_df$cell_type != "Mast", ]
result_df_sig <- result_df[result_df$PValue < 0.05, ]
result_df_sig_epi <- result_df_sig[result_df_sig$cell_type == "Epithelial", ]


epi_idx <- str_detect(rownames(cpm_all_select_signal), "Epithelial")
signal_all_pca_result <- run_epi_pca(
  cpm_all_select_signal[epi_idx, ], n_x = 10,"All signaling genes", standardize = TRUE)

library(cluster)
library(dplyr)


pca_mat <- signal_all_pca_result$sample_coord[, paste0("Dim.",1:10)]
meta <- signal_all_pca_result$sample_coord[, c("group", "dataset")]


df <- cbind(meta, pca_mat)


batch_sil_res <- df %>%
    group_by(group) %>%
    group_modify(~{
        sub_df <- .x

        tab <- table(sub_df$dataset)
        if(nrow(sub_df) < 3 || length(tab) <2 || any(tab <2)){
            return(tibble(ASW_batch = NA, n_sample = nrow(sub_df), n_batch = length(tab)))
        }


        dist_pca <- dist(sub_df[, paste0("Dim.",1:10)], method = "euclidean")
        # silhouette：cluster label = dataset
        sil_obj <- silhouette(as.integer(factor(sub_df$dataset)), dist_pca)
        asw <- mean(sil_obj[,3])
        tibble(ASW_batch = asw, n_sample = nrow(sub_df), n_batch = length(tab))
    }) %>% ungroup()

batch_sil_res


# 4. Fig 3a - flowchart

library("tiff")
img <- rasterGrob(readTIFF(file.path(fig_dir, "figure4a.tif"),
                           native = TRUE), interpolate = TRUE)
p4_flowchart <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
  geom_blank() +
  annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

# 5. Epithelial pathway network and Louvain modules

p_all <- draw_pathway_network_gg(
  expr_matrix = sp_result$cell_fc_list[["Epithelial"]],
  anno_df = all_gene_pathway,
  fc_df = sp_result$sig_data_all %>% filter(cell == "Epithelial"),
  cell_type = "Epithelial",
  r_threshold = 0.5,
  seed = 115,
  filter_kegg_class = NULL,
  show_label = TRUE,
  label_size = 3.2,
  show_kegg_legend = TRUE,
  p_cutoff = 0.05,
  remove_isolated = TRUE
)
print(p_all$p)

cor_path <- p_all$cor_long
g <- p_all$g

# Louvain community detection on the pathway network
com <- igraph::cluster_louvain(g)
V(g)$louvain_module <- igraph::membership(com)
cat("Network modularity (Louvain):", igraph::modularity(com), "\n")

module_path_list <- split(V(g)$name, V(g)$louvain_module)
module_anno_df <- tibble(
  pathway = unlist(module_path_list),
  module_id = rep(names(module_path_list), lengths(module_path_list))
) %>%
  mutate(module_label = case_when(
    module_id == "1" ~ "Glycan remodeling(Mixed)",
    module_id == "2" ~ "Tissue-specialized catabolism(Down)",
    module_id == "3" ~ "Central carbon allocation(Up)",
    module_id %in% c("4", "5") ~ "Other & Heme & transport"
  ))

row_anno <- module_anno_df %>%
  select(pathways = pathway, module_label) %>%
  column_to_rownames("pathways")

# 6. Fig 3c - gene x pathway correlation heatmap


epi_gene_df <- as.data.frame(
  cpm_all_select_signal_fcAN[grepl("Epi", rownames(cpm_all_select_signal_fcAN)),]
)

# Keep only epithelial correlations for pathways present in the network
select_data <- temp_corSCTnormalized %>%
  filter(cell == "Epithelial",
         pathways %in% unique(c(cor_path$from, cor_path$to)))

# Gene filter: significant strong correlation (|rho| > 0.5, P < 0.05) with >= 2 pathways

keep_gene_list_candidate <- select_data %>%
  mutate(hit = as.integer(abs(rho.rho) > 0.5 & pvalue < 0.05)) %>%
  group_by(gene) %>%
  summarise(total_hit = sum(hit, na.rm = TRUE), .groups = "drop") %>%
  filter(total_hit > 1) %>%
  pull(gene)

epi_fc <- sp_result$cell_fc_list$Epithelial

# run
meta_result <- gene_pathway_meta_cor(
  epi_gene_df = as.data.frame(epi_gene_df[,keep_gene_list_candidate]),
  epi_fc      = as.data.frame(epi_fc),
  COR_METHOD  = "spearman",
  MIN_N       = 4,
  MIN_K       = 3,
  merge_groups = list(c("KUL3","KUL5"))
)


cmp <- merge(
  meta_result[, c("gene", "pathway", "rho_pool", "p", "fdr")],
  select_data[, c("gene", "pathways", "rho.rho", "pvalue")],
  by.x = c("gene", "pathway"), by.y = c("gene", "pathways")
)
names(cmp)[names(cmp) == "rho.rho"] <- "rho_merged"


per_gene_sigA <- do.call(rbind, lapply(split(cmp, cmp$gene), function(x) {
  xs <- x[x$pvalue < 0.05, ]
  if (nrow(xs) < 3) return(NULL)
  data.frame(
    gene       = x$gene[1],
    n_sig_pw   = nrow(xs),
    cor_rho    = cor(xs$rho_pool, xs$rho_merged,method = "spearman"),
    sign_agree = mean(sign(xs$rho_pool) == sign(xs$rho_merged))
  )
}))
per_gene_sigA <- per_gene_sigA[order(-per_gene_sigA$cor_rho), ]

cat(sprintf("A: %d genes with >=3 merged-significant pathways\n", nrow(per_gene_sigA)))
cat(sprintf("median cor: %.3f | cor>0.7: %d | cor<0: %d\n",
            median(per_gene_sigA$cor_rho),
            sum(per_gene_sigA$cor_rho > 0.7),
            sum(per_gene_sigA$cor_rho < 0)))
head(per_gene_sigA, 20)
keep_gene_list<-keep_gene_list_candidate[keep_gene_list_candidate %in% per_gene_sigA$gene[per_gene_sigA$cor_rho>0.7 & per_gene_sigA$sign_agree>0.7]]


mat_wide_filtered <- select_data %>%
  filter(gene %in% keep_gene_list) %>%
  pivot_wider(id_cols = pathways, names_from = gene, values_from = rho.rho) %>%
  column_to_rownames("pathways")

cat("Total genes:", length(unique(select_data$gene)), "\n",
    "Genes passing filter (|rho| > 0.5 & P < 0.05 in >= 2 pathways):", length(keep_gene_list), "\n",
    "Remaining pathways:", nrow(mat_wide_filtered), "\n")

# Heatmap (row annotation only) - also used to derive the gene dendrogram
ht_genes <- pheatmap::pheatmap(
  mat_wide_filtered,
  scale = "none",
  color = cor_palette,
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 12,
  treeheight_col = 12,
  annotation_colors = list(module_label = module_colors),
  main = "Genes with strong significant correlation (|rho| > 0.5, P < 0.05) with >= 2 metabolic pathways"
)

# Cut gene dendrogram into 7 subclusters (MASM1-7)
gene_cluster <- cutree(ht_genes$tree_col, k = 7)
saveRDS(gene_cluster, file = "gene_cluster.rds")

col_anno <- data.frame(Gene_subcluster = paste0("MASM", gene_cluster),
                       row.names = names(gene_cluster))[colnames(mat_wide_filtered), , drop = FALSE]


#

# sig_pathways <- unique(select_data$pathways[select_data$pvalue < 0.05])
# cat(sprintf("merged-significant pathways (any gene): %d\n", length(sig_pathways)))
#
# cmp_sig <- cmp[cmp$pathway %in% sig_pathways, ]
#
# per_gene_sigB <- do.call(rbind, lapply(split(cmp_sig, cmp_sig$gene), function(x) {
# if (nrow(x) < 5) return(NULL)
# data.frame(
# gene       = x$gene[1],
# n_pathways = nrow(x),
# cor_rho    = cor(x$rho_pool, x$rho_merged,method = "spearman"),
# sign_agree = mean(sign(x$rho_pool) == sign(x$rho_merged))
# )
# }))
# per_gene_sigB <- per_gene_sigB[order(-per_gene_sigB$cor_rho), ]
#

# nrow(per_gene_sigB)))
# cat(sprintf("median cor: %.3f | cor>0.7: %d | cor<0: %d\n",
# median(per_gene_sigB$cor_rho),
# sum(per_gene_sigB$cor_rho > 0.7),
# sum(per_gene_sigB$cor_rho < 0)))
# head(per_gene_sigB, 20)
#

# cat(sprintf("\noverall on merged-sig pairs: Pearson r = %.3f, sign agree = %.1f%%\n",
# cor(cmp_sig$rho_pool, cmp_sig$rho_merged),
# 100 * mean(sign(cmp_sig$rho_pool) == sign(cmp_sig$rho_merged))))
#
# library(ggplot2)
# ggplot(cmp_sig, aes(rho_pool, rho_merged)) +
# geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
# geom_point(alpha = 0.2, size = 0.8) +
# geom_abline(slope = 1, intercept = 0, linetype = 2) +
# geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
# labs(title = "Meta vs merged rho, restricted to merged-significant pathways",
# x = "rho (REML meta)", y = "rho (merged cohorts)") +
# theme_bw()

# 7. Fig 3b - gene number per MASM subcluster

plot_data <- data.frame(gene = names(gene_cluster),
                        cluster = as.factor(gene_cluster)) %>%
  dplyr::count(cluster, name = "n")

# Representative genes shown under each subcluster label
masm_label_genes <- list(
  "1" = c("MTOR", "HIF1A", "TGFA"),
  "2" = c("GLI1", "SMAD4", "EGLN3"),
  "3" = c("TFEB", "PARD6B", "KRT23"),
  "4" = c("AXIN2", "TCF7L1", "NOX1"),
  "5" = c("MYC", "RPS6", "E2F5"),
  "6" = c("AREG", "CDKN2B", "MMP2"),
  "7" = c("CCNB2", "CCNA2", "SKP2")
)

# Build rich-text axis labels: bold MASM title + gene names in cluster color
cluster_levels <- levels(plot_data$cluster)
axis_labels <- sapply(cluster_levels, function(cl) {
  paste0("<span style='color:black; font-weight:bold;'>MASM", cl, "</span><br>",
         "<span style='color:", masm_colors_2[cl], ";'>",
         paste(masm_label_genes[[cl]], collapse = "<br>"), "</span>")
})
names(axis_labels) <- cluster_levels

p_bar <- ggplot(plot_data, aes(x = cluster, y = n, fill = cluster)) +
  geom_col(width = 0.6, color = "white", alpha = 0.8) +
  geom_text(aes(label = n), vjust = -0.5, size = 4) +
  scale_fill_manual(values = masm_colors_2) +
  scale_x_discrete(labels = axis_labels) +
  theme_hd_plain() +
  labs(x = NULL, y = "Gene number per MASM") +
  theme(axis.text.x = element_markdown(angle = 0, hjust = 0.5, lineheight = 1.2, size = 10),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "lines"))
p_bar

# 8. Fig 3c - annotated correlation heatmap (via plot_gene_path_cor_heatmap)

genes_to_mark <- unlist(masm_label_genes, use.names = FALSE)

res_cor_heatmap <- plot_gene_path_cor_heatmap(
  mat_input = mat_wide_filtered,
  module_anno_df = module_anno_df,
  col_anno = col_anno,
  mark_gene_vec = genes_to_mark,
  plot_title = "Correlation matrix\nSignaling gene changes vs metabolic pathway reallocation"
)
heat_grob <- grid.grabExpr(draw(res_cor_heatmap$draw_plot))


run_epi_pca <- function(mat, title, standardize = FALSE, ...) {
  if (standardize) mat <- standardize_rowname(mat)
  plot_epi_pca(epi_all_mat = mat, patient_infor = patient_infor,
               ncp = 4, top_var_n = 12, plot_title = title, ...)
}

# a. Co-reallocation pathways (network modules)
pathway_pca_result <- run_epi_pca(
  build_epi_pathway_mat(fc_pathwayMeanSCT, as.character(module_anno_df$pathway)),
  "Co-reallocation pathways")

# b. All signaling genes
epi_idx <- str_detect(rownames(cpm_all_select_signal), "Epithelial")
signal_all_pca_result <- run_epi_pca(
  cpm_all_select_signal[epi_idx, ], "All signaling genes", standardize = TRUE)

# c. MASM-related signaling genes
signal_masm_pca_result <- run_epi_pca(
  cpm_all_select_signal[epi_idx, names(gene_cluster)],
  "MASM related signaling genes", standardize = TRUE)

# d. Projected MASM scores: PC1 loadings from tumor-vs-normal FC, projected onto
# tumor expression
epi_gene_df_all <- as.data.frame(cpm_all_select_signal[epi_idx, names(gene_cluster)])

all_loading_list <- readRDS( "./all_loading_list.rds")
all_pc1_df <- readRDS( "./all_pc1_df.rds")

masm_score_res <- calc_masm_weight_score(
  expr_all_mat = epi_gene_df_all,
  loading_list = all_loading_list,
  pc1_df = all_pc1_df,
  patient_infor = patient_infor
)
all_data_tumor <- masm_score_res$tumor_combined_df

masm_pca_result <- run_epi_pca(
  masm_score_res$pca_input_mat, "Projected MASM score",
  standardize = TRUE, show_var_arrow = TRUE)

print(pathway_pca_result$plot)
print(signal_all_pca_result$plot)
print(signal_masm_pca_result$plot)
print(masm_pca_result$plot)


pall4_1 <- ggarrange(p4_flowchart, p_bar,
                     labels = letters[1:2], nrow = 1, ncol = 2,
                     font.label = panel_label_style, widths = c(3, 1.5))
pall4_2 <- ggarrange(pathway_pca_result$plot, signal_all_pca_result$plot, signal_masm_pca_result$plot, masm_pca_result$plot, labels = letters[4:7], nrow = 1,
                     font.label = panel_label_style,common.legend = T,legend = "bottom")
pall4_3 <- ggarrange(pall4_1, heat_grob, pall4_2,
                     labels = c("", letters[3], ""), nrow = 3,
                     font.label = panel_label_style, heights = c(1, 2, 0.8))

save_tiff(pall4_3, "figure3.tiff", width = 18, height = 20)


# 9. Fig 3d-g - MASM PC1 scores and correlation with pathway reallocation

# Epithelial expression of the clustered signaling genes (patient level)

masm_load_res_fc <- calc_masm_loading(
  expr_fc_mat = epi_gene_df,
  gene_cluster = gene_cluster,
  flip_check = "mean_loading",
  export_csv = TRUE
)

# Extract intermediate results
all_pc1_df <- masm_load_res_fc$pc1_df
all_loading_list <- masm_load_res_fc$loading_list

saveRDS(all_loading_list,file = "./all_loading_list.rds")


epi_fc <- sp_result$cell_fc_list$Epithelial

# Patient x MASM PC1 matrix joined with pathway reallocation FC
all_pc_path <- cbind(all_pc1_df, epi_fc[all_pc1_df$patient, ])
all_pc_path$Dataset <- str_split_fixed(rownames(all_pc_path), "_", 3)[, 1]

# Scatter of MASM score vs pathway reallocation (Spearman)
plot_masm_cor <- function(masm, pathway) {
  ggplot(all_pc_path, aes(x = .data[[masm]], y = .data[[pathway]])) +
    geom_point(aes(color = Dataset), size = 2.5) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(values = dataset_colors) +
    theme_hd_minimal_plain() +
    theme(legend.position = "bottom") +
    stat_cor(method = "spearman") +
    xlab(paste(masm, "score"))
}

pcor1 <- plot_masm_cor("MASM7", "Glycolysis/gluconeogenesis")
pcor2 <- plot_masm_cor("MASM5", "Folate metabolism")
pcor3 <- plot_masm_cor("MASM2", "Leukotriene metabolism")
pcor4 <- plot_masm_cor("MASM4", "Leukotriene metabolism")
pcor6 <- plot_masm_cor("MASM7", "Purine synthesis")
pcor7 <- plot_masm_cor("MASM5", "Fatty acid synthesis")

all_pc1_df_data<-all_pc1_df[,-1]
rownames(all_pc1_df_data)<-all_pc1_df$patient
meta_res <- gene_pathway_meta_cor(
    epi_gene_df = as.data.frame(all_pc1_df_data),
    epi_fc      = as.data.frame(epi_fc),
    COR_METHOD  = "spearman",
    MIN_N       = 4,
    MIN_K       = 3,
    merge_groups = list(c("KUL3","KUL5"))
)

res_cmp <- compare_meta_global_cor(
    all_pc1_df = all_pc1_df_data,
    epi_fc = epi_fc,
    meta_res = meta_res,
    COR_METHOD = "spearman",
    merged_p_cutoff = 0.05,
    use_cortest = FALSE
)


cor_res <- res_cmp$cmp_sig %>%
    group_by(gene) %>%
    summarise(
        rho = cor(rho_pool, rho_merged, method = "spearman"),
        p_value = suppressWarnings(cor.test(rho_pool, rho_merged, method = "spearman")$p.value),
        n_pathway = n_distinct(pathway),
        .groups = "drop"
    ) %>%
    mutate(gene = forcats::fct_reorder(gene, rho))

p_bar1=ggplot(cor_res, aes(x = gene, y = rho)) +
    geom_bar(stat = "identity", width = 0.65,
             fill= "#225477", linewidth = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.6, color = "grey20") +
    scale_fill_manual(values = "#6FA8DC") +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    labs(
        title = "Concordance of gene–pathway correlations\nMeta-analysis vs merged dataset",
        x = "MASM", y = expression("Spearman "*rho)
    ) +
    theme_hd_plain() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )  +scale_y_continuous(breaks = seq(0, 1, 0.2),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.05)))


cmp <- res_cmp$cmp_all %>%
    mutate(sig = p_merged < 0.05)


cor_lab <- cmp %>%
    group_by(gene) %>%
    summarise(
        R_all = cor(rho_pool, rho_merged),
        n_sig = sum(sig),
        R_sig = if (n_sig >= 5) cor(rho_pool[sig], rho_merged[sig]) else NA_real_,
        .groups = "drop"
    ) %>%
    mutate(label = ifelse(is.na(R_sig),
                          sprintf("R[all] == %.2f", R_all),
                          sprintf("R[all] == %.2f * ';' ~ R[sig] == %.2f", R_all, R_sig)))

p_masm<-ggplot(cmp, aes(rho_pool, rho_merged)) +
    geom_point(aes(color = sig), size = 1.3, alpha = 0.8) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "black"),
                       labels = c("ns", "p(merged) < 0.05")) +
    geom_text(data = cor_lab, aes(label = label),
              x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
              parse = TRUE, size = 3.2) +
    facet_wrap(~gene,nrow=1) +
    labs(color = NULL)+theme_hd_minimal_plain()+theme(legend.position = "bottom")


score_from_loadings <- function(mat, genes, w) {
  genes <- intersect(genes, colnames(mat))
  w <- w[genes]
  x <- mat[, genes, drop = FALSE]
  numer <- rowSums(sweep(x, 2, w, `*`), na.rm = TRUE)
  denom <- as.vector((!is.na(x)) %*% abs(w))
  score <- numer / denom
  score[denom / sum(abs(w)) < 0.7] <- NA
  score
}

loo_fixed_loading <- function(cl_name, loading_list, expr_fc_mat, metab_mat, pc1_df) {
  ld <- loading_list[[cl_name]]
  w  <- setNames(ld$loading, ld$gene)
  genes <- intersect(ld$gene, colnames(expr_fc_mat))

  ds <- vapply(strsplit(rownames(expr_fc_mat), "_", fixed = TRUE), `[`, character(1), 1)


  s_full <- pc1_df[[cl_name]]
  prof_full <- as.numeric(cor(s_full, metab_mat, use = "pairwise.complete.obs"))


  s_full_w <- score_from_loadings(expr_fc_mat, genes, w)
  cor_check <- cor(s_full, s_full_w, use = "pairwise.complete.obs")


  prof_full_ds <- sapply(unique(ds), function(d) {
    idx <- ds == d
    if (sum(idx) < 5) return(rep(NA_real_, ncol(metab_mat)))
    as.numeric(cor(s_full[idx], metab_mat[idx, ], use = "pairwise.complete.obs"))
  })   # pathways x datasets

  res <- lapply(genes, function(g) {
    rest <- setdiff(genes, g)
    if (length(rest) < 2) return(NULL)


    s_loo <- score_from_loadings(expr_fc_mat, rest, w[rest])
    prof_loo <- as.numeric(cor(s_loo, metab_mat, use = "pairwise.complete.obs"))
    r_merged <- cor(prof_full, prof_loo, use = "pairwise.complete.obs")


    r_ds <- sapply(seq_len(ncol(prof_full_ds)), function(j) {
      idx <- ds == colnames(prof_full_ds)[j]
      if (sum(idx) < 5) return(NA_real_)
      p_loo_ds <- as.numeric(cor(s_loo[idx], metab_mat[idx, ],
                                 use = "pairwise.complete.obs"))
      cor(prof_full_ds[, j], p_loo_ds, use = "pairwise.complete.obs")
    })
    names(r_ds) <- colnames(prof_full_ds)

    data.frame(module = cl_name, dropped = g, cor_check = cor_check,
               r_merged = r_merged,
               r_ds_median = median(r_ds, na.rm = TRUE),
               r_ds_min = min(r_ds, na.rm = TRUE),
               worst_ds = names(which.min(r_ds)))
  })
  do.call(rbind, res)
}

loo_all <- do.call(rbind, lapply(names(all_loading_list), loo_fixed_loading,
                                 all_loading_list, epi_gene_df, epi_fc, all_pc1_df_data))


module_summary <- do.call(rbind, lapply(split(loo_all, loo_all$module), function(x) {
  w <- x[which.min(x$r_merged), ]
  data.frame(module = x$module[1], pivotal_gene = w$dropped,
             min_r_merged = round(w$r_merged, 3),
             r_ds_median = round(w$r_ds_median, 3),
             worst_ds = w$worst_ds,
             cor_check = round(w$cor_check, 4))
}))
print(module_summary)
write.csv(loo_all, "loo_fixed_loading_all.csv", row.names = FALSE)
write.csv(module_summary, "loo_fixed_loading_summary.csv", row.names = FALSE)
plot_df <- module_summary %>%
  mutate(module = forcats::fct_reorder(module, min_r_merged))
pbar_loo<-ggplot(plot_df, aes(x = module, y = min_r_merged)) +
    geom_bar(stat = "identity", width = 0.65,
             fill = "#6FA8DC", color = "#225477", linewidth = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.6, color = "grey20") +
    geom_text(aes(label = pivotal_gene),
              hjust = -0.2, size = 3, angle=45,color = "#225477") +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       limits = c(0, 1.25),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
        title = "Robustness of MASM–pathway correlations",
        x = "MASM",
        y = "Minimum Spearman \u03c1\n after gene removal"
    ) +
    theme_hd_plain() +
    theme(

        axis.text.x = element_text(angle = 45, hjust = 1)
    )

library(tidytext)

p_loo<-ggplot(loo_all, aes(x = reorder_within(dropped, r_merged, module), y = r_merged)) +
    geom_col(width = 0.7, fill = "grey65") +
    geom_point(aes(y = r_ds_min), shape = 21, size = 2,
               fill = "white", color = "#0072B2", stroke = 1) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50") +
    scale_x_reordered() +
    facet_grid(~module, scales = "free_x", space = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_hd_minimal_plain() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.spacing = unit(0.5, "lines")) +
    labs(x = NULL, y = "Spearman correlation coefficient")


library(dplyr); library(tibble); library(metafor)
library(forestploter); library(grid); library(cowplot)


EX_MOD <- "MASM7"; EX_PW <- "Glycolysis/gluconeogenesis"

common_ids <- intersect(rownames(all_pc1_df_data), rownames(epi_fc))
x_all <- all_pc1_df_data[common_ids, EX_MOD]
y_all <- epi_fc[common_ids, EX_PW]
ds <- vapply(strsplit(common_ids, "_", fixed = TRUE), `[`, character(1), 1)


kul_idx <- ds %in% c("KUL3", "KUL5")
for (g in c("KUL3", "KUL5")) {
  gi <- which(ds == g)
  x_all[gi] <- x_all[gi] - mean(x_all[gi], na.rm = TRUE)
  y_all[gi] <- y_all[gi] - mean(y_all[gi], na.rm = TRUE)
}
ds[kul_idx] <- "KUL"
ds <- factor(ds)

meta_df <- do.call(rbind, lapply(levels(ds), function(d) {
  idx <- which(ds == d)
  ct <- cor.test(x_all[idx], y_all[idx], method = "spearman", exact = FALSE)
  data.frame(study = d, n = length(idx), rho = unname(ct$estimate), p_ds = ct$p.value)
}))


fit <- rma(yi = atanh(pmax(pmin(meta_df$rho, 0.9999), -0.9999)),
           vi = 1 / (meta_df$n - 3), method = "REML", test = "knha")
pool_info <- list(pooled_rho = tanh(fit$beta[1]),
                  pooled_rho_lower = tanh(fit$ci.lb),
                  pooled_rho_upper = tanh(fit$ci.ub),
                  pooled_pvalue = fit$pval,
                  I2 = fit$I2, tau2 = fit$tau2, Q_pvalue = fit$QEp)


fp_res <- draw_masm_meta_forest(meta_df, pool_info)
p_b <- fp_res$plot
p_b


mod_coh<-run_module_coherence(t(epi_all_tumor), all_loading_list,n_perm      = 1000)
plot_coh<-mod_coh$plot+ylim(-0.2,0.5)


# 10. Assemble Figure 3

pall4_s1 <- ggarrange(pcor3, pcor4,pcor6,pcor7,plot_coh,
                     labels = letters[1:5], nrow = 1, ncol = 5,
                     font.label = panel_label_style,common.legend = T,legend="bottom")
pall4_s2 <- ggarrange(pall4_s1,p_masm,p_loo,labels = c("",letters[6:7]), nrow = 3, ncol = 1,font.label = panel_label_style)


save_tiff(pall4_s2, "supplementaryfigure9.tiff", width = 18, height = 20)

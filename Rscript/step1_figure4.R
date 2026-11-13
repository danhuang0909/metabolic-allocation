# Figure 4 and associated supplementary figures.
# Run this script from the project root with Rscript.

project_dir <- normalizePath(
  Sys.getenv("METABOLISM_PROJECT_DIR", unset = getwd()),
  mustWork = TRUE
)
source(file.path(project_dir, "Rscript", "step0_config.R"))
setwd(data_dir)

# ==============================================================================
# MASM module scores: PCA, patient heatmap, pseudotime trajectory
# ==============================================================================

# ---- 1. Setup ----------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(ggpubr)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(monocle3)
library(grid)

all_meta <- readRDS("meta_all_sample.rds")
panel_label_style <- list(size = 18, color = "black", face = "bold")
class_palette <- c("Normal" = "#3B9AB2", "Tumor" = "#F21A00")
msi_color <- c("MSS" = "#3498db", "MSI" = "#e67e22")
full_masm <- paste0("Tumor_MASM", 1:7)
masm_cols <- paste0("Tumor_MASM", 1:7)
cell_types <- c("B cell", "Endothelial", "Fibroblast", "Myeloid", "Plasma", "T cell")
ct_pattern <- c("B cell", "Endo", "Fibro", "Myeloid", "Plasma", "T cell")

dataset_colors <- c("CRC-SG1" = "#7FDBB3", "CRC-SG2" = "#77FFDD", "KUL3" = "#F5AD75",
                    "KUL5" = "#85A4E6", "Pelka" = "#F593CE", "SMC" = "#CAF77B","KUL3_5"="#F5AD75")


# Move Normal/Tumor suffix to the end of rownames ("<id>_<group>")
standardize_rowname <- function(mat) {
  rn_tbl <- tibble(old = rownames(mat)) %>%
    mutate(group = str_extract(old, "Normal|Tumor"),
           id_clean = str_remove(old, "_Normal|_Tumor"),
           new = paste0(id_clean, "_", group))
  `rownames<-`(mat, rn_tbl$new)
}

# Epithelial pathway matrix (Normal + Tumor rows) for a given pathway set
build_epi_pathway_mat <- function(fc_result, pathways) {
  epi_rows <- function(mat) mat[str_detect(rownames(mat), "Epithelial"), pathways, drop = FALSE]
  epi_normal_mat <- epi_rows(fc_result$normal)
  epi_tumor_mat <- epi_rows(fc_result$tumor)
  rownames(epi_normal_mat) <- paste0(rownames(epi_normal_mat), "_Normal")
  rownames(epi_tumor_mat) <- paste0(rownames(epi_tumor_mat), "_Tumor")
  rbind(epi_normal_mat, epi_tumor_mat)
}

# GAM smooth of MASM scores along pseudotime for a given module set
plot_masm_pseudotime <- function(meta_long, modules) {
  ggplot(meta_long %>% filter(MASM %in% modules),
         aes(x = pseudotime, y = score)) +
    geom_smooth(method = "gam", linewidth = 1.2, se = TRUE) +
    facet_wrap(~MASM, scales = "free_y", ncol = 1) +
    labs(x = "Pseudotime", y = "MASM Module Score") +
    theme_hd_minimal_plain()
}

# FeaturePlot of MASM scores with consistent styling
plot_masm_features <- function(obj, features) {
  FeaturePlot(obj, features = features, pt.size = 0.25, order = TRUE,
              min.cutoff = "q1", max.cutoff = "q99", ncol = 1, raster = FALSE) &
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "plain"))
}


run_forest <- function(prop_df, cell_type, masm_module, output_tsv) {
  draw_meta_cor_forest(
    celltype_prop_df = prop_df,
    all_data_tumor = all_data_tumor_all_sample ,
    cell_type = cell_type,
    masm_module = masm_module,
    title = paste(cell_type, "proportion vs", sub("Tumor_", "", masm_module)),
    output_tsv = output_tsv
  )
}


# ---- 2. Data loading -----------------------------------------------------------
all_gene_sct_count_matrix <- readRDS(file.path(data_dir, "all_gene_sct_count_matrix.rds"))
expr_signal<-readRDS("./pexpr_signal_signalinggenes_all_samples.rds")
all_loading_list <- readRDS("all_loading_list.rds")

 single_cell_masm <- as.data.frame(calc_masm_weight_score_simple(t( expr_signal), all_loading_list))
meta_data <- readRDS("all_meta_data.rds")
# Per-patient pseudobulk CPM of all genes, and KEGG pathway gene sets
cpm_all <- readRDS("sum_all_genes_count_rawRNA_cpm.rds")
kegg_path_genes <- readRDS("all_kegg_pathway_genes.rds")
all_gene_pathway <- readRDS("all_genes_pathway_infor_select.rds")
allcells_pathwayMeanSCT <- readRDS("allcells_pathwayMeanSCTNormalizedByMetacount.Rds")
patient_infor <- readRDS("single_cell_patient_infor_all_samples.rds")
gene_cluster <- readRDS("gene_cluster.rds")
pseudo_res <- readRDS("pseudo_seurat_all_sample.rds")
pathway_mean_allsct_epi_path_all <- readRDS("pathway_meanallsct_seurat_epi_path_all_samples.rds")

# Signaling-pathway genes (KEGG "signal" entries), excluding metabolic genes
kegg_signal <- kegg_path_genes[grepl("signal", kegg_path_genes$description), ]
cpm_all_select_signal <- cpm_all[, colnames(cpm_all) %in% kegg_signal$symbol &
                                   !colnames(cpm_all) %in% unique(all_gene_pathway$gene)]
cpm_all_select_signal_fcAN<-getFCdataSinglelogAllnormaleach(cpm_all_select_signal,str_split_fixed(rownames(cpm_all_select_signal),"_",4)[,3])

fc_pathwayMeanSCT <- getFCresult(allcells_pathwayMeanSCT, meta_data)
sp_result <- prepare_pathway_stats(fc_data = fc_pathwayMeanSCT,
                                   classification = all_gene_pathway)


result_df <- readRDS("sum_all_genes_count_rawRNA_test.rds")
result_df <- result_df[result_df$cell_type != "Mast", ]
result_df_sig <- result_df[result_df$PValue < 0.05, ]
result_df_sig_epi <- result_df_sig[result_df_sig$cell_type == "Epithelial", ]


non_epi_cpm_mat<-readRDS( "./pseudo_bulk_output/non_epi_merged_cpm.rds")
no_epi_all_tumor_pseudo <- as.data.frame(non_epi_cpm_mat[grepl("Tumor_", rownames(non_epi_cpm_mat)) , ])
rownames(no_epi_all_tumor_pseudo)<-gsub("Tumor_","",rownames(no_epi_all_tumor_pseudo))

all_data_tumor_pseudono_epi <- calc_masm_weight_score_simple(no_epi_all_tumor_pseudo, all_loading_list)

# ---- 3. Epithelial pathway network and Louvain modules -------------------------

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

g <- p_all$g
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

# ---- 4. Fig 4a-d: PCA panels ---------------------------------------------------

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

# ---- 5. Local-pathway PCA on MASM-correlated pathways (exploratory) ------------


rownames(all_pc1_df) <- all_pc1_df$patient
all_pc1_df <- all_pc1_df[, -1]
epi_fc <- sp_result$cell_fc_list$Epithelial
temp_cor <- getcorSpearmantest(all_pc1_df, epi_fc[rownames(all_pc1_df), ])
masm_cor_pathways <- unique(temp_cor$pathway[abs(as.numeric(temp_cor$rho.rho)) > 0.3])

pathway_pca_resultMASM <- run_epi_pca(
  build_epi_pathway_mat(fc_pathwayMeanSCT, masm_cor_pathways),
  "Local pathway allocation\nin all MASM related pathways")
print(pathway_pca_resultMASM$plot)

# ---- 6. Fig 4e: consistency of MASM scores -------------------------------------

corr_table <- calc_masm_corr_from_df(all_data_tumor)
corr_bar_plot <- plot_masm_corr_bar(corr_table) +
  ggtitle("Consistency of MASM scores") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(corr_bar_plot)

# ---- 7. Fig 4f: tumor-projected MASM score heatmap (all patients) --------------

# Project MASM scores onto all-sample pseudobulk tumor epithelial expression
epi_all_tumor <- as.data.frame(pseudo_res$cpm_matrix[grepl("Epi", rownames(pseudo_res$cpm_matrix)), ])

masm_score_res_all <- calc_masm_weight_score(
  expr_all_mat = epi_all_tumor,
  loading_list = all_loading_list,
  pc1_df = all_pc1_df,
  patient_infor = patient_infor
)
all_data_tumor_all_sample <- masm_score_res_all$tumor_combined_df

# Patients x MASM score matrix, short patient ids as rownames
mat_plot <- as.matrix(all_data_tumor_all_sample[, 1:7])
pid_short <- str_split_fixed(rownames(mat_plot), "_", 3)[, 3]
rownames(mat_plot) <- pid_short
colnames(mat_plot) <- gsub("Tumor_", "", colnames(mat_plot))
pcs <- prcomp(mat_plot, center = TRUE, scale. = TRUE)
pc_df <- data.frame(pcs$x[, 1:2], dataset = patient_infor[pid_short,"dataset"], MSI = patient_infor[pid_short,"MSI"])

ggplot(pc_df, aes(PC1, PC2, color = dataset, shape = MSI)) +
geom_point(size = 2.5, alpha = 0.85) +
stat_ellipse(aes(group = MSI), type = "norm", linetype = 2) +
theme_bw()


masm_pca_ggplot(all_data_tumor_all_sample, color_var = "MSI", point_size = 2.6)+ggtitle("scRNA")


# Clinical annotation aligned to matrix rows
row_anno_df <- patient_infor[rownames(mat_plot), ] %>%
    select(dataset, MSI, sideness, stage)


all(rownames(row_anno_df) == rownames(mat_plot))

# Annotation colors: gray fallback for any level missing from the base list
base_color <- list(
    dataset  = dataset_colors,
    MSI      = msi_color,
    sideness = c("Left" = "#73b873", "Right" = "#e86868"),
    stage    = c("Advanced" = "#c82423", "Early" = "#2874A6")
)
fill_anno_color <- function(anno_df, base_color_list) {
    lapply(names(base_color_list), function(var) {
        col_vec <- base_color_list[[var]]
        miss_vals <- setdiff(unique(as.character(anno_df[[var]])), names(col_vec))
        c(col_vec, setNames(rep("#999999", length(miss_vals)), miss_vals))
    }) %>% setNames(names(base_color_list))
}
anno_color_list <- fill_anno_color(row_anno_df, base_color)


# incorrectly said "row"
mat_scaled <- scale(mat_plot)


row_anno <- rowAnnotation(
    df = row_anno_df,
    col = anno_color_list,
    annotation_name_side = "top",
    annotation_legend_param = list(dataset = list(ncol = 1, title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)), MSI = list(ncol = 1, title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)),
                                   sideness = list(ncol = 1, title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)), stage = list(ncol = 1, title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 14)))
)

ht <- Heatmap(
    mat_scaled,
    name = "z-score",
    left_annotation = row_anno,     # top_annotation → left_annotation
    clustering_method_rows = "ward.D",
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,column_names_gp = gpar(fontsize = 16),
    row_dend_width = unit(10, "mm"),
    column_dend_height = unit(15, "mm"),
    col = colorRamp2(seq(-2, 2, length.out = 100),
                     colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)),
    row_title = "Tumor-projected MASM score across patients",
    rect_gp = gpar(col = "grey80", lwd = 0.8),
    heatmap_legend_param = list(direction = "vertical",title_gp = gpar(fontsize = 14),
    labels_gp = gpar(fontsize = 14) ) ,
    column_names_rot = 45
)
heat_grob <- grid.grabExpr(draw(ht, heatmap_legend_side = "right",
                                annotation_legend_side = "right", merge_legends = TRUE))


# ---- 8. Fig 4g: MASM composition in tumor epithelial cells ---------------------


# then MASM scores projected per cell
res_tumor <- plot_masm_proportion(
  seurat_obj = pathway_mean_allsct_epi_path_all,
  filter_group = "Tumor",
  plot_title_prefix = "Tumor Epithelial",
  facet_ncol_box = 7,
  box_output = "Tumor_MASM2_5_7_boxplot.pdf"
)
print(res_tumor$box_plot)
p_composition <- res_tumor$stack_plot2

# ---- 9. Fig 4h-j: monocle3 pseudotime trajectory --------------------------------
library(monocle3)
cds <- new_cell_data_set(
  expression_data = pathway_mean_allsct_epi_path_all@assays$RNA@counts,
  cell_metadata = pathway_mean_allsct_epi_path_all@meta.data,
  gene_metadata = data.frame(
    gene_id = rownames(pathway_mean_allsct_epi_path_all),
    gene_short_name = rownames(pathway_mean_allsct_epi_path_all),
    row.names = rownames(pathway_mean_allsct_epi_path_all)
  )
)

# Reuse the Seurat UMAP embedding (cell order must match)
reducedDims(cds)$UMAP <- pathway_mean_allsct_epi_path_all@reductions$umap@cell.embeddings[colnames(cds), ]

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds, root_pr_nodes = "Y_557")

colData(cds)$pseudotime <- pseudotime(cds)

saveRDS(cds, file = "./pathway_monocle3_select_epi_cds.rds")
cds<-readRDS(file = "./pathway_monocle3_select_epi_cds.rds")
meta_long <- as.data.frame(colData(cds)) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = starts_with("MASM"), names_to = "MASM", values_to = "score")
meta_long$dataset2=ifelse(meta_long$dataset=="Pelka","Pelka","non-Pelka")
p_pseudo <- plot_masm_pseudotime(meta_long, c("MASM2", "MASM5", "MASM6", "MASM7"))
p_feature <- plot_masm_features(pathway_mean_allsct_epi_path_all,
                                c("MASM2", "MASM5", "MASM6", "MASM7"))


mods <- paste0("MASM",1:7)
pair_df <- combn(mods,2)
select_pairs <- paste(pair_df[1,], pair_df[2,], sep = "-")
select_pairs
out5 <- plot_masm_pair_high_barplot(
    pathway_mean_allsct_epi_path_all, quant_cut = 0.75,
    select_pairs = select_pairs,add_sig = F);out5$plot
out0<- plot_masm_pair_high_barplot(
    pathway_mean_allsct_epi_path_all, value_cut = 0,
    select_pairs = c("MASM1-MASM2", "MASM3-MASM4", "MASM1-MASM4",
                     "MASM5-MASM7", "MASM6-MASM7"),add_sig = F);out0$plot


p_masm67 <- plot_umap_two_module_groups_highlight(
obj = pathway_mean_allsct_epi_path_all,
mod1 = "MASM6",
mod2 = "MASM7",col_cohigh = "#FFE02EDB",
quant_cut = 0.75
)+guides(color = guide_legend(nrow = 2,byrow = TRUE,
    override.aes = list(size = 4, stroke = 1),
    keywidth = unit(4,"mm"),keyheight = unit(4,"mm")))

p_masm14 <- plot_umap_two_module_groups_highlight(
obj = pathway_mean_allsct_epi_path_all,
mod1 = "MASM1",
mod2 = "MASM4",col_cohigh = "#FFE02EDB",
quant_cut = 0.75
)+guides(color = guide_legend(nrow = 2,byrow = TRUE,
    override.aes = list(size = 4, stroke = 1),
    keywidth = unit(4,"mm"),keyheight = unit(4,"mm")))


p_masm57 <- plot_umap_two_module_groups_highlight(
obj = pathway_mean_allsct_epi_path_all,
mod1 = "MASM5",
mod2 = "MASM7",col_cohigh = "#FFE02EDB",
quant_cut = 0.75
)+guides(color = guide_legend(nrow = 2,byrow = TRUE,
    override.aes = list(size = 4, stroke = 1),
    keywidth = unit(4,"mm"),keyheight = unit(4,"mm")))


p_masm34 <- plot_umap_two_module_groups_highlight(
obj = pathway_mean_allsct_epi_path_all,
mod1 = "MASM3",
mod2 = "MASM4",col_cohigh = "#FFE02EDB",
quant_cut = 0.75
)+guides(color = guide_legend(nrow = 2,byrow = TRUE,
    override.aes = list(size = 4, stroke = 1),
    keywidth = unit(4,"mm"),keyheight = unit(4,"mm")))


expr_df <- as.matrix(counts(cds)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  left_join(., colData(cds) %>% as.data.frame() %>% rownames_to_column("cell_id") %>% select(cell_id, pseudotime)) %>% left_join(meta_long[,c("cell.ID","dataset2")],by=c("cell_id"="cell.ID"))

p7=ggplot(expr_df, aes(x = pseudotime, y = .data[["Oxidative phosphorylation"]],col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+
    theme_bw()

p8=ggplot(expr_df, aes(x = pseudotime, y = .data[["Fatty acid synthesis"]],col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()

p6=ggplot(expr_df, aes(x = pseudotime, y = .data[["Purine synthesis"]],col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()

p5=ggplot(expr_df, aes(x = pseudotime, y = .data[["Glutathione metabolism"]],col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()

p11=ggplot(meta_long[meta_long$MASM=="MASM5",], aes(x = pseudotime, y = score,col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+ylab("MASM5")+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()
p12=ggplot(meta_long[meta_long$MASM=="MASM7",], aes(x = pseudotime, y = score,col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+
    ylab("MASM7")+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()

p13=ggplot(meta_long[meta_long$MASM=="MASM1",], aes(x = pseudotime, y = score,col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+ylab("MASM1")+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()
p14=ggplot(meta_long[meta_long$MASM=="MASM2",], aes(x = pseudotime, y = score,col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+
    ylab("MASM2")+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()

p15=ggplot(meta_long[meta_long$MASM=="MASM3",], aes(x = pseudotime, y = score,col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+ylab("MASM3")+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()
p16=ggplot(meta_long[meta_long$MASM=="MASM4",], aes(x = pseudotime, y = score,col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+
    ylab("MASM4")+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()
p17=ggplot(meta_long[meta_long$MASM=="MASM6",], aes(x = pseudotime, y = score,col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+
    ylab("MASM6")+ xlab("Pseudotime based on metabolic pathways")+
    labs(color = "Dataset") +
    theme_hd_plain()


pPyrimi6=ggplot(expr_df, aes(x = pseudotime, y = .data[["Purine synthesis"]],col=dataset2))+
    geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6)) +
    theme_hd_plain()
pPyrimi6
ggarrange(p12,p11,p5,p7,p6,p8,nrow=3,ncol=2)

p=FeaturePlot(pathway_mean_allsct_epi_path_all,
features = c("MASM5", "MASM7"),
pt.size = 0.25,
            order = TRUE,
            min.cutoff = "q1",
            max.cutoff = "q99",
            ncol = 1,
            raster = FALSE,
            blend = TRUE,
            cols = c("lightgray","#d62728","#1f77b4") # gray, MYC-red, FAS-blue
)
df <- pathway_mean_allsct_epi_path_all@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  left_join(pathway_mean_allsct_epi_path_all@meta.data %>% rownames_to_column("cell") %>% select(cell, MASM5, MASM7))

p_masm57<-ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = MASM5), alpha = 0.4, size=0.2) +
  geom_point(aes(color = MASM7), alpha = 0.3, size=0.2) +
  theme_bw()


colData(cds)$Group <- factor(colData(cds)$class_new, levels = c("Normal", "Tumor"))
p_umap <- plot_cells(
  cds,
  color_cells_by = "Group",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  cell_size = 0.5,
  alpha = 0.7,
  label_roots = TRUE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  label_cell_groups = FALSE,
  graph_label_size = 3.5,
  reduction_method = "UMAP"
) +
  scale_color_manual(values = class_palette) +
  ggtitle("Epithelial pseudotime trajectory based on co-reallocation pathways") +
  theme_hd_plain()
print(p_umap)


p_umap2 <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    trajectory_graph_color = "black",
    cell_size = 0.5,
    alpha = 0.7,
    label_roots = TRUE,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    label_cell_groups = FALSE,
    graph_label_size = 3.5,
    reduction_method = "UMAP"
) + ggtitle("Epithelial pseudotime trajectory\nbased on co-reallocation pathways") +
  theme_hd_plain()


  filter_epi<-all_meta[all_meta$Cell_type=="Epithelial",]
  filter_epi<- as.data.frame(colData(cds))
# inputs <- list(
# expr         = all_gene_sct_count_matrix[,filter_epi$cell.ID],          # genes x cells, log-normalized

# loading_list = all_loading_list,         # list(MASM1..MASM7 = named numeric)

# )
# res_single <- run_sc_diagnosis(inputs, out_dir = "results/sc_diag",
# pt_col = "pseudotime", patient_col = "patient")
# new_masm<-res_single$new_scores
# colnames(new_masm)<-paste("Rank",colnames(new_masm),sep="_")
# filter_epi<-cbind(filter_epi ,new_masm[filter_epi$cell.ID,])
# filter_epi$dataset2=ifelse(filter_epi$dataset=="Pelka","Pelka","non-Pelka")
# p6=ggplot(filter_epi ,aes(x = pseudotime, y = Rank_MASM5,col=dataset2))+
# geom_smooth(method = "gam",  se = TRUE, formula = y ~ s(x, k = 6))+
# theme_bw()
# p6
non_epi_meta<-readRDS(file = "meta_all_sample_non_epi.rds")
celltype_prop_df_detail <- calc_celltype_prop(non_epi_meta,
                                              "cell_type_new")
celltype_prop_df_detail$Cell_type <- gsub("\n", " ", celltype_prop_df_detail$Cell_type)


meta_results_detail <- run_all_meta_cor_pairs(
  celltype_prop_df = celltype_prop_df_detail,
  all_data_tumor = all_data_tumor_all_sample,
  dataset_col = "dataset",
  patient_col = "patient",
  p_threshold = 0.05
)
p_masm_cell <- plot_masm_cell_heatmap(
  sig_results_df = meta_results_detail$sig_results,
  full_masm_list = full_masm,
  full_cell_list = unique(celltype_prop_df_detail$Cell_type),
  title = "Pooled correlation\nMASM modules vs cell-type proportions"
)+theme(plot.margin = margin(12, 15, 22, 12))
print(p_masm_cell)

# ---- 5. Forest plots for selected associations ------------------------------------

result1 <- run_forest(celltype_prop_df_detail, "CD8+ Effector", "Tumor_MASM6", "Meta_Cor_Epithelial_MASM6.tsv")
result2 <- run_forest(celltype_prop_df_detail, "CD8+ Effector", "Tumor_MASM7", "Meta_Cor_Epithelial_MASM7.tsv")
result3 <- run_forest(celltype_prop_df_detail, "Treg", "Tumor_MASM2", "Meta_Cor_Treg_MASM2.tsv")
result4 <- run_forest(celltype_prop_df_detail, "CD8+ T", "Tumor_MASM7", "Meta_Cor_Effector_MASM7.tsv")


cellchat_dir <- file.path(data_dir, "cellchat_patient_result")
file_list <- list.files(cellchat_dir, pattern = "\\TRUE.rds$", full.names = TRUE)
names(file_list) <- basename(file_list) %>% sub(pattern = "\\.rds$", replacement = "")
cellchat_list <- lapply(file_list, readRDS)
cat("Loaded", length(cellchat_list), "CellChat objects\n")

# Extract pathway-level communication probabilities to long format
extract_cc_pathway <- function(cc_obj, patient_id) {
  as.data.frame.table(cc_obj@netP$prob) %>%
    set_names(c("source_cell", "target_cell", "pathway", "prob")) %>%
    mutate(patient = patient_id)
}

patient_names <- names(cellchat_list) %>%
  str_remove("cellchat_") %>%
  str_remove("_popTRUE$")
cc_all_long <- map2_dfr(cellchat_list, patient_names, extract_cc_pathway)

# Join with MASM scores (patient key matches after suffix removal)
masm_long <- all_data_tumor_all_sample %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = starts_with("Tumor_MASM"),
               names_to = "masm_name", values_to = "masm_score")
cc_masm <- cc_all_long %>% inner_join(masm_long, by = "patient")

# Keep communications active in >= 50% of patients within >= 3 datasets
comm_pass <- cc_masm %>%
  group_by(masm_name, pathway, source_cell, target_cell, dataset) %>%
  summarise(frac_active = mean(prob > 1e-7), .groups = "drop") %>%
  filter(frac_active >= 0.5) %>%
  count(masm_name, pathway, source_cell, target_cell, name = "n_ds_pass") %>%
  filter(n_ds_pass >= 3)
cc_masm_filtered <- cc_masm %>%
  inner_join(comm_pass, by = c("masm_name", "pathway", "source_cell", "target_cell"))

# Spearman correlation per communication pair, requiring >= 8 valid patients
safe_spearman <- function(x, y) {
  tryCatch({
    ct <- cor.test(x, y, method = "spearman")
    tibble(rho = unname(ct$estimate), p_raw = ct$p.value)
  }, error = function(e) tibble(rho = NA_real_, p_raw = NA_real_))
}

cor_res <- cc_masm_filtered %>%
  group_by(masm_name, pathway, source_cell, target_cell) %>%
  group_modify(~{
    n_valid <- sum(!is.na(.x$masm_score) & !is.na(.x$prob))
    if (n_valid >= 8) safe_spearman(.x$masm_score, .x$prob) %>% mutate(n_valid = n_valid)
    else tibble(rho = NA_real_, p_raw = NA_real_, n_valid = n_valid)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_raw, method = "BH"))

candidate_hits <- cor_res %>%
  filter(n_valid >= 8) %>%
  drop_na(rho, p_raw) %>%
  filter(p_adj < 0.05, abs(rho) >= 0.2)
cat("Candidate hits:", nrow(candidate_hits), "\n")

# ---- 7. Per-cohort correlations and cross-cohort consistency -----------------------

# Merge KUL3 and KUL5 into one cohort group
cc_masm_filtered <- cc_masm_filtered %>%
  mutate(ds_group = ifelse(dataset %in% c("KUL3", "KUL5"), "KUL3_5", dataset))

safe_ds_spearman <- function(sub_df) {
  n_pat <- n_distinct(sub_df$patient)
  n_valid <- sum(!is.na(sub_df$masm_score) & !is.na(sub_df$prob))
  out <- tibble(ds_group = unique(sub_df$ds_group), n_pat_ds = n_pat, n_valid_ds = n_valid)
  if (n_pat < 5 || n_valid < 5) return(out %>% mutate(rho_ds = NA_real_, p_ds = NA_real_))
  tryCatch({
    ct <- cor.test(sub_df$masm_score, sub_df$prob, method = "spearman")
    out %>% mutate(rho_ds = unname(ct$estimate), p_ds = ct$p.value)
  }, error = function(e) out %>% mutate(rho_ds = NA_real_, p_ds = NA_real_))
}

per_cohort_cor <- cc_masm_filtered %>%
  group_by(masm_name, pathway, source_cell, target_cell, ds_group) %>%
  group_modify(~safe_ds_spearman(.x)) %>%
  ungroup() %>%
  mutate(source_cell = gsub("\n", " ", source_cell),
         target_cell = gsub("\n", " ", target_cell),
         pair_id = str_c(masm_name, pathway, source_cell, target_cell, sep = "_"))
write_tsv(per_cohort_cor, "cellchat_per_cohort_cor.tsv")

# Direction consistency across cohorts: |rho| <= 0.1 counts as neutral;
# conflicting if >= 2 cohorts positive AND >= 2 negative
consist_summary <- per_cohort_cor %>%
  drop_na(rho_ds) %>%
  group_by(pair_id, masm_name, pathway, source_cell, target_cell) %>%
  summarise(
    n_effective_cohort = n(),
    n_pos = sum(rho_ds > 0.1),
    n_neg = sum(rho_ds < -0.1),
    n_neutral = sum(abs(rho_ds) <= 0.1),
    n_cohort_sig = sum(p_ds < 0.05),
    rho_mean = mean(rho_ds),
    rho_sd = sd(rho_ds),
    .groups = "drop"
  ) %>%
  mutate(conflict = (n_pos >= 2) & (n_neg >= 2))

consistent_all <- consist_summary %>%
  filter(n_effective_cohort >= 4, !conflict, n_cohort_sig >=1, abs(rho_mean) > 0.3)

select_patient <- per_cohort_cor %>% filter(pair_id %in% consistent_all$pair_id)

# ---- 8. Random-effects meta-analysis per consistent pair -------------------------

# Fisher-z pooled correlation (REML + Hartung-Knapp) for one CCI pair
pool_pair_meta <- function(meta_df, p_threshold = 0.05) {
  meta_obj <- metagen(
    TE = atanh(meta_df$rho),
    seTE = 1 / sqrt(meta_df$n - 3),
    studlab = meta_df$study,
    sm = "ZCOR",
    hakn = TRUE,
    method.tau = "REML"
  )
  info <- meta_df %>% slice(1)
  status <- ifelse(meta_obj$pval.random < p_threshold, "*** SIGNIFICANT", "Not significant")

  message(sprintf("Pair %s | Overall: rho = %.3f (%.3f to %.3f), p = %.4f, I2 = %.1f%%",
                  info$pair_id, tanh(meta_obj$TE.random), tanh(meta_obj$lower.random),
                  tanh(meta_obj$upper.random), meta_obj$pval.random, meta_obj$I2 * 100))
  message(paste0("  Status: ", status))
  tibble(
    pair_id = info$pair_id, masm_name = info$masm_name, pathway = info$pathway,
    source_cell = info$source_cell, target_cell = info$target_cell,
    n_studies = nrow(meta_df), n_total = sum(meta_df$n),
    pooled_rho = tanh(meta_obj$TE.random),
    pooled_rho_lower = tanh(meta_obj$lower.random),
    pooled_rho_upper = tanh(meta_obj$upper.random),
    pooled_pvalue = meta_obj$pval.random,
    I2 = meta_obj$I2, tau2 = meta_obj$tau2, Q = meta_obj$Q, Q_pvalue = meta_obj$p.Q,
    model_type = "random (REML + hakn)",
    is_significant = meta_obj$pval.random < p_threshold,
    status = status
  )
}

meta_final <- bind_rows(lapply(unique(select_patient$pair_id), function(pid) {
  meta_df <- select_patient %>%
    filter(pair_id == pid) %>%
    rename(study = ds_group, n = n_pat_ds, rho = rho_ds) %>%
    select(study, n, rho, pair_id, masm_name, pathway, source_cell, target_cell)
  if (nrow(meta_df) < 2) {
    message("Pair: ", pid, " | insufficient studies (<2), skipping")
    return(NULL)
  }
  pool_pair_meta(meta_df)
})) %>%
  mutate(padj_pooled = p.adjust(pooled_pvalue, method = "fdr"))

# Keep pairs that are also significant in the meta-analysis

consistent_all <- consistent_all %>%
  filter(pair_id %in% meta_final$pair_id[meta_final$is_significant]) %>%
  mutate(masm_name_new = gsub("Tumor_", "", masm_name))

# ---- 9. Fig 4m: flowchart ----------------------------------------------------------

img <- rasterGrob(readTIFF(file.path(fig_dir, "figure6k.tif"),
                           native = TRUE), interpolate = TRUE)
p6_flowchart <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
  geom_blank() +
  annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

# ---- 10. Fig 4h: number of conserved CCIs per MASM (diverging tile plot) -----------

source_heat_df <- consistent_all %>%
  mutate(sign_group = ifelse(rho_mean > 0, "Positive", "Negative")) %>%
  count(masm_name_new, source_cell, sign_group, name = "cci_count") %>%
  mutate(
    sign_group = factor(sign_group, levels = c("Positive", "Negative")),
    masm_name_new = factor(masm_name_new, levels = paste0("MASM", 1:7)),
    fill_value = ifelse(sign_group == "Positive", cci_count, -cci_count),
    label_num = abs(fill_value)
  )
max_range <- max(abs(source_heat_df$fill_value))

p_count <- ggplot(source_heat_df, aes(x = source_cell, y = sign_group)) +
  geom_tile(aes(fill = fill_value), color = "white", linewidth = 0.3) +
  geom_text(data = . %>% filter(label_num > 0), aes(label = label_num),
            size = 2.7, color = "black") +
  scale_fill_gradient2(low = "#003380", mid = "#f8f8f8", high = "#990000",
                       midpoint = 0, limits = c(-max_range, max_range), guide = "none") +
  facet_wrap(~masm_name_new, ncol = 1, strip.position = "right") +
  labs(title = "Number of conserved cell-cell interactions\nassociated with MASM states",
       x = "Source Cell (Sender)", y = "Association with MASM score") +
  theme_hd_minimal_plain() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1, size = 13),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 12, angle = 90),
        plot.margin = margin(12, 15, 22, 12))

# ---- 11. Fig 4i-j: CCI circle plots --------------------------------------------------


cell_col_full <- c(
  # T cell family (green)
  "T cell" = "#B3DE69", "CD8+ T" = "#99C84A", "CD8+ Effector" = "#7FB030", "Treg" = "#D0E899",
  # Fibroblast family (coral red)
  "Fibroblast" = "#FB8072", "myCAF" = "#F86454", "ecmCAF" = "#FC9C8C",
  # Myeloid family (blue)
  "Myeloid" = "#80B1D3", "Inflammatory macrophage" = "#5A98C2",
  "TAM/M2-like myeloid" = "#A3C4E0", "cDC2" = "#3E7FBF",
  # Other main types
  "B cell" = "#8DD3C7", "Endothelial" = "#FFFFB3",
  "Epithelial" = "#BEBADA", "Plasma" = "#FDB462"
)

# Assign contrast colors to any cell type missing from the palette

all_cells <- sort(unique(c(consistent_all$source_cell, consistent_all$target_cell)))
miss_cells <- setdiff(all_cells, names(cell_col_full))
conflict_extra_cols <- c("#E41A1C", "#377EB8", "#984EA3", "#FF7F00", "#4DAF4A")
if (length(miss_cells) > 0) {
  cell_col_full <- c(cell_col_full,
                     setNames(conflict_extra_cols[seq_along(miss_cells)], miss_cells))
}
cell_colors <- cell_col_full[all_cells]

p_list <- list()
for (m in paste0("MASM", 1:7)) {
  df_sub <- consistent_all %>% filter(masm_name_new == m)
  if (nrow(df_sub) == 0) {
    message(m, ": no qualified pairs, skip")
    next
  }
  p <- plot_cci_circle(df_sub, m, all_cells, cell_colors)+theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  if (!is.null(p)) p_list[[m]] <- p
}


p_combined <- wrap_plots(p_list, ncol = 4, nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
p_combined

# ---- 12. Fig 4k-l: example pair forest plot and scatter ------------------------------


pid_example <- "Tumor_MASM6_TENASCIN_ecmCAF_Epithelial"
res_fp <- draw_cci_meta_forest(
  target_pair_id = pid_example,
  select_patient = select_patient,
  meta_final = meta_final
)
print(res_fp$plot)

pid_example <- "Tumor_MASM7_IFN-II_CD8+ Effector_Epithelial"
res_fp2 <- draw_cci_meta_forest(
  target_pair_id = pid_example,
  select_patient = select_patient,
  meta_final = meta_final
)
print(res_fp2$plot)


cc_target <- cc_masm_filtered %>%
  filter(source_cell == "Epithelial", target_cell == "CD8+ T",
         pathway == "MIF", masm_name == "Tumor_MASM7")

p_scatter <- ggplot(cc_target, aes(x = masm_score, y = prob)) +
  geom_point(aes(color = ds_group), size = 2.8, alpha = 0.75) +
  stat_cor(method = "spearman") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(name = "Dataset", values = dataset_colors) +
  labs(x = "MASM7 score",
       y = "Cell communication probability") + ggtitle("Epithelial → CD8+ T (MIF)")+ theme_hd_minimal_plain()
print(p_scatter)


# ---- 10. Assemble Figure 4 -------------------------------------------------------

pall4_1 <- ggarrange(p11,p12,p8,pPyrimi6,
                     labels = letters[2:5], nrow = 2, ncol = 2,
                     font.label = panel_label_style, common.legend = TRUE, legend = "bottom")


pall4_2 <- ggarrange(p_umap, pall4_1 , nrow = 1, ncol = 2,labels = c(letters[1],""),
                     font.label = panel_label_style,widths = c(1,1.5))
pall4_3 <- ggarrange(result1$plot,result2$plot, nrow = 2, ncol = 1,labels = c(letters[7:8]),
                     font.label = panel_label_style)

pall4_4 <- ggarrange(p_masm_cell, pall4_3,p_count  , nrow = 1, ncol = 3,labels = c(letters[6],"",letters[9]),
                     font.label = panel_label_style,widths = c(1.5,1,1.5))
pall4_5 <- ggarrange(res_fp$plot,res_fp2$plot, nrow = 2, ncol = 1,labels = c(letters[12:13]),
                     font.label = panel_label_style)
pall4_6 <- ggarrange(p_list$MASM6, p_list$MASM7,nrow = 1, ncol = 2,labels = c(letters[10:11]),
                     font.label = panel_label_style,widths = c(1,1),common.legend = T,legend = "bottom")
pall4_7 <- ggarrange(pall4_6,pall4_5 , nrow = 1, ncol = 2,widths = c(2,1))

pall4_8<- ggarrange(pall4_2, pall4_4,pall4_7 , nrow = 3, ncol = 1,
                     font.label = panel_label_style)


save_tiff(pall4_8, "figure4.tiff", width = 20, height = 20, path = fig_dir)

# ==============================================================================
# Supplementary Figure 10
# ==============================================================================
palls0<-ggarrange(p_umap2,p14,p17,p13,nrow = 1, ncol = 4,labels = c(letters[1:4]),
                     font.label = panel_label_style)
palls1 <- ggarrange(p15,p16,out5$plot, nrow = 1, ncol = 3,labels = c(letters[5:7]),
                     font.label = panel_label_style,widths = c(1,1,2))
palls2 <- ggarrange(p_masm57 , p_masm67,p_masm34, p_masm14 , nrow = 1, ncol = 4,labels = c(letters[8:11]),font.label = panel_label_style)

palls3<- ggarrange(palls0,palls1, palls2  , nrow = 3, ncol = 1)


save_tiff(palls3, "supplementaryfigure10.tiff", width = 20, height = 20, path = fig_dir)


# ==============================================================================
# Supplementary Figure 11
# ==============================================================================

temp=all_meta[all_meta$class_new=="Tumor",]
x=aggregate(temp[,c("nCount_RNA","nFeature_RNA","percent.mt")],list(patient=temp$patient,cell_types=temp$Cell_type),mean)
x_epi=x[x$cell_types=="Epithelial",]
 x_epi<- x_epi %>% left_join(all_data_tumor_all_sample,by="patient")
  x_epi$dataset[x_epi$dataset %in% c("KUL3","KUL5")]="KUL3_5"
p1qc=ggplot(x_epi,aes(x=nCount_RNA,y=Tumor_MASM6,col=dataset))+geom_point()+stat_cor()+theme_hd_plain()+
  scale_color_manual(name = "Dataset", values = dataset_colors) +xlab("UMI counts") +ylab("MASM6 score")
p2qc=ggplot(x_epi,aes(x=percent.mt,y=Tumor_MASM6,col=dataset))+geom_point()+stat_cor()+theme_hd_plain()+
  scale_color_manual(name = "Dataset", values = dataset_colors) +xlab("Mitochondrial fraction") +ylab("MASM6 score")
p3qc=ggplot(x_epi,aes(x=nFeature_RNA,y=Tumor_MASM6,col=dataset))+geom_point()+stat_cor()+theme_hd_plain()+
  scale_color_manual(name = "Dataset", values = dataset_colors) +xlab("Number of Detected genes") +ylab("MASM6 score")
p4qc=ggplot(x_epi,aes(x=nFeature_RNA,y=Tumor_MASM7,col=dataset))+geom_point()+stat_cor()+theme_hd_plain()+
  scale_color_manual(name = "Dataset", values = dataset_colors) +xlab("Number of Detected genes") +ylab("MASM7 score")


masm_cols <- paste0("Tumor_MASM", 1:7)
res_all <- run_masm_meta(
  epi_mat = all_data_tumor_all_sample, cpm_mat = pseudo_res$cpm_matrix,
  all_loading_list = all_loading_list, patient_infor = patient_infor,
  cell_types = cell_types, ct_pattern = ct_pattern, masm_cols = masm_cols,
  samples_keep = NULL, pelka_alone = F, min_n = 5)

masm6_data<-res_all$meta_all[res_all$meta_all$fi=="Tumor_MASM6" & res_all$meta_all$ej=="Tumor_MASM6",]

library(tidyverse)

df_masm6 <- masm6_data %>%
  mutate(celltype = factor(celltype, levels = celltype[order(-rho_pool)]))

star_of <- function(p) ifelse(is.na(p), "",
                              ifelse(p < 0.001, "***",
                                     ifelse(p < 0.01,  "**",
                                            ifelse(p < 0.05,  "*", "ns"))))

p_masm6_cor<-ggplot(df_masm6, aes(x = celltype, y = rho_pool)) +

  geom_errorbar(aes(ymin = ci_l, ymax = ci_u),
                width = 0.22, linewidth = 0.6, color = "grey40") +
  geom_col(width = 0.6, fill = "#C4698A", color = "grey25", linewidth = 0.4) +

  geom_text(aes(label = star_of(pval), y = ci_u),
            vjust = -0.6, size = 4, fontface = "bold",
            color = "grey25") +

  geom_text(aes(label = sprintf("%.2f", rho_pool)),
            vjust = 1.4, size = 3, color = "white", fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Cell type",
       y = expression("Pooled Spearman's " * rho),
       title = "MASM6 in epithelial correlation\nwith MASM6 in other cell types") +
  theme_hd_plain()+theme(axis.text.x = element_text(angle = 30,hjust=1))


all_meta_main <- all_meta %>%
  filter(!Cell_type %in% c("Entericglial", "Peri", "Mast", "Schwann", "SmoothMuscle"),
         class_new == "Tumor")
celltype_prop_df <- calc_celltype_prop(all_meta_main, "Cell_type")

p_stack_celltypes <- plot_patient_cellprop_stack(df = celltype_prop_df) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"))
print(p_stack_celltypes)

library(tidyverse)
library(ggpubr)   # stat_compare_means
epi_all_tumor <- as.data.frame(pseudo_res$cpm_matrix[grepl("Tumor_Epi", rownames(pseudo_res$cpm_matrix)) , ])
rownames(epi_all_tumor)<-gsub("Tumor_","",rownames(epi_all_tumor))
all_data_tumor <- calc_masm_weight_score_simple(epi_all_tumor, all_loading_list)


masm_cols <- intersect(names(all_data_tumor), paste0("MASM", 1:7))

df_wide <- bind_rows(
  all_data_tumor %>%
    mutate(core_id = str_split_fixed(patient, "_",3)[,3],
           celltype  = "Epithelial"),
  all_data_tumor_pseudono_epi %>%
    mutate(core_id =str_split_fixed(patient, "_",3)[,3],
           celltype  = "Non-epithelial")
) %>%
  select(core_id, celltype, all_of(masm_cols)) %>%
  pivot_longer(all_of(masm_cols), names_to = "MASM", values_to = "value") %>%

  group_by(core_id, MASM, celltype) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = celltype, values_from = value) %>%

  filter(!is.na(Epithelial), !is.na(`Non-epithelial`)) %>%
  mutate(MASM = factor(MASM, levels = paste0("MASM", 1:7)))


df_long <- df_wide %>%
  pivot_longer(c(Epithelial, `Non-epithelial`),
               names_to = "celltype", values_to = "value") %>%
  mutate(celltype = factor(celltype,
                           levels = c("Non-epithelial", "Epithelial")))


p_exp=ggplot(df_long, aes(x = celltype, y = value, fill = celltype)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  geom_point(position = position_jitter(width = 0.12),
             size = 0.8, alpha = 0.4, color = "grey40") +
  scale_fill_manual(values = c("Non-epithelial" = "#6699CC",
                               "Epithelial"     = "#E67E22")) +
  facet_wrap(~MASM, nrow = 1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    label.x = 1.5,
    size = 4,
    hide.ns = TRUE
  ) +
  labs(x = "", y = "MASM score", fill = "Cell type") +
  theme_hd_minimal_plain() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )


res_full_cor <- run_full_cross_cor_matrix(
  celltype_prop_df = celltype_prop_df,
  all_data_tumor = all_data_tumor_all_sample,
  output_path = "Full_Cell_MASM_Cross_Correlation.pdf"
)

# Pooled correlations: main cell types (CRC-SG1 excluded)
meta_results_main <- run_all_meta_cor_pairs(
  celltype_prop_df = celltype_prop_df,
  all_data_tumor = all_data_tumor_all_sample,
  dataset_col = "dataset",
  patient_col = "patient",
  p_threshold = 0.05
)
p_all <- plot_block_heatmap_gg(
  sig_results_df = meta_results_main$sig_results,
  full_masm_list = full_masm,
  full_cell_list = unique(celltype_prop_df$Cell_type),
  title = "Meta-analysis pooled correlation"
)
print(p_all)


df_epi <- meta_results_main$all_results %>%
    filter(Var1 == "Epithelial", Var2Type == "MASM") %>%
    mutate(MASM = sub(".*(MASM[1-7])", "\\1", Var2),
           MASM = factor(MASM, levels = paste0("MASM", 1:7)))

p_epi<-ggplot(df_epi, aes(x = MASM, y = pooled_rho)) +

    geom_errorbar(aes(ymin = pooled_rho_lower, ymax = pooled_rho_upper),
                  width = 0.22, linewidth = 0.6, color = "grey40") +
    geom_col(width = 0.6, fill = "#7FA6C9", color = "grey25", linewidth = 0.4) +

    geom_text(aes(label = "ns", y = pooled_rho_upper),
              vjust = -0.7, size = 3.8, color = "grey30", fontface = "bold")  +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = "MASM module",
         y = expression("Pooled Spearman's " * rho),
         title = "Correlation between MASM modules\nand epithelial proportion") +
    theme_hd_plain() +theme(axis.text.x = element_text(angle = 30,hjust=1))

palls11 <- ggarrange(p2qc,p3qc,p4qc,p_epi,nrow=1,ncol=4,labels = c(letters[4:7]),
                     font.label = panel_label_style,common.legend = T,legend="bottom")

palls12 <- ggarrange(p_exp,p_masm6_cor,p1qc,nrow = 1, ncol = 3,labels = c(letters[1:3]),
                     font.label = panel_label_style,widths = c(2,1,1))

palls13 <- ggarrange(p_stack_celltypes, p_list$MASM1,p_list$MASM2,p_list$MASM3,p_list$MASM4,p_list$MASM5 , nrow = 2, ncol = 3,labels = c(letters[8:13]),font.label = panel_label_style,common.legend = T,legend="bottom")

palls14<- ggarrange(palls12, palls11 ,palls13 , nrow = 3, ncol = 1,heights = c(1,1,2))


save_tiff(palls14, "supplementaryfigure11.tiff", width = 20, height = 20, path = fig_dir)




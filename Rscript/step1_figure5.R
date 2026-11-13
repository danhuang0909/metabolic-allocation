# Figure 5 and cross-cohort subtype, clinical, and drug analyses.
# Run this script from the project root with Rscript.

project_dir <- normalizePath(
  Sys.getenv("METABOLISM_PROJECT_DIR", unset = getwd()),
  mustWork = TRUE
)
source(file.path(project_dir, "Rscript", "step0_config.R"))
source(file.path(project_dir, "Rscript", "analysis_functions.R"))
setwd(data_dir)


# ==========================================================================

# Fig 7A: MASM states vs established CRC subtype systems (CMS / iCMS / PDS / SMI)
# Fig 7B: clinical associations (cross-cohort meta-Cox) + drug vulnerabilities
# Supplementary: MASM correlation matrices and supporting analyses
#

# - iCMS is now called with the iCMS.SSC package (DQ criterion, matching

# continuous scores are available as i2/i3/i2i3
# - MASM <-> subtype correlation matrices moved to the supplement
# - additional bulk cohorts: GSE39582, GSE248381, GSE156915, metaCRC,
# AC-ICAM (coad_silu_2022), GSE17538, GSE38832
# ==========================================================================

# ---- 1. Setup ---------------------------------------------------------------


# meta-analysis in section 6 needs lower_multi/upper_multi

library(tidyverse)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)
library(coxphf)
library(readxl)
library(meta)
library(PDSclassifier)
library(iCMS.SSC)   # devtools::install_github("CRCrepository/iCMS.SSC")
panel_label_style <- list(size = 18, color = "black", face = "bold")
full_masm <- paste0("Tumor_MASM", 1:7)

cms_colors   <- c(CMS1 = "#D55E00", CMS2 = "#0072B2", CMS3 = "#009E73", CMS4 = "#CC79A7")
surv_colors  <- c("#EE7621", "#00BFFF")
msi_colors   <- c("#00BFFF", "#EE7621")
cor_palette  <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
                                   "#F4A582", "#D6604D", "#B2182B"))(100)
masm_cols7   <- paste0("Tumor_MASM", 1:7)
# continuous subtype scores available in every run_gse_array_workflow plot_df:
# d.CMS1-4 (sign-flipped CMS distances), i2/i3 (iCMS.SSC DQ scores), PDS1-3, SMI
subtype_cols <- c(paste0("d.CMS", 1:4), "i2", "i3", paste0("PDS", 1:3), "SMI")


# iCMS classification with iCMS.SSC - DQ criterion, matching run_gse_array_workflow
# (which already calls iCMS.DQ internally). Use this only for matrices that do NOT
# go through the workflow (scRNA epithelial pseudobulk, DepMap cell lines).
# expr_log: genes x samples, log-transformed, gene symbols as rownames
# output columns: i2, i3, i2i3 (continuous scores), nearest.icms, confident.icms

run_icms <- function(expr_log) {
  set.seed(123)
  res <- iCMS.SSC::iCMS.DQ(as.matrix(expr_log))
  data.frame(patient = colnames(expr_log),
             i2 = res$i2, i3 = res$i3, i2i3 = res$i2i3,
             iCMS_nearest   = res$nearest.icms,
             iCMS_confident = res$confident.icms,
             stringsAsFactors = FALSE)
}

# iCMS.KNN expects log-transformed values; log2(x + 1) only if clearly linear
ensure_log <- function(mat) {

  if (max(mat, na.rm = TRUE) > 50) log2(mat + 1) else mat
}

# Collapse a probe-level matrix to gene level (mean over probes per gene)
collapse_probes_to_genes <- function(expr_mat, probe2gene) {
  # probe2gene: data.frame with columns ID and symbol (first two columns used)
  colnames(probe2gene)[1:2] <- c("ID", "symbol")
  probe2gene <- probe2gene %>%
    filter(symbol != "" & !is.na(symbol) & symbol != "---") %>%
    mutate(symbol = sub("\\s*///.*", "", symbol)) %>%
    distinct()
  expr_mat$ID <- rownames(expr_mat)
  anno <- merge(probe2gene, expr_mat, by = "ID")
  anno <- anno[!duplicated(anno[, c("ID", "symbol")]), ]
  collapsed <- aggregate(anno[, !colnames(anno) %in% c("ID", "symbol")],
                         by = list(symbol = anno$symbol), FUN = mean, na.rm = TRUE)
  mat <- as.data.frame(collapsed[, -1])
  rownames(mat) <- collapsed$symbol
  mat
}

# PCA biplot of MASM scores colored by CMS (shared by TCGA / metaCRC panels)
plot_masm_pca_biplot <- function(score_df, class_var, masm_cols = paste0("MASM", 1:7),
                                 shape_var = NULL, scale_factor = 5,
                                 icms_arrows = FALSE, title = NULL) {
  pca_res <- prcomp(score_df[, masm_cols], scale. = TRUE)
  sample_coord <- as.data.frame(pca_res$x[, 1:2]) %>%
    cbind(score_df) %>%
    filter(!is.na(.data[[class_var]]), .data[[class_var]] != "NOLBL")
  loading <- as.data.frame(pca_res$rotation[, 1:2])
  loading$feature <- gsub("Tumor_", "", rownames(loading))
  loading <- loading %>% mutate(PC1 = PC1 * scale_factor, PC2 = PC2 * scale_factor)

  point_aes <- if (!is.null(shape_var)) {
    aes(x = PC1, y = PC2, color = .data[[class_var]], shape = .data[[shape_var]])
  } else {
    aes(x = PC1, y = PC2, color = .data[[class_var]])
  }
  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.8) +
    geom_point(data = sample_coord, point_aes, size = 3) +
    geom_segment(data = loading,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_text(data = loading, aes(x = PC1, y = PC2, label = feature),
              size = 4, fontface = 2, nudge_y = 0.1, nudge_x = 0.04) +
    scale_color_manual(name = "CMS types", values = cms_colors) +
    labs(x = paste0("PC1 (", round(pca_res$sdev[1]^2 / sum(pca_res$sdev^2) * 100, 1), "%)"),
         y = paste0("PC2 (", round(pca_res$sdev[2]^2 / sum(pca_res$sdev^2) * 100, 1), "%)")) +
    theme_hd_plain() +
    theme(panel.grid = element_blank(), legend.position = "right")
  if (!is.null(shape_var)) p <- p + scale_shape_manual(values = c(MSI = 21, MSS = 19))
  if (icms_arrows) {
    p <- p +
      geom_segment(aes(x = -2.2, y = 4.5, xend = -0.1, yend = 4.5),
                   arrow = arrow(length = unit(0.15, "cm"), ends = "first"),
                   color = "black", linewidth = 0.3) +
      geom_segment(aes(x = 0.1, y = 4.5, xend = 2.2, yend = 4.5),
                   arrow = arrow(length = unit(0.15, "cm"), ends = "last"),
                   color = "black", linewidth = 0.3) +
      annotate("text", x = -2.4, y = 4.5, label = "iCMS3", size = 3.5, hjust = 1) +
      annotate("text", x = 2.4, y = 4.5, label = "iCMS2", size = 3.5, hjust = 0)
  }
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

# ---- 2. Single-cell data and MASM scores ------------------------------------
patient_infor <- readRDS("./single_cell_patient_infor_all_samples.rds")
meta_data <- readRDS(file.path(data_dir, "all_meta_data.rds"))
all_gene_pathway <- readRDS("./all_genes_pathway_infor_select.rds")
cpm_all <- readRDS("./sum_all_genes_count_rawRNA_cpm.rds")
kegg_path_genes <- readRDS("all_kegg_pathway_genes.rds")
gene_cluster <- readRDS("gene_cluster.rds")
all_loading_list <- readRDS("all_loading_list.rds")
pseudo_res <- readRDS(file.path(data_dir, "pseudo_seurat_all_sample.rds"))

# signalling-pathway genes excluding genes already assigned to metabolic pathways
kegg_signal <- kegg_path_genes[grepl("signal", kegg_path_genes$description), ]
cpm_all_select_signal <- cpm_all[, colnames(cpm_all) %in% kegg_signal$symbol &
                                   !colnames(cpm_all) %in% unique(all_gene_pathway$gene)]

pseudo_all_df_count <- t(pseudo_res$pseudo_all_raw_count)
groups <- str_split_fixed(rownames(pseudo_all_df_count), "_", 4)[, 2]
temp_pseudo <- getEdgeRdata(t(pseudo_all_df_count), groups)
cpm_pseudo <- temp_pseudo$count

cpm_pseudo_signal <- as.data.frame(
  cpm_pseudo[, colnames(cpm_pseudo) %in% kegg_signal$symbol &
               !colnames(cpm_pseudo) %in% unique(all_gene_pathway$gene)])
rn <- str_split_fixed(rownames(cpm_pseudo_signal), "_", 3)
rownames(cpm_pseudo_signal) <- paste(rn[, 1], rn[, 2], "pseudo", rn[, 3], sep = "_")

epi_all_tumor <- as.data.frame(pseudo_res$cpm_matrix[grepl("Tumor_Epi", rownames(pseudo_res$cpm_matrix)) , ])
rownames(epi_all_tumor)<-gsub("Tumor_","",rownames(epi_all_tumor))
fibro_all_tumor <- as.data.frame(pseudo_res$cpm_matrix[grepl("Tumor_Fib", rownames(pseudo_res$cpm_matrix)) , ])
rownames(fibro_all_tumor)<-gsub("Tumor_","",rownames(fibro_all_tumor))

epi_all_tumor_pseudo <- as.data.frame(cpm_pseudo_signal[grepl("Tumor_", rownames(cpm_pseudo_signal)) , ])
rownames(epi_all_tumor_pseudo)<-gsub("Tumor_","",rownames(epi_all_tumor_pseudo))


cell_types <- c("B cell", "Endothelial", "Myeloid", "Plasma", "T cell","Fibroblas")

# store results
masm_list <- list()
mat_list<-list()
for(ct in cell_types){
  # match rownames starting with cell type label, pseudo_res$cpm_matrix
  pp="Tumor_"
  pattern <- paste0(pp, ct, "")
  ct_mat <- as.data.frame(pseudo_res$cpm_matrix[grepl(pattern, rownames(pseudo_res$cpm_matrix)), ])
  # remove celltype prefix from rownames
  rownames(ct_mat) <- gsub(pp, "", rownames(ct_mat))
  # calculate MASM weight score
  ct_masm <- calc_masm_weight_score_simple(ct_mat, all_loading_list)
  masm_list[[ct]] <- ct_masm
  mat_list[[ct]]<-ct_mat
}

all_data_tumor <- calc_masm_weight_score_simple(epi_all_tumor, all_loading_list)
all_data_tumor_pseudo <- calc_masm_weight_score_simple(epi_all_tumor_pseudo, all_loading_list)
non_epi_cpm_mat<-readRDS( "./pseudo_bulk_output/non_epi_merged_cpm.rds")
no_epi_all_tumor_pseudo <- as.data.frame(non_epi_cpm_mat[grepl("Tumor_", rownames(non_epi_cpm_mat)) , ])
rownames(no_epi_all_tumor_pseudo)<-gsub("Tumor_","",rownames(no_epi_all_tumor_pseudo))

all_data_tumor_pseudono_epi <- calc_masm_weight_score_simple(no_epi_all_tumor_pseudo, all_loading_list)


# ---- 3. Subtype calls on the scRNA cohort (dedicated packages) --------------

# pseudo_res$cpm_matrix is log-CPM before running
epi_tumor_mat <- t(as.matrix(epi_all_tumor))
icms_sc <- run_icms(ensure_log(epi_tumor_mat))

icms_sc$patient <- str_split_fixed(icms_sc$patient, "_", 3)[, 3]

pds_calls_sc <- PDSpredict(data.frame(gene = rownames(epi_tumor_mat), epi_tumor_mat,
                                      check.names = FALSE),
                           species = "human", threshold = 0.6)
smi_data_sc <- calculateSMI(epi_tumor_mat, datatype = "bulk", species = "human")

pds_calls_sc$patient <- str_split_fixed(pds_calls_sc$Sample_ID, "_", 3)[, 3]
pds_calls_sc$SMI <- smi_data_sc[pds_calls_sc$Sample_ID, "SMI"]

subtype_sc <- icms_sc %>%
  full_join(pds_calls_sc %>% dplyr::select(patient, PDS = PDS_call, SMI), by = "patient")
# PDSpredict output columns verified: Sample_ID, PDS1, PDS2, PDS3, PDS_call

# ---- 4. TCGA-COREAD bulk cohort ---------------------------------------------
tcga_coread <- readRDS("TCGA_COREAD_STAR_tumor_CPM_genesymbol.rds")

tcga_re <- run_gse_array_workflow(expr_mat = tcga_coread, loading_list = all_loading_list,
                                  scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)
saveRDS(tcga_re, "TCGA_COREAD_STAR_classtypes.rds")
# iCMS already computed inside run_gse_array_workflow (iCMS.DQ):
# tcga_re$plot_df carries i2, i3, i2i3, nearest.icms, confident.icms
tcga_re<-readRDS("TCGA_COREAD_STAR_classtypes.rds")
pan_cancer_clinical <- read_excel("pan_cancer_TCGA_clinical_data_supplementary.xlsx")
colon_clinical <- as.data.frame(
  pan_cancer_clinical[pan_cancer_clinical$bcr_patient_barcode %in% colnames(tcga_coread), ])
clinical_molecular_public_all <- read.delim("./clinical_molecular_public_all.txt")
clinical_tcga <- clinical_molecular_public_all[
  clinical_molecular_public_all$sample %in% colnames(tcga_coread), ]
rownames(clinical_tcga) <- clinical_tcga$sample
clinical_tcga$patient <- clinical_tcga$sample

clinical_all <- tcga_re$plot_df %>%
  left_join(clinical_tcga, by = "patient") %>%
  left_join(colon_clinical, by = c("patient" = "bcr_patient_barcode"))


clinical_all$`TCGA iCMS2` <- clinical_all$i2
clinical_all$`TCGA iCMS3` <- clinical_all$i3

clinical_all$bcr_patient_barcode <- clinical_all$patient
clinical_all$gender <- clinical_all$gender.x
clinical_all$stage <- clinical_all$ajcc_pathologic_tumor_stage
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA")] <- "Stage I"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB")] <- "Stage IV"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("[Not Available]", "[Discrepancy]")] <- NA
clinical_all$stage2 <- ifelse(clinical_all$stage %in% c("Stage I", "Stage II"), "Early",
                              ifelse(clinical_all$stage %in% c("Stage III", "Stage IV"), "Advance", NA))
clinical_all$MSI <- toupper(clinical_all$msi)
clinical_all <- clinical_all %>%
  mutate(kras_mut_label = ifelse(kras_mut == 1, "Yes", "No"),
         braf_mut_label = ifelse(braf_mut == 1, "Yes", "No"),
         pt_lable = paste("pT", pt), pm_lable = paste("pM", pm), pn_lable = paste("pN", pn),
         iCMS_class = confident.icms)

# switch to nearest.icms if you want every sample classified

# ---- 5. Additional bulk cohorts ----------------------------------------------
# --- 5.1 GSE39582 (GPL570) ---
GSE39582_annotation_only <- read.delim(file.path(data_dir, "GSE39582_annotation_only.txt"), header = FALSE)
clinical_gse39582 <- parse_geo_anno(GSE39582_annotation_only)
GSE39582_series_matrix <- read.delim(file.path(data_dir, "GSE39582_series_matrix.txt"),
                                     row.names = 1, comment.char = "!")
GPL570_clean <- read.delim(file.path(data_dir, "GPL570_clean.txt"))
GSE39582_mat <- collapse_probes_to_genes(GSE39582_series_matrix,
                                         GPL570_clean %>% dplyr::select(ID, Gene.symbol))

GSE39582_re <- run_gse_array_workflow(GSE39582_mat, all_loading_list,
                                      scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)
saveRDS(GSE39582_re, "GSE39582_classtypes.rds")
GSE39582_re<-readRDS("GSE39582_classtypes.rds")
cms_labels_public_all <- read.delim(file.path(data_dir, "cms_labels_public_all.txt"))
gse39582_cms <- cms_labels_public_all[cms_labels_public_all$sample %in% clinical_gse39582$GSM, ]
rownames(gse39582_cms) <- gse39582_cms$sample
all_GSE39582 <- GSE39582_re$plot_df %>%
  left_join(clinical_gse39582, by = c("patient" = "GSM"))
all_GSE39582$cms_label <- gse39582_cms[all_GSE39582$patient, "CMS_final_network_plus_RFclassifier_in_nonconsensus_samples"]
all_GSE39582$stage <- all_GSE39582$tnm_stage
all_GSE39582$stage[all_GSE39582$tnm_stage %in% c(0, 1, 2)] <- "early"
all_GSE39582$stage[all_GSE39582$tnm_stage %in% c(3, 4)] <- "advance"

# --- 5.2 GSE248381 (GPL23985) ---
gpl_raw <- read.delim("GPL23985-23107.txt", skip = 7, sep = "\t", header = TRUE)
GSE248381_series_matrix <- read.delim(file.path(data_dir, "GSE248381_series_matrix.txt"),
                                      row.names = 1, comment.char = "!")
GSE248381_mat <- collapse_probes_to_genes(GSE248381_series_matrix,
                                          gpl_raw %>% dplyr::select(ID, Gene.Symbol))


GSE248381_re <- run_gse_array_workflow(GSE248381_mat, all_loading_list,
                                       scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)
saveRDS(GSE248381_re, "GSE248381_classtypes.rds")

# --- 5.3 GSE156915 (gene-level rows with "///" separators in rownames) ---
GSE156915_series_matrix <- read.delim(file.path(data_dir, "GSE156915_series_matrix.txt"),
                                      row.names = 1, comment.char = "!")
expr_gse_df <- as.data.frame(GSE156915_series_matrix)
expr_gse_df$gene <- sub("\\s*///.*", "", rownames(GSE156915_series_matrix))
expr_matGSE156915 <- expr_gse_df %>%
  filter(gene != "") %>%
  group_by(gene) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  as.data.frame()
rownames(expr_matGSE156915) <- expr_matGSE156915$gene
expr_matGSE156915 <- expr_matGSE156915[, -1]

GSE156915_re <- run_gse_array_workflow(expr_matGSE156915, all_loading_list,
                                       scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)
saveRDS(GSE156915_re, "GSE156915_classtypes.rds")
GSE156915_re<-readRDS("GSE156915_classtypes.rds")
# GSE156915 annotation from the GEO series matrix header
# grep '^!' GSE156915_series_matrix.txt > GSE156915_annotation_only.txt
anno <- read.delim("GSE156915_annotation_only.txt", FALSE)
colnames(anno) <- c("tag", paste0("S", seq(ncol(anno) - 1)))
df_anno <- anno %>%
  filter(tag %in% c("!Sample_geo_accession", "!Sample_characteristics_ch1")) %>%
  pivot_longer(-tag, names_to = "sid", values_to = "raw") %>%
  mutate(raw = str_remove_all(raw, '^"|"$')) %>% filter(raw != "")
gsm <- df_anno %>% filter(tag == "!Sample_geo_accession") %>% dplyr::select(sid, GSM = raw)
clin <- df_anno %>% filter(tag == "!Sample_characteristics_ch1") %>%
  separate(raw, c("v", "val"), sep = ":\\s*", extra = "merge") %>%
  mutate(across(c(v, val), str_squish)) %>% filter(v != "")
clin_gse156915 <- inner_join(gsm, clin, by = "sid") %>%
  group_by(GSM, v) %>% summarise(val = unique(val), .groups = "drop") %>%
  pivot_wider(names_from = v, values_from = val) %>%
  rename(DDIR_call = `dna damage immune response call`,
         DDR_score = `dna damage repair deficient score`,
         CMS = `consensus molecular subtype`,
         CRIS_subtype = `colorectal cancer intrinsic sub-type`,
         MSI = msi, Tissue_type = tissue, KRAS = kras, NRAS = nras,
         BRAF = braf, PIK3CA = pik3ca, TP53 = tp53)
all_GSE156915 <- inner_join(clin_gse156915, GSE156915_re$plot_df, by = c("GSM" = "patient"))

saveRDS(all_GSE156915 ,file="all_GSE156915.Rds")
# --- 5.4 metaCRC ---
load("CRC_metacohort_July_2024.RData")  # provides exprs.data
gene_by_all <- aggregate(x = exprs.data[, -c(1:3)], by = list(Gene.Symbol = exprs.data$symbol),
                         FUN = mean, na.rm = TRUE)
rownames(gene_by_all) <- gene_by_all$Gene.Symbol
gene_all_meta <- as.data.frame(gene_by_all[, -1])

metaCRC_re <- run_gse_array_workflow(gene_all_meta, all_loading_list,
                                     scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)
saveRDS(metaCRC_re, "metaCRC_classtypes.rds")
all_metaCRC <- metaCRC_re$plot_df

# --- 5.5 AC-ICAM (coad_silu_2022) ---
silu_colon <- readRDS("coad_silu_2022_data.rds")
silu_colon_re <- run_gse_array_workflow(silu_colon$data, all_loading_list,
                                        scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)
saveRDS(silu_colon_re, "silu_colon_classtypes.rds")

all_acicam <- silu_colon_re$plot_df %>%
  left_join(silu_colon$clinical, by = c("patient" = "SAMPLE_ID")) %>%
  transmute(patient, across(all_of(masm_cols7)),
            across(all_of(subtype_cols)),
            iCMS_nearest = nearest.icms, iCMS_confident = confident.icms,
            cms_label = ifelse(CMS %in% c("CMS1", "CMS2", "CMS3", "CMS4"), CMS, "NOLBL"),
            OS_status  = ifelse(OS_STATUS == "1:DECEASED", 1, 0),
            OS_month   = as.numeric(OS_MONTHS),
            PFS_status = ifelse(PFS_STATUS == "1:Recurred", 1, 0),
            PFS_month  = as.numeric(PFS_MONTHS),
            age_at_diagnosis = as.numeric(AGE_AT_DX),
            Sex = SEX, stage = factor(AJCC_PATH_STAGE)) %>%
  filter(!is.na(PFS_month), !is.na(PFS_status))

# --- 5.6 GSE17538 / GSE38832 ---
gse17538_colon <- readRDS("gse17538_data.rds")
gse38832_colon <- readRDS("gse38832_data.rds")
re_17538 <- run_gse_array_workflow(gse17538_colon$data, all_loading_list,
                                   scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)
re_38832 <- run_gse_array_workflow(gse38832_colon$data, all_loading_list,
                                   scale_gene = TRUE, scale_score = TRUE, pds_threshold = 0.6)


build_cohort_df <- function(re, clinical) {
  df <- as.data.frame(re$plot_df %>%
                        left_join(clinical, by = c("patient" = "gsm")))
  df$cms_label <- as.character(df$prediction)
  df$cms_label[is.na(df$cms_label)] <- "NOLBL"
  df
}
all_GSE17538 <- build_cohort_df(re_17538, gse17538_colon$clinical) %>%
  filter(!is.na(ajcc_stage)) %>%  # drop 6 adenomas
  mutate(OS_status = case_when(overall_event..death.from.any.cause. == "death" ~ 1,
                               overall_event..death.from.any.cause. == "no death" ~ 0,
                               TRUE ~ NA_real_),
         OS_month = as.numeric(overall.survival.follow.up.time),
         DSS_status = case_when(dss_event..disease.specific.survival..death.from.cancer. == "death" ~ 1,
                                dss_event..disease.specific.survival..death.from.cancer. == "no death" ~ 0,
                                TRUE ~ NA_real_),
         DSS_month = as.numeric(dss_time),
         DFS_status = case_when(dfs_event..disease.free.survival..cancer.recurrence. == "recurrence" ~ 1,
                                dfs_event..disease.free.survival..cancer.recurrence. == "no recurrence" ~ 0,
                                TRUE ~ NA_real_),
         DFS_month = as.numeric(dfs_time),
         Age = as.numeric(age), Sex = gender, stage = factor(ajcc_stage))
all_GSE38832 <- build_cohort_df(re_38832, gse38832_colon$clinical) %>%
  mutate(DSS_status = case_when(grepl("^1", dss_event..disease.specific.survival.) ~ 1,
                                grepl("^0", dss_event..disease.specific.survival.) ~ 0,
                                TRUE ~ NA_real_),
         DSS_month = as.numeric(dss_time..disease.specific.survival.time..months.),
         DFS_status = case_when(grepl("^1", dfs_event..disease.free.survival.) ~ 1,
                                grepl("^0", dfs_event..disease.free.survival.) ~ 0,
                                TRUE ~ NA_real_),
         DFS_month = as.numeric(dfs_time..disease.free.survival.time..months.),
         stage = factor(ajcc_stage))


#all_data_bulk<-rbind(cbind(dataset="TCGA",tcga_re$plot_df),cbind(dataset="GSE39582",GSE39582_re$plot_df)
#,cbind(dataset="GSE17538",re_17538$plot_df),cbind(dataset="GSE38832",re_38832$plot_df)
#,cbind(dataset="metaCRC",metaCRC_re$plot_df[!metaCRC_re$plot_df$patient #%in%c(GSE39582_re$plot_df$patient,re_17538$plot_df$patient,re_38832$plot_df$patient), #]),cbind(dataset="GSE156915",GSE156915_re$plot_df)
#,cbind(dataset="GSE248381",GSE248381_re$plot_df),cbind(dataset="silu_colon",silu_colon_re$plot_df))

all_data_bulk<-rbind(cbind(dataset="TCGA",tcga_re$plot_df)
                     ,cbind(dataset=pheno.data[metaCRC_re$plot_df$patient,"dataset"],metaCRC_re$plot_df),cbind(dataset="GSE156915",GSE156915_re$plot_df)
                     ,cbind(dataset="GSE248381",GSE248381_re$plot_df),cbind(dataset="silu_colon",silu_colon_re$plot_df))

all_data_bulk$iCMS_call=gsub("i","iCMS",all_data_bulk$confident.icms)
all_data_bulk$CMS=all_data_bulk$prediction
all_data_bulk$cms_label<-all_data_bulk$CMS

saveRDS(all_data_bulk,file="all_data_bulk.rds")
all_data_bulk<-readRDS(file="all_data_bulk.rds")

my_subtype_colors <- c(

    "CMS1" ="#E28700",
    "CMS2" = "#003380",
    "CMS3" = "#782F80",
    "CMS4" = "#59C141",


    "iCMS2"   = "#0072B2",
    "iCMS3"   = "#D55E00",
       "i2"   = "#0072B2",
    "i3"   = "#D55E00",
           "MSI"   = "#0072B2",
    "MSS"   = "#D55E00",

    "PDS1" = "#4088D0",
    "PDS2" = "#0F6B44",
    "PDS3" = "#E69F00",
    "Unclassified" = "gray80"
)


p_cms <- plot_pca_biplot(all_data_bulk,
                         color_col = "CMS",
                         title = "MASM distribution across CMS subtypes",
                         custom_colors = my_subtype_colors)


p_pds <- plot_pca_biplot(all_data_bulk,
                         color_col = "PDS_call",
                         title = "MASM distribution across PDS subtypes",
                         custom_colors = my_subtype_colors)

p_icms <- plot_pca_biplot(all_data_bulk,
                         color_col = "iCMS_call",
                         title = "MASM distribution across iCMS subtypes",
                         custom_colors = my_subtype_colors)


ggarrange(p_cms ,p_pds ,p_icms,nrow=1,ncol=3)


df <- masm_prepare(all_data_bulk)
ht_sys <- masm_system_heatmap(df)

ht_grob <- grid.grabExpr(
  draw(ht_sys$heatmap,
       column_title = "MASM dimensions are associated with\nmultiple molecular classification systems",
       column_title_gp = gpar(fontsize = 14), heatmap_legend_side = "bottom"),
  wrap = TRUE)

cms <- masm_scan(df, "CMS")
masm_report(cms)
 masm_plot_heatmap(cms, type = "pair")
  masm_plot_scatter(cms)


  pds <- masm_scan(df, "PDS_call")
  masm_report(pds)
  masm_plot_heatmap(pds, type = "both")

  icms <- masm_scan(df, "iCMS")
  masm_report(icms)
  masm_plot_scatter(icms)


p_scatter1<- plot_masm_scatter(df = df,
                         x_var = "Tumor_MASM4",
                         y_var = "Tumor_MASM7",
                          group_type = "CMS",
                             group_select = c("CMS1","CMS3"),auc_label_x =-1)
p_scatter1$plot


p_scatter2<- plot_masm_scatter(df = df,
                               x_var = "Tumor_MASM5",
                               y_var = "Tumor_MASM6",
                               group_type = "CMS",
                             group_select = c("CMS2","CMS4"),auc_label_x =-1)
p_scatter2$plot

p_scatter3<- plot_masm_scatter(df = df,
                               x_var = "Tumor_MASM1",
                               y_var = "Tumor_MASM6",
                               group_type = "CMS",
                             group_select = c("CMS4","CMS3"),auc_label_x =-1)
p_scatter3$plot


p_scatter4<- plot_masm_scatter(df = df,
                               x_var = "Tumor_MASM4",
                               y_var = "Tumor_MASM6",
                               group_type = "PDS_call",
                               group_select = c("PDS2","PDS3"),auc_label_x =-2)
p_scatter4$plot


res_lodo <- run_final_model_lodo(all_data_bulk)
df_plot <- res_lodo$per_cohort
df_summary <- res_lodo$summary


df_summary <- df_summary %>%
    arrange(pooled_acc) %>%
    mutate(model = fct_inorder(model))


df_plot <- df_plot %>%
    mutate(model = factor(model, levels = levels(df_summary$model)))


model_colors <- c(
    "iCMS call"                   = "#BDBDBD",
    "PDS call"                    = "#969696",
    "SMI"                         = "#737373",
    "SMI + PDS + iCMS"            = "#525252",
    "all 7 MASM"                  = "#2166AC",
    "iCMS + all MASM"             = "#FDB863",
    "PDS + all MASM"              = "#F4A261",
    "SMI + all MASM"              = "#E08214",
    "SMI + PDS + iCMS + all MASM" = "#B2182B"
)


p_lodo_high_end <- ggplot() +

    geom_line(data = df_plot,
              aes(x = model, y = accuracy, group = heldout),
              color = "grey90", alpha = 0.6, linewidth = 0.7) +


    geom_point(data = df_plot,
               aes(x = model, y = accuracy, color = model),
               alpha = 0.5, size = 2.5) +
    geom_point(data = df_summary,
               aes(x = model, y = pooled_acc, fill = model),
               shape = 23, size = 5.5, color = "white", stroke = 1.2) +

    geom_text(data = df_summary,
              aes(x = model, y = pooled_acc, label = sprintf("%.3f", pooled_acc)),
              vjust = -2.5, size = 3.5, fontface = "bold", color = "black") +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors) +
    scale_y_continuous(limits = c(0.35, 0.85), breaks = seq(0.4, 0.8, 0.1)) +
    labs(
        x = NULL,
        y = "Prediction accuracy of CMS subtypes",
        title = "MASM states improve CMS classification\nacross independent CRC cohorts" ) +
    theme_hd_plain() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
        legend.position = "none",

        panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),

    )


print(p_lodo_high_end)


p_cor<-ggplot(all_data_bulk,aes(x=Tumor_MASM7,y=SMI))+geom_point()+stat_cor(method = "spearman")+theme_hd_plain()

rho_df <- cor_one(all_data_bulk) %>%
  mutate(star = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                          p < 0.05 ~ "*", TRUE ~ ""),
         i = as.integer(gsub("MASM", "", MASM)))
rho_plot_df <- rho_df %>%
    mutate(

        direction = ifelse(rho > 0, "Positive", "Negative"),

        MASM = reorder(MASM, rho)
    )


p_diverging <- ggplot(rho_plot_df, aes(x = rho, y = MASM, fill = direction)) +

    geom_col(width = 0.65, alpha = 0.9) +


    geom_text(aes(label = sprintf("%.2f %s", rho, star),
                  hjust = ifelse(rho > 0, -0.15, 1.15)),
              size = 4.5,  color = "black") +


    scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +


    geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +


    labs(
        x = "Spearman Correlation (rho) with SMI",
        y = NULL,
        title = "MASM modules align with the\nSMI proliferation-maturation axis"
    ) +xlim(-0.65,0.4)+


    theme_hd_plain() +
    theme(
        legend.position = "none",

        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
    )


print(p_diverging)


g_within <- masm_gradient_heatmap(df)

expr_list_full<-readRDS("all_data_bulk_exp_data_expr_list_full.rds")
#
# TCGA_cell <- read_excel("NIHMS958212-supplement-2.xlsx")
# TCGA_cell_colon <- TCGA_cell[TCGA_cell$`TCGA Participant Barcode` %in% clinical_all$patient, ]
# cell_data <- TCGA_cell_colon %>%
# rename(patient = `TCGA Participant Barcode`) %>%
# mutate(across(-c(patient, `TCGA Study`, `Immune Subtype`, `TCGA Subtype`),
# ~ as.numeric(as.character(.))))
#
# merged_data <- clinical_all %>%
# rename_with(~ str_replace(., "Tumor_MASM", "MASM"), starts_with("Tumor_MASM")) %>%
# inner_join(cell_data %>% select(patient, Stromal_Fraction = `Stromal Fraction`,
# Lymphocytes, Macrophages), by = "patient") %>%
# mutate(score_CMS4 =  d.CMS4,
# is_CMS4 = case_when(cms_label == "CMS4" ~ 1,
# cms_label %in% c("CMS1", "CMS2", "CMS3") ~ 0,
# TRUE ~ NA_real_))
# merged_data_clean <- merged_data %>%
# filter(complete.cases(score_CMS4,  Stromal_Fraction, MSI,
# SMI,
# Macrophages, Lymphocytes, MASM6, MASM7))
# merged_data_clean$PDS_call[merged_data_clean$PDS_call=="Mixed"]=NA
# # continuous iCMS covariates (`TCGA iCMS2`/`TCGA iCMS3`) now come from iCMS.SSC
# # DQ scores i2/i3 (renamed above) - run_nested_cms4_lm works unchanged
# lm_result <- run_nested_cms4_lm(merged_data_clean)
# plot_lm_res <- plot_lm_coefficients(lm_result)
#
# res_fig <- plot_regression_forest_r2(
# lm_base = lm_result$Model4,
# lm_full = lm_result$Model5,
# var_levels = c("Lymphocytes", "Macrophages", "Stromal Fraction", "MASM6"),
# var_rename = c("Stromal_Fraction" = "Stromal Fraction")
# )
#
# res_fig$plot
# res_fig$anova_p
# res_fig$forest_df
#
#
# lm_result_sub <- run_nested_cms4_lm(merged_data_clean[which(merged_data_clean$prediction=="CMS4"),])
# plot_lm_res_sub <- plot_lm_coefficients(lm_result_sub)
#
# res_fig_sub <- plot_regression_forest_r2(
# lm_base = lm_result_sub$Model4,
# lm_full = lm_result_sub$Model5,
# var_levels = c("Lymphocytes", "Macrophages", "Stromal Fraction", "MASM6"),
# var_rename = c("Stromal_Fraction" = "Stromal Fraction")
# )
#
# res_fig_sub$plot
# res_fig_sub$anova_p
# res_fig_sub$forest_df
#
# p_lm<-ggarrange(res_fig_sub$p_forest,res_fig_sub$p_r2,nrow=1,ncol=2,widths = c(2,1))


 datasets<-unique(all_data_bulk$dataset)

 affy_cohorts<- datasets[ datasets %in% c("TCGA" ,"silu_colon")]
 all_data_bulk$score_CMS4<-all_data_bulk$d.CMS4

 all_data_bulk$dataset[all_data_bulk$dataset=="silu_colon"]="AC-ICAM"
 expr_list_full$`AC-ICAM`=expr_list_full$silu_colon


 res_cms4 <- run_cms4_stromal_meta(
 all_data_bulk, expr_list_full,
   id_col   = "patient",
   affy_cohorts = affy_cohorts,
   filter_col = "prediction", filter_val = "CMS4",
   min_n = 0,
   cache_file = "results/cms4within_stromal_meta_estimate.rds"
 )


 figs <- plot_cms4_stromal_forest( res_cms4 )
  p_lm_all <- ggarrange(figs$p_forest+theme(plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm")), figs$p_r2+theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.1), "cm"))+ggtitle(" \n   "),                            nrow = 1, ncol = 2, widths = c(2.5, 1))
 p_lm_all

 res_me <- run_masm_estimate_assoc(
 all_data_bulk,  expr_list_full,
   id_col = "patient",
  affy_cohorts =  affy_cohorts,  min_n = 0,
  cache_file = "results/masm_estimate_assoc.rds",
 per_sample_csv = "results/MASM_estimate_per_sample.csv"
 )
 p_me <- plot_masm_estimate_bar(res_me)


 d <- res_me$cor_by_cohort %>%
  filter(target == "stromal", MASM == "Tumor_MASM6") %>%
  arrange(rho) %>%
  mutate(dataset = factor(dataset, levels = dataset))


p_col <- intersect(c("padj", "p", "p.value", "pval"), colnames(d))
if (length(p_col)) {
  pv <- d[[p_col[1]]]
  d$star <- ifelse(is.na(pv), "",
            ifelse(pv < 1e-4, "****",
            ifelse(pv < 1e-3, "***",
            ifelse(pv < 0.01, "**",
            ifelse(pv < 0.05, "*", "ns")))))
} else {
  d$star <- "****"
  message("p , ****")
}


d$lab_y <- d$rho + 0.05 * sign(d$rho) * max(abs(d$rho), na.rm = TRUE)

p_cor_stromal <- ggplot(d, aes(dataset, rho, fill = rho > 0)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  geom_text(aes(y = lab_y, label = paste0(sprintf("%.2f", rho), "\n", star)),
            size = 3.4, lineheight = 0.9) +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.22))) +
  labs(title = "MASM6 correlation with stromal content across cohorts",
       x = NULL, y = "Spearman's rho") +
  theme_hd_minimal_plain() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

print(p_cor_stromal)


p_masm6 <- ggplot(all_data_bulk[all_data_bulk$prediction=="CMS4",], aes(x = Tumor_MASM6, y = d.CMS4)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.25) +
  stat_cor(vjust = 0.5,method = "spearman") +
  labs(x = "MASM6", y = "CMS4 proximity score")+ggtitle("Within CMS4 subtype")+ theme_hd_plain()

pall60 <- ggarrange(p_cms ,p_pds ,p_icms, ht_grob,p_scatter1$plot,p_scatter2$plot,p_scatter4$plot, p_diverging  , p_lodo_high_end,
                    labels = letters[1:9], nrow = 3, ncol = 3,
                    font.label = panel_label_style)
pall61 <- ggarrange(g_within$grob, p_lm_all,
                    labels = letters[10:11], nrow = 1, ncol = 2,
                    font.label = panel_label_style,widths = c(1,2))
pall6 <- ggarrange(pall60 ,pall61, nrow = 2, ncol = 1,
                    font.label = panel_label_style,heights = c(3,1))

save_tiff(pall6, "figure5.tiff", width = 19, height = 22, path = fig_dir)


all_data_combined <- bind_rows(all_data_tumor %>% mutate(src="Epi"),
                                           all_data_tumor_pseudo %>% mutate(src="pseudo"))
res_epi <- run_masm_meta_pseudo(all_data_combined)
all_data_combined_non_epi <- bind_rows(all_data_tumor_pseudono_epi %>% mutate(src="NonEpi"),
all_data_tumor_pseudo %>% mutate(src="pseudo"))


resnon_epi <- run_masm_meta_pseudo(all_data_combined_non_epi )


p_tumor_bulk<-ggplot(res_epi$meta, aes(x = module, y = rho_pool)) +

    geom_errorbar(aes(ymin = ci_l, ymax = ci_u),
                  width = 0.22, linewidth = 0.6, color = "grey40") +
    geom_col(width = 0.6, fill = "#4C8C6A", color = "grey25", linewidth = 0.4) +

    geom_text(aes(label = star, y = ci_u),
              vjust = -0.6, size = 4.2, fontface = "bold", color = "grey25") +

    geom_text(aes(label = sprintf("%.2f", rho_pool)),
              vjust = 1.4, size = 3.1, color = "white", fontface = "bold") +
    scale_y_continuous(limits = c(0, 1.05),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "MASM module",
         y = expression("Pooled Spearman's " * rho),
         title = "Epithelial vs pseudo-bulk MASM score") +
  theme_hd_plain()+theme(axis.text.x = element_text(angle = 30,hjust=1))

p_tumornon_epi_bulk<-ggplot(resnon_epi$meta, aes(x = module, y = rho_pool)) +

    geom_errorbar(aes(ymin = ci_l, ymax = ci_u),
                  width = 0.22, linewidth = 0.6, color = "grey40") +
    geom_col(width = 0.6, fill = "#4C8C6A", color = "grey25", linewidth = 0.4) +

    geom_text(aes(label = star, y = ci_u),
              vjust = -0.6, size = 4.2, fontface = "bold", color = "grey25") +

    geom_text(aes(label = sprintf("%.2f", rho_pool)),
              vjust = 1.4, size = 3.1, color = "white", fontface = "bold") +
    scale_y_continuous(limits = c(0, 1.05),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "MASM module",
         y = expression("Pooled Spearman's " * rho),
         title = "Non-epithelial vs pseudo-bulk MASM score") +
  theme_hd_plain()+theme(axis.text.x = element_text(angle = 30,hjust=1))


p_cms_tcga <- plot_pca_biplot(all_data_bulk[all_data_bulk$dataset=="TCGA",],
                          color_col = "CMS",
                          title = "MASM distribution across CMS subtypes\n(RNA-seq bulk)",
                          custom_colors = my_subtype_colors,feature_cols = paste0("Tumor_MASM", c(1:7)))

p_cms_gse <- plot_pca_biplot(all_data_bulk[all_data_bulk$dataset=="GSE39582",],
                          color_col = "CMS",
                          title = "MASM distribution across CMS subtypes\n(Microarray bulk)",
                          custom_colors = my_subtype_colors,feature_cols = paste0("Tumor_MASM", c(1:7)), flip_pc_x =T, flip_pc_y =T)

p_masm4 <- ggplot(all_data_bulk[all_data_bulk$prediction=="CMS2",], aes(x = Tumor_MASM4, y = d.CMS2)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.25) +
  stat_cor(vjust = 0.5,method = "spearman") +
  labs(x = "MASM4", y = "CMS2 proximity score")+ggtitle("Within CMS2 subtype")+ theme_hd_plain()

p_masm1 <- ggplot(all_data_bulk[all_data_bulk$prediction=="CMS3",], aes(x = Tumor_MASM1, y = d.CMS3)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.25) +
  stat_cor(vjust = 0.5,method = "spearman") +
  labs(x = "MASM1", y = "CMS3 proximity score")+ggtitle("Within CMS3 subtype")+ theme_hd_plain()

p_masm5 <- ggplot(all_data_bulk[all_data_bulk$prediction=="CMS4",], aes(x = Tumor_MASM5, y = d.CMS4)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.25) +
  stat_cor(vjust = 0.5,method = "spearman") +
  labs(x = "MASM5", y = "CMS4 proximity score")+ggtitle("Within CMS4 subtype")+ theme_hd_plain()

p_masm7 <- ggplot(all_data_bulk[all_data_bulk$prediction=="CMS4",], aes(x = Tumor_MASM7, y = d.CMS4)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.25) +
  stat_cor(vjust = 0.5,method = "spearman") +
  labs(x = "MASM7", y = "CMS4 proximity score")+ggtitle("Within CMS4 subtype")+ theme_hd_plain()


 #coh_all <- run_module_coherence_all_cohorts(expr_list_full, all_loading_list)
coh_all<-readRDS(file="coh_all.rds")
plot_cohort_all<-coh_all$plot+theme_hd_plain()+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position="right")

plot_cohort_all

pall60s <- ggarrange(p_tumor_bulk,plot_cohort_all,
                    labels = letters[c(1,3)], nrow = 1, ncol = 2,
                    font.label = panel_label_style,widths = c(1,2))
pall61s <- ggarrange(p_tumornon_epi_bulk,p_cms_tcga ,p_cms_gse ,
                    labels = letters[c(2,4:5)], nrow = 1, ncol = 3,
                    font.label = panel_label_style,common.legend = T,legend="bottom")
pall62s <- ggarrange(p_scatter3$plot ,p_masm1,p_masm4,p_masm5,p_masm7,p_cor_stromal,
                    labels = letters[c(6:11)], nrow = 2, ncol = 3,
                    font.label = panel_label_style)
pall63s<-ggarrange(pall60s,pall61s,pall62s,nrow=3,ncol=1,heights = c(1,1,2),labels = c(letters[1],"",""), font.label = panel_label_style)

save_tiff(pall63s, "supplementaryfigure12.tiff", width = 18, height = 18, path = fig_dir)


library(tidyverse)


sample_by_dataset <- split(rownames(pheno.data), pheno.data$dataset)


expr_list <- lapply(sample_by_dataset, function(samp_vec){

  keep_samp <- intersect(samp_vec, colnames(gene_all_meta))
  gene_all_meta[, keep_samp, drop = FALSE]
})


names(expr_list)

expr_list_full <- c(list(
    "TCGA"     = tcga_coread,
    "silu_colon"  = silu_colon$data,
    "GSE248381" = GSE248381_mat,
    "GSE156915" = expr_matGSE156915
),expr_list)
cellaudit_full <- run_celltype_audit(expr_list_full, mcp_genes)

saveRDS(expr_list_full,file="all_data_bulk_exp_data_expr_list_full.rds")
expr_list_full<-readRDS("all_data_bulk_exp_data_expr_list_full.rds")

est_audit <- run_estimate_audit(expr_list_full,
                                affy_cohorts = c("GSE39582","GSE17538","GSE38832",
                                                 "metaCRC","GSE156915","GSE248381"))


# ---- 6. Per-cohort Cox models + cross-cohort meta-analysis -------------------
cols_tcga     <- c(Age = "age_at_initial_pathologic_diagnosis", Sex = "gender", Stage = "stage")
cols_gse39582 <- c(Age = "age_at_diagnosis", Sex = "Sex", Stage = "stage")
cols_acicam   <- c(Age = "age_at_diagnosis", Sex = "Sex", Stage = "stage")
cols_17538    <- c(Age = "Age", Sex = "Sex", Stage = "stage")
cols_38832    <- c(Stage = "stage")  # no age/sex in this cohort

cox_os <- list(
  "TCGA-COREAD" = batch_masm_cms_cox(clinical_all, "OS.time", "OS", covar_list = cols_tcga),
  "AC-ICAM"     = batch_masm_cms_cox(all_acicam, "OS_month", "OS_status", covar_list = cols_acicam),
  "GSE17538"    = batch_masm_cms_cox(all_GSE17538 %>% filter(!is.na(OS_status), !is.na(OS_month)),
                                     "OS_month", "OS_status", covar_list = cols_17538),
  "GSE39582"    = batch_masm_cms_cox(all_GSE39582, "OS_month", "OS_status", covar_list = cols_gse39582)
)
cox_dss <- list(
  "TCGA-COREAD" = batch_masm_cms_cox(clinical_all, "DSS.time", "DSS", covar_list = cols_tcga),
  "GSE17538"    = batch_masm_cms_cox(all_GSE17538 %>% filter(!is.na(DSS_status), !is.na(DSS_month)),
                                     "DSS_month", "DSS_status", covar_list = cols_17538),
  "GSE38832"    = batch_masm_cms_cox(all_GSE38832 %>% filter(!is.na(DSS_status), !is.na(DSS_month)),
                                     "DSS_month", "DSS_status", covar_list = cols_38832)
)


within_all_os <- map_dfr(paste0("CMS", 1:4), function(st) {
map_dfr(names(cohort_cfgs_os), function(nm) {
cfg <- cohort_cfgs_os[[nm]]
cox_within_subtype(cfg$df, subtype = st,
time_col = cfg$time, status_col = cfg$status,
cohort = nm, endpoint = paste0("OS|", st),
covars = cfg$covars)
})
})

within_all_dss <- map_dfr(paste0("CMS", 1:4), function(st) {
map_dfr(names(cohort_cfgs_dss), function(nm) {
cfg <- cohort_cfgs_dss[[nm]]
cox_within_subtype(cfg$df, subtype = st,
time_col = cfg$time, status_col = cfg$status,
cohort = nm, endpoint = paste0("DSS|", st),
covars = cfg$covars)
})
})

within_dss$res_long %>% filter(MASM == "Tumor_MASM1") %>%
select(cohort, endpoint, HR_uni, p_uni, HR_full, p_full, n, events) %>%
as.data.frame()


p_d1 <- panel_forest_within(within_os, within_os$meta,
                           "Tumor_MASM1", "CMS4", "full")

p_d2 <- panel_forest_within(within_os, within_os$meta,
                           "Tumor_MASM6", "CMS2", "full")
p_d3 <- panel_forest_within(within_os, within_os$meta,
                           "Tumor_MASM7", "CMS4", "full")
p_d4 <- panel_forest_within(within_os, within_os$meta,
                           "Tumor_MASM4", "CMS3", "full")

ggarrange(grid.draw(ht_within_os),p_d1,p_d2,p_d3,p_d4,nrow=2,ncol=3)


   within_pds  <- run_within_class(cohort_cfgs_os, "OS",
                                   class_col = "PDS_call",
                                   subtypes  = c("PDS1","PDS2","PDS3"))
   within_icms <- run_within_class(cohort_cfgs_os, "OS",
                                   class_col = "iCMS_call",
                                   subtypes  = c("iCMS2","iCMS3"))


# add its PFS model as a separate endpoint list if needed


# Pool multivariable log-HRs per MASM x CMS stratum across cohorts (REML + HK)

# functions file (additive change) - without them the SE cannot be derived
pool_cox_meta <- function(res_list) {
  bind_rows(res_list, .id = "cohort") %>%
    mutate(
      logHR = log(HR_multi),
      se    = (log(upper_multi) - log(lower_multi)) / (2 * 1.96)
    ) %>%
    filter(is.finite(logHR), is.finite(se), se > 0) %>%
    group_by(MASM_module, CMS_subtype) %>%
    filter(n() >= 2) %>%
    summarise(
      k = n(),
      meta_res = list(metagen(logHR, se, method.tau = "REML", hakn = TRUE, data = cur_data())),
      HR_pool  = exp(map_dbl(meta_res, ~ .x$TE.random)),
      P_pool   = map_dbl(meta_res, ~ .x$pval.random),
      .groups = "drop"
    ) %>%
    select(-meta_res)
}

meta_os  <- pool_cox_meta(cox_os)
meta_dss <- pool_cox_meta(cox_dss)

# Forest plot of the pooled HRs per MASM (one panel per endpoint)
plot_meta_cox_forest <- function(meta_df, title) {
  ggplot(meta_df, aes(x = HR_pool, y = MASM_module, color = CMS_subtype)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    scale_x_log10() +
    scale_color_manual(values = cms_colors) +
    labs(x = "Pooled multivariable HR (log scale)", y = NULL, title = title) +
    theme_hd_plain()
}
p_meta_os  <- plot_meta_cox_forest(meta_os,  "OS - cross-cohort meta-analysis")
p_meta_dss <- plot_meta_cox_forest(meta_dss, "DSS - cross-cohort meta-analysis")

# ---- 7. FIGURE 7A panels: MASM vs established subtypes -----------------------
# a) PCA biplot of MASM scores (TCGA-COREAD), CMS color + MSI shape
pca_biplot <- plot_masm_pca_biplot(
  clinical_all %>% mutate(MASM1 = Tumor_MASM1, MASM2 = Tumor_MASM2, MASM3 = Tumor_MASM3,
                          MASM4 = Tumor_MASM4, MASM5 = Tumor_MASM5, MASM6 = Tumor_MASM6,
                          MASM7 = Tumor_MASM7),
  class_var = "cms_label", shape_var = "MSI", icms_arrows = TRUE,
  title = "TCGA-COREAD cohort")

# b) cross-cohort meta-correlation: MASM1-7 vs subtype scores (CMS/PDS/SMI)
cor_one_cohort <- function(plot_df) {
  map_dfr(masm_cols7, function(m) {
    map_dfr(subtype_cols, function(s) {
      ok <- complete.cases(plot_df[[m]], plot_df[[s]])
      tibble(MASM = m, score = s, n = sum(ok),
             rho = if (sum(ok) >= 10) cor(plot_df[[m]][ok], plot_df[[s]][ok],
                                          method = "spearman") else NA_real_)
    })
  })
}
cor_by_cohort <- list(
  "TCGA-COREAD" = cor_one_cohort(clinical_all),
  "GSE39582"    = cor_one_cohort(all_GSE39582),
  "GSE248381"   = cor_one_cohort(GSE248381_re$plot_df),
  "GSE156915"   = cor_one_cohort(GSE156915_re$plot_df),
  "metaCRC"     = cor_one_cohort(all_metaCRC),
  "AC-ICAM"     = cor_one_cohort(silu_colon_re$plot_df),
  "GSE17538"    = cor_one_cohort(all_GSE17538),
  "GSE38832"    = cor_one_cohort(all_GSE38832)
)

meta_cor <- bind_rows(cor_by_cohort, .id = "cohort") %>%
  filter(!is.na(rho)) %>%
  group_by(MASM, score) %>%
  filter(n() >= 2) %>%
  summarise(mc = list(metacor(rho, n, sm = "ZCOR")), .groups = "drop") %>%
  mutate(rho_pool = map_dbl(mc, "cor.random"),
         p_pool   = map_dbl(mc, "pval.random"),
         k        = map_int(mc, "k"),
         sig      = case_when(p.adjust(p_pool, "BH") < 0.001 ~ "***",
                              p.adjust(p_pool, "BH") < 0.01  ~ "**",
                              p.adjust(p_pool, "BH") < 0.05  ~ "*",
                              TRUE ~ "")) %>%
  select(-mc)

cor_mat_meta <- meta_cor %>% select(MASM, score, rho_pool) %>%
  pivot_wider(names_from = score, values_from = rho_pool) %>%
  column_to_rownames("MASM") %>% as.matrix()
sig_mat_meta <- meta_cor %>% select(MASM, score, sig) %>%
  pivot_wider(names_from = score, values_from = sig) %>%
  column_to_rownames("MASM") %>% as.matrix()

ht_meta_cor <- Heatmap(
  cor_mat_meta[, subtype_cols], name = "Pooled rho",
  col = colorRamp2(c(-0.6, 0, 0.6), c("#2166AC", "white", "#B2182B")),
  cell_fun = function(j, i, x, y, w, h, f) {
    if (sig_mat_meta[i, subtype_cols][i, j] != "")
      grid.text(sig_mat_meta[i, subtype_cols][i, j], x, y,
                gp = gpar(fontsize = 9, fontface = "bold"))
  },
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_title = "MASM", column_title = "Subtype scores (meta across bulk cohorts)",
  rect_gp = gpar(col = "white", lwd = 0.5))
ht_meta_cor_g <- grid.grabExpr(draw(ht_meta_cor))

# c) MASM scores by package-based iCMS class (TCGA-COREAD)
p_icms_box <- plot_masm_box(clinical_all %>% filter(!is.na(iCMS_class)),
                            group_var = "iCMS_class", masm_cols = masm_cols7) +
  ggtitle("TCGA-COREAD cohort") + xlab("iCMS (iCMS.SSC, confident calls)")

# d) LASSO multinomial classification of CMS from MASM scores

lasso_result <- run_multinomial_lasso_cv(clinical_all)

# e/f) CMS4 score associations (TCGA cell-deconvolution covariates)
TCGA_cell <- read_excel("NIHMS958212-supplement-2.xlsx")
TCGA_cell_colon <- TCGA_cell[TCGA_cell$`TCGA Participant Barcode` %in% clinical_all$patient, ]
cell_data <- TCGA_cell_colon %>%
  rename(patient = `TCGA Participant Barcode`) %>%
  mutate(across(-c(patient, `TCGA Study`, `Immune Subtype`, `TCGA Subtype`),
                ~ as.numeric(as.character(.))))

merged_data <- clinical_all %>%
  rename_with(~ str_replace(., "Tumor_MASM", "MASM"), starts_with("Tumor_MASM")) %>%
  inner_join(cell_data %>% select(patient, Stromal_Fraction = `Stromal Fraction`,
                                  Lymphocytes, Macrophages), by = "patient") %>%
  mutate(score_CMS4 = -1 * d.CMS4,
         is_CMS4 = case_when(cms_label == "CMS4" ~ 1,
                             cms_label %in% c("CMS1", "CMS2", "CMS3") ~ 0,
                             TRUE ~ NA_real_))
merged_data_clean <- merged_data %>%
  filter(complete.cases(score_CMS4,  Stromal_Fraction, MSI,
                        SMI,
                        Macrophages, Lymphocytes, MASM6, MASM7))
merged_data_clean$PDS_call[merged_data_clean$PDS_call=="Mixed"]=NA
# continuous iCMS covariates (`TCGA iCMS2`/`TCGA iCMS3`) now come from iCMS.SSC
# DQ scores i2/i3 (renamed above) - run_nested_cms4_lm works unchanged
lm_result <- run_nested_cms4_lm(merged_data_clean)
plot_lm_res <- plot_lm_coefficients(lm_result)

p_masm6 <- ggplot(merged_data_clean, aes(x = MASM6, y = score_CMS4)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.25) +
  stat_cor(vjust = 0.5) +
  labs(x = "MASM6", y = "CMS4 score") + theme_hd_plain()

p_icms_cms4 <- ggplot(merged_data_clean, aes(x = `TCGA iCMS3`, y = score_CMS4)) +
  geom_point(alpha = 0.55, size = 1.7, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#e74c3c", alpha = 0.25) +
  stat_cor(vjust = 0.5) +
  labs(x = "iCMS3", y = "CMS4 score") + theme_hd_plain()

# g) CMS x iCMS concordance (package calls vs CMS labels)
concord <- clinical_all %>%
  filter(!is.na(iCMS_class), !is.na(cms_label), cms_label != "NOLBL") %>%
  count(cms_label, iCMS_class) %>%
  group_by(cms_label) %>% mutate(frac = n / sum(n)) %>% ungroup()
p_concord <- ggplot(concord, aes(x = cms_label, y = frac, fill = iCMS_class)) +
  geom_col(position = "fill", color = "white") +
  scale_fill_manual(values = c(i2 = "#D55E00", i3 = "#0072B2"), na.value = "grey80") +
  labs(x = "CMS", y = "Fraction", fill = "iCMS") + theme_hd_plain()
# iCMS.SSC class labels are "i2"/"i3" (verified against the package source)

# assemble Figure 7A
p7a_row1 <- ggarrange(pca_biplot, ht_meta_cor_g, nrow = 1, ncol = 2,
                      widths = c(1.2, 1), labels = letters[1:2],
                      font.label = list(size = 18, color = "black", face = "bold"))
p7a_row2 <- ggarrange(p_icms_box, p_concord, lasso_result$plot_coefficient,
                      nrow = 1, ncol = 3, labels = letters[3:5],
                      font.label = list(size = 18, color = "black", face = "bold"))
p7a_row3 <- ggarrange(p_masm6, p_icms_cms4, plot_lm_res$plot, nrow = 1, ncol = 3,
                      labels = letters[6:8],
                      font.label = list(size = 18, color = "black", face = "bold"))
figure7A <- ggarrange(p7a_row1, p7a_row2, p7a_row3, nrow = 3, ncol = 1)
save_tiff(figure7A, "figure7A_subtypes.tiff", width = 20, height = 18)

# ---- 8. FIGURE 7B panels: clinical associations ------------------------------
# a) KM curve, MASM3 DFI (median split, TCGA-COREAD)
plot_km <- function(df, feature, time_col, status_col, title,
                    lower_q = NULL, upper_q = NULL) {
  if (is.null(lower_q)) lower_q <- median(df[[feature]], na.rm = TRUE)
  if (is.null(upper_q)) upper_q <- lower_q
  d <- get_surveRaw2(df, feature, lower_q, upper_q, status_col, time_col)$df
  d$Log2FC <- d$FC
  fit <- survfit(Surv(time, vital_status) ~ Log2FC, data = d)
  ggsurvplot(fit, pval = TRUE)$plot +
    ggtitle(title) +
    scale_color_manual(name = gsub("Tumor_", "", feature), values = surv_colors,
                       labels = c("High", "Low")) +
    theme_hd_plain() + theme_hd_plain2()
}
p_surv_masm3 <- plot_km(clinical_all, "Tumor_MASM3", "DFI.time", "DFI",
                        "DFI (disease-free interval)")

# b) multivariate Cox forest, TCGA-COREAD DFI (Firth correction)

# plot_cox_gradient_surv(full_covariate = TRUE)
covars_tcga_full <- c(setNames(masm_cols7, gsub("Tumor_", "", masm_cols7)),
                      iCMS2 = "TCGA iCMS2", iCMS3 = "TCGA iCMS3",
                      Age = "age_at_initial_pathologic_diagnosis", Sex = "gender",
                      Stage = "stage", MSI = "msi", CMS = "cms_label", CIMP = "cimp")
p_dfi <- PlotCox(clinical_all, "DFI.time", "DFI", covariates = covars_tcga_full, use_firth = TRUE)
p_pfi <- PlotCox(clinical_all, "PFI.time", "PFI", covariates = covars_tcga_full, use_firth = FALSE)
p_os  <- PlotCox(clinical_all, "OS.time",  "OS",  covariates = covars_tcga_full, use_firth = FALSE)


res_grad_cms4 <- plot_cox_gradient_surv(clinical_all[clinical_all$cms_label == "CMS4", ],
                                        "Tumor_MASM1", "PFI.time", "PFI", n_level = 8)

res_grad_cms4_os <- plot_cox_gradient_surv(clinical_all[clinical_all$cms_label == "CMS4", ],
                                           "Tumor_MASM1", "OS.time", "OS", n_level = 8)
res_grad_cms1 <- plot_cox_gradient_surv(clinical_all[clinical_all$cms_label == "CMS1", ],
                                        "Tumor_MASM6", "PFI.time", "PFI", n_level = 8)
res_grad_cms1_os <- plot_cox_gradient_surv(clinical_all[clinical_all$cms_label == "CMS1", ],
                                           "Tumor_MASM6", "OS.time", "OS", n_level = 8)


# ---- 9. FIGURE 7B panels: drug vulnerabilities (DepMap cell lines) -----------
Model <- read.csv("./drugData/Model.csv")
select_colon_sample <- Model[Model$OncotreePrimaryDisease == "Colorectal Adenocarcinoma" &
                               Model$PrimaryOrMetastasis == "Primary", ]
exp <- read.csv("./drugData/Expression_Public_24Q2_subsetted.csv")
rownames(exp) <- exp[, 1]
exp_matrix <- apply(exp[, -1], 1, as.numeric)
rownames(exp_matrix) <- colnames(exp[, -1])
# orientation verified against getDrugResultnew/getDrugResultnewGene:
# rows = genes, columns = cell lines (ModelID)
exp_matrix <- exp_matrix[, colnames(exp_matrix) %in% select_colon_sample$ModelID]

drugall_prism <- read.csv("./drugData/Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted_NAsdropped.csv")
drugall_ctd   <- read.csv("./drugData/Drug_sensitivity_AUC_(CTD^2)_subsetted.csv")
drug1 <- read.csv("./drugData/Drug_sensitivity_AUC_(Sanger_GDSC1)_subsetted.csv")
drug2 <- read.csv("./drugData/Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted_NAsdropped.csv")
drugall_gdsc <- merge(drug1, drug2, by = "X")

Model_select <- Model[Model$ModelID %in% colnames(exp_matrix), ]
rownames(Model_select) <- Model_select$ModelID
Model_select$patient <- Model_select$ModelID

masm_score_res_ccl <- as.data.frame(calc_masm_weight_score_simple(
  as.data.frame(t(exp_matrix)), all_loading_list, scale_gene = T, scale_score = T))
temp_ccl <- as.data.frame(masm_score_res_ccl[, -1])
rownames(temp_ccl) <- as.character(masm_score_res_ccl[, 1])
sc_ccle=temp_ccl
temp_ccl$MSI <- Model_select[rownames(temp_ccl), "LegacyMolecularSubtype"]

temp_ccl$MSI[temp_ccl$MSI == ""] <- "MSS"

# DepMap expression is already log2(TPM+1); iCMS.SSC provides continuous i2/i3

icms_ccl <- run_icms(exp_matrix)
idx_ccl <- match(rownames(temp_ccl), icms_ccl$patient)
temp_ccl$iCMS2 <- icms_ccl$i2[idx_ccl]
temp_ccl$iCMS3 <- icms_ccl$i3[idx_ccl]
temp_ccl$iCMS_call<-gsub("i","iCMS",icms_ccl$iCMS_confident[idx_ccl])
#temp_ccl$iCMS_class <- icms_ccl$iCMS_confident[idx_ccl]

gdsc_drug  <- getDrugResultnew(drugall_gdsc, exp_matrix, temp_ccl)
ctd_drug   <- getDrugResultnew(drugall_ctd, exp_matrix, temp_ccl)
prism_drug <- getDrugResultnew(drugall_prism, exp_matrix, temp_ccl)

all_drug_cor <- rbind(cbind(gdsc_drug$all_cor, database = "GDSC"),
                      cbind(ctd_drug$all_cor, database = "CTD^2"),
                      cbind(prism_drug$all_cor, database = "PRISM"))
all_drug_cor$drug_names <- str_split_fixed(all_drug_cor$drug, "\\.\\.", 3)[, 1]
merged_df <- preprocessDrugCor(all_drug_cor)
plot_df <- filterSigPairs(merged_df, min_db = 2)

target_drugs <- c("SEPANTRONIUM.BROMIDE", "PEVONEDISTAT", "AZ.960", "RUXOLITINIB",
                  "EXELBINE", "DOCETAXEL", "DINACICLIB", "APO866")
main_anno <- c("SEPANTRONIUM.BROMIDE" = "Sepantronium\n(Survivin/BIRC5)",
               "PEVONEDISTAT"  = "Pevonedistat\n(NAE/NEDD8i)",
               "AZ.960"        = "AZ.960\n(JAK2i)",
               "RUXOLITINIB"   = "Ruxolitinib\n(JAK1/2i)",
               "EXELBINE"      = "Exelbine\n(Microtubule)",
               "DOCETAXEL"     = "Docetaxel\n(Microtubule)",
               "DINACICLIB"    = "Dinaciclib\n(CDKi)",
               "APO866"        = "APO866\n(NAMPTi)")
p_main <- plot_drug_dot(plot_df, target_drugs, main_anno,
                        title = "MASM-specific drug vulnerabilities",
                        filename = "Main_Drug_Dotplot.pdf", width = 6, height = 4.5)

# gene-level drug correlations (supplement)
# samples x genes, as getDrugResultnewGene expects cell lines in rows
exp_sub <- t(exp_matrix[rownames(exp_matrix) %in% names(gene_cluster), ])
res_gene  <- getDrugResultnewGene(drugall = drugall_prism, exp_matrix = exp_sub)
res_gene2 <- getDrugResultnewGene(drugall = drugall_ctd, exp_matrix = exp_sub)
res_gene3 <- getDrugResultnewGene(drugall = drugall_gdsc, exp_matrix = exp_sub)
all_drug_cor_gene <- rbind(cbind(res_gene3$all_cor, database = "GDSC"),
                           cbind(res_gene2$all_cor, database = "CTD^2"),
                           cbind(res_gene$all_cor, database = "PRISM"))
all_drug_cor_gene$drug_names <- str_split_fixed(all_drug_cor_gene$drug, "\\.\\.", 3)[, 1]
all_drug_cor_gene$type <- all_drug_cor_gene$gene
plot_df_gene <- filterSigPairs(preprocessDrugCor(all_drug_cor_gene), min_db = 2)
plot_df_gene$module <- paste("MASM", gene_cluster[plot_df_gene$type], sep = "")

# drug meta forests (main figure)
plot_drug_forest <- function(df, drug, module, title = NULL) {
  draw_drug_single_meta_forest(df, target_drug = drug, target_module = module,
                               title = ifelse(is.null(title), paste(drug, "vs", module), title))$plot
}
pdrug1 <- plot_drug_forest(plot_df, "PEVONEDISTAT", "MASM2")
pdrug4=plot_drug_forest(plot_df, "EXELBINE", "MASM4")

pdrug2 <- plot_drug_forest(plot_df, "PEVONEDISTAT", "MASM4")
pdrug3 <- plot_drug_forest(plot_df, "SEPANTRONIUM.BROMIDE", "MASM4")

# assemble Figure 7B
p7b_row1 <- ggarrange(p_surv_masm3, p_dfi$p_multi, p_meta_os, nrow = 1, ncol = 3,
                      labels = letters[1:3],
                      font.label = list(size = 18, color = "black", face = "bold"),
                      widths = c(1, 1.2, 1.2))
p7b_row2 <- ggarrange(p_main, pdrug1, pdrug4, nrow = 1, ncol = 3,
                      labels = letters[4:6],
                      font.label = list(size = 18, color = "black", face = "bold"),
                      widths = c(1.4, 1, 1))
figure7B <- ggarrange(p7b_row1, p7b_row2, nrow = 2, ncol = 1, heights = c(1, 1))
save_tiff(figure7B, "figure7B_clinical_drug.tiff", width = 20, height = 14)

# ---- 10. Supplement: MASM correlation matrices --------------------------------
# 10.1 scRNA epithelial MASM correlation, annotated with iCMS/PDS calls
sc_anno <- subtype_sc %>% column_to_rownames("patient")
masm_mat_sc <- all_data_tumor %>%
  select(patient, all_of(masm_cols7)) %>%
  column_to_rownames("patient")
cor_sc <- cor(masm_mat_sc, use = "complete.obs")
sc_anno <- sc_anno[rownames(cor_sc), c("iCMS_confident", "PDS"), drop = FALSE]
ha_sc <- HeatmapAnnotation(iCMS = sc_anno$iCMS_confident, PDS = sc_anno$PDS)
ht_sc <- grid.grabExpr(draw(Heatmap(
  cor_sc, name = "rho", col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
  top_annotation = ha_sc, rect_gp = gpar(col = "white", lwd = 0.5),
  column_title = "Epithelial MASM correlation (scRNA)")))

# 10.2 TCGA-COREAD: MASM + continuous subtype scores (iCMS shown as annotation)
corr_data <- clinical_all[, c(masm_cols7, subtype_cols)]
M <- cor(corr_data, method = "pearson", use = "complete.obs")

rownames(M) <- gsub("Tumor_", "", rownames(M))
colnames(M) <- gsub("Tumor_", "", colnames(M))
heat_corr_final <- grid.grabExpr({
  pheatmap(M, color = cor_palette, breaks = seq(-1, 1, length.out = 101),
           cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE,
           number_color = "black", fontsize_number = 10, border_color = "white",
           cellwidth = 25, cellheight = 25, main = "", angle_col = 45,
           fontsize = 11, treeheight_row = 20, treeheight_col = 20)
})

# 10.3 cell lines: MASM + package iCMS score correlation
cor_mat_ccl <- cor(temp_ccl[, c(paste0("MASM", 1:7), "iCMS2", "iCMS3")],
                   use = "complete.obs")
heat_corr_ccl <- grid.grabExpr({
  pheatmap(cor_mat_ccl, color = cor_palette, breaks = seq(-1, 1, length.out = 101),
           cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE,
           number_color = "black", fontsize_number = 8, border_color = "white",
           cellwidth = 22, cellheight = 22, main = "", angle_col = 45,
           fontsize = 9, treeheight_row = 20, treeheight_col = 20)
})
heat_corr_ccl_final <- annotate_figure(
  heat_corr_ccl, top = text_grob("Correlation matrix (colon cancer cell lines)",
                                 face = "plain", size = 14, vjust = 0))

# 10.4 per-cohort correlation heatmaps (MASM + PDS + SMI + CMS distances)
cor_heat_cols <- c(masm_cols7, paste0("PDS", 1:3), "SMI", paste0("d.CMS", 1:4),
                   "i2", "i3", "i2i3")
# i2/i3/i2i3 now come from iCMS.SSC inside run_gse_array_workflow
heat_gse39582  <- plot_cor_heatmap(cor(GSE39582_re$plot_df[, cor_heat_cols],  use = "complete.obs"),
                                   title_text = "Correlation matrix (GSE39582)")
heat_gse248381 <- plot_cor_heatmap(cor(GSE248381_re$plot_df[, cor_heat_cols], use = "complete.obs"),
                                   title_text = "Correlation matrix (GSE248381)")
heat_gse156915 <- plot_cor_heatmap(cor(GSE156915_re$plot_df[, cor_heat_cols], use = "complete.obs"),
                                   title_text = "Correlation matrix (GSE156915)")
heat_metaCRC   <- plot_cor_heatmap(cor(metaCRC_re$plot_df[, cor_heat_cols],   use = "complete.obs"),
                                   title_text = "Correlation matrix (metaCRC)")
heat_coread    <- plot_cor_heatmap(cor(tcga_re$plot_df[, cor_heat_cols],      use = "complete.obs"),
                                   title_text = "Correlation matrix (TCGA-COREAD)")


# 10.5 quadrant plots (exploratory)
plot_quadrant <- function(df, class_var, title) {
  ggplot(df, aes(x = scale(Tumor_MASM4 - Tumor_MASM2),
                 y = scale(Tumor_MASM6 - Tumor_MASM7), col = .data[[class_var]])) +
    geom_point() + theme_hd_minimal() + ggtitle(title) +
    scale_color_manual(name = "CMS types", values = cms_colors) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed")
}
pp_quadrant_tcga    <- plot_quadrant(clinical_all, "cms_label", "TCGA-COREAD")
pp_quadrant_metacrc <- plot_quadrant(metaCRC_re$plot_df, "prediction", "CRC_metacohort_July_2024")


p_allsup_corr <- ggarrange(
  ht_sc, heat_corr_final, heat_corr_ccl_final,
  heat_coread, heat_gse39582, heat_metaCRC,
  heat_gse248381, heat_gse156915, ggarrange(pp_quadrant_tcga, pp_quadrant_metacrc, nrow = 1),
  nrow = 3, ncol = 3, labels = letters[1:9],
  font.label = list(size = 18, color = "black", face = "bold"))
save_tiff(p_allsup_corr, "supplementaryfigure12_correlations.tiff", width = 24, height = 21)

# ---- 11. Supplement: subtype support ------------------------------------------
# 11.1 overlap between iCMS marker genes and MASM gene clusters
ICMS_genelist <- as.data.frame(read_excel("ICMS_genelist.xlsx"))
rownames(ICMS_genelist) <- ICMS_genelist$GENE
ICMS_genelist$ICMS_state <- ifelse(ICMS_genelist$ICMS %in% c("iCMS2_Up", "iCMS3_Down"),
                                   "iCMS2", "iCMS3")
ICMS_genelist <- as.matrix(ICMS_genelist)
all_gene_num <- table(unique(kegg_signal$symbol) %in% all_gene_pathway$gene)[1]
res_gene_overlap <- icms_masm_overlap_test(ICMS_genelist, gene_cluster, universe = all_gene_num)
p_bar_overlap <- plot_overlap_bar(res_gene_overlap$hyper_result) +
  ylim(0, 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 11.2 CMS enrichment per MASM gene cluster
res_cms_enrich <- run_MASM_CMS_pairwise_enrich(gene_cluster = gene_cluster,
                                               universe_gene = all_gene_num)
cms_num_plot <- res_cms_enrich$bubble_plot +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

# 11.3 PCA of MASM scores colored by MSI / iCMS class
pca_sc  <- masm_pca_ggplot(all_data_tumor, color_var = "MSI", point_size = 2.6) +
  ggtitle("Epithelial cells (scRNA)")
pca_icms <- masm_pca_ggplot(clinical_all %>% filter(!is.na(iCMS_class)),
                            color_var = "iCMS_class", point_size = 2.6) +
  ggtitle("TCGA-COREAD cohort")
pca_ccl_df <- temp_ccl
colnames(pca_ccl_df) <- gsub("MASM", "Tumor_MASM", colnames(pca_ccl_df))
pca_ccl <- masm_pca_ggplot(pca_ccl_df, color_var = "MSI", point_size = 2.6) +
  ggtitle("Colon cancer cell lines")

# 11.4 ROC for MSI prediction (MASM + package iCMS scores)
# roc_objects keys verified against batch_roc_analysis: "<outcome>|<predictor>"
roc_df_sc <- all_data_tumor %>%
  left_join(icms_sc %>% select(patient, iCMS2 = i2, iCMS3 = i3), by = "patient") %>%
  mutate(msi_bin = case_when(MSI == "MSI" ~ 1, MSI == "MSS" ~ 0, TRUE ~ NA_real_))
roc_result_sc <- batch_roc_analysis(roc_df_sc, pred_vars = c("iCMS2", "iCMS3", masm_cols7),
                                    outcome_dict = list(MSI = "msi_bin"), min_per_group = 8)
roc_msi_sc <- PlotROCclinical_summary(roc_result_sc, target_outcome = "MSI")

roc_df_tcga <- clinical_all %>%
  mutate(msi_bin = case_when(MSI == "MSI" ~ 1, MSI == "MSS" ~ 0, TRUE ~ NA_real_)) %>%
  rename_with(~ str_replace(., "^Tumor_MASM", "TCGA_MASM"), starts_with("Tumor_MASM"))
roc_result_tcga <- batch_roc_analysis(
  roc_df_tcga, pred_vars = c("TCGA iCMS2", "TCGA iCMS3", paste0("TCGA_MASM", 1:7)),
  outcome_dict = list(MSI = "msi_bin"), min_per_group = 8)
roc_msi_tcga <- PlotROCclinical_summary(roc_result_tcga, target_outcome = "MSI")

p_roc <- plot_two_roc(roc_obj1 = roc_result_sc$roc_objects$`MSI|iCMS2`,
                      roc_obj2 = roc_result_sc$roc_objects$`MSI|Tumor_MASM3`,
                      name1 = "iCMS2", name2 = "MASM3",
                      plot_title = "MSI classification (epithelial features)",
                      anno_x = 0.55, anno_y = 0.1)

# 11.5 stage association forest (scRNA)
res_nonmeta <- masm_stage_forest_facet_masm_nonmeta(
  df = all_data_tumor, min_per_group = 2, masm_cols = masm_cols7)

p_forest <- res_nonmeta$plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_allsup_subtype <- ggarrange(
  res_gene_overlap$venn_plot, p_bar_overlap, cms_num_plot,
  pca_sc, pca_icms, pca_ccl,
  roc_msi_sc$p_bar + ggtitle("MSI prediction (scRNA)"),
  roc_msi_tcga$p_bar + ggtitle("MSI prediction (TCGA-COREAD)"),
  p_roc,
  nrow = 3, ncol = 3, labels = letters[1:9],
  font.label = list(size = 18, color = "black", face = "bold"))
save_tiff(p_allsup_subtype, "supplementaryfigure13_subtypes.tiff", width = 22, height = 20)

# ---- 12. Supplement: clinical associations -------------------------------------
box_tcga <- function(group_var, xlab, df = clinical_all) {
  plot_masm_box(df, group_var = group_var, masm_cols = masm_cols7) +
    ggtitle("TCGA-COAD cohort") + xlab(xlab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size = 12))
}
p_cms_tcga  <- box_tcga("cms_label", "CMS labels", clinical_all %>% filter(cms_label != "NOLBL"))
p_cimp_tcga <- box_tcga("cimp", "CIMP")

p_stage_tcga <- box_tcga("stage2", "Stages", clinical_all %>% filter(!is.na(stage2)))
p_kras_tcga  <- box_tcga("kras_mut_label", "KRAS mutation")
p_braf_tcga  <- box_tcga("braf_mut_label", "BRAF mutation")
p_pt_tcga <- box_tcga("pt_lable", "Pathological tumour (pT) status",
                      clinical_all %>% filter(!is.na(pt)))
p_pm_tcga <- box_tcga("pm_lable", "Distant metastasis (pM) status",
                      clinical_all %>% filter(!is.na(pm)))
p_pn_tcga <- box_tcga("pn_lable", "Lymph-node (pN) status",
                      clinical_all %>% filter(!is.na(pn)))

p_sc_stage <- plot_masm_box(all_data_tumor, group_var = "Stage_Binary_final",
                            masm_cols = masm_cols7) +
  ggtitle("Epithelial cells (scRNA)") + xlab("Stages") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size = 12))

# MASM vs TCGA deconvolution / genomic features (ComplexHeatmap)
feature_groups <- list(
  B_cell = c("B Cells Memory", "B Cells Naive"),
  T_cell = c("T Cells CD4 Memory Activated", "T Cells CD4 Memory Resting",
             "T Cells CD4 Naive", "T Cells CD8", "T Cells Follicular Helper",
             "T Cells gamma delta", "T Cells Regulatory Tregs"),
  Macrophage = c("Macrophages M0", "Macrophages M1", "Macrophages M2"),
  DC = c("Dendritic Cells Activated", "Dendritic Cells Resting"),
  Other_immune = c("Monocytes", "NK Cells Activated", "NK Cells Resting",
                   "Plasma Cells", "Mast Cells Activated", "Mast Cells Resting"),
  Immune_summary = c("Lymphocytes", "Neutrophils...60", "Eosinophils...61",
                     "Mast Cells", "Dendritic Cells", "Macrophages"),
  Stroma = c("Stromal Fraction", "Leukocyte Fraction", "TIL Regional Fraction"),
  Function = c("Proliferation", "Wound Healing", "Macrophage Regulation",
               "Lymphocyte Infiltration Signature Score",
               "IFN-gamma Response", "TGF-beta Response"),
  Genomic = c("SNV Neoantigens", "Indel Neoantigens", "Silent Mutation Rate",
              "Nonsilent Mutation Rate", "Number of Segments",
              "Fraction Altered", "Aneuploidy Score",
              "Homologous Recombination Defects", "Intratumor Heterogeneity"),
  Receptor = c("BCR Evenness", "BCR Shannon", "BCR Richness",
               "TCR Shannon", "TCR Richness", "TCR Evenness"),
  Other = c("CTA Score", "Th1 Cells", "Th2 Cells", "Th17 Cells"))
feature_df <- map_dfr(names(feature_groups),
                      ~ tibble(feature = feature_groups[[.x]], group = .x)) %>%
  filter(feature %in% colnames(cell_data)) %>%
  mutate(feature_clean = case_when(feature == "Neutrophils...60" ~ "Neutrophils",
                                   feature == "Eosinophils...61" ~ "Eosinophils",
                                   TRUE ~ str_replace(feature, " ", "_")))
cell_data_clean <- cell_data %>%
  select(patient, all_of(feature_df$feature)) %>%
  rename(Neutrophils = `Neutrophils...60`, Eosinophils = `Eosinophils...61`) %>%
  rename_with(~ str_replace(., " ", "_"), .cols = -c(patient, Neutrophils, Eosinophils))

merged_cor_data <- cell_data_clean %>%
  inner_join(clinical_all %>%
               rename_with(~ str_replace(., "Tumor_MASM", "MASM"), starts_with("Tumor_MASM")),
             by = "patient")
message("Matched samples: ", nrow(merged_cor_data))

cor_results <- map_dfr(paste0("MASM", 1:7), function(m) {
  map_dfr(feature_df$feature_clean, function(f) {
    x <- as.numeric(merged_cor_data[[m]]); y <- as.numeric(merged_cor_data[[f]])
    valid <- !is.na(x) & !is.na(y); n <- sum(valid)
    ct <- if (n > 5) cor.test(x[valid], y[valid], method = "spearman") else NULL
    tibble(MASM = m, Feature = f,
           FeatureGroup = feature_df$group[feature_df$feature_clean == f][1],
           n = n, rho = if (is.null(ct)) NA else ct$estimate,
           p_value = if (is.null(ct)) NA else ct$p.value)
  })
}) %>% mutate(fdr = p.adjust(p_value, method = "BH"),
              sig = case_when(fdr < 0.001 ~ "***", fdr < 0.01 ~ "**",
                              fdr < 0.05 ~ "*", TRUE ~ ""))
write_csv(cor_results, "TCGA_MASM_all_features_correlation.csv")

cor_matrix <- cor_results %>% select(MASM, Feature, rho) %>%
  pivot_wider(names_from = Feature, values_from = rho) %>%
  column_to_rownames("MASM") %>% as.matrix()
sig_matrix <- cor_results %>% select(MASM, Feature, sig) %>%
  pivot_wider(names_from = Feature, values_from = sig) %>%
  column_to_rownames("MASM") %>% as.matrix()
feature_order <- feature_df %>% arrange(group, feature_clean) %>% pull(feature_clean)
cor_matrix <- cor_matrix[, feature_order]
sig_matrix <- sig_matrix[, feature_order]
col_groups <- feature_df %>% arrange(group, feature_clean) %>% pull(group)
group_colors <- c(B_cell = "#E69F00", T_cell = "#56B4E9", Macrophage = "#009E73",
                  DC = "#F0E442", Other_immune = "#0072B2", Immune_summary = "#D55E00",
                  Stroma = "#CC79A7", Function = "#999999", Genomic = "#E64B35",
                  Receptor = "#4DBBD5", Other = "#00A087")
ht <- Heatmap(
  cor_matrix, name = "Spearman rho",
  col = colorRamp2(c(-0.6, 0, 0.6), c("#2166AC", "white", "#B2182B")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(sig_matrix[i, j]) && sig_matrix[i, j] != "")
      grid.text(sig_matrix[i, j], x = x, y = y, gp = gpar(fontsize = 8, fontface = "bold"))
  },
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_groups, column_gap = unit(2, "mm"),
  row_title = "MASM Module", column_title = "TCGA Features",
  column_names_rot = 135, row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 9), rect_gp = gpar(col = "white", lwd = 0.5),
  top_annotation = HeatmapAnnotation(Group = col_groups, col = list(Group = group_colors),
                                     show_annotation_name = FALSE,
                                     annotation_legend_param = list(title = "Feature Group", nrow = 6)),
  heatmap_legend_param = list(title = "Spearman rho", at = c(-0.6, -0.3, 0, 0.3, 0.6)))
ht_cell <- grid.grabExpr(draw(ht, merge_legend = TRUE))

p_allsup_clinical <- ggarrange(
  p_cms_tcga, p_cimp_tcga, p_kras_tcga, p_braf_tcga,
  p_stage_tcga, p_sc_stage, p_pt_tcga, p_pn_tcga, p_pm_tcga,
  nrow = 3, ncol = 3, labels = letters[1:9],
  font.label = list(size = 18, color = "black", face = "bold"))
p_allsup_clinical2 <- ggarrange(
  p_allsup_clinical,
  ggarrange(p_forest + ggtitle("Epithelial cells (scRNA)"), ht_cell, nrow = 1,
            labels = letters[10:11],
            font.label = list(size = 18, color = "black", face = "bold")),
  nrow = 2, ncol = 1, heights = c(3, 2))
save_tiff(p_allsup_clinical2, "supplementaryfigure14_clinical.tiff", width = 23, height = 21)

# ---- 13. Supplement: drugs ------------------------------------------------------
supp_drugs <- setdiff(unique(plot_df$drug_names), target_drugs)
supp_anno <- c("A.154475" = "Lestaurtinib (FLT3i)", "AK.55409" = "AK-55409 (HSP90i)",
               "AR.12" = "AR-12 (ER stress/chaperone)", "ASP.7487" = "ASP-7487 (ATRi)",
               "AT.7519" = "AT-7519 (CDKi)", "AT7867" = "AT7867 (AKT/p70S6Ki)",
               "BMS.345541" = "BMS-345541 (IKK/NF-kBi)", "GSK.2141795" = "Uprosertib (AKTi)",
               "JNJ.26854165" = "JNJ-26854165 (MDM2i)", "PI.103" = "PI-103 (PI3K/mTORi)",
               "TW.37" = "TW-37 (Bcl-2i)", "PARBENDAZOLE" = "Parbendazole (Microtubule)",
               "SB.225002" = "SB-225002 (CXCR2i)", "KX2.391" = "KX2-391 (Src/tubulin)",
               "BREFELDIN.A" = "Brefeldin-A (Golgi transport)", "AT13387" = "AT13387 (HSP90i)")
p_supp <- plot_drug_dot(plot_df, supp_drugs, supp_anno,
                        title = "Additional significant drugs (>=2 databases)",
                        filename = "Supplementary_Drug_Dotplot_Final.pdf",
                        width = 7, height = 11)

bubble_data <- plot_df_gene %>%
  filter(n_sig_db == 3) %>%
  mutate(neg_log10P = -log10(pvalue),
         pair_label = paste0(drug_names, " \n ", type, "|", module),
         pair_label = factor(pair_label, levels = unique(pair_label)),
         database = str_replace(database, "\\^2", "2"))
p_bubble_gene <- ggplot(bubble_data, aes(x = database, y = pair_label)) +
  geom_point(aes(size = neg_log10P, fill = rho), shape = 21, stroke = 0.3, color = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       limits = c(-0.8, 0.8), midpoint = 0, name = "rho") +
  scale_size_continuous(range = c(3, 10), name = "-log10(P)") +
  labs(x = "Database", y = "Drug-gene pair",
       title = "Consistent significant correlations across 3 datasets") +
  theme_hd_plain() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1))

p_msi_ccl <- plot_masm_box(temp_ccl, group_var = "MSI",
                           masm_cols = c(paste0("MASM", 1:7))) +
  ggtitle("Colon cancer cell lines") + xlab("MSI")


pdrug_gene1 <- plot_drug_forest(plot_df_gene, "PEVONEDISTAT", "PPP1R3D",
                                "Pevonedistat vs PPP1R3D (MASM3 genes)")
pdrug_gene2 <- plot_drug_forest(plot_df_gene, "SEPANTRONIUM.BROMIDE", "CDKN2B",
                                "Sepantronium bromide vs CDKN2B (MASM6 genes)")
pdrug_gene4 <- plot_drug_forest(plot_df_gene, "SEPANTRONIUM.BROMIDE", "MOSMO",
                                "Sepantronium bromide vs MOSMO (MASM6 genes)")
pdrug_gene5 <- plot_drug_forest(plot_df_gene, "SEPANTRONIUM.BROMIDE", "RPRM",
                                "Sepantronium bromide vs RPRM (MASM4 genes)")

p_allsup_drug <- ggarrange(
  p_supp, p_msi_ccl, p_bubble_gene,
  ggarrange(pdrug_gene1, pdrug_gene2, pdrug_gene4, pdrug_gene5, pdrug2, pdrug3,
            nrow = 2, ncol = 3),
  nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"),
  font.label = list(size = 18, color = "black", face = "bold"))
save_tiff(p_allsup_drug, "supplementaryfigure15_drugs.tiff", width = 22, height = 18)


res <- run_model_system_pipeline(all_loading_list)
icms_sys <- run_icms_across_systems(all_data_bulk, res, exp_matrix, all_loading_list)
 cross <- run_mmcms_cross(res)
gg    <- run_guest_gradient_gemm(res, cross)

emt_sys <- run_emt_signature_analysis(res, exp_matrix, sc_ccle,
                                                                            all_data_bulk = all_data_bulk,
                                                                                bulk_expr_list = expr_list_full)
pw_sys <- run_pathway_signature_analysis(res, exp_matrix, sc_ccle,
                                        all_data_bulk = all_data_bulk,
                                          bulk_expr_list = expr_list_full)


scepi <- compute_scepi_all(epi_all_tumor, all_loading_list)
cor_plus <- bind_rows(pw_sys$cor, scepi$pathway_cor)
PP<-plot_masm_pathway_heatmap(cor_plus, thresh = 0.4,
out_file = "results/Fig_Hallmark_4systems.pdf", width = 10)


scfibro <- compute_scepi_all(fibro_all_tumor, all_loading_list)

est_audit <- run_estimate_audit(expr_list_full)


temp_res_GSE200997_epi<-readRDS("GSE200997_cpm.rds")
GSE200997_epi_tumor<-temp_res_GSE200997_epi[grepl("^T",rownames(temp_res_GSE200997_epi)),]

sc_all <- compute_sc_all_celltypes(list(scEpi = epi_all_tumor, scFibro = mat_list$Fibroblas, scTcell =mat_list$`T cell`, scBcell =mat_list$`B cell`, scEndothelial = mat_list$Endothelial, scMyeloid = mat_list$Myeloid, scPlasma = mat_list$Plasma,scEpi2=GSE200997_epi_tumor), all_loading_list)


sc_all <- compute_sc_all_celltypes(list(scEpi = epi_all_tumor, CCLE=t(exp_matrix),  scFibro = mat_list$Fibroblas, scTcell =mat_list$`T cell`, scBcell =mat_list$`B cell`, scEndothelial = mat_list$Endothelial, scMyeloid = mat_list$Myeloid, scPlasma = mat_list$Plasma,scnon_epi=no_epi_all_tumor_pseudo), all_loading_list)
GSE35602_exp<-readRDS( "GSE35602_exp_matrix.rds")
temp_epi_data<-GSE35602_exp$expr_SxG[rownames(GSE35602_exp$expr_SxG) %in% GSE35602_exp$pheno[GSE35602_exp$pheno$`tissue:ch1` =="epithelium","geo_accession"],]
rownames(temp_epi_data)<-paste("GSE35602","Epi",rownames(temp_epi_data),sep="_")
sc_all <- compute_sc_all_celltypes(list(scEpi = epi_all_tumor, scEpi2=GSE200997_epi_tumor,scnon_epi=no_epi_all_tumor_pseudo,epi=temp_epi_data), all_loading_list)

sc_all_mss <- compute_sc_all_celltypes(list(scEpi = epi_all_tumor, scEpi2=GSE200997_epi_tumor,scnon_epi=no_epi_all_tumor_pseudo,epi=temp_epi_data, CCLE=t(exp_matrix)), all_loading_list)


mod_res <- run_module_cor(all_data_bulk, other_systems = c(list(CCLE = sc_ccle, GEMM6363 = res$sc18, GEMM218776 = res$sc39), lapply(sc_all, "[[", "masm_scores")))


epi_bulk_masm<-GSE35602_re$plot_df[GSE35602_re$plot_df$patient %in% GSE35602_exp$pheno[GSE35602_exp$pheno$`tissue:ch1` =="epithelium","geo_accession"],]
colnames(epi_bulk_masm)<-gsub("Tumor_","",colnames(epi_bulk_masm))

all_list_masm<-list(CCLE = sc_ccle, scEpi = all_data_tumor,scNonEpi=all_data_tumor_pseudono_epi, GEMM6363 = res$sc18, GEMM218776 = res$sc39)
saveRDS(all_list_masm,file="all_list_masm.rds")
mod_res <- run_module_cor(all_data_bulk, other_systems = list(CCLE = sc_ccle, scEpi = all_data_tumor,scNonEpi=all_data_tumor_pseudono_epi, GEMM6363 = res$sc18, GEMM218776 = res$sc39))
#scEpi2 = GSE200997_epi_masm,Epi=epi_bulk_masm,
 ht_sum <- plot_module_cor_summary(
   mod_res,out_file = "MASM_module_correlations_display" )
ht_sum$heatmap

sign_ok <- apply(sign( ht_sum$rho[,1:4]), 1, function(x) all(x==1)||all(x==-1))
#, GEMM6363 = res$sc18, GEMM218776 = res$sc39

cor_plus <- bind_rows(pw_sys$cor, sc_all$scEpi$pathway_cor, sc_all$scFibro$pathway_cor, sc_all$scTcell$pathway_cor, sc_all$scBcell$pathway_cor, sc_all$scEndothelial$pathway_cor, sc_all$scMyeloid$pathway_cor, sc_all$scPlasma$pathway_cor)
PP<-plot_masm_pathway_heatmap(cor_plus, thresh = 0.4,
out_file = "results/Fig_Hallmark_4systems.pdf", width = 10)

cor_plus <- bind_rows(pw_sys$cor, sc_all$scEpi$pathway_cor, sc_all$scEpi2$pathway_cor, sc_all$scnon_epi$pathway_cor,sc_all$epi$pathway_cor)
PP<-plot_masm_pathway_heatmap(cor_plus, thresh = 0.4,
out_file = "results/Fig_Hallmark_4systems.pdf", width = 10)


PP<- plot_masm_pathway_core(cor_plus,   sys_map = c(scEpi     = "scEpi_discovery",
                                               scEpi2    = "scEpi_verification",
                                               Epi       = "Bulk_epi",
                         human     = "human_bulk",CCLE      = "CCLE",
                                               scnon_epi = "sc_nonEpi"),thresh = 0.4,
                            exclude_paths = c("Custom_EMT_TF"),
                            out_file = "Fig_masm_pathway_core.pdf");PP$ht

PP$mat


 mut_res <- run_mutation_anchor(datasets = list(GSE156915 = all_GSE156915,
                   GSE39582  = all_GSE39582,
                   TCGA      = clinical_all),
   out_prefix = "results/MASM_mutation_anchor")

  filter_epi<-all_meta[all_meta$Cell_type=="Epithelial",]
  filter_epi<- as.data.frame(colData(cds))
   inputs <- list(
 expr         = all_gene_sct_count_matrix[,filter_epi$cell.ID],          # genes x cells, log-normalized
   meta         = filter_epi$patient,
   loading_list = all_loading_list,         # list(MASM1..MASM7 = named numeric)
   scores       = single_cell_masm
 )
 res_single <- run_sc_diagnosis(inputs, out_dir = "results/sc_diag",
                         pt_col = "pseudotime", patient_col = "patient")


 ph <- setNames(as.character(res_single$cc$cells$phase), res_single$cc$cells$cell)
g2m_sens <- run_g2m_sensitivity(all_gene_sct_count_matrix[, filter_epi$cell.ID],
                                 all_loading_list, ph,
                                 modules = c("MASM7","MASM5"),
                                 rm_sets = list(cell_cycle_core = c("CCNB2","CCNA2","SKP2"),
                                                cyclins = c("CCNB2","CCNA2"),
                                                PRMT1_only = c("PRMT1")),
                                 score_method = "rank")
 print(g2m_sens)


  base_scores <- res_single$new_scores
 topo_sens <- run_target_sensitivity(all_gene_sct_count_matrix[, filter_epi$cell.ID],
                                   all_loading_list, target_mat = as.matrix(base_scores),
                                   modules = c("MASM5","MASM7"),
                                   rm_sets = list(cell_cycle_core = c("CCNB2","CCNA2","SKP2"),
                                                  ribosomal_RPS6 = c("RPS6"),
                                                  MYC_axis = c("MYC","RHEB","RUVBL1")),
                                   also_single_gene = FALSE, score_method = "rank")


  grp <- setNames(as.numeric(clinical_all$kras_mut), clinical_all$sample)
 anch_sens <- run_group_sensitivity(tcga_coread, all_loading_list, grp,
                                  modules = c("MASM2","MASM3","MASM4","MASM5","MASM7"),
                                  rm_sets = list(cell_cycle_core = c("CCNB2","CCNA2","SKP2"),
                                                 ribosomal_RPS6 = c("RPS6")),
                                  score_method = "weighted")


 m6_sens <- run_target_sensitivity(epi_all_tumor, all_loading_list,
                                 target_mat = as.matrix(sc_all$scEpi$emt_rows[,2:5]),
                                 modules = "MASM6",
                                 rm_sets = list(AREG_only = c("AREG"),
                                                secreted_signal = c("AREG","CCL4L2"),
                                                tgfb_matrix = c("FHL1","CDKN2B")),
                                 score_method = "weighted")
  GSE200997_epi_masm<-calc_masm_weight_score_simple(GSE200997_epi_tumor, all_loading_list)
 mod_cohsig<-run_module_coherence(GSE200997_epi_tumor, all_loading_list)
 mod_coh<-run_module_coherence(epi_all_tumor, all_loading_list)
 mod_coh_tcga<-run_module_coherence(t(tcga_coread), all_loading_list)
 mod_coh_ccle<-run_module_coherence(t(exp_matrix), all_loading_list)
 mod_coh_bulkmetacrc<-run_module_coherence(t(gene_all_meta), all_loading_list)
 mod_coh_hs39<-run_module_coherence(t(res$hs39), all_loading_list)
  mod_coh_hs18<-run_module_coherence(t(res$hs18), all_loading_list)

  coh_list <- list(discovery =  mod_coh, bulk_meta = mod_coh_bulkmetacrc,
                   tcga = mod_coh_tcga, ccle = mod_coh_ccle)
   res_wk <- run_weak_gene_audit(coh_list, as.matrix(epi_all_tumor), all_loading_list,
                                 bulk_layers = c("discovery", "bulk_meta", "tcga"))
   res_wk2 <- run_weak_gene_audit(coh_list, epi_all_tumor, all_loading_list,
                                  group_vec = group_vec_aligned)
  sys_list <- list(discovery = epi_all_tumor, tcga = t(tcga_coread),
                  ccle = t(exp_matrix))
   path_sets <- list(MYC_V2 =gsets$HALLMARK_MYC_TARGETS_V2, E2F = gsets$HALLMARK_E2F_TARGETS,
                     G2M = gsets$HALLMARK_G2M_CHECKPOINT, EMT = gsets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)
   res_xs <- run_cross_system(sys_list, all_loading_list, path_sets,
                              core_genes = setdiff(names(gene_cluster), res_wk$drop_genes))

    gsets <- get_hallmark_sets("human")

    gt18 <- plot_masm_genotype_box(res$sc18, out_file = "MASM_GEMM6363_gene")
    gt39 <- plot_masm_genotype_box(res$sc39, out_file = "MASM_GEMM218776_gene")


    cohort_cfgs_os <- list(
  "TCGA-COREAD" = list(df = clinical_all, time = "OS.time", status = "OS",
                       covars = c(Age = "age_at_initial_pathologic_diagnosis",
                                  Sex = "gender", Stage = "stage")),
  "AC-ICAM"     = list(df = all_acicam, time = "OS_month", status = "OS_status",
                       covars = c(Age = "age_at_diagnosis", Sex = "Sex", Stage = "stage"))
)

cohort_cfgs_dss <- list(
  "TCGA-COREAD" = list(df = clinical_all, time = "DSS.time", status = "DSS",
                       covars = c(Age = "age_at_initial_pathologic_diagnosis",
                                  Sex = "gender", Stage = "stage"))
)

compendium <- build_compendium_cfgs(pheno.data, all_data_bulk)
cohort_cfgs_dfs <- compendium$DFS
cohort_cfgs_os  <- c(cohort_cfgs_os,
                     compendium$OS [setdiff(names(compendium$OS),  names(cohort_cfgs_os))])
cohort_cfgs_dss <- c(cohort_cfgs_dss,
                     compendium$DSS[setdiff(names(compendium$DSS), names(cohort_cfgs_dss))])
names(cohort_cfgs_os); names(cohort_cfgs_dss); names(cohort_cfgs_dfs)

within_os   <- run_within_all(cohort_cfgs_os, "OS")
 within_dfs  <- run_within_all(cohort_cfgs_dfs, "DFS")
 os_stage    <- run_surv_by_stage(cohort_cfgs_os, "OS")
  os_res      <- run_surv_pipeline(cohort_cfgs_os, "OS")
 dfs_res      <- run_surv_pipeline(cohort_cfgs_dfs, "DFS")

 p_d<-panel_forest_pretty(within_os, within_os$meta, "Tumor_MASM6", "CMS2", "full");p_d
p_d2<-panel_forest_pretty(within_os, within_os$meta, "Tumor_MASM7", "CMS4", "full");p_d2
p_d3<-panel_forest_pretty(within_os, within_os$meta, "Tumor_MASM7", "CMS4", "uni");p_d3
p_d4<-panel_forest_pretty(within_os, within_os$meta, "Tumor_MASM6", "CMS2", "uni");p_d4

 p_da1<-panel_forest_os(os_res, meta_df = os_res$meta, masm = "Tumor_MASM7",
                layer = "full", endpoint = "OS")
  p_da2<-panel_forest_os(os_res, meta_df = os_res$meta, masm = "Tumor_MASM6",
                layer = "full", endpoint = "OS")
    p_da3<-panel_forest_os(dfs_res, meta_df = dfs_res$meta, masm = "Tumor_MASM7",
                layer = "full", endpoint = "DFS")
  p_df<-panel_forest_pretty(within_dfs, within_dfs$meta, "Tumor_MASM7", "CMS4", "full",endpoint = "DFS");p_df

  ggarrange(p_d,p_d2,p_da1,p_da2,p_df,p_da3,nrow=2,ncol=3)


 grp <- list(Apc_driven  = c("A", "AP", "APN", "AK"),
            Kras_driven = c("KP", "KPN","BP", "BPN"),
            Braf_driven = c("s"))
gt39g <- plot_masm_genotype_box(res$sc39, geno_groups = grp,
                                geno_order = c("Apc_driven", "Kras_driven", "Braf_driven"),
                                out_file = "MASM_GEMM218776_")
gt39g$pairwise

hres <- run_human_mss_mirror(expr_list_full, list(TCGA=clinical_all,GSE39582=all_GSE39582),
                             cohorts = c("TCGA", "GSE39582"),
                             out_file = "MASM_MSS", strict_wt=F)
g <- plot_mutation_mirror_modules(gres, hres, out_file = "MASM__module")

pal = c(Apc_driven = "#56B4E9", MAPK_WT = "#56B4E9",
APC_only   = "#0072B2",
Kras_driven = "#009E73", KRAS_mut = "#009E73",
Braf_driven = "#E69F00", BRAF_mut = "#E69F00")


  p <- ggplot(g$plot$data[g$plot$data$module=="MASM1",], aes(grp, z, fill = grp)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey75", linewidth = 0.35) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.62,
                 color = "grey20", linewidth = 0.35) +
    geom_jitter(width = 0.16, size = 0.7, alpha = 0.35, color = "grey25") +
    geom_text(data = kw_df, aes(x = -Inf, y = Inf, label = lab),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.4,
              size = base_size / 3.6, color = "grey30") +
    (if (is.null(star_df)) NULL else
      geom_text(data = star_df, aes(x = grp, y = y, label = star),
                inherit.aes = FALSE, vjust = 0, size = base_size / 3.2)) +
    facet_grid(module ~ cohort, scales = "free") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = x_lab) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.22))) +
    labs(title = "Driver genotype shifts MASM scores: concordant across modules, mirrored between mouse and human",
         x = NULL, y = "Module score (z within cohort)") +
    th +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text.y = element_text(face = "bold"))


all_mss<-hres$data
all_mss$dataset=hres$data$cohort
pw_sys_MSS <- run_pathway_signature_analysis(res, exp_matrix, sc_ccle,
                                        all_data_bulk = all_mss,
                                          bulk_expr_list = list(TCGA=expr_list_full$TCGA[,colnames(expr_list_full$TCGA) %in% all_mss$patient],GSE39582=expr_list_full$GSE39582[,colnames(expr_list_full$GSE39582) %in% all_mss$patient]))


PP<- plot_masm_pathway_heatmap(pw_sys_MSS$cor, thresh = 0.4,
                            out_file = "Fig_masm_pathway_core.pdf");PP$ht

PP$mat

mod_res_human_mouse <- run_module_cor(all_data_bulk, other_systems = list(TCGA_MSS = hres$per_cohort$TCGA$scores,GSE39582_MSS = hres$per_cohort$GSE39582$scores, GEMM6363 = res$sc18, GEMM218776 = res$sc39))

mss_gse <- colnames(expr_list_full$GSE39582)[colnames(expr_list_full$GSE39582) %in%
                hres$per_cohort$GSE39582$scores$patient]
mss_tcga <- colnames(expr_list_full$TCGA)[colnames(expr_list_full$TCGA) %in%
                hres$per_cohort$TCGA$scores$patient]
   pw_res <- run_masm_pathway_mss_mouse(
     expr_list  = list(TCGA = expr_list_full$TCGA[, mss_tcga],
                       GSE39582 = expr_list_full$GSE39582[, mss_gse],
                       GEMM6363 = res$vst18, GEMM218776 = res$vst39),
     score_list = list(TCGA = hres$per_cohort$TCGA$scores,
                       GSE39582 = hres$per_cohort$GSE39582$scores,
                       GEMM6363 = res$sc18, GEMM218776 = res$sc39),
  mouse_systems = c("GEMM6363", "GEMM218776"),
     meta_groups = list(human_MSS = c("TCGA", "GSE39582"),
                        mouse     = c("GEMM6363", "GEMM218776")),
     out_file = "MASM_pathway_MSS")


topo <- plot_masm_topo_compare(
    mod_res_human_mouse$long,
    col_groups = list(
        human_meta   =unique(all_data_bulk$dataset),
        TCGA_MSS     = "TCGA_MSS",
        GSE39582_MSS = "GSE39582_MSS",
        GEMM6363     = "GEMM6363",
        GEMM218776   = "GEMM218776"),
    ref_col = "human_meta",
    out_file = "MASM_heatmap")

exclude_pairs <- names(sign_ok[sign_ok==F])
 ht_sumHM <- plot_module_cor_summary(
   mod_res_human_mouse,sys_map = c(GSE39582_MSS    = "GSE39582_MSS",
                                                TCGA_MSS   = "TCGA_MSS",
                                                GEMM6363      = "GEMM6363",
                                                GEMM218776 = "GEMM218776"),mouse_systems = "",exclude_pairs=exclude_pairs,types="mouse",out_file = "human_mouseMASM_module_correlations_display" )
ht_sumHM$heatmap


ol <- run_masm_overlap_sensitivity(
  expr_list      = expr_list_full,
  loading_list   = all_loading_list,
  species        = "human",
  focus_modules  = paste0("MASM", 5:7),
  focus_pathways = c("HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2",
                     "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"))

pw_m39 <- masm_pw_cor(res$hs39,  res$sc39,  system = "sc39")
pw_m18 <- masm_pw_cor(res$hs18,  res$sc18,  system = "sc18")
tcga_mss <- masm_pw_cor(expr_list_full$TCGA[, mss_tcga],  hres$per_cohort$TCGA$scores,   system = "TCGA_MSS")
GSE3958_mss <- masm_pw_cor(expr_list_full$GSE39582[, mss_gse],  hres$per_cohort$GSE39582$scores,   system = "GSE39582_MSS")


long <- bind_rows(pw_m39, pw_m18, tcga_mss,GSE3958_mss)
long <- masm_pw_meta(long, groups = list(human_bulk = c("tcga_mssC", "GSE3958_mss")))

plot_masm_pw_heatmap(
  long[gsub("HALLMARK_", "", long$pathway) %in% rownames(PP$mat), ],
  systems    = c("TCGA_MSS", "GSE39582_MSS", "sc39", "sc18"),
  sys_rename = c("sc39" = "GEMM218776",
                 "sc18" = "GEMM6363"),
  thresh     = 0.3)


p3 <- pca_biplot_masm(rbind(temp39,temp18), colour_by = "mCMS",pal = my_subtype_colors);p3
p3c <- pca_biplot_masm(rbind(temp39,temp18), colour_by = "mCMS_C");p3c


gres <- run_gemm_gradient(list(sc18 = res$sc18, sc39 = res$sc39),
cohort_labels = c(sc18 = "GSE6363",
sc39 = "GSE218776"),
out_file = "MASM_GEMMscore")

cross_B <- run_mmcms_cross(res, "B")
cross_C <- run_mmcms_cross(res, "C")

temp18<-cbind(res$sc18,mCMS=cross_B$`18`[res$sc18$patient,"mCMS"],mCMS_C=cross_C$`18`[res$sc18$patient,"mCMS"])
temp39<-cbind(res$sc39,mCMS=cross_B$`39`[res$sc39$patient,"mCMS"],mCMS_C=cross_C$`39`[res$sc39$patient,"mCMS"])
all_temp=rbind(temp39,temp18)
x1_mouse <- aggregate(all_temp[2:8], list(CMS = all_temp$mCMS), mean)        # mouse
x2_human <- aggregate(all_data_bulk[3:9], list(CMS = all_data_bulk$CMS), mean) # human


x18=aggregate(temp18[2:8],list(CMS=temp18$mCMS),mean)
x39=aggregate(temp39[2:8],list(CMS=temp39$mCMS),mean)

res_cor_HM <- plot_cms_transfer(x2_human,x1_mouse ,cor_method="spearman")
res_cor_MM <- plot_cms_transfer(x18,x39,label_x ="GEMM6363",label_y = "GEMM218776"  ,cor_method="spearman")

p <- plot_cms_cor_bar(res_cor_HM, res_cor_MM)
p

ggplot(res_cor_HM$data[res_cor_HM$data$CMS=="CMS4",],
       aes(x=vx,y=vy,label=module))+
    geom_point(size=4)+
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color="#2166AC", linewidth=1)+
    ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20)+
    theme_hd_minimal_plain()+
    stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95)+
    xlab("Mean MASM score in Human")+
    ylab("Mean MASM score in Mouse")+ggtitle("CMS4 subtype")


p_cms_sub <- plot_pca_biplot(all_data_bulk,
                          color_col = "CMS",
                          title = "MASM distribution across CMS subtypes in human bulk",
                          custom_colors = my_subtype_colors,feature_cols = paste0("Tumor_MASM", c(1,5,7)))

p_cms_mouse <- plot_pca_biplot(rbind(temp39,temp18),
                          color_col = "mCMS", flip_pc_x = T,
                          title = "MASM distribution across CMS subtypes in GEMM mouse",
                          custom_colors = my_subtype_colors,feature_cols = paste0("MASM", c(1,5,7)))

ggarrange(p_cms_sub ,p_cms_mouse)


temp_x=all_GSE39582[,c("patient",paste0("Tumor_MASM", c(1:7)),"MMR_status","stage")]
colnames(temp_x)=c("patient",paste0("Tumor_MASM", c(1:7)),"MSI","stage")
temp_x$MSI[temp_x$MSI=="pMMR"]="MSS"
temp_x$MSI[temp_x$MSI=="dMMR"]="MSI"
temp_bulk=rbind(clinical_all[,c("patient",paste0("Tumor_MASM", c(1:7)),"MSI","stage")],temp_x)
p_icms_sub <- plot_pca_biplot(temp_bulk,
                          color_col = "MSI",
                          title = "MASM distribution  in human bulk",
                          custom_colors = my_subtype_colors,feature_cols = paste0("Tumor_MASM", c(1,2,4)))


p_icms_sub_ccle <- plot_pca_biplot(temp_ccl, flip_pc_x = T,
                          color_col = "MSI",
                          title = "MASM distribution  in CCLE",
                          custom_colors = my_subtype_colors,feature_cols = paste0("MASM", c(1,2,4)))
ggarrange(p_icms_sub,p_icms_sub_ccle)


pp=plot_state_trajectory(
    human_mean = human_mean,
    mouse_list = list(GSE18 = mouse_mean18, GSE39 = mouse_mean39),
    label_ys   = c("GEMM6363", "GEMM218776")
)

human_mean_cms<-x2_human[,-1]
rownames(human_mean_cms)<-x2_human[,1]
x18_mouse<-x18[,-1]
rownames(x18_mouse)<-x18[,1]

x39_mouse<-x39[,-1]
rownames(x39_mouse)<-x39[,1]

p_cms_HM <- plot_state_trajectory( human_mean = human_mean_cms,
    mouse_list = list(GEMM6363 = x18_mouse,GEMM218776 = x39_mouse),
                                  group_levels = paste0("CMS", 1:4),
                                  group_name   = "CMS",
                                  group_label  = "Consensus Molecular Subtype (CMS)",
                                  out_file = "Fig_CMS_trajectory")
p_cms_HM$plot

# Figure 7 and Supplementary Figure 14.
# Run this script from the project root with Rscript.

project_dir <- normalizePath(
  Sys.getenv("METABOLISM_PROJECT_DIR", unset = getwd()),
  mustWork = TRUE
)
source(file.path(project_dir, "Rscript", "step0_config.R"))
source(file.path(project_dir, "Rscript", "analysis_functions.R"))
setwd(data_dir)


subtype_cols <- c(paste0("d.CMS", 1:4), "i2", "i3", paste0("PDS", 1:3), "SMI")
panel_label_style <- list(size = 18, color = "black", face = "bold")

masm_cols7   <- paste0("Tumor_MASM", 1:7)
img <- rasterGrob(readTIFF(file.path(fig_dir, "concept2_01.tif"), native = TRUE), interpolate = TRUE)

p_concept <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
  geom_blank() +
  annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()
load("CRC_metacohort_July_2024.RData")
silu_colon <- readRDS("coad_silu_2022_data.rds")
silu_colon_re<-readRDS( "silu_colon_classtypes.rds")
all_loading_list <- readRDS("all_loading_list.rds")
gene_cluster <- readRDS("gene_cluster.rds")

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


GSE39582_annotation_only <- read.delim(file.path(data_dir, "GSE39582_annotation_only.txt"), header = FALSE)
clinical_gse39582 <- parse_geo_anno(GSE39582_annotation_only)
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

all_data_bulk$iCMS_call=gsub("i","iCMS",all_data_bulk$confident.icms)
all_data_bulk$CMS=all_data_bulk$prediction
all_data_bulk$cms_label<-all_data_bulk$CMS

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


res_grad_cms6_os <- plot_cox_gradient_surv(cohort_cfgs_os$`AC-ICAM`$df,
                                           "Tumor_MASM6", "OS_month", "OS_status", n_level = 8);res_grad_cms6_os$plot+ggtitle("AC-ICAM (adjusted)")+xlab("Time (Months)")
plot_os1<-res_grad_cms6_os$plot+ggtitle("AC-ICAM (adjusted)")+xlab("Time (Months)")
res_grad_cms7_os <- plot_cox_gradient_surv(cohort_cfgs_os$`AC-ICAM`$df,
                                           "Tumor_MASM7", "OS_month", "OS_status", n_level = 8);
plot_os2<-res_grad_cms7_os$plot+ggtitle("AC-ICAM (adjusted for stage & CMS)")+xlab("Time (Months)")


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
#temp_ccl$iCMS_class <- icms_ccl$iCMS_confident[idx_ccl]

gdsc_drug  <- getDrugResultnew(drugall_gdsc, exp_matrix, temp_ccl)
ctd_drug   <- getDrugResultnew(drugall_ctd, exp_matrix, temp_ccl)
prism_drug <- getDrugResultnew(drugall_prism, exp_matrix, temp_ccl)

all_drug_cor <- rbind(cbind(gdsc_drug$all_cor, database = "GDSC"),
                      cbind(ctd_drug$all_cor, database = "CTD^2"),
                      cbind(prism_drug$all_cor, database = "PRISM"))
all_drug_cor$drug_names <- str_split_fixed(all_drug_cor$drug, "\\.\\.", 3)[, 1]
all_drug_cor<-all_drug_cor[!grepl("..GDSC1.",all_drug_cor$drug),]
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

p1_drug1<-ggplot(gdsc_drug$ccl_signal_drug, aes(x =MASM2, y = `PEVONEDISTAT..GDSC2.1529.`)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,col="black") +stat_cor(method = "spearman",size=5, r.digits = 1)+theme_hd_plain()+ylab("AUC (drug response)")+ggtitle("PEVONEDISTAT")+xlab("MASM2 score in CCLE of GDSC cohort")+ylim(0.25,1.2)
p1_drug2<-ggplot(gdsc_drug$ccl_signal_drug, aes(x =MASM4, y = `PEVONEDISTAT..GDSC2.1529.`)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,col="black") +stat_cor(method = "spearman",size=5, r.digits = 1)+theme_hd_plain()+ylab("AUC (drug response)")+ggtitle("PEVONEDISTAT")+xlab("MASM4 score in CCLE of GDSC cohort")+ylim(0.25,1.2)

pall7_1<-ggarrange(plot_os1,plot_os2,p_d,p_d2,nrow=1,ncol=4,widths = c(1,1,1.2,1.2),labels = c(letters[1:4]), font.label = panel_label_style)

pall7_2<-ggarrange(pdrug1,p1_drug1,pdrug2,p1_drug2,nrow=1,ncol=4,labels = c(letters[5:8]), font.label = panel_label_style)


pall7_3<-ggarrange(pall7_1,pall7_2, p_concept ,nrow=3,labels = c("","",letters[9]), font.label = panel_label_style,heights = c(1,1,2.2))


save_tiff(pall7_3, "figure7.tiff", width = 19, height = 22, path = fig_dir)


pdrug_all<-plotVennResult(ctd_drug$all_cor,gdsc_drug$all_cor,prism_drug$all_cor,"MASM4",c("CTD^2","GDSC","PRISM"),size_x = 0)

pdrug_overlap<- pdrug_all$p1+theme(plot.title = element_text(hjust = 0.5),plot.margin = unit(c(1,1,1,1), "cm"))+annotate("text", x = 0, y = 0, label = expression(CTD^2), size = 4, angle = 60, hjust = -0.7, vjust = -6)+annotate("text", x = 0, y = 0, label = "GDSC", size = 4, hjust = -0.5, vjust = 20)+annotate("text", x = 0, y = 0, label = "PRISM", size = 4, angle = -60, hjust = 0.1, vjust = -17)

pall141 <- ggarrange(p_da1,p_da2,p_da3,p_df,pdrug_overlap,p_main,
                    labels = letters[1:6], nrow = 3 ,ncol = 2,
                    font.label = panel_label_style)


save_tiff(pall141, "supplementaryfigure14.tiff", width = 14, height = 15, path = fig_dir)

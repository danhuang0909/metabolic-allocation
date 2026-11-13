# Figure 6 and Supplementary Figure 13.
# Run this script from the project root with Rscript.

project_dir <- normalizePath(
  Sys.getenv("METABOLISM_PROJECT_DIR", unset = getwd()),
  mustWork = TRUE
)
source(file.path(project_dir, "Rscript", "step0_config.R"))
source(file.path(project_dir, "Rscript", "analysis_functions.R"))
setwd(data_dir)


panel_label_style <- list(size = 18, color = "black", face = "bold")

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

all_data_bulk<-readRDS(file="all_data_bulk.rds")
all_list_masm<-readRDS(file="all_list_masm.rds")
gene_cluster <- readRDS("gene_cluster.rds")
all_loading_list <- readRDS("all_loading_list.rds")
res <- run_model_system_pipeline(all_loading_list)
RNAseq<-c("TCGA","silu_colon")

temp_res_GSE200997_epi<-readRDS("GSE200997_cpm.rds")
GSE200997_epi_tumor<-temp_res_GSE200997_epi[grepl("^T",rownames(temp_res_GSE200997_epi)),]
 GSE200997_epi_masm<-calc_masm_weight_score_simple(GSE200997_epi_tumor, all_loading_list)
 GSE35602_re<-readRDS("GSE35602_classtypes.rds")
 GSE35602_exp<-readRDS( "GSE35602_exp_matrix.rds")
temp_epi_data<-GSE35602_exp$expr_SxG[rownames(GSE35602_exp$expr_SxG) %in% GSE35602_exp$pheno[GSE35602_exp$pheno$`tissue:ch1` =="epithelium","geo_accession"],]
rownames(temp_epi_data)<-paste("GSE35602","Epi",rownames(temp_epi_data),sep="_")
 epi_bulk_masm<-GSE35602_re$plot_df[GSE35602_re$plot_df$patient %in% GSE35602_exp$pheno[GSE35602_exp$pheno$`tissue:ch1` =="epithelium","geo_accession"],]
colnames(epi_bulk_masm)<-gsub("Tumor_","",colnames(epi_bulk_masm))

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


pseudo_res <- readRDS(file.path(data_dir, "pseudo_seurat_all_sample.rds"))
epi_all_tumor <- as.data.frame(pseudo_res$cpm_matrix[grepl("Tumor_Epi", rownames(pseudo_res$cpm_matrix)) , ])
rownames(epi_all_tumor)<-gsub("Tumor_","",rownames(epi_all_tumor))
all_data_tumor <- calc_masm_weight_score_simple(epi_all_tumor, all_loading_list)

non_epi_cpm_mat<-readRDS( "./pseudo_bulk_output/non_epi_merged_cpm.rds")
no_epi_all_tumor_pseudo <- as.data.frame(non_epi_cpm_mat[grepl("Tumor_", rownames(non_epi_cpm_mat)) , ])
rownames(no_epi_all_tumor_pseudo)<-gsub("Tumor_","",rownames(no_epi_all_tumor_pseudo))

all_data_tumor_pseudono_epi <- calc_masm_weight_score_simple(no_epi_all_tumor_pseudo, all_loading_list)

mod_res <- run_module_cor(all_data_bulk[all_data_bulk$dataset %in%RNAseq ,], other_systems = list(scEpi1 = GSE200997_epi_masm,Epi=epi_bulk_masm,CCLE = sc_ccle, scEpi = all_data_tumor[grepl("Pelka",all_data_tumor$patient),],scEpi2 = all_data_tumor[!grepl("Pelka",all_data_tumor$patient),],scNonEpi=all_data_tumor_pseudono_epi, GEMM6363 = res$sc18, GEMM218776 = res$sc39))
mod_res2<- run_module_cor(all_data_bulk[!all_data_bulk$dataset %in%RNAseq ,])
mod_res$human_metaMicro=mod_res2$human_meta
#scEpi2 = GSE200997_epi_masm,Epi=epi_bulk_masm,
ht_sum <- plot_module_cor_summary(sys_map = c(scEpi    = "scEpi_Pelka\n(n=56)",
                                              scEpi2   = "scEpi_nonPelka\n(n=52)",
                                              scEpi1 ="scEpi_verification\n(n=12)",
                                              Epi ="Epi_bulk\n(n=17)",
                                              scNonEpi = "sc_nonEpi\n(n=108)",
                                              CCLE     = "CCLE\n(n=58)"),human_label = "Bulk_RNAseq\n(meta,2 cohorts)",mouse_label = "Mouse\n(meta,2cohorts)",
                                  mod_res,out_file = "MASM_module_correlations_display",gap_after = 6  , orange_pairs = c("MASM1-MASM2", "MASM2-MASM4", "MASM1-MASM4"),
                                  blue_pairs   = c("MASM1-MASM5", "MASM5-MASM7", "MASM1-MASM7"),)
ht_sum$heatmap


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


tcga_coread <- readRDS("TCGA_COREAD_STAR_tumor_CPM_genesymbol.rds")
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

 all_GSE156915 <- readRDS(file.path(data_dir, "all_GSE156915.Rds"))

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
p_pca_ccle<-ggarrange(p_icms_sub,p_icms_sub_ccle,nrow=2,common.legend = T)


 mut_res <- run_mutation_anchor(datasets = list(GSE156915 = all_GSE156915,
                   GSE39582  = all_GSE39582,
                   TCGA      = clinical_all),
   out_prefix = "results/MASM_mutation_anchor")

 mut_res$heatmaps$meta

 ggarrange(ht_sum$heatmap ,p_pca_ccle, mut_res$heatmaps$meta,nrow=1,ncol=3)
library(grid)

ht1 <- grid.grabExpr(ComplexHeatmap::draw(ht_sum$heatmap))
ht3 <- grid.grabExpr(ComplexHeatmap::draw(mut_res$heatmaps$meta, heatmap_legend_side = "bottom",padding = unit(c(2, 8, 2, 2), "mm")), wrap = TRUE)


cross_B <- run_mmcms_cross(res, "B")
cross_C <- run_mmcms_cross(res, "C")

temp18<-cbind(res$sc18,mCMS=cross_B$`18`[res$sc18$patient,"mCMS"],mCMS_C=cross_C$`18`[res$sc18$patient,"mCMS"])
temp39<-cbind(res$sc39,mCMS=cross_B$`39`[res$sc39$patient,"mCMS"],mCMS_C=cross_C$`39`[res$sc39$patient,"mCMS"])
all_temp=rbind(cbind(Dataset="GEMM GSE218776",temp39),cbind(Dataset="GEMM E-MTAB-6363",temp18))
x1_mouse <- aggregate(all_temp[3:9], list(CMS = all_temp$mCMS), mean)        # mouse
x2_human <- aggregate(all_data_bulk[3:9], list(CMS = all_data_bulk$CMS), mean) # human
colnames(x2_human )<-gsub("Tumor_","",colnames(x2_human ))

x18=aggregate(temp18[2:8],list(CMS=temp18$mCMS),mean)
x39=aggregate(temp39[2:8],list(CMS=temp39$mCMS),mean)

res_cor_HM <- plot_cms_transfer(x1_mouse ,x2_human,cor_method="spearman")
res_cor_MM <- plot_cms_transfer(x18,x39,label_x ="GEMM E-MTAB-6363",label_y = "GEMM GSE218776"  ,cor_method="spearman")

p_bar <- plot_class_cor_bar(res_cor_HM, res_cor_MM)
p_bar

expr_list_full<-readRDS("all_data_bulk_exp_data_expr_list_full.rds")
gres <- run_gemm_gradient(list(sc18 = res$sc18, sc39 = res$sc39),
cohort_labels = c(sc18 = "E-MTAB-6363",
sc39 = "GSE218776"),
out_file = "MASM_GEMMscore")
hres <- run_human_all_mirror(expr_list_full,
                             list(TCGA = clinical_all, GSE39582 = all_GSE39582,GSE156915=all_GSE156915),
                             out_file = "MASM_sample")
g <- plot_mutation_mirror_modules(gres, hres, out_file = "MASM__module")

all_mss<-hres$data

p_cms_sub <- plot_pca_biplot(all_data_bulk,
                          color_col = "CMS",
                          title = "Human samples",size_large = 5.5,
                          custom_colors = my_subtype_colors,feature_cols = paste0("Tumor_MASM", c(1,5,7)))

p_cms_mouse <- plot_pca_biplot(all_temp,
                          color_col = "mCMS", flip_pc_x = T,
                          title = "GEMM mouse samples",size_large = 4.5,
                          custom_colors = my_subtype_colors,feature_cols = paste0("MASM", c(1,5,7)))

ggarrange(p_cms_sub ,p_cms_mouse)

p_scatter_hm=ggplot(res_cor_HM$data[res_cor_HM$data$CMS=="CMS4",],
       aes(x=vx,y=vy,label=module))+
    geom_point(size=4)+
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color="#2166AC", linewidth=1)+
    ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20)+
    theme_hd_minimal_plain()+
    stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95)+
    xlab("Mean MASM score in Human")+
    ylab("Mean MASM score in Mouse")+ggtitle("CMS4 subtype")

plot_data <- g$plot$data %>% filter(module == "MASM1")
main_cohorts <- c("Mouse GEMM GSE218776", "Human GSE39582")
pal = c(Apc_driven = "#56B4E9", MAPK_WT = "#56B4E9",
APC_only   = "#0072B2",
Kras_driven = "#009E73", KRAS_mut = "#009E73",
Braf_driven = "#E69F00", BRAF_mut = "#E69F00")
plist_main <- base::lapply(base::seq_along(main_cohorts), function(i) {
    ct <- main_cohorts[i]
    make_cohort_plot(
        d = plot_data %>% dplyr::filter(cohort == ct),
        pal = pal,
        show_ylab = (i == 1)
    )
})
fig6f_main <- patchwork::wrap_plots(plist_main, nrow = 1)

sc_all <- compute_sc_all_celltypes(list(scEpi = epi_all_tumor, scEpi2=GSE200997_epi_tumor,scnon_epi=no_epi_all_tumor_pseudo,epi=temp_epi_data), all_loading_list)

cor_plus <- bind_rows(pw_sys$cor, sc_all$scEpi$pathway_cor, sc_all$scEpi2$pathway_cor, sc_all$scnon_epi$pathway_cor,sc_all$epi$pathway_cor)
PP<-plot_masm_pathway_heatmap(cor_plus, thresh = 0.4,
out_file = "results/Fig_Hallmark_4systems.pdf", width = 10)


PP<- plot_masm_pathway_core(cor_plus,   sys_map = c(scEpi     = "scEpi_discovery",
                                               scEpi2    = "scEpi_verification",
                                               Epi       = "Bulk_epi",
                         human     = "human_bulk",mouse     = "Mouse",CCLE      = "CCLE",
                                               scnon_epi = "sc_nonEpi"),thresh = 0.4,
                            exclude_paths = c("Custom_EMT_TF"),
                            out_file = "Fig_masm_pathway_core.pdf");PP$ht


p_pca_ccle<-ggarrange(p_icms_sub,p_icms_sub_ccle,nrow=1,ncol = 2,legend = "bottom",common.legend = T,labels = c(letters[2:3]), font.label = panel_label_style)

pw_sys <- run_pathway_signature_analysis(res, exp_matrix, sc_ccle,
                                        all_data_bulk = all_data_bulk,
                                          bulk_expr_list = expr_list_full)


p_pca_mouse<-ggarrange(p_cms_sub,p_cms_mouse,nrow=1,ncol = 2,legend = "bottom",common.legend = T,labels = c(letters[4:5]), font.label = panel_label_style)
p_pca<-ggarrange(p_pca_ccle,p_pca_mouse,nrow=2)
pall6_1<-ggarrange(ht1, p_pca, nrow = 1, ncol = 2,labels = c(letters[1],""), font.label = panel_label_style,widths = c(1.5,2))
pall6_2<-ggarrange(p_scatter_hm,ht3,plist_main[[1]],plist_main[[2]],nrow=1,ncol=4,labels = c(letters[6:9]), font.label = panel_label_style,widths = c(1,1.5,0.75,0.75))
ht_pathway <- grid.grabExpr(ComplexHeatmap::draw(PP$ht))
pall6_3<-ggarrange(pall6_1,pall6_2,ht_pathway,nrow=3,labels = c("","",letters[10]), font.label = panel_label_style,heights = c(1.5,1,1.2))


save_tiff(pall6_3, "figure6.tiff", width = 19, height = 20, path = fig_dir)


 smi_sys <- run_smi_across_systems(all_data_bulk, res, exp_matrix, all_loading_list)
 pds_sys <- run_pds_across_systems(all_data_bulk, res, exp_matrix, all_loading_list)

PDS_call18 = get_pds(res$vst18, "mouse", "GEMM E-MTAB-6363")
PDS_call39 = get_pds(res$vst39, "mouse", "GEMM GSE218776")

all_temp_pds<-rbind(cbind(dataset="GEMM218776",PDS_call39 ),cbind(dataset="GEMM6363",PDS_call18))


p_pds_mouse <- plot_pca_biplot(all_temp_pds,
                          color_col = "PDS_call", flip_pc_x = T,
                          title = "GEMM mouse samples",size_large = 4.5,
                          custom_colors = my_subtype_colors,feature_cols = paste0("MASM", c(1,5,7)))

p_pds_sub <- plot_pca_biplot(all_data_bulk,
                          color_col = "PDS_call",
                          title = "Human samples",size_large = 5.5,
                          custom_colors = my_subtype_colors,feature_cols = paste0("Tumor_MASM", c(1,5,7)))


ggarrange(p_pds_sub ,p_pds_mouse)


res_h <- run_masm_pds(all_data_bulk)


res_m <- run_masm_pds(all_temp_pds,min_group_n = 1)

human_mean <- masm_group_means(all_data_bulk, "PDS_call")
mouse_mean <- masm_group_means(all_temp_pds,     "PDS_call")

res_cor_pds <- plot_pds_transfer(mouse_mean, human_mean,
                         out_file = "Fig_PDS_transfer.pdf")
res_cor_pds$cor
res_cor_pds$r_all

p_scatterPDS=ggplot(res_cor_pds$data[res_cor_pds$data$PDS=="PDS1" & res_cor_pds$data$module %in%c("MASM1","MASM5","MASM7","MASM6"),],
       aes(x=vx,y=vy,label=module))+
    geom_point(size=4)+
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color="#2166AC", linewidth=1)+
    ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20)+
    theme_hd_minimal_plain()+
    stat_cor(method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95)+
    xlab("Mean MASM score in Human")+
    ylab("Mean MASM score in Mouse")+ggtitle("PDS subtype")


mouse_mean18 <- masm_group_means(PDS_call18,     "PDS_call",min_group_n = 1)
mouse_mean39 <- masm_group_means(PDS_call39,     "PDS_call",min_group_n = 1)

res_cor_pdsMM <- plot_pds_transfer(mouse_mean18, mouse_mean39,
                         out_file = "Fig_PDS_transfer.pdf")


p_bar_pds<-plot_class_cor_bar(res_cor_pds, res_cor_pdsMM,
                   x_col = "PDS", x_levels = paste0("PDS", 1:3))


pp=plot_state_trajectory(
    human_mean = human_mean,
    mouse_list = list(GSE18 = mouse_mean18, GSE39 = mouse_mean39),
    label_ys   = c("GEMM E-MTAB-6363", "GEMM GSE218776")
)

pp$plot

mouse_mean_cms<-x1_mouse[,-1]
rownames(mouse_mean_cms)<-x1_mouse[,1]

human_mean_cms<-x2_human[,-1]
rownames(human_mean_cms)<-x2_human[,1]


x18_mouse<-x18[,-1]
rownames(x18_mouse)<-x18[,1]
x39_mouse<-x39[,-1]
rownames(x39_mouse)<-x39[,1]
p_cms_HM <- plot_state_trajectory( human_mean = human_mean_cms,
mouse_list = list(`GEMM E-MTAB-6363` = x18_mouse,`GEMM GSE218776`= x39_mouse),
group_levels = paste0("CMS", 1:4),
group_name   = "CMS",
group_label  = "Consensus Molecular Subtype (CMS)",
out_file = "Fig_CMS_trajectory")
p_cms_HM$plot


all_data_bulk$dataset[all_data_bulk$dataset=="silu_colon"]="AC-ICAM"
metab <- run_metab_across_systems(all_data_bulk, res, exp_matrix, all_loading_list)


metab$plot


sup_cohorts <- c("Mouse GEMM GSE6363", "Human TCGA","Human GSE156915"  )
plot_data$cohort<-as.character(plot_data$cohort)
plist_sup <- base::lapply(base::seq_along(sup_cohorts), function(i) {
    ct <- sup_cohorts[i]
    make_cohort_plot(
        d = plot_data %>% dplyr::filter(cohort == ct),
        pal = pal,
        show_ylab = (i == 1)
    )
})
plist_sup[[1]]=plist_sup[[1]] +ggtitle("Mouse GEMM E-MTAB-6363")
fig6f_sup <- patchwork::wrap_plots(plist_sup, nrow = 1)


pall131 <- ggarrange(p_pds_sub ,p_pds_mouse,p_bar,
                    labels = letters[1:3], nrow = 1, ncol = 3,
                    font.label = panel_label_style)
pall132 <- ggarrange(p_cms_HM$plot,pp$plot,nrow = 2, ncol = 1,labels = c(letters[4:5]),
                    font.label = panel_label_style,common.legend = T,legend = "bottom")

pall133 <- ggarrange(pall131,pall132, fig6f_sup,metab$plot_full,nrow = 4, ncol = 1,labels = c("","",letters[6:7]),
                    font.label = panel_label_style,heights = c(1,2,1,1))

save_tiff(pall133, "supplementaryfigure13.tiff", width = 18, height = 22, path = fig_dir)


mod_coh<-run_weighted_module_audit(epi_all_tumor, all_loading_list)

mod_coh_ccle<-run_weighted_module_audit(t(exp_matrix), all_loading_list)
mod_coh_tcga<-run_weighted_module_audit(t(tcga_coread), all_loading_list)

mod_coh_hs39<-run_weighted_module_audit(t(res$hs39), all_loading_list)
mod_coh_hs18<-run_weighted_module_audit(t(res$hs18), all_loading_list)
for( x in expr_list_full){
   audtc <- run_module_coherence(x, all_loading_list)
}
mod_coh<-run_module_coherence(t(epi_all_tumor), all_loading_list)

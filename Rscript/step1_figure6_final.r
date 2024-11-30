# Load necessary data files
load("./data/load_all_files2.rdata")
setwd("/data/Project/Project1/result_new")
source("./Rscript/functions_plot.r")
patient_infor <- readRDS(file="./data/single_cell_patient_infor.rds")
kegg_path_genes <- readRDS(file="./data/all_kegg_pathway_genes.rds")

# Filter KEGG pathways related to signaling
kegg_signal = kegg_path_genes[grepl("signal", kegg_path_genes$description),]
pathwayRatioNormalized <- readRDS(file="./pathwayRatioNormalized.rds")
fc_pathwayRatioNormalized <- getFCresult(pathwayRatioNormalized, meta_data)

# Load combined single-cell data
combined_data <- readRDS(file="./data/combined_data_single_cell.rds")

# Load CPM data for epithelial and pseudo tumor samples
cpm_epi_tumor_icms_data <- readRDS(file="./data/icms_SCT_count.rds")
cpm_pseudo_tumor_icms_data <- readRDS(file="./data/cpm_pseudo_tumor_icms_data.rds")

# Rename columns for epithelial data
colnames(cpm_epi_tumor_icms_data) <- paste("Epithelial", colnames(cpm_epi_tumor_icms_data))

# Combine epithelial and pseudo tumor data
combined_data_icms <- as.data.frame(cbind(cpm_epi_tumor_icms_data, cpm_pseudo_tumor_icms_data$mean_x[rownames(cpm_epi_tumor_icms_data),]))

# Combine with existing data
combined_data_all <- cbind(combined_data, combined_data_icms[rownames(combined_data),])
combined_data_all[is.na(combined_data_all)] = 0

# Define signaling module labels
labels_x = c("CGS", "TED", "POS", "ACP", "MIC")

# Select relevant columns for analysis
combined_data_all1 <- combined_data_all[, c("Pseudo iCMS2", "Pseudo iCMS3", "Epithelial iCMS2", "Epithelial iCMS3", 
                                            paste("Pseudo", labels_x, sep=" "), paste("Epithelial", labels_x, sep=" "))]

# Initialize result storage for correlation analysis
result_cor_icms <- c()

# Calculate Spearman correlation between signaling modules and iCMS scores
for (i in labels_x) {
  for (j in c("iCMS2", "iCMS3")) {
    temp_r <- cor.test(combined_data_all[, paste("Epithelial", i)], combined_data_all[, paste("Epithelial", j)], method="spearman")
    result_cor_icms <- rbind(result_cor_icms, c(type=i, type2=j, rho=round(temp_r$estimate, 2), pvalue=temp_r$p.value))
  }
}

# Convert results to a data frame and format
result_cor_icms <- as.data.frame(result_cor_icms)
result_cor_icms$rho.rho <- as.numeric(result_cor_icms$rho.rho)
result_cor_icms$pvalue <- as.numeric(result_cor_icms$pvalue)
result_cor_icms$type <- factor(result_cor_icms$type, levels=result_cor_icms[result_cor_icms$type2=="iCMS2",]$type[order(result_cor_icms[result_cor_icms$type2=="iCMS2",]$rho.rho)])

# Create bar plot for correlation results
p_cor_icms = ggplot(result_cor_icms, aes(x=type, y=rho.rho)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_hd_plain() +
  ylab("Spearman's rho") +
  geom_text(aes(label=rho.rho), size=3.5) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  facet_wrap(~type2) +
  xlab("Metabolism-regulating signaling modules")

# Create scatter plots for TED signaling module against iCMS scores
p1_ted_icms <- ggplot(combined_data_all1, aes(x=`Epithelial TED`, y=`Epithelial iCMS2`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  xlab("Average log2FC of\nTED signaling module") +
  ggtitle("Epithelial cells") +
  ylab("iCMS2 score in epithelial\ncells of tumor samples")

p1_ted_icms2 <- ggplot(combined_data_all1, aes(x=`Epithelial TED`, y=`Epithelial iCMS3`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  xlab("Average log2FC of\nTED signaling module") +
  ggtitle("Epithelial cells") +
  ylab("iCMS3 score in epithelial\ncells of tumor samples")

# ROC analysis for iCMS2
roc_obj = roc(patient_infor[rownames(combined_data_all), "MSI"], combined_data_all[, "Epithelial iCMS2"], plot=FALSE, print.auc=FALSE)
proc_icms2 = ggroc(roc_obj) +
  annotate(geom="segment", x=1, y=0, xend=0, yend=1, linetype="dashed", color="grey") +
  annotate("text", x=0.5, y=0.15, label=paste("AUC =", round(roc_obj$auc, 3)), size=5) +
  ggtitle(paste("iCMS3 score in epithelial cells of tumor samples\n", "(MSI VS MSS)")) +
  theme_hd_minimal_plain() +
  xlab("Specificity") +
  ylab("Sensitivity") +
  theme(plot.margin=unit(c(1, 1, 1, 1), "cm"))

# Plot ROC for MSI
p_msi = PlotROCclinical(combined_data_all1, patient_infor, "MSI")
ggarrange(p_msi$p_bar, p_msi$p_roc)

# Plot ROC for sideness
p_sideness = PlotROCclinical(combined_data_all1, patient_infor, "sideness")
ggarrange(p_sideness$p_bar, p_sideness$p_roc)

# Combine patient information with combined data
combined_data_infor <- cbind(combined_data, patient_infor[rownames(combined_data),])

# Scatter plot for TED scores in epithelial vs pseudo data
p_ted_msi <- ggscatter(combined_data_infor, x="Epithelial TED", y="Pseudo TED", col="MSI") +
  stat_cor(method="spearman", size=5) +
  geom_abline(slope=1, linetype="dashed") +
  scale_color_manual(name="Type", values=c("#00BFFF", "#EE7621"))

# Prepare data for box plots comparing TED scores by MSI
combined_data_infor2 = rbind(cbind(combined_data_infor[, c("MSI", "sideness")], ted=combined_data_infor[, c("Epithelial TED")], type="Epithelial cells"),
                             cbind(combined_data_infor[, c("MSI", "sideness")], ted=combined_data_infor[, c("Pseudo TED")], type="Pseudo-bulk data"))

p_box_msi1 <- ggplot(combined_data_infor2, aes(x=MSI, y=ted, fill=MSI)) +
  geom_boxplot(width=0.3) +
  stat_compare_means(method="wilcox") +
  theme_hd_plain() +
  scale_fill_manual(name="Side", values=c("#00BFFF", "#EE7621")) +
  xlab("MSI") +
  theme(legend.position="none") +
  ylab("Average log2FC of\nTED signaling module") +
  facet_wrap(~type)

# Prepare data for box plots comparing MIC scores by sideness
combined_data_infor2 = rbind(cbind(combined_data_infor[, c("MSI", "sideness")], mic=combined_data_infor[, c("Epithelial MIC")], type="Epithelial cells"),
                             cbind(combined_data_infor[, c("MSI", "sideness")], mic=combined_data_infor[, c("Pseudo MIC")], type="Pseudo-bulk data"))

p_box_side1 <- ggplot(combined_data_infor2, aes(x=sideness, y=mic, fill=sideness)) +
  geom_boxplot(width=0.3) +
  stat_compare_means(method="wilcox") +
  theme_hd_plain() +
  scale_fill_manual(name="Side", values=c("#FFAEB9", "#B0E2FF")) +
  xlab("Side") +
  theme(legend.position="none") +
  ylab("Average log2FC of\nMIC signaling module") +
  facet_wrap(~type)

# Scatter plot for TED scores based on sideness
p_ted_side <- ggscatter(combined_data_infor, x="Epithelial TED", y="Pseudo TED", col="sideness") +
  stat_cor(method="spearman", size=5) +
  geom_abline(slope=1, linetype="dashed") +
  scale_color_manual(name="Type", values=c("#00BFFF", "#EE7621"))

# Load signal pathway data
signal_path <- readRDS(file="./data/signal_path.rds")

# Create heatmap for correlation coefficients
bk <- c(seq(-0.8, 0.8, by=0.01))
a = pheatmap(signal_path$all_select$Epithelial, 
             col=c(colorRampPalette(colors=c("#00BFFF", "white"))(length(bk)/2), 
                   colorRampPalette(colors=c("white", "#EE7621"))(length(bk)/2)),
             legend_breaks=seq(-0.8, 0.8, 0.1), 
             breaks=bk, 
             main="Spearman's correlation coefficient\nLog2FC(pathway) vs Log2FC(signal genes)", 
             fontsize_row=6, 
             fontsize_col=4)

# Cluster assignments for heatmap
labels_x = c("CGS", "TED", "POS", "ACP", "MIC")
cluster_assignments <- cutree(a$tree_col, k=5)

# Hierarchical clustering for signal pathways
data_x = signal_path$all_select$Epithelial
row_km <- cutree(hclust(dist(data_x)), k=5)

# Load TCGA results
tcga_result <- readRDS(file="./data/tcga_result.Rds")
cpm_tcga_all <- tcga_result$count

# Process sample identifiers
temp_all = as.data.frame(str_split_fixed(rownames(cpm_tcga_all), "-", 5))
temp_all$sample = paste(temp_all[,1], temp_all[,2], temp_all[,3], substr(temp_all[,4], 1, 2), sep=".")
temp_all$patient = paste(temp_all[,1], temp_all[,2], temp_all[,3], sep="-")
temp_all$type = ifelse(temp_all[,4] == "11A", "Normal", ifelse(temp_all[,4] %in% c("01A", "01B", "01C"), "Tumor", "others"))
rownames(cpm_tcga_all) <- paste(temp_all$type, temp_all$patient, sep="_")

# Filter signal genes from TCGA data
cpm_tcga_signal_all = cpm_tcga_all[, colnames(cpm_tcga_all) %in% kegg_signal$symbol & !colnames(cpm_tcga_all) %in% unique(all_gene_pathway$gene)]
temp <- aggregate(cpm_tcga_signal_all, list(samples=rownames(cpm_tcga_signal_all)), mean)
rownames(temp) <- temp[,1]
cpm_tcga_signal_all <- temp[,-1]

# Calculate fold change for TCGA signal data
cpm_tcga_signal_fcAN_all <- getFCdataSinglelogAllnormal(cpm_tcga_signal_all)
cpm_tcga_signal_fcAN_each <- getFCdataSinglelog(cpm_tcga_signal_all)

# Get fold change results for TCGA
cpm_TCGA_fc_allnorm <- GetSignalResultFCBulk(cpm_tcga_signal_fcAN_all, cluster_assignments, paste("TCGA", labels_x))
cpm_TCGA_fc_alleach <- GetSignalResultFCBulk(cpm_tcga_signal_fcAN_each, cluster_assignments, paste("TCGA", labels_x))

# Select DPM for TCGA
tcga_dpm <- cpm_tcga_signal_fcAN_all[, colnames(cpm_tcga_signal_fcAN_all) %in% names(cluster_assignments[cluster_assignments == 2])]

# Transpose fold change data for TCGA
temp_x = as.data.frame(t(cpm_TCGA_fc_alleach$mean_x))

# Create scatter plots for TCGA TED and MIC modules
p1_tcga_ted_mic <- ggplot(temp_x, aes(x=`TCGA TED`, y=`TCGA MIC`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5, label.y=1) +
  theme_hd_plain() +
  xlab("Average log2FC of\nTED signaling module") +
  ylab("Average log2FC of\nMIC signaling module") +
  ggtitle("TCGA cohort")

p1_tcga_pos_acp <- ggplot(temp_x, aes(x=`TCGA ACP`, y=`TCGA POS`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5, label.y=1) +
  theme_hd_plain() +
  xlab("Average log2FC of\nACP signaling module") +
  ylab("Average log2FC of\nPOS signaling module") +
  ggtitle("TCGA cohort")

# Create correlation heatmap for TCGA
cor_heatmap <- pheatmap(cor(temp_x), col=c(colorRampPalette(colors=c("#00BFFF", "white", "#EE7621"))(100)))

# Load pseudo data
pseudo_all_df_count <- readRDS(file="./data/pseudo_all_df_count.Rds")
groups = str_split_fixed(rownames(pseudo_all_df_count), "_", 4)[,2]
temp_pseudo <- getEdgeRdata(t(pseudo_all_df_count), groups)
cpm_pseudo <- temp_pseudo$count

# Load pathway information
all_gene_pathway <- readRDS(file="./data/all_genes_pathway_infor_select.rds")
tcga_pathway <- getBulkdataRatioPathway(cpm_tcga_all, all_gene_pathway)
pseudo_pathway <- getBulkdataRatioPathway(cpm_pseudo, all_gene_pathway)

# Combine TCGA and pseudo data for analysis
combined_tcga_signal_path_each <- cbind(tcga_pathway$cpm_Bulk_pathwayratio_fcAN_each, t(cpm_TCGA_fc_alleach$mean_x[, rownames(tcga_pathway$cpm_Bulk_pathwayratio_fcAN_each)]))
cpm_pseudo_pathwayratio_fcAN_each <- pseudo_pathway$cpm_Bulk_pathwayratio_fcAN_each
rownames(cpm_pseudo_pathwayratio_fcAN_each) <- str_split_fixed(rownames(cpm_pseudo_pathwayratio_fcAN_each), "_", 2)[,1]
pseudo_signal_path_each <- cbind(cpm_pseudo_pathwayratio_fcAN_each[, rownames(signal_path$all_select$Epithelial)], combined_data[rownames(cpm_pseudo_pathwayratio_fcAN_each), grepl("Pseudo", colnames(combined_data))])

# Calculate correlations between TCGA and pseudo data
x = cor(combined_tcga_signal_path_each, method="spearman")
x = x[!grepl("TCGA", rownames(x)), grepl("TCGA", rownames(x))]

# Calculate correlations for pseudo data
x_epi_all <- cor(pseudo_signal_path_each, method="spearman")
x_epi_all <- x_epi_all[!grepl("Pseudo", rownames(x_epi_all)), grepl("Pseudo", rownames(x_epi_all))]

# Combine correlation results for plotting
x_cor_all <- as.data.frame(cbind(x[rownames(x_epi_all),], x_epi_all))

# Create scatter plots for correlations in TCGA
p1_tcga_epi_ted <- ggplot(x_cor_all, aes(x=`Pseudo TED`, y=`TCGA TED`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  xlab("Local metabolic allocation\nVS TED module(pseudo-bulk data)") +
  ylab("Local metabolic allocation\nVS TED module(TCGA cohort)") +
  ggtitle("Spearman's correlation of log2FC(T vs N)")

p1_tcga_epi_mic <- ggplot(x_cor_all, aes(x=`Pseudo MIC`, y=`TCGA MIC`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  xlab("Local metabolic allocation\nVS MIC module(pseudo-bulk data)") +
  ylab("Local metabolic allocation\nVS MIC module(TCGA cohort)") +
  ggtitle("Spearman's correlation of log2FC(T vs N)")

# Initialize storage for correlation results between TCGA and pseudo data
result_cor_tcga <- c()

# Calculate Spearman correlation for each signaling module
for(i in labels_x) {
  temp_r <- cor.test(x_cor_all[, paste("Pseudo", i)], x_cor_all[, paste("TCGA", i)], method="spearman")
  result_cor_tcga <- rbind(result_cor_tcga, c(type=i, rho=round(temp_r$estimate, 2), pvalue=temp_r$p.value))
}

# Convert results to a data frame and format
result_cor_tcga <- as.data.frame(result_cor_tcga)
result_cor_tcga$rho.rho <- as.numeric(result_cor_tcga$rho.rho)
result_cor_tcga$pvalue <- as.numeric(result_cor_tcga$pvalue)
result_cor_tcga$type <- factor(result_cor_tcga$type, levels=result_cor_tcga$type[order(result_cor_tcga$rho.rho)])

# Create bar plot for TCGA correlation results
p_cor_tcga = ggplot(result_cor_tcga, aes(x=type, y=rho.rho)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_hd_plain() +
  ggtitle("Metabolic changes under signaling\nmodule regulation (TCGA cohort vs pseudo-bulk data)") +
  ylab("Spearman's rho") +
  geom_text(aes(label=rho.rho), vjust=0.1, size=3.5) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Metabolism-regulating signaling modules")

# Scatter plot for TED signaling module vs leukotriene metabolism
p_tcga_ted_path <- ggplot(as.data.frame(combined_tcga_signal_path_each), aes(x=`TCGA TED`, y=`Leukotriene metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\nleukotriene metabolism") +
  xlab("TED signaling module") +
  ggtitle("Log2FC (T VS N) in TCGA cohort")

# Scatter plot for MIC signaling module vs N-glycan metabolism
p_tcga_ted_path2 <- ggplot(as.data.frame(combined_tcga_signal_path_each), aes(x=`TCGA MIC`, y=`N-glycan metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\nN-glycan metabolism") +
  xlab("MIC signaling module") +
  ggtitle("Log2FC (T VS N) in TCGA cohort")

# Select epithelial-related fold changes
epi_select = fc_pathwayRatioNormalized$FC[grepl("Epi", rownames(fc_pathwayRatioNormalized$FC)), rownames(signal_path$all_select$Epithelial)]
rownames(epi_select) <- str_split_fixed(rownames(epi_select), "_", 3)[, 3]

# Combine data for analysis
combined_data_path <- cbind(combined_data, epi_select[rownames(combined_data),])

# Scatter plot for Epithelial MIC vs N-glycan metabolism
p1_signal_path_epi <- ggplot(combined_data_path, aes(x=`Epithelial MIC`, y=`N-glycan metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\nN-glycan metabolism") +
  xlab("MIC signaling module") +
  ggtitle("Log2FC (T VS N) in epithelial cells")

# Scatter plot for Pseudo MIC vs N-glycan metabolism
p1_signal_path_pseudo = ggplot(pseudo_signal_path_each, aes(x=`Pseudo MIC`, y=`N-glycan metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\nN-glycan metabolism") +
  xlab("MIC signaling module") +
  ggtitle("Log2FC (T VS N) in pseudo-bulk data")

# Load clinical data
clinical_molecular_public_all <- read.delim("./data/clinical_molecular_public_all.txt")
clinical_tcga = clinical_molecular_public_all[clinical_molecular_public_all$sample %in% as.character(cpm_TCGA_fc_alleach$df_x$Var1), ]
rownames(clinical_tcga) <- clinical_tcga$sample

# Combine clinical data with TCGA fold change data
tcga_df <- cpm_TCGA_fc_alleach$df_x[as.character(cpm_TCGA_fc_alleach$df_x$Var1) %in% rownames(clinical_tcga),]
tcga_df_infor <- cbind(tcga_df, clinical_tcga[as.character(tcga_df$Var1),])
tcga_df_infor$msi = toupper(tcga_df_infor$msi)

# Boxplot for TCGA ICMS scores by MSI status
p_tcga_icms <- ggplot(tcga_df_infor[!is.na(tcga_df_infor$msi),], aes(fill=msi, y=value, x=msi)) +
  geom_boxplot() +
  facet_wrap(~Var2, nrow=1) +
  stat_compare_means(label="p.signif", method="wilcox.test") +
  scale_fill_manual(name="Type", values=c("#00BFFF", "#EE7621")) +
  ylab("Average of Log2FC in TCGA") +
  theme_hd_plain() +
  theme(axis.title.x=element_blank(), legend.position="none", strip.background=element_blank(), axis.text.x=element_text(angle=90, hjust=1))

# Boxplot for TCGA ICMS scores by BRAF mutation status
tcga_df_infor$braf_mut2 = ifelse(tcga_df_infor$braf_mut == 0, "Without", "With")
p_tcga_braf <- ggplot(tcga_df_infor[!is.na(tcga_df_infor$braf_mut2),], aes(fill=as.factor(braf_mut2), y=value, x=as.factor(braf_mut2))) +
  geom_boxplot() +
  facet_wrap(~Var2, nrow=1) +
  stat_compare_means(label="p.signif", method="wilcox.test") +
  scale_fill_manual(name="Type", values=c("#00BFFF", "gray")) +
  ylab("Average of Log2FC in TCGA") +
  theme_hd_plain() +
  theme(legend.position="none", strip.background=element_blank(), axis.text.x=element_text(angle=90, hjust=1)) +
  xlab("BRAF mutation")

# Boxplot for TCGA scores by CMS label
p_tcga_cms <- ggplot(tcga_df_infor[tcga_df_infor$cms_label != "NOLBL",], aes(fill=Var2, y=value, x=Var2)) +
  geom_boxplot() +
  facet_wrap(~cms_label, nrow=1) +
  ylab("Average of Log2FC in TCGA") +
  theme_hd_plain() +
  theme(axis.title.x=element_blank(), legend.position="none", strip.background=element_blank(), axis.text.x=element_text(angle=90, hjust=1)) +
  scale_fill_manual(name="Class", values=c("#9AFF9A", "#AEEEEE", "#EECFA1", "#FFBBFF", "pink")) +
  geom_hline(yintercept=0, linetype="dashed") +
  stat_compare_means(method="kruskal.test", label="p.signif", position=position_nudge(x=1))

# Boxplot for TCGA scores by CIMP status
p_tcga_cimp <- ggplot(tcga_df_infor[!is.na(tcga_df_infor$cimp),], aes(fill=cimp, y=value, x=cimp)) +
  geom_boxplot() +
  facet_wrap(~Var2, nrow=1) +
  ylab("Average of Log2FC in TCGA") +
  theme_hd_plain() +
  theme(axis.title.x=element_blank(), legend.position="none", strip.background=element_blank(), axis.text.x=element_text(angle=90, hjust=1)) +
  scale_fill_manual(name="Class", values=c("steelblue", "#AEEEEE", "gray")) +
  geom_hline(yintercept=0, linetype="dashed") +
  stat_compare_means(method="kruskal.test", label="p.signif", position=position_nudge(x=1))

# Combine plots for final figure
proc_ted <- p_msi$p_roc + ggtitle(paste("Average log2FC of TED signaling module\n", "(MSI VS MSS)"))
p6all <- ggarrange(p_cor_icms, p1_ted_icms, p_box_msi1, p_msi$p_bar, proc_ted, p_tcga_icms, p_tcga_braf, p_tcga_cimp, p_tcga_cms, nrow=3, ncol=3, labels=letters[c(1:9)], font.label=list(size=18, color="black", face="bold"))

# Save combined figure as TIFF
tiff(paste("./figures/", "figure6new.tiff", sep=""), height=20, width=20, res=300, units="in", compression="lzw")
print(p6all)
dev.off()

# Save combined figure as PDF
pdf(paste("./figures/", "figure6new.pdf", sep=""), height=20, width=20)
print(p6all)
dev.off()

# Prepare additional plots for supplementary figure
proc_mic <- p_sideness$p_roc + ggtitle(paste("Average log2FC of MIC signaling module\n", "(Sideness: Right VS Left)"))
p6all_sub <- ggarrange(p1_ted_icms2, proc_icms2, p_box_side1, p_sideness$p_bar, proc_mic, as_ggplot(cor_heatmap$gtable), p1_tcga_ted_mic, p1_tcga_pos_acp, p_cor_tcga, p1_tcga_epi_ted, p1_tcga_epi_mic, p_tcga_ted_path, p_tcga_ted_path2, p1_signal_path_pseudo, p1_signal_path_epi, nrow=5, ncol=3, labels=letters[c(1:15)], font.label=list(size=18, color="black", face="bold"))

# Save supplementary figure as TIFF
tiff(paste("./figures/", "supfigure9.tiff", sep=""), height=20, width=20, res=300, units="in", compression="lzw")
print(p6all_sub)
dev.off()


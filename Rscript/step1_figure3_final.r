# Load necessary data files
load("./data/load_all_files2.rdata")
setwd("/data/Project/Project1/result_new")
source("./Rscript/functions_plot.r")

# Extract metadata from the filtered dataset
meta_data <- all_dataset_SCT_filter@meta.data

# Load pathway mean data
allcells_pathwayMeanSCTMainClass <- readRDS(file="./data/allcells_pathwayMeanSCTMainClass.Rds")
allcells_pathwayMeanSCTMainClassNormalized <- readRDS(file="./data/allcells_pathwayMeanSCTMainClassNormalizedByMetacount.Rds")
pathwayRatioNormalized <- readRDS(file="./pathwayRatioNormalized.rds")

# Calculate fold change results for various pathway datasets
fc_pathwayRatioNormalized <- getFCresult(pathwayRatioNormalized, meta_data)
fc_mainclass <- getFCresult(allcells_pathwayMeanSCTMainClass, meta_data)
fc_mainclassNormalized <- getFCresult(allcells_pathwayMeanSCTMainClassNormalized, meta_data)

# Plot heatmap for normalized fold change results
gb_heatmapmainclassNormalized <- plotFChetamap(t(fc_mainclassNormalized$FC_sum_sig), "Log2FC of global metabolic allocation")
ggarrange(gb_heatmapmainclassNormalized)

# Prepare data for boxplot visualization
data_x = fc_mainclass$FC_sum_sig
df_data_x <- reshape2::melt(as.matrix(data_x))  
data_xNormalized = fc_mainclassNormalized$FC_sum_sig
df_data_xNormalized <- reshape2::melt(as.matrix(data_xNormalized))  
sort_name = colnames(fc_mainclassNormalized$FC_sum_sig[order(fc_mainclassNormalized$FC_sum_sig["Epithelial",])])

# Combine datasets for plotting
data_all <- rbind(cbind(df_data_x, type="Average pathway\nactivity"), cbind(df_data_xNormalized, type="Global metabolic\nallocation"))
data_all$Var2 <- factor(data_all$Var2, levels=sort_name)

# Create boxplot for fold changes
p_fc <- ggplot(data=data_all, aes(x=Var1, y=value, fill=type)) +
  geom_boxplot() + 
  theme_hd_plain() + 
  theme(axis.text.x = element_text(angle=60, hjust=1), axis.title.x = element_blank()) + 
  ylab("Log2FC(TvsN)") + 
  geom_hline(yintercept=0, linetype="dashed") + 
  scale_fill_manual(name="Type", values=c("#00BFFF", "#B2DFF0")) + 
  guides(fill=guide_legend(nrow=2))

# Load prediction data
pre_class <- readRDS("./prediction_data.rds")
all_reglobal = pre_class$all_pre[!grepl("_", pre_class$all_pre$path),]

# Plot global metabolic allocation predictions
p_NT_global <- plotPredictionResultNT(all_reglobal) + ylab("Global metabolic allocation")

# Load normalized prediction data
ratioNormalize_pre <- readRDS("ratioNormalize_prediction_all.rds")
all_re = ratioNormalize_pre$all_pre

# Plot local metabolic allocation predictions
p_NT_local <- plotPredictionResultNT(all_re) + ylab("Local metabolic allocation")
temp_local <- all_re[all_re$type == "N_Tall",]
temp_local$path = factor(temp_local$path, levels=temp_local$path[order(temp_local$accuracy)])

# Create bar plot for local metabolic allocation
pNTall_local = ggplot(temp_local, aes(x=path, y=accuracy)) +
  geom_bar(stat="identity", position="dodge", fill="steelblue") + 
  theme_hd_plain() + 
  ggtitle("Tumor vs. normal comparison in all cell types") + 
  geom_hline(yintercept=0.6, linetype="dashed") + 
  ylab("10-fold cross-validation \nclassifier accuracy") + 
  xlab("Local metabolic allocation") + 
  coord_flip()

# Aggregate mean AUC results by patient and cell type
mean_auc_rec2_class <- aggregate(t(allcells_pathwayMeanSCTMainClassNormalized[, rownames(meta_data)]), 
                                 list(dataset=meta_data$dataset, type=meta_data$class_new, 
                                      cell=meta_data$Cell_type, patient=meta_data$patient), 
                                 mean)

# Remove Mast cells from the dataset
mean_auc_rec2_class <- mean_auc_rec2_class[mean_auc_rec2_class$cell != "Mast",]
mean_auc_rec2_class$cell <- factor(mean_auc_rec2_class$cell, 
                                   levels=c("Plasma", "Epithelial", "Myeloid", "Fibroblast", "B cell", "T cell", "Endothelial"))

# Plot various metabolic pathways
p_metac1 = plotpathwayNT(mean_auc_rec2_class, "Nucleotide metabolism") + scale_y_log10() + ylab("Global metabolic allocation")
p_metac2 = plotpathwayNT(mean_auc_rec2_class, "Energy metabolism") + scale_y_log10() + ylab("Global metabolic allocation")
p_metac3 = plotpathwayNT(mean_auc_rec2_class[mean_auc_rec2_class$cell != "Mast", ], "Metabolism of cofactors and vitamins") + 
  scale_y_log10() + ylab("Global metabolic allocation") + 
  ggtitle("Metabolism of cofactors\nand vitamins")
p_metac4 = plotpathwayNT(mean_auc_rec2_class[mean_auc_rec2_class$cell != "Mast", ], "Glycan metabolism") + 
  scale_y_log10() + ylab("Global metabolic allocation")

# Arrange metabolic pathway plots into a single figure
ggarrange(p_metac1, p_metac2, p_metac4, p_metac3, nrow=1, ncol=4, common.legend=T, legend="bottom")

# Log transformation for normalization
fc_mainclassNormalized$N_FC_all$Normallog = log2(fc_mainclassNormalized$N_FC_all$Normal)

# Scatter plot for Myeloid cells
p2_scatter = ggscatter(fc_mainclassNormalized$N_FC_all[fc_mainclassNormalized$N_FC_all$cell == "Myeloid", ], 
                       x="Normallog", y="FC", add="reg.line", conf.int=TRUE, 
                       add.params=list(color="#EE7621", fill="lightgray")) +
  stat_cor(method="spearman") + 
  theme_hd_plain() + 
  ylab("Log2FC of global metabolic allocation\n(T vs N)") + 
  xlab("Global metabolic allocation in normal samples") + 
  ggtitle("Myeloid") + 
  geom_text(aes(label=pathway), angle=0) + 
  xlim(-14, -9)

# Calculate Spearman correlation for each cell type
mean_df_all <- fc_mainclassNormalized$N_FC_all
all_sig <- c()
for(cc in unique(mean_df_all$cell)) {
  temp_data_x = mean_df_all[mean_df_all$cell == cc, ]
  p_sig <- cor.test(temp_data_x[,"FC"], temp_data_x[,"Normal"], method="spearman")
  all_sig <- rbind(all_sig, cbind(cell=cc, sp_p=p_sig$p.value, sp_cor=p_sig$estimate))
}
all_sig <- as.data.frame(all_sig)
all_sig$sp_p = as.numeric(all_sig$sp_p)
all_sig$sp_cor = round(as.numeric(all_sig$sp_cor), 3)
all_sig$sp_sig = ifelse(all_sig$sp_p < 0.05, "Pvalue<0.05", "NS")

# Create bar plot for Spearman correlation significance
p_sig <- ggplot(all_sig, aes(x=cell, y=sp_cor, fill=sp_sig)) +
  geom_bar(stat="identity") + 
  theme_hd_plain() + 
  scale_fill_manual(name="Sig", values=c("gray", "#EE7621")) + 
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank()) + 
  ylab("Spearman's rho of global metabolic allocation\nbetween normal samples and log2FC (T vs N)") 

# Reshape and prepare data for average absolute log2 fold change
mean_auc_dfClass_FC_sum_dfNormalized <- reshape2::melt(as.matrix(abs(fc_mainclassNormalized$FC_sum)))
colnames(mean_auc_dfClass_FC_sum_dfNormalized) <- c("cell", "variable", "value")
mean_auc_dfClass_FC_sum_dfNormalized <- mean_auc_dfClass_FC_sum_dfNormalized[mean_auc_dfClass_FC_sum_dfNormalized$cell != "Mast",]

# Create boxplot for average absolute log2 fold change
p_average = ggplot(mean_auc_dfClass_FC_sum_dfNormalized, aes(x=cell, y=value, fill=cell)) +
  geom_boxplot() + 
  theme_hd_plain() + 
  scale_fill_manual(name="", values=brewer.pal(n=7, name="Set3")) + 
  theme(legend.position="none", axis.text.x=element_blank()) + 
  ylab("Average abs(log2(FC)) of global\nallocation (tumor vs normal)") + 
  xlab("Cell type")
print(p_average)

# Aggregate distance data for local allocation
diss_nt_sum <- aggregate(diss_nt, list(cell=str_split_fixed(rownames(diss_nt), "_", 3)[, 2]), mean)
diss_nt_sum_df <- reshape2::melt(diss_nt_sum, id.vars=c("cell"))

# Create boxplot for average distance of local allocation
p_distance = ggplot(diss_nt_sum_df, aes(x=cell, y=value, fill=cell)) +
  geom_boxplot() + 
  theme_hd_plain() + 
  scale_fill_manual(name="", values=brewer.pal(n=7, name="Set3")) + 
  theme(axis.text.x=element_blank()) + 
  ylab("Average distance of the local\nallocation (tumor Vs normal)") + 
  xlab("Cell type")
print(p_distance)

# Arrange average and distance plots into a single figure
ggarrange(p_average, p_distance, nrow=1, ncol=2, common.legend=T, legend="bottom")

# Create boxplot for specific metabolic pathways
p_distance2 <- PlotboxRaw2(c("Glycan metabolism", "Canonical amino acid metabolism", "Carbohydrate metabolism", "Lipid metabolism", "Metabolism of cofactors and vitamins"), 
                           cell_order=c("B cell", "T cell", "Myeloid", "Plasma", "Endothelial", "Fibroblast", "Epithelial"), 
                           t(diss_nt)) + 
  guides(fill=guide_legend(title="Cell type", nrow=2)) + 
  ylab("Distance of the local\nallocation (tumor Vs normal)")

# Aggregate pathway means by patient
pathway_mean_patient <- aggregate(t(as.matrix(pathway_mean_allsct_seu@assays$RNA@counts)), 
                                  list(dataset=pathway_mean_allsct_seu@meta.data$dataset, 
                                       type=pathway_mean_allsct_seu$class_new, 
                                       cell=pathway_mean_allsct_seu@meta.data$Cell_type, 
                                       patient=pathway_mean_allsct_seu$patient), 
                                  mean)

# Prepare dataframe for patient pathways
pathway_mean_df_patient <- pathway_mean_patient[, -c(1:4)]
rownames(pathway_mean_df_patient) <- paste(pathway_mean_patient[, 1], pathway_mean_patient[, 2], 
                                           pathway_mean_patient[, 3], pathway_mean_patient[, 4], sep="_")

# Create PCA for carbohydrate metabolism
pca_car <- PlotpcaRatio(c("Hyaluronan metabolism", "Glycolysis/gluconeogenesis", "Citric acid cycle", 
                          "Pentose phosphate pathway", "Butanoate metabolism"), 
                        t(pathway_mean_df_patient[grepl("Epithelial", rownames(pathway_mean_df_patient)), ]), 
                        "Carbohydrate metabolism") + 
  ggtitle("Local metabolic allocation\nCarbohydrate metabolism in normal tissue") + 
  guides(color=guide_legend(override.aes=list(size=3)))

# Prepare carbohydrate ratio data for PCA
carbo_ratio = pathwayRatioNormalizedpatient[grepl("Epithelial|Endothelial|Myeloid|Fibroblast", rownames(pathwayRatioNormalizedpatient)), 
                                            colnames(pathwayRatioNormalizedpatient) %in% unique(all_gene_pathway[all_gene_pathway$kegg_class == "Carbohydrate metabolism", 2])]
prc <- prcomp(carbo_ratio)

# Create biplot for PCA results
temp_x <- str_split_fixed(rownames(carbo_ratio), "_", 4)
colnames(temp_x) <- c("dataset", "class", "cell", "patient")
pc_result = plotPcaBiplot(carbo_ratio, temp_x, "Local metabolic allocation\ncarbohydrate metabolism", 
                          labels=c("Hyaluronan metabolism", "Glycolysis/gluconeogenesis", 
                                   "Citric acid cycle", "Pentose phosphate pathway", "Butanoate metabolism"))
col_x = brewer.pal(n=7, name="Set3")
col_x[2] = "#EEE685"
bi_plot_carbo = pc_result$bi_plot + scale_color_manual(name="Cell type", values=col_x[-c(1, 6, 7)])

# Prepare data for specific carbohydrate pathways
temp = str_split_fixed(rownames(carbo_ratio), "_", 4)
colnames(temp) <- c("dataset", "type", "cell", "patient")
carbo_ratio_all <- cbind(temp, carbo_ratio)

# Plot specific carbohydrate pathways
p_carbo1 <- plotpathwayNT(carbo_ratio_all, "Glycolysis/gluconeogenesis") + scale_y_log10() + ylab("Local metabolic allocation")
p_carbo2 <- plotpathwayNT(carbo_ratio_all, "Butanoate metabolism") + scale_y_log10() + ylab("Local metabolic allocation")

# Create boxplot for carbohydrate metabolism
p_box_carbo <- PlotboxFC2(c("Hyaluronan metabolism", "Glycolysis/gluconeogenesis", "Citric acid cycle", 
                            "Pentose phosphate pathway", "Butanoate metabolism"), 
                          t(fc_pathwayRatioNormalized$FC[grepl("Epithelial|Endothelial|Myeloid|Fibroblast", rownames(fc_pathwayRatioNormalized$FC)), ]), 
                          fc_pathwayRatioNormalized$all_path_sigClass[fc_pathwayRatioNormalized$all_path_sigClass$wcx_p_pair < 0.05, ], 
                          "Log2FC of local metabolic allocation") + 
  geom_hline(yintercept=0, linetype="dashed") + 
  guides(fill=guide_legend(nrow=2), color=guide_legend(title="Sig", nrow=2)) + 
  ggtitle("Carbohydrate metabolism") + 
  coord_flip() + 
  scale_fill_manual(name="Cell type", values=col_x[-c(1, 6, 7)]) + 
  xlab("Pathways")

# Prepare lipid ratio data for PCA
lipid_ratio = pathwayRatioNormalizedpatient[grepl("Plasma|Endothelial|Myeloid|Fibroblast", rownames(pathwayRatioNormalizedpatient)), 
                                            colnames(pathwayRatioNormalizedpatient) %in% unique(all_gene_pathway[all_gene_pathway$kegg_class == "Lipid metabolism", 2])]
prc <- prcomp(lipid_ratio)
temp_x <- str_split_fixed(rownames(lipid_ratio), "_", 4)
colnames(temp_x) <- c("dataset", "class", "cell", "patient")
pc_result = plotPcaBiplot(lipid_ratio, temp_x, "Local allocation in lipid metabolism")
col_x = brewer.pal(n=7, name="Set3")
col_x[2] = "#EEE685"
bi_plot_lipid = pc_result$bi_plot + scale_color_manual(name="Cell type", values=col_x[-c(1, 3, 7)])

# Create boxplot for lipid metabolism
p_box_lipid <- PlotboxFC2(c("Fatty acid synthesis", "Triacylglycerol synthesis", "Arachidonic acid metabolism", 
                            "Fatty acid degradation", "Sphingolipid metabolism", "Ether lipid metabolism"), 
                          t(fc_pathwayRatioNormalized$FC[grepl("Plasma|Endothelial|Myeloid|Fibroblast", rownames(fc_pathwayRatioNormalized$FC)), ]), 
                          fc_pathwayRatioNormalized$all_path_sigClass[fc_pathwayRatioNormalized$all_path_sigClass$wcx_p_pair < 0.05, ], 
                          "Log2FC of local metabolic allocation") + 
  xlab("Pathway") + 
  geom_hline(yintercept=0, linetype="dashed") + 
  guides(fill=guide_legend(nrow=2), color=guide_legend(title="Sig", nrow=2)) + 
  ggtitle("Lipid metabolism") + 
  coord_flip() + 
  scale_fill_manual(name="Cell type", values=col_x[-c(1, 3, 7)])

# Arrange multiple plots into a single figure
p1all <- ggarrange(p_fc, gb_heatmapmainclassNormalized, p_NT_global, p2_scatter, 
                   nrow=1, ncol=4, 
                   labels=c(letters[c(1:4)]), 
                   font.label=list(size=18, color="black", face="bold"), 
                   widths=c(1, 2.3, 2.5, 1.5))

# Arrange metabolic pathway plots into a single figure
p2all <- ggarrange(p_metac3, p_metac1, p_metac2, 
                   nrow=1, ncol=3, 
                   labels=c(letters[c(6:8)]), 
                   font.label=list(size=18, color="black", face="bold"), 
                   common.legend=T, legend="bottom")

# Arrange average and distance plots into a single figure
p3all <- ggarrange(p_average, p_distance, 
                   nrow=1, ncol=2, 
                   labels=c(letters[9:10]), 
                   font.label=list(size=18, color="black", face="bold"), 
                   common.legend=T, legend="bottom")

# Combine all previous plots into one larger figure
p4all <- ggarrange(p_sig, p2all, p3all, 
                   nrow=1, ncol=3, 
                   labels=c(letters[c(5)], "", ""), 
                   font.label=list(size=18, color="black", face="bold"), 
                   widths=c(1, 3, 1))

# Arrange distance and local allocation plots into a single figure
p5all <- ggarrange(p_distance2, p_NT_local, pNTall_local, 
                   nrow=1, ncol=3, 
                   labels=c(letters[c(11:13)]), 
                   font.label=list(size=18, color="black", face="bold"), 
                   widths=c(1, 2, 1.5))

# Arrange carbohydrate plots into a single figure
p6all <- ggarrange(bi_plot_carbo, p_box_carbo, bi_plot_lipid, p_box_lipid, 
                   nrow=1, ncol=4, 
                   labels=c(letters[c(14:17)]), 
                   font.label=list(size=18, color="black", face="bold"))

# Combine all figures into the final arrangement
p7all <- ggarrange(p1all, p4all, p5all, p6all, 
                   nrow=4, ncol=1)

# Save the final combined figure as a TIFF file
tiff(paste("./figures/", "figure3_new.tiff", sep=""), 
     height=27, width=26, res=300, units="in", compression="lzw")
print(p7all)
dev.off()

# Save the final combined figure as a PDF file
pdf(paste("./figures/", "figure3_new.pdf", sep=""), height=27, width=26)
print(p7all)
dev.off()

# Optionally save the current workspace
# save.image("/data/Project/Project1/result_new/rdata/step_figure3.RData")

# Create and arrange a heatmap for metabolic activity
gb_heatmapmainclass <- plotFChetamap(t(fc_mainclass$FC_sum_sig), "Log2FC of metabolic activity")
ggarrange(gb_heatmapmainclass)

# Plot for Glycan metabolism
p_metac4 <- plotpathwayNT(mean_auc_rec2_class[mean_auc_rec2_class$cell != "Mast", ], "Glycan metabolism") + 
  scale_y_log10() + ylab("Global metabolic allocation")

# Prepare nucleotide ratio data for PCA
nuc_ratio <- pathwayRatioNormalizedpatient[grepl("Epithelial|Myeloid|Fibroblast", rownames(pathwayRatioNormalizedpatient)), 
                                           colnames(pathwayRatioNormalizedpatient) %in% unique(all_gene_pathway[all_gene_pathway$kegg_class == "Nucleotide metabolism", 2])]
prc <- prcomp(nuc_ratio)

# Prepare data for PCA biplot
temp_x <- str_split_fixed(rownames(nuc_ratio), "_", 4)
colnames(temp_x) <- c("dataset", "class", "cell", "patient")
pc_result <- plotPcaBiplot(nuc_ratio, temp_x, "Local metabolic allocation\nin nucleotide metabolism", scale=2)

# Set colors for the plot
col_x <- brewer.pal(n=7, name="Set3")
col_x[2] <- "#EEE685"
bi_plot_nuc <- pc_result$bi_plot + scale_color_manual(name="Cell type", values=col_x[-c(1, 2, 6, 7)])

# Create boxplot for nucleotide metabolism pathways
p_box_nuc <- PlotboxFC2(c("Purine synthesis", "Pyrimidine catabolism", "Pyrimidine synthesis", "Purine catabolism"), 
                        t(fc_pathwayRatioNormalized$FC[grepl("Epithelial|Myeloid|Fibroblast", rownames(fc_pathwayRatioNormalized$FC)), ]), 
                        fc_pathwayRatioNormalized$all_path_sigClass[fc_pathwayRatioNormalized$all_path_sigClass$wcx_p_pair < 0.05, ], 
                        "Log2FC of local metabolic allocation") + 
  geom_hline(yintercept=0, linetype="dashed") + 
  guides(fill=guide_legend(nrow=2), color=guide_legend(title="Sig", nrow=2)) + 
  ggtitle("Nucleotide metabolism") + 
  coord_flip() + 
  scale_fill_manual(name="Cell type", values=col_x[-c(1, 2, 6, 7)]) + 
  xlab("Pathways")

# Arrange all supplementary plots into a single figure
psupall <- ggarrange(gb_heatmapmainclass, p_metac4, p_carbo1, p_carbo2, p_nuc1, p_nuc2, bi_plot_nuc, p_box_nuc, 
                     nrow=4, ncol=2, 
                     labels=c(letters[c(1:8)]), 
                     font.label=list(size=18, color="black", face="bold"))

# Save the supplementary figure as a TIFF file
tiff(paste("./figures/", "supfigure5.tiff", sep=""), height=20, width=15, res=300, units="in", compression="lzw")
print(psupall)
dev.off()

                   


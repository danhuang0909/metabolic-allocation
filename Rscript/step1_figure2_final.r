# Load the required data
load("./load_all_files2.rdata")

# Set the working directory for results
setwd("./result_new")
# Source the plotting functions
source("./Rscript/functions_plot.r")

##############################Plot flowchart
# Define x and y positions for the flowchart
xpos <- 1:5
ypos <- xpos**2
data_frame = data.frame(xpos = xpos,
                        ypos = ypos)

# Print the data points
print("Data points")
print(data_frame)

# Load the TIFF library for image handling
library(tiff)

# Specify the path to the TIFF image
path <- "./figures/figure2a.tif"
# Read the TIFF image
img <- readTIFF(path, native = TRUE)
# Convert the image to a raster graphic
img <- rasterGrob(img, interpolate=TRUE)

# Plotting the data with the image as background
p2_flowchart = qplot(xpos, ypos, geom="blank") +
  annotation_custom(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(color=I("White")) + theme_void() + 
  theme(plot.margin = unit(c(-2,0.2,0.2,0.2), "cm"))

# Aggregate data for pathway mean normalized by patient class
pathway_mean_patientClassNormalized <- getAggregateDate(allcells_pathwayMeanSCTMainClassNormalized, pathway_mean_allsctClassNormalized@meta.data)
# Create a heatmap for normalized metabolic allocation in normal samples
gb_heatmapNormalized = plotHeatmapPathway(pathway_mean_patientClassNormalized, type_x="Normal", title_x="Global metabolic allocation in normal samples")
# Arrange the heatmap
ggarrange(gb_heatmapNormalized)

# Aggregate data for class normalized pathways
pathway_classNormalized = getAggregateDate2(allcells_pathwayMeanSCTMainClassNormalized, pathway_mean_allsctClassNormalized@meta.data)
# Transpose the data for normal samples
pathway_classnormal_Normalized = t(pathway_classNormalized[grepl("Normal", rownames(pathway_classNormalized)), ])
# Clean up column names
colnames(pathway_classnormal_Normalized) <- gsub("_Normal", "", colnames(pathway_classnormal_Normalized))

# Aggregate data for pathways
pathway_class = getAggregateDate2(allcells_pathwayMeanSCTMainClass, pathway_mean_allsct_seuClass@meta.data)
# Transpose data for normal samples
pathway_classnormal = t(pathway_class[grepl("Normal", rownames(pathway_class)), ])
# Clean up column names
colnames(pathway_classnormal) <- gsub("_Normal", "", colnames(pathway_classnormal))

# Plot line pathways for SCT and normalized data
plineSCT <- plotLinePathway(pathway_classnormal, "Normal sample SCT")
plinenormalized <- plotLinePathway(pathway_classnormal_Normalized, "Normal sample Normalized")

# Define color palette
col_x = brewer.pal(n = 7, name = "Set3")
col_x[2] = "#EEE685"

# Combine rank data for plotting
rank_data_all <- rbind(cbind(plineSCT$rank_data, type="Average pathway activity"), cbind(plinenormalized$rank_data, type="Global metabolic allocation"))
rank_data_all$type = factor(rank_data_all$type, levels=c("Average pathway activity", "Global metabolic allocation"))
rank_data_all$Var2 <- factor(rank_data_all$Var2, levels=levels(rank_data_all$Var2)[c(1:4, 6:7, 5, 8:11, 13, 12)])

# Create a rank plot
p_rank0 <- ggplot(data=rank_data_all, aes(x=Var2, y=value, col=Var1, group=Var1)) +
  geom_point() + geom_line() + theme_hd_plain2() + 
  theme(axis.title.y = element_blank()) + 
  scale_color_manual(name="", values=col_x) + 
  ylab("Rank") + 
  facet_wrap(~type, nrow=2, scales="free_x") + 
  coord_flip() + 
  theme(legend.justification = "right") + 
  guides(color=guide_legend(override.aes=list(size=3)))

# Aggregate metabolic activity data
metabolic_activity = aggregate(combined.sct_filter@meta.data$metabolic_activity, 
                               list(dataset=combined.sct_filter@meta.data$dataset, 
                                    cell=combined.sct_filter@meta.data$Cell_type, 
                                    patient=combined.sct_filter@meta.data$patient, 
                                    type=combined.sct_filter@meta.data$class_new), 
                               mean)

# Filter for normal samples excluding Mast cells
metabolic_activity = metabolic_activity[metabolic_activity[,4] == "Normal" & metabolic_activity[,2] != "Mast", ]
rownames(metabolic_activity) <- paste(metabolic_activity[,1], metabolic_activity[,2], metabolic_activity[,3], sep="_")

# Calculate Spearman correlation for metabolic activity
xx = getcorSpearmanForall(metabolic_activity[,5, drop=F], t(pathway_classnormal))
xx_mean <- aggregate(xx$rho.rho, list(path=xx$pathways), mean)
xx$pathways <- factor(xx$pathways, levels=xx_mean$path[order(xx_mean$x)])

# Create a correlation plot
p_cor_meta1 <- ggplot(xx, aes(x=pathways, y=rho.rho)) +
  geom_boxplot() + geom_point(aes(color=cell)) + 
  theme_hd_plain2() + 
  xlab("Average pathway \nactivity") + 
  ylab("Spearman's rho") + 
  coord_flip() + 
  ylim(0, 1) + 
  theme(axis.text.x = element_text(angle=30, hjust=1), 
        plot.title = element_text(hjust=1)) + 
  ggtitle("Correlation with total metabolic activity") + 
  scale_color_manual(name="", values=col_x)

# Add pathway activity to metabolic activity data
x = metabolic_activity
x$path <- pathway_classnormal["Canonical amino acid metabolism", rownames(x)]
cols_x = brewer.pal(n = 7, name = "Set3")[c(3, 4, 2, 5, 7, 6, 1)]
cols_x[3] = "#EEE685"

# Create another correlation plot
p_cor_meta2 <- ggplot(x, aes(x=path, y=x, fill=I("gray"))) +
  geom_point(aes(color=cell)) + 
  scale_color_manual(name="", values=brewer.pal(n = 7, name = "Set3")) + 
  theme_hd_plain() + 
  guides(color=guide_legend(override.aes=list(size=3))) + 
  xlab("Average pathway activity of canonical amino acid metabolism") + 
  ylab("Total metabolic activity")

# Process normalized pathway data and create boxplots
temp = as.data.frame(t(pathway_classnormal_Normalized))
temp$cells = str_split_fixed(rownames(temp), "_", 3)[,2]
temp$cells <- factor(temp$cells, levels=c("Epithelial", "Fibroblast", "Endothelial", "Myeloid", "T cell", "Plasma", "B cell"))

# Combine data for boxplots
temp_combine <- rbind(data.frame(path="Nucleotide metabolism", value=temp[, c("Nucleotide metabolism")], Cell_types=temp$cells),
                      data.frame(path="Xenobiotics metabolism", value=temp[, c("Xenobiotics metabolism")], Cell_types=temp$cells))
temp_combine$path <- factor(temp_combine$path, levels=c("Xenobiotics metabolism", "Nucleotide metabolism"))

# Create a boxplot for nucleotides and xenobiotics
p_nuc_xeno <- ggplot(temp_combine, aes(x=Cell_types, y=value, fill=Cell_types)) +
  geom_boxplot() + 
  scale_fill_manual(values=brewer.pal(n = 7, name = "Set3")[c(3, 4, 2, 5, 7, 6, 1)]) + 
  theme_hd_minimal_plain() + 
  theme(axis.text.x = element_blank()) + 
  ylab("Global metabolic allocation") + 
  facet_wrap(~path, nrow=2) + 
  xlab("Cell type")

# Create another combined dataset for carbohydrate and amino acid metabolism
temp_combine2 <- rbind(data.frame(path="Carbohydrate metabolism", value=temp[, c("Carbohydrate metabolism")], Cell_types=temp$cells),
                       data.frame(path="Canonical amino acid metabolism", value=temp[, c("Canonical amino acid metabolism")], Cell_types=temp$cells))
temp_combine2$path <- factor(temp_combine2$path, levels=c("Canonical amino acid metabolism", "Carbohydrate metabolism"))

# Create boxplots for carbohydrate and amino acid metabolism
p_car_cano <- ggplot(temp_combine2, aes(x=Cell_types, y=value, fill=Cell_types)) +
  geom_boxplot() + 
  scale_fill_manual(values=brewer.pal(n = 7, name = "Set3")[c(3, 4, 2, 5, 7, 6, 1)]) + 
  theme_hd_minimal_plain() + 
  theme(axis.text.x = element_blank()) + 
  ylab("Global metabolic allocation") + 
  facet_wrap(~path, nrow=2)

# Create a boxplot for protein metabolism
ggplot(temp, aes(x=cells, y=`Protein metabolism`, fill=cells)) +
  geom_boxplot() + 
  scale_fill_manual(values=brewer.pal(n = 7, name = "Set3")[c(3, 4, 2, 5, 7, 6, 1)]) + 
  theme_hd_minimal_plain() + 
  theme(axis.text.x = element_blank()) + 
  ylab("Global metabolic allocation")

# Scatter plot for nucleotide and xenobiotics metabolism
cols_x = brewer.pal(n = 7, name = "Set3")[c(3, 4, 2, 5, 7, 6, 1)]
cols_x[3] = "#EEE685"
nucle_scatter <- ggplot(temp, aes(x=`Nucleotide metabolism`, y=`Xenobiotics metabolism`, color=cells, fill=I("gray"))) +
  geom_point() + 
  scale_color_manual(values=cols_x) + 
  theme_hd_plain() + 
  ggtitle("Global metabolic allocation") + 
  guides(color=guide_legend(override.aes=list(size=3)))

# Load prediction data
pre_class <- readRDS("./prediction_data.rds")
all_reglobal = pre_class$all_pre[!grepl("_", pre_class$all_pre$path), ]

# Plot global metabolic allocation predictions
p_celltype_global <- plotPredictionResult(all_reglobal) + 
  ylab("Global metabolic allocation") + 
  theme(legend.justification = "right") + 
  xlab("Cell type")

# Load normalized prediction data
ratioNormalize_pre <- readRDS("ratioNormalize_prediction_all.rds")
all_re = ratioNormalize_pre$all_pre

# Plot local metabolic allocation predictions
p_celltype_local <- plotPredictionResult(all_re) + 
  ylab("Local metabolic allocation") + 
  theme(legend.justification = "right") + 
  xlab("Cell type")

# Aggregate pathway means by patient
pathway_mean_patient <- aggregate(t(as.matrix(pathway_mean_allsct_seu@assays$RNA@counts)), 
                                  list(dataset=pathway_mean_allsct_seu@meta.data$dataset, 
                                       type=pathway_mean_allsct_seu$class_new, 
                                       cell=pathway_mean_allsct_seu@meta.data$Cell_type, 
                                       patient=pathway_mean_allsct_seu$patient), 
                                  mean)

# Prepare dataframe for normal patients
pathway_mean_df_patient <- pathway_mean_patient[,-c(1:4)]
rownames(pathway_mean_df_patient) <- paste(pathway_mean_patient[,1], pathway_mean_patient[,2], pathway_mean_patient[,3], pathway_mean_patient[,4], sep="_")
pathway_normal_patient = t(pathway_mean_df_patient[grepl("Normal", rownames(pathway_mean_df_patient)), ])
colnames(pathway_normal_patient) <- gsub("_Normal", "", colnames(pathway_normal_patient))

# Create boxplot for carbohydrate metabolism
pbox_car = Plotbox2(c("Citric acid cycle", "Hyaluronan metabolism", "Glycolysis/gluconeogenesis", "Pentose phosphate pathway", "Butanoate metabolism"), 
                    pathway_normal_patient, "Carbohydrate metabolism")
pbox_car <- pbox_car + guides(fill=guide_legend(title="Cell type", nrow=1)) + 
  ggtitle("Carbohydrate metabolism in normal tissue") + 
  ylab("Local metabolic allocation")

# Create PCA for carbohydrate metabolism
pca_car = PlotpcaRatio(c("Hyaluronan metabolism", "Glycolysis/gluconeogenesis", "Citric acid cycle", "Pentose phosphate pathway", "Butanoate metabolism"), 
                       pathway_normal_patient, "Carbohydrate metabolism") + 
  ggtitle("Local metabolic allocation\ncarbohydrate metabolism in normal tissue") + 
  guides(color=guide_legend(override.aes=list(size=3)))

# Create boxplot for amino acid metabolism
pbox_ami = Plotbox2(c("Lysine metabolism", "Arginine and proline metabolism", "Methionine and cysteine metabolism", "Alanine, aspartate and glutamate metabolism"), 
                    pathway_normal_patient, "Canonical amino acid metabolism", cell_order=c("B cell", "T cell", "Myeloid", "Plasma", "Endothelial", "Fibroblast", "Epithelial"))
pbox_ami = pbox_ami + guides(fill=guide_legend(title="Cell type", nrow=1)) + 
  ggtitle("Canonical amino acid metabolism in normal tissue") + 
  ylab("Local metabolic allocation")

# Create PCA for amino acid metabolism
pca_ami = PlotpcaRatio(c("Lysine metabolism", "Arginine and proline metabolism", "Methionine and cysteine metabolism", "Alanine, aspartate and glutamate metabolism"), 
                       pathway_normal_patient, "Canonical amino acid metabolism")
pca_ami = pca_ami + guides(fill=guide_legend(title="Cell type", nrow=2)) + 
  ggtitle("Local metabolic allocation\ncanonical amino acid metabolism in normal tissue")

# Arrange plots
p0all = ggarrange(p_cor_meta1, p_cor_meta2, nrow=2, ncol=1, 
                  labels=c(letters[c(1, 2)]), 
                  font.label=list(size=18, color="black", face="bold"), 
                  common.legend=T, legend="bottom") 

p1all = ggarrange(p0all, p2_flowchart, gb_heatmapNormalized, nrow=1, ncol=3, 
                  labels=c("", letters[c(3, 4)]), 
                  font.label=list(size=18, color="black", face="bold"), 
                  widths=c(1.2, 2, 2))

p2all = ggarrange(p_rank0, p_nuc_xeno + theme(legend.position="none"), nrow=2, ncol=1, 
                  labels=c(letters[c(5, 8)]), 
                  font.label=list(size=18, color="black", face="bold"))

p3all = ggarrange(p_celltype_global, p_celltype_local, nrow=2, ncol=1, 
                  labels=c(letters[c(6, 9)]), 
                  font.label=list(size=18, color="black", face="bold"), 
                  legend="bottom", common.legend=T)

p4all = ggarrange(pbox_car, pbox_ami, nrow=1, ncol=2, 
                  labels=c(letters[c(11, 12)]), 
                  font.label=list(size=18, color="black", face="bold"), 
                  common.legend=T, legend="bottom")

p5all = ggarrange(nucle_scatter, pca_car, pca_ami, nrow=3, ncol=1, 
                  labels=c(letters[c(7, 10, 13)]), 
                  font.label=list(size=18, color="black", face="bold"), 
                  common.legend=T, legend="bottom")

p6all = ggarrange(p2all, p3all, nrow=1, ncol=2, widths=c(2, 3))
p7all = ggarrange(p6all, p4all, nrow=2, ncol=1, heights=c(2, 1))
p8all = ggarrange(p7all, p5all, nrow=1, ncol=2, widths=c(2.5, 1))
p9all = ggarrange(p1all, p8all, nrow=2, ncol=1, heights=c(1, 3))

# Save the final figure as TIFF
tiff(paste("./figures/", "figure2_new.tiff", sep=""), height=27, width=23, res=300, units="in", compression="lzw")
print(p9all)
dev.off()

# Save the final figure as PDF
pdf(paste("./", "figure2_new.pdf", sep=""), height=27, width=23)
print(p9all)
dev.off()




# Aggregate pathway data using RNA counts and metadata
pathway_sct = getAggregateDate2(pathway_mean_allsct_seu@assays$RNA@counts, pathway_mean_allsct_seu@meta.data)

# Transpose the data for normal samples
pathway_sctnormal = t(pathway_sct[grepl("Normal", rownames(pathway_sct)), ])
# Clean up column names by removing "_Normal"
colnames(pathway_sctnormal) <- gsub("_Normal", "", colnames(pathway_sctnormal))

# Calculate Spearman correlation for metabolic activity using pathway data
xx_sct = getcorSpearmanForall(metabolic_activity[, 5, drop=F], t(pathway_sctnormal))
# Aggregate mean Spearman correlation values by pathway
xx_mean_sct <- aggregate(xx_sct$rho.rho, list(path=xx_sct$pathways), mean)
# Set pathway levels for plotting
xx_sct$pathways <- factor(xx_sct$pathways, levels = xx_mean_sct$path[order(xx_mean_sct$x)])

# Create a boxplot for Spearman correlation with total metabolic activity
p_cor_meta_sct <- ggplot(xx_sct, aes(x=pathways, y=rho.rho)) +
  geom_boxplot() + 
  geom_point(aes(color=cell)) + 
  theme_hd_plain2() + 
  xlab("Average pathway activity") + 
  ylab("Spearman's rho") + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust=1)) + 
  ggtitle("Correlation with total metabolic activity") + 
  scale_color_manual(name="", values=col_x) + 
  geom_vline(xintercept=0.5, linetype="dashed")

# Calculate Spearman correlation for normalized pathway data
xx_Normalized = getcorSpearmanForall(metabolic_activity[, 5, drop=F], t(pathway_classnormal_Normalized))
# Aggregate mean Spearman correlation values for normalized pathways
xx_mean_Normalized <- aggregate(xx_Normalized$rho.rho, list(path=xx_Normalized$pathways), mean)
# Set pathway levels for normalized data
xx_Normalized$pathways <- factor(xx_Normalized$pathways, levels = xx_mean_Normalized$path[order(xx_mean_Normalized$x)])

# Create a boxplot for normalized Spearman correlation
p_cor_meta_Normalized <- ggplot(xx_Normalized, aes(x=pathways, y=rho.rho)) +
  geom_boxplot() + 
  geom_point(aes(color=cell)) + 
  theme_hd_plain2() + 
  xlab("Global metabolic allocation") + 
  ylab("Spearman's rho") + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust=1)) + 
  ggtitle("Correlation with total metabolic activity") + 
  scale_color_manual(name="", values=col_x)

# Prepare metabolic activity data for plotting
x = metabolic_activity
x$path <- pathway_classnormal_Normalized["Canonical amino acid metabolism", rownames(x)]
cols_x = brewer.pal(n = 7, name = "Set3")[c(3, 4, 2, 5, 7, 6, 1)]
cols_x[3] = "#EEE685"

# Create a scatter plot for metabolic allocation
p_cor_meta2_normalized <- ggplot(x, aes(x=path, y=x, fill=I("gray"))) +
  geom_point(aes(color=cell)) + 
  scale_color_manual(name="", values=brewer.pal(n = 7, name = "Set3")) + 
  theme_hd_plain() + 
  guides(color=guide_legend(override.aes=list(size=3))) + 
  xlab("Global metabolic allocation\n(Canonical amino acid metabolism)") + 
  ylab("Total metabolic activity")

# Create a global metabolic allocation bar plot
p1_global <- ggplot(all_reglobal[all_reglobal$type == "cell_cell", ], aes(x=cell, y=accuracy, fill=cell2)) +
  geom_bar(stat="identity", position='dodge') + 
  facet_wrap(~path, nrow=4) + 
  theme_hd_minimal_plain() + 
  theme(axis.text.x = element_text(angle=90, hjust=1), legend.position="bottom", axis.title.x=element_blank()) + 
  guides(fill=guide_legend(nrow=1)) + 
  scale_fill_manual(name="Cell type", values=brewer.pal(n = 7, name = "Set3")) + 
  ylab("10-fold cross-validation classifier accuracy for pairwise cell type comparison") + 
  ggtitle("Global metabolic allocation") + 
  geom_hline(yintercept=0.8, linetype="dashed")

# Save the global metabolic allocation plot as TIFF
tiff(paste("./figures/", "supfigure3.tiff", sep=""), height=16, width=16, res=300, units="in", compression="lzw")
print(p1_global)
dev.off()

# Create a local metabolic allocation bar plot
p1_local <- ggplot(all_re[all_re$type == "cell_cell", ], aes(x=cell, y=accuracy, fill=cell2)) +
  geom_bar(stat="identity", position='dodge') + 
  facet_wrap(~path, nrow=4) + 
  theme_hd_minimal_plain() + 
  theme(axis.text.x = element_text(angle=90, hjust=1), legend.position="bottom", axis.title.x=element_blank()) + 
  guides(fill=guide_legend(nrow=1)) + 
  scale_fill_manual(name="Cell type", values=brewer.pal(n = 7, name = "Set3")) + 
  ylab("10-fold cross-validation classifier accuracy for pairwise cell type comparison") + 
  ggtitle("Local metabolic allocation") + 
  geom_hline(yintercept=0.9, linetype="dashed")

# Save the local metabolic allocation plot as TIFF
tiff(paste("./figures/", "supfigure4.tiff", sep=""), height=16, width=16, res=300, units="in", compression="lzw")
print(p1_local)
dev.off()

# Process normalized pathway data for scatter plot
temp = as.data.frame(t(pathway_classnormal_Normalized))
temp$cells = str_split_fixed(rownames(temp), "_", 3)[, 2]
temp$cells <- factor(temp$cells, levels=c("Epithelial", "Fibroblast", "Endothelial", "Myeloid", "T cell", "Plasma", "B cell"))

# Create scatter plot for canonical amino acid and carbohydrate metabolism
cols_x = brewer.pal(n = 7, name = "Set3")[c(3, 4, 2, 5, 7, 6, 1)]
cols_x[3] = "#EEE685"
carcno_scatter <- ggplot(temp, aes(x=`Canonical amino acid metabolism`, y=`Carbohydrate metabolism`, color=cells, fill=I("gray"))) +
  geom_point() + 
  scale_color_manual(values=cols_x) + 
  theme_hd_plain() + 
  ggtitle("Global metabolic allocation") + 
  guides(color=guide_legend(override.aes=list(size=3)))

# Create boxplot for xenobiotics metabolism
pbox_Xeno = Plotbox2(unique(all_gene_pathway[all_gene_pathway$kegg_class == "Xenobiotics metabolism", "final_name"]), pathway_normal_patient, "Xenobiotics metabolism", cell_order=c("B cell", "T cell", "Myeloid", "Plasma", "Endothelial", "Fibroblast", "Epithelial"))
pbox_Xeno = pbox_Xeno + guides(fill=guide_legend(title="Cell types", nrow=1)) + 
  ggtitle("Xenobiotics metabolism in normal tissue") + 
  ylab("Local metabolic allocation")

# Create PCA for xenobiotics metabolism
pca_Xeno <- PlotpcaRatio(unique(all_gene_pathway[all_gene_pathway$kegg_class == "Xenobiotics metabolism", "final_name"]), pathway_normal_patient, "Xenobiotics metabolism")
pca_Xeno <- pca_Xeno + guides(fill=guide_legend(title="Metabolism", nrow=2)) + 
  ggtitle("Local metabolic allocation\nxenobiotics metabolism in normal tissue") + 
  guides(color=guide_legend(override.aes=list(size=3)))

# Create PCA for nucleotide metabolism
pca_Nuc <- PlotpcaRatio(unique(all_gene_pathway[all_gene_pathway$kegg_class == "Nucleotide metabolism", "final_name"]), pathway_normal_patient, "Nucleotide metabolism")
pca_Nuc <- pca_Nuc + guides(fill=guide_legend(title="Metabolism", nrow=2)) + 
  ggtitle("Local metabolic allocation\nnucleotide metabolism in normal tissue") + 
  guides(color=guide_legend(override.aes=list(size=3)))

# Create boxplot for nucleotide metabolism
pbox_Nuc = Plotbox2(unique(all_gene_pathway[all_gene_pathway$kegg_class == "Nucleotide metabolism", "final_name"]), pathway_normal_patient, "Nucleotide metabolism", cell_order=c("B cell", "T cell", "Myeloid", "Plasma", "Endothelial", "Fibroblast", "Epithelial"))
pbox_Nuc = pbox_Nuc + guides(fill=guide_legend(title="Cell types", nrow=1)) + 
  ggtitle("Nucleotide metabolism in normal tissue") + 
  ylab("Local metabolic allocation")

# Create scatter plot for citric acid cycle vs glycolysis
i = "Citric acid cycle"
j = "Glycolysis/gluconeogenesis"
temp = data.frame(i=pathway_normal_patient[i, ], j=pathway_normal_patient[j, ], cell=str_split_fixed(colnames(pathway_normal_patient), "_", 3)[, 2], dataset=str_split_fixed(colnames(pathway_normal_patient), "_", 3)[, 1])
ggplot(temp, aes(x=i, y=j)) +
  geom_point(aes(col=cell)) + 
  theme_hd_plain() + 
  scale_color_manual(values=brewer.pal(n = 7, name = "Set3")) + 
  xlab(i) + 
  ylab(j) + 
  theme(legend.position="left") + 
  ggtitle("Normal")

# Arrange subplots for final presentation
subpall1 <- ggarrange(p_cor_meta_Normalized, p_cor_meta2_normalized, pca_Nuc, pca_Xeno, carcno_scatter, nrow=5, ncol=1, labels=c(letters[c(2:6)]), font.label=list(size=18, color="black", face="bold"), common.legend=T, legend="none")
subpall2 <- ggarrange(p_cor_meta_sct, subpall1, nrow=1, ncol=2, labels=c(letters[c(1)], ""), font.label=list(size=18, color="black", face="bold"), widths=c(1.5, 1), common.legend=T, legend="bottom")

# Save the final combined figure as TIFF
tiff(paste("./figures/", "supfigure2.tiff", sep=""), height=18, width=16, res=300, units="in", compression="lzw")
print(subpall2)
dev.off()



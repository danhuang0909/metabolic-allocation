# Load necessary data files
load("./data/load_all_files2.rdata")
setwd("/data/Project/Project1/result_new")
source("./Rscript/functions_plot.r")

# Extract metadata from the filtered dataset
meta_data <- all_dataset_SCT_filter@meta.data

# Load required libraries
library(circlize)
library(tiff)

# Load KEGG pathway genes and filter for signaling pathways
kegg_path_genes <- readRDS(file="./data/all_kegg_pathway_genes.rds")
kegg_signal <- kegg_path_genes[grepl("signal", kegg_path_genes$description), ]

# Load CPM (Counts Per Million) data and test results
cpm_all <- readRDS(file="./data/sum_all_genes_count_rawRNA_cpm.rds")
result_df <- readRDS("./data/sum_all_genes_count_rawRNA_test.rds")
result_df <- result_df[result_df$cell_type != "Mast", ]
result_df_sig <- result_df[result_df$PValue < 0.05, ]
cpm_all <- cpm_all[!grepl("Mast", rownames(cpm_all)), ]

# Select signaling pathway related genes
cpm_all_select_signal <- cpm_all[, colnames(cpm_all) %in% kegg_signal$symbol & !colnames(cpm_all) %in% unique(all_gene_pathway$gene)]

# Calculate fold change data for selected signaling genes
cpm_all_select_signal_fc <- getFCdataSinglelog(cpm_all_select_signal)
cpm_all_select_signal_fcAN <- getFCdataSinglelogAllnormaleach(cpm_all_select_signal, str_split_fixed(rownames(cpm_all_select_signal), "_", 4)[, 3])

# Load normalized pathway mean data
allcells_pathwayMeanSCTMainClassNormalized <- readRDS(file="./data/allcells_pathwayMeanSCTMainClassNormalizedByMetacount.Rds")
meta_data <- all_dataset_SCT_filter@meta.data
fc_mainclassNormalized <- getFCresult(allcells_pathwayMeanSCTMainClassNormalized, meta_data)
pathwayRatioNormalized <- readRDS(file="./pathwayRatioNormalized.rds")
fc_pathwayRatioNormalized <- getFCresult(pathwayRatioNormalized, meta_data)

# Set colors for plots
col_x <- brewer.pal(n=7, name="Set3")
col_x[2] <- "#EEE685"

# UMAP visualization for local metabolic allocation
u1 <- umap(fc_pathwayRatioNormalized$FC)
df <- data.frame(Dim1 = u1$layout[, 1], Dim2 = u1$layout[, 2], CellType = str_split_fixed(rownames(u1$layout), "_", 3)[, 2])
umap_local <- ggplot(df, aes(x=Dim1, y=Dim2, col=CellType)) + 
  geom_point() + 
  scale_color_manual(name="Cell type", values=col_x) + 
  theme_hd_plain() + 
  ggtitle("Log2FC of local\nmetabolic allocation")

# UMAP visualization for global metabolic allocation
u1 <- umap(fc_mainclassNormalized$FC)
df <- data.frame(Dim1 = u1$layout[, 1], Dim2 = u1$layout[, 2], CellType = str_split_fixed(rownames(u1$layout), "_", 3)[, 2])
umap_global <- ggplot(df, aes(x=Dim1, y=Dim2, col=CellType)) + 
  geom_point() + 
  scale_color_manual(name="Cell type", values=col_x) + 
  theme_hd_plain() + 
  ggtitle("Log2FC of global\nmetabolic allocation")

# UMAP visualization for signaling pathway related genes
u1 <- umap(cpm_all_select_signal_fc)
df <- data.frame(Dim1 = u1$layout[, 1], Dim2 = u1$layout[, 2], CellType = str_split_fixed(rownames(u1$layout), "_", 3)[, 2])
umap_signal <- ggplot(df, aes(x=Dim1, y=Dim2, col=CellType)) + 
  geom_point() + 
  scale_color_manual(name="Cell type", values=col_x) + 
  theme_hd_plain() + 
  ggtitle("Log2FC of signaling\npathway related genes")

# Calculate Spearman correlation for signaling pathway related genes
temp_corsignalRatio <- getcorSpearmanForall(cpm_all_select_signal_fc, fc_pathwayRatioNormalized$FC)
saveRDS(temp_corsignalRatio, file="./data/signalgenes_correlation_pathwayRatio.Rds")
temp_corsignalRatio <- readRDS(file="./data/signalgenes_correlation_pathwayRatio.Rds")
temp_corsignalRatio_sig <- temp_corsignalRatio[temp_corsignalRatio$pvalue < 0.05, ]
path_sig <- fc_pathwayRatioNormalized$all_path_sigClass[fc_pathwayRatioNormalized$all_path_sigClass$wcx_p < 0.05, ]
temp_corsignalRatio_sig_sig <- temp_corsignalRatio_sig[paste(temp_corsignalRatio_sig$cell, temp_corsignalRatio_sig$gene) %in% paste(result_df_sig$cell_type, result_df_sig$genes) & 
                                                         paste(temp_corsignalRatio_sig$cell, temp_corsignalRatio_sig$pathways) %in% paste(path_sig$cell, path_sig$path), ]

# Plot heatmap for significant correlations
pdf("./figures/signal_pathway_correlation_sigAllNormal.pdf", width=28, height=20)
all_sig <- plotHeatmapSignalPath(temp_corsignalRatio_sig_sig, "all sig")
dev.off()

# Plot heatmap for selected signal genes
pdf("./figures/selected_signal_genes_pheatmap_pathwayAllNormal.pdf", width=20, height=20)
signal_path <- GetSigheatmapCor(all_sig)
dev.off()

# Save the signal path data
saveRDS(signal_path, file="./data/signal_path.rds")

# Prepare data for bar plot
pathways <- signal_path$paths_all
genes_all <- as.data.frame(signal_path$genes_all)

pbar_plot <- PlotBarCount(temp_corsignalRatio_sig_sig[paste(temp_corsignalRatio_sig_sig$cell, temp_corsignalRatio_sig_sig$gene) %in% paste(genes_all[, 1], genes_all[, 2]) & 
                                                        paste(temp_corsignalRatio_sig_sig$cell, temp_corsignalRatio_sig_sig$pathways) %in% paste(pathways[, 1], pathways[, 2]), ], "")

# Count opposing pathways
corr_count <- pbar_plot$corr_counts
opposing_count <- corr_count %>%
  group_by(cell) %>%
  summarize(all_count = sum(pos_count != 0 | neg_count != 0),
            opp_count_high = sum(pos_count > 3 & abs(neg_count) > 3)) %>% 
  as.data.frame()

# Prepare data for plotting opposing counts
opposing_count_df <- reshape2::melt(opposing_count, id.vars="cell")
opposing_count_df$cell <- factor(opposing_count_df$cell, levels=c("Epithelial", "T cell", "B cell", "Fibroblast", "Myeloid", "Endothelial", "Plasma"))

# Create bar plot for opposing counts
p_count <- ggplot(opposing_count_df, aes(x=cell, y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(labels=c("All metabolic pathway under\nsignal gene regulation", "Metabolic pathway under\ndiverse signal gene regulation"), 
                    values=c("#528B8B", "#98F5FF")) +
  theme_hd_plain() +
  theme(axis.text.x=element_text(angle=90)) +
  theme(legend.position="bottom", legend.direction="vertical") +
  ylab("Count")

# Create a cell-cell interaction matrix
cell_types <- unique(genes_all$cell)
cell_cell_matrix <- expand.grid(cell1=cell_types, cell2=cell_types) %>%
  mutate(n_genes = apply(., 1, function(x) {
    length(intersect(
      genes_all$path[genes_all$cell %in% x[1]],
      genes_all$path[genes_all$cell %in% x[2]]
    ))
  })) %>%
  mutate(n_genes_cell2 = apply(., 1, function(x) {
    n_distinct(genes_all$path[genes_all$cell %in% x[2]])
  })) %>%
  mutate(n_genes_ratio = apply(., 1, function(x) {
    length(intersect(
      genes_all$path[genes_all$cell %in% x[1]],
      genes_all$path[genes_all$cell %in% x[2]]
    )) / length(genes_all$path[genes_all$cell %in% x[2]])
  }))

# Output results
print(cell_cell_matrix)

# Create statistics for cell types
cell_type_stats <- data.frame(
  cell = unique(genes_all$cell),
  n_total_genes = sapply(unique(genes_all$cell), function(c) n_distinct(genes_all$path[genes_all$cell == c])),
  n_unique_genes = sapply(unique(genes_all$cell), function(c) n_distinct(genes_all$path[genes_all$cell == c & !genes_all$path %in% genes_all$path[genes_all$cell != c]])),
  pct_unique_genes = sapply(unique(genes_all$cell), function(c) {
    n_distinct(genes_all$path[genes_all$cell == c & !genes_all$path %in% genes_all$path[genes_all$cell != c]]) / 
      n_distinct(genes_all$path[genes_all$cell == c])
  })
)

# Update the cell-cell matrix with gene ratios
cell_cell_matrix$n_genes_ratio2 = cell_type_stats[cell_cell_matrix$cell2, "pct_unique_genes"]
cell_cell_matrix[cell_cell_matrix$cell1 == cell_cell_matrix$cell2, "n_genes_ratio"] = cell_cell_matrix[cell_cell_matrix$cell1 == cell_cell_matrix$cell2, "n_genes_ratio2"]
cell_cell_matrix$cell1 <- paste(as.character(cell_cell_matrix$cell1), cell_type_stats[as.character(cell_cell_matrix$cell1), "n_total_genes"], sep=" ")
cell_cell_matrix$cell2 <- paste(as.character(cell_cell_matrix$cell2), cell_type_stats[as.character(cell_cell_matrix$cell2), "n_total_genes"], sep=" ")

# Load circlize library for visualization
library(circlize)
library(patchwork)

# Create a chord diagram for cell-cell interactions
names(col_x) <- sort(as.character(unique(cell_cell_matrix$cell1)))
p_circos <- wrap_elements(chordDiagram(cell_cell_matrix[, 1:3], transparency=0.5, grid.col=col_x))
circos.clear()

# Save the chord diagram as a TIFF file
tiff(filename="./circos.tiff", width=2080, height=2080)
par(cex=7, mar=c(0, 0, 0, 0))
chordDiagram(cell_cell_matrix[, 1:3], transparency=0.5, grid.col=col_x, 
             annotationTrack=c("grid"), annotationTrackHeight=mm_h(55))

# Add cell names to the diagram
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index=si, track.index=1)
  ylim = get.cell.meta.data("ylim", sector.index=si, track.index=1)
  circos.text(mean(xlim), mean(ylim), si, sector.index=si, track.index=1, 
              facing="bending.inside", niceFacing=TRUE, col="black")
}
dev.off()

# Function to plot TIFF figures
plottiffFigures <- function(path_x) {
  xpos <- 1:5
  ypos <- xpos^2
  data_frame <- data.frame(xpos=xpos, ypos=ypos)
  img <- readTIFF(path_x, native=TRUE)
  img <- rasterGrob(img, interpolate=TRUE)
  
  # Plotting the data
  p2 <- qplot(xpos, ypos, geom="blank") +
    annotation_custom(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(color=I("White")) + 
    theme_void()
  return(p2)
}

# Plot the circos diagram
p_circos <- plottiffFigures("./circos.tiff")

# Prepare data for Venn diagram
venn_data <- genes_all %>%
  group_by(path) %>%
  summarize(n_cell_types = n_distinct(cell)) %>%
  group_by(n_cell_types) %>%
  summarize(n_genes = n())

# Plot the bar chart for shared signaling module genes
p_bar <- ggplot(venn_data, aes(x=n_cell_types, y=n_genes)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=n_genes), vjust=0.1, color="black", size=5) + 
  ylim(0, 500) +  # Add text labels
  labs(x="Number of cell types", y="Number of shared signaling\nmodule genes") +
  theme_hd_plain()

# Load schematic figure
p_schematic <- plottiffFigures("./figures/Figure 4_sub.tif")

# Prepare data for heatmap of signaling module gene expression
data_x = signal_path$all_select$Epithelial
rownames(data_x) <- gsub(" - lacto and neolacto series", "", rownames(data_x))  # Shorten long names
row_km <- cutree(hclust(dist(data_x)), k=5)
col_km <- cutree(hclust(dist(t(data_x))), k=5)
col_fun = colorRamp2(c(min(data_x), median(as.matrix(data_x)), max(data_x)), c("#00BFFF", "white", "#EE7621"))

# Create and draw the heatmap
gb_heatmap <- grid.grabExpr(draw(Heatmap(data_x, 
                                         row_split=row_km,
                                         col=col_fun,
                                         show_row_names=TRUE,
                                         clustering_method_columns="ward.D2",
                                         row_names_max_width=unit(11.5, 'cm'),
                                         clustering_method_rows="ward.D2",
                                         row_gap=unit(1, "mm"),
                                         row_title_rot=0,
                                         name="Spearman's Rho",
                                         row_title=NULL,
                                         column_split=col_km,
                                         row_names_gp=gpar(fontsize=14, fontface="plain"),
                                         top_annotation=HeatmapAnnotation(foo=anno_block(gp=gpar(fill=c("#9AFF9A", "#AEEEEE", "#EECFA1", "#FFBBFF", "pink"), col="white"),
                                                                                         labels=c("Module 1", "Module 2", "Module 3", "Module 4", "Module 5"), 
                                                                                         labels_gp=gpar(col="black", fontsize=14, fontface="plain"))),
                                         column_names_gp=gpar(fontsize=8, fontface="plain"),
                                         column_title="Spearman’s correlation of log2 fold change between signaling module gene expression and local metabolic allocation in epithelial cells",
                                         column_title_gp=gpar(fontsize=14, fontface="plain"),
                                         heatmap_legend_param=list(direction="horizontal", title_position="leftcenter")),
                                 heatmap_legend_side="top", 
                                 padding=unit(c(10, 10, 30, 0), "points")))

# Display the heatmap
ggarrange(gb_heatmap)





# Arrange UMAP plots for global, local, and signaling pathways
p4all1 <- ggarrange(umap_global, umap_local, umap_signal, 
                    nrow=1, ncol=3, 
                    labels=c(letters[c(1:3)]), 
                    font.label=list(size=18, color="black", face="bold"), 
                    common.legend=TRUE, legend="bottom")

# Arrange count, circos, and bar plots
p4all2 <- ggarrange(p_count, p_circos, p_bar, 
                    nrow=3, ncol=1, 
                    labels=c(letters[c(6:8)]), 
                    font.label=list(size=18, color="black", face="bold"), 
                    heights=c(1.5, 2, 1))

# Arrange schematic and bar plots
p4all3 <- ggarrange(p_schematic, pbar_plot$p1, 
                    nrow=1, ncol=2, 
                    labels=c(letters[c(4:5)]), 
                    font.label=list(size=18, color="black", face="bold"), 
                    widths=c(1.5, 2))

# Combine previous arrangements into a larger plot
p4all4 <- ggarrange(p4all1, p4all3, 
                    nrow=2, ncol=1, 
                    font.label=list(size=18, color="black", face="bold"), 
                    heights=c(1, 2))

# Combine the larger plot with count plots
p4all5 <- ggarrange(p4all4, p4all2, 
                    nrow=1, ncol=2, 
                    font.label=list(size=18, color="black", face="bold"), 
                    widths=c(4, 1))

# Combine everything with the heatmap
p4all6 <- ggarrange(p4all5, gb_heatmap, 
                    nrow=2, 
                    labels=c("", letters[c(9)]), 
                    font.label=list(size=18, color="black", face="bold"))

# Save the combined figure as a TIFF file
tiff(paste("./figures/", "figure4_new.tiff", sep=""), 
     height=26, width=20, res=300, units="in", compression="lzw")
print(p4all6)
dev.off()

# Save the combined figure as a PDF file
pdf(paste("./figures/", "figure4_new.pdf", sep=""), 
    height=26, width=20)
print(p4all6)
dev.off()

# Uncomment to save the current workspace
# save.image("/data/Project/Project1/result_new/rdata/step_figure4.RData")

# Calculate Spearman correlation for signaling pathway classes
temp_corsignalClass <- getcorSpearmanForall(cpm_all_select_signal_fc, fc_mainclassNormalized$FC)
saveRDS(temp_corsignalClass, file="./data/signalgenes_correlation_pathwayClass.Rds")
temp_corsignalClass <- readRDS(file="./data/signalgenes_correlation_pathwayClass.Rds")
temp_corsignalClass_sig <- temp_corsignalClass[temp_corsignalClass$pvalue < 0.05, ]

# Filter significant pathways
pathClass_sig <- fc_mainclassNormalized$all_path_sigClass[fc_mainclassNormalized$all_path_sigClass$wcx_p < 0.05, ]
temp_corsignalClass_sig_sig <- temp_corsignalClass_sig[paste(temp_corsignalClass_sig$cell, temp_corsignalClass_sig$gene) %in% 
                                                         paste(result_df_sig$cell_type, result_df_sig$genes) & 
                                                         paste(temp_corsignalClass_sig$cell, temp_corsignalClass_sig$pathways) %in% 
                                                         paste(pathClass_sig$cell, pathClass_sig$path), ]

# Plot heatmap for significant correlations
all_sigClass <- plotHeatmapSignalPath(temp_corsignalClass_sig, "all sig")
signal_pathClass <- GetSigheatmapCor(all_sigClass)
pathwaysClass <- signal_pathClass$paths_all
genes_allClass <- as.data.frame(signal_pathClass$genes_all)

# Create bar plot for global signaling module genes
p_barcount_global <- PlotBarCount(temp_corsignalClass_sig_sig[paste(temp_corsignalClass_sig_sig$cell, 
                                                                    temp_corsignalClass_sig_sig$gene) %in% 
                                                                paste(genes_allClass[, 1], genes_allClass[, 2]) & 
                                                                paste(temp_corsignalClass_sig_sig$cell, 
                                                                      temp_corsignalClass_sig_sig$pathways) %in% 
                                                                paste(pathwaysClass[, 1], pathwaysClass[, 2]), ], "")
p_barcount_global_p1 <- p_barcount_global$p1 + 
  ylab("Number of signaling module genes correlated with changes in global metabolic allocation") + 
  theme(axis.text.y=element_text(size=14, hjust=1))

# Load pathway mean data
fc_pathwayMeanSCT <- readRDS(file="fc_pathwayMeanSCT.rds")

# Load signaling correlation data
temp_corsignalSCT <- readRDS(file="./signalgenes_correlation_pathwaySCT.Rds")
temp_corsignalSCT_sig <- temp_corsignalSCT[temp_corsignalSCT$pvalue < 0.05, ]

# Filter significant pathways for SCT
pathSCT_sig <- fc_pathwayMeanSCT$all_path_sigClass[fc_pathwayMeanSCT$all_path_sigClass$wcx_p < 0.05, ]
temp_corsignalSCT_sig_sig <- temp_corsignalSCT_sig[paste(temp_corsignalSCT_sig$cell, temp_corsignalSCT_sig$gene) %in% 
                                                     paste(result_df_sig$cell_type, result_df_sig$genes) & 
                                                     paste(temp_corsignalSCT_sig$cell, temp_corsignalSCT_sig$pathways) %in% 
                                                     paste(pathSCT_sig$cell, pathSCT_sig$path), ]

# Plot heatmap for significant correlations in SCT
all_sigSCT <- plotHeatmapSignalPath(temp_corsignalSCT_sig, "all sig")
signal_pathSCT <- GetSigheatmapCor(all_sigSCT)
pathwaysSCT <- signal_pathSCT$paths_all
genes_allSCT <- as.data.frame(signal_pathSCT$genes_all)

# Create bar plot for SCT signaling module genes
p_barcount_SCT <- PlotBarCount(temp_corsignalSCT_sig_sig[paste(temp_corsignalSCT_sig_sig$cell, 
                                                               temp_corsignalSCT_sig_sig$gene) %in% 
                                                           paste(genes_allSCT[, 1], genes_allSCT[, 2]) & 
                                                           paste(temp_corsignalSCT_sig_sig$cell, 
                                                                 temp_corsignalSCT_sig_sig$pathways) %in% 
                                                           paste(pathwaysSCT[, 1], pathwaysSCT[, 2]), ], "")
p_barcount_SCT_p1 <- p_barcount_SCT$p1 + 
  ylab("Number of signaling module genes correlated with changes in pathway metabolic activity") + 
  theme(axis.text.y=element_text(size=14, hjust=1))

# Arrange global and SCT bar plots
p4all_sub1 <- ggarrange(p_barcount_global_p1, p_barcount_SCT_p1, 
                        nrow=2, 
                        labels=c(letters[c(1:2)]), 
                        font.label=list(size=18, color="black", face="bold"), 
                        heights=c(1, 4), 
                        common.legend=TRUE, legend="bottom")

# Save the combined bar plot as a TIFF file
tiff(paste("./figures/", "supfigure6.tiff", sep=""), 
     height=23, width=20, res=300, units="in", compression="lzw")
print(p4all_sub1)
dev.off()

# Create heatmap for T cell signaling module gene expression
data_x = signal_path$all_select$`T cell`
rownames(data_x) <- gsub(" - lacto and neolacto series", "", rownames(data_x))
row_km <- cutree(hclust(dist(data_x)), k=2)
col_km <- cutree(hclust(dist(t(data_x))), k=2)
col_fun = colorRamp2(c(min(data_x), median(as.matrix(data_x)), max(data_x)), c("#00BFFF", "white", "#EE7621"))

# Draw the heatmap for T cells
gb_heatmap_t <- grid.grabExpr(draw(Heatmap(data_x, 
                                           row_split=row_km, 
                                           col=col_fun, 
                                           show_row_names=TRUE, 
                                           clustering_method_columns="ward.D2", 
                                           row_names_max_width=unit(11.5, 'cm'), 
                                           clustering_method_rows="ward.D2", 
                                           row_gap=unit(1, "mm"), 
                                           row_title_rot=0, 
                                           name="Spearman's Rho", 
                                           row_title=NULL, 
                                           column_split=col_km, 
                                           row_names_gp=gpar(fontsize=14, fontface="plain"), 
                                           column_names_gp=gpar(fontsize=8, fontface="plain"), 
                                           column_title="Spearman’s correlation of log2 fold change between signaling module gene expression and local metabolic allocation in T cells", 
                                           column_title_gp=gpar(fontsize=14, fontface="plain"), 
                                           heatmap_legend_param=list(direction="horizontal", title_position="leftcenter")), 
                                   heatmap_legend_side="top", 
                                   padding=unit(c(10, 10, 30, 0), "points")))

# Display the heatmap for T cells
ggarrange(gb_heatmap_t)

# Create heatmap for plasma signaling module gene expression
data_x = signal_path$all_select$Plasma
rownames(data_x) <- gsub(" - lacto and neolacto series", "", rownames(data_x))
row_km <- cutree(hclust(dist(data_x)), k=2)
col_km <- cutree(hclust(dist(t(data_x))), k=2)
col_fun = colorRamp2(c(min(data_x), median(as.matrix(data_x)), max(data_x)), c("#00BFFF", "white", "#EE7621"))

# Draw the heatmap for plasma
gb_heatmap_Plasma <- grid.grabExpr(draw(Heatmap(data_x, 
                                                row_split=row_km, 
                                                col=col_fun, 
                                                show_row_names=TRUE, 
                                                clustering_method_columns="ward.D2", 
                                                row_names_max_width=unit(11.5, 'cm'), 
                                                clustering_method_rows="ward.D2", 
                                                row_gap=unit(1, "mm"), 
                                                row_title_rot=0, 
                                                name="Spearman's Rho", 
                                                row_title=NULL, 
                                                column_split=col_km, 
                                                row_names_gp=gpar(fontsize=14, fontface="plain"), 
                                                column_names_gp=gpar(fontsize=8, fontface="plain"), 
                                                column_title="Spearman’s correlation of log2 fold change between signaling module gene expression and local metabolic allocation in Plasma", 
                                                column_title_gp=gpar(fontsize=14, fontface="plain"), 
                                                heatmap_legend_param=list(direction="horizontal", title_position="leftcenter")), 
                                        heatmap_legend_side="top", 
                                        padding=unit(c(10, 10, 30, 0), "points")))

# Display the heatmap for plasma
ggarrange(gb_heatmap_Plasma)

# Arrange T cell and plasma heatmaps
p4all_sub2 <- ggarrange(gb_heatmap_t, gb_heatmap_Plasma, 
                        nrow=2, 
                        labels=c(letters[c(1:2)]), 
                        font.label=list(size=18, color="black", face="bold"), 
                        heights=c(1, 1))

# Save the combined heatmap as a TIFF file
tiff(paste("./figures/", "supfigure7.tiff", sep=""), 
     height=23, width=20, res=300, units="in", compression="lzw")
print(p4all_sub2)
dev.off()




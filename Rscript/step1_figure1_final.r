



# Load the required libraries
library("jpeg")
library("ggplot2")
library("grid")
library("tiff")
library("reshape2")
library("dplyr")
library("ggpubr")  # For combining plots
library("RColorBrewer")  # For color palettes

# Load the combined filtered data
combined.sct_filter <- readRDS(file = "./data/all_meta_6cohorts_NGpaper_integrated_result_meta_genes_filtered.rds")
all_dataset_SCT_filter <- readRDS(file = "./data/all_dataset_SCT_filter.rds")

# Standardize cell type names
all_dataset_SCT_filter$Cell_type[all_dataset_SCT_filter$Cell_type == "Myloid"] <- "Myeloid"
all_dataset_SCT_filter$cell_clustered[all_dataset_SCT_filter$cell_clustered == "Myloid"] <- "Myeloid"
all_dataset_SCT_filter$Cell_type[all_dataset_SCT_filter$Cell_type == "Bcell"] <- "B cell"
all_dataset_SCT_filter$Cell_type[all_dataset_SCT_filter$Cell_type == "Tcell"] <- "T cell"
all_dataset_SCT_filter$cell_clustered[all_dataset_SCT_filter$cell_clustered == "Bcell"] <- "B cell"
all_dataset_SCT_filter$cell_clustered[all_dataset_SCT_filter$cell_clustered == "Tcell"] <- "T cell"
combined.sct_filter$Cell_type[combined.sct_filter$Cell_type == "Bcell"] <- "B cell"
combined.sct_filter$Cell_type[combined.sct_filter$Cell_type == "Tcell"] <- "T cell"
combined.sct_filter$Cell_type[combined.sct_filter$Cell_type == "Myloid"] <- "Myeloid"

# Load gene pathway information
all_gene_pathway <- readRDS(file = "./data/all_genes_pathway_infor_select.rds")
# Filter pathways to include only those present in the dataset
all_gene_pathway <- all_gene_pathway[all_gene_pathway$gene %in% rownames(all_dataset_SCT_filter@assays$SCT@counts), ]

# Load custom plotting functions
source("./Rscript/functions_plot.r")

############################## Plot flowchart
# Set positions for flowchart points
xpos <- 1:5
ypos <- xpos ** 2
data_frame <- data.frame(xpos = xpos, ypos = ypos)

print("Data points")
print(data_frame)

# Load the image for the flowchart
path <- "./figures/figure1a.tif"
img <- readTIFF(path, native = TRUE)
img <- rasterGrob(img, interpolate = TRUE)

# Create the flowchart plot
p1_flowchart <- qplot(xpos, ypos, geom = "blank") +
  annotation_custom(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(color = I("White")) + theme_void()

############################## Plot UMAP
# Extract UMAP embeddings and combine with cell type metadata
umap <- combined.sct_filter@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(cell_type = combined.sct_filter@meta.data$Cell_type, class_new = combined.sct_filter@meta.data$class_new)

# Calculate median UMAP coordinates for each cell type
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

# Create UMAP plot
p2_umap <- DimPlot(combined.sct_filter, label = FALSE, group.by = "class_new") +
  ggtitle("Metabolic genes") +
  scale_color_manual(values = c("#00BFFF", "#EE7621")) +
  theme_void(base_size = 14, base_family = "Helvetica")

# Add arrows and annotations to the UMAP plot
p2_umap <- p2_umap + 
  geom_segment(aes(x = min(umap$UMAP_1), y = min(umap$UMAP_2), 
                   xend = min(umap$UMAP_1) + 5, yend = min(umap$UMAP_2)), 
               colour = "black", size = 0.1, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = min(umap$UMAP_1), y = min(umap$UMAP_2), 
                   xend = min(umap$UMAP_1), yend = min(umap$UMAP_2) + 5), 
               colour = "black", size = 0.1, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x = min(umap$UMAP_1) + 4, y = min(umap$UMAP_2) - 1.5, 
           label = "UMAP 1", color = "black", size = 4.7) + 
  annotate("text", x = min(umap$UMAP_1) - 1.5, y = min(umap$UMAP_2) + 4, 
           label = "UMAP 2", color = "black", size = 4.7, angle = 90) +
  geom_text(aes(label = cell_type, x = UMAP_1, y = UMAP_2, size = I(4.7), fontface = "plain"), 
            data = cell_type_med) +
  theme(plot.title = element_text(lineheight = .8, size = 14, hjust = 0.5), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        legend.text = element_text(size = 14))

######################################## Plot marker gene dotplot
# Define marker genes for each cell type
markers <- as.data.frame(rbind(
  cbind(cell = "B cell", gene = c("MS4A1", "CD19")),
  cbind(cell = "Plasma", gene = c("CD79A", "MZB1", "JCHAIN")),
  cbind(cell = "T cell", gene = c("CD3E", "CD3D", "CD3G")),
  cbind(cell = "Epithelial", gene = c("EPCAM", "MUC1", "CDH1")),
  cbind(cell = "Mast", gene = c("HDC", "CPA3", "TPSAB1")),
  cbind(cell = "Myeloid", gene = c("FCGR3A", "CD14", "CD163")),
  cbind(cell = "Endothelial", gene = c("PECAM1", "VWF", "CLDN5")),
  cbind(cell = "Fibroblast", gene = c("FBLN1", "PDGFRA", "FGF7"))
))

# Split markers by cell type
x <- split(markers$gene, markers$cell)
x <- x[c(1, 7, 8, 3, 5, 6, 2, 4)]

# Set factor levels for cell clustering
all_dataset_SCT_filter$cell_clustered <- factor(all_dataset_SCT_filter$cell_clustered, 
                                                levels = c("B cell", "Plasma", "T cell", "Epithelial", "Mast", "Myeloid", "Endothelial", "Fibroblast"))

# Create a dot plot for marker genes
p1_dot <- DotPlot(object = all_dataset_SCT_filter, features = x, group.by = "cell_clustered") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  ggtitle("")

# Set color palette for the dot plot
cc1 <- brewer.pal(n = 8, name = "Set3")
cc1 <- cc1[c(1:4, 8, 5:7)]

# Create a customized dot plot with specified colors
p2_dot <- DotPlot(object = all_dataset_SCT_filter, features = x, group.by = "cell_clustered", 
                  cols = c("lightgrey", "#EE7621")) +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), 
        strip.text = element_text(angle = 60), 
        axis.title = element_blank(), 
        legend.position = "right") + 
  ggtitle("") + 
  guides(size = guide_legend(title = "Percent\nexpressed"))
p2_dot$guides$colour$title <- "Average\nexpression"

############################# Plot cell type percentage barplot
# Create a table of cell type percentages
x <- as.matrix(table(paste(combined.sct_filter$class_new, combined.sct_filter$dataset, sep = "_"), 
                     combined.sct_filter$Cell_type))
x1 <- apply(x, 1, function(x) x / sum(x))
x2 <- reshape2::melt(x1)
x2$type <- x2$Var2

# Create a bar plot for cell type composition
p3_cell <- ggplot(x2, aes(x = Var2, fill = Var1, y = value)) +
  geom_bar(stat = "identity") +
  theme_hd_plain() +
  theme(legend.position = "right") +
  xlab("Sample type") +
  ylab("Cell type composition") +
  labs(fill = "Cell Type") +
  coord_flip() +
  scale_fill_manual(values = brewer.pal(n = 8, name = "Set3")[c(1:4, 8, 5:7)])

############################# Plot cell type percentage barplot (by patient and dataset)
combined.sct_filter$Class <- combined.sct_filter$class_new
x <- as.matrix(table(paste(combined.sct_filter$Class, combined.sct_filter$patient, combined.sct_filter$dataset, sep = "_"), 
                     combined.sct_filter$Cell_type))
x1 <- apply(x, 1, function(x) x / sum(x))
x2 <- reshape2::melt(x1)
x2$type <- x2$Var2

# Create a bar plot for frequency of cell types
p3 <- ggplot(x2, aes(x = Var2, fill = Var1, y = value)) +
  geom_bar(stat = "identity") +
  theme_hd() +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1)) +
  ggtitle("Cell type composition") +
  xlab("Sample type") +
  ylab("Frequency") +
  labs(fill = "Cell Type")

############################# Plot cell type percentage boxplot
# Extract dataset and type from Var2
x2$dataset <- str_split_fixed(x2$Var2, "_", 3)[, 3]
x2$type <- str_split_fixed(x2$Var2, "_", 3)[, 1]

# Create a box plot for cell type composition
p2_cell <- ggplot(x2[x2$Var1 != "Mast", ], aes(x = type, y = value, color = dataset)) +
  geom_boxplot() +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        strip.text = element_text(angle = 60), 
        plot.margin = unit(c(1, 2, 1, 2), "cm"), 
        axis.title.x = element_blank()) + 
  ylab("Cell type composition") + 
  facet_wrap(~Var1, nrow = 1) + 
  geom_signif(comparisons = list(c("Normal", "Tumor")), map_signif_level = TRUE, color = "black") + 
  scale_color_manual(name = "Dataset", values = c("#7FDBB3", "#F5AD75", "#85A4E6", "#F593CE", "#CAF77B"))

################################## Plot correlation heatmap
# Load gene pathway data again
all_gene_pathway <- readRDS(file = "/data/Project/Project1/data/all_genes_pathway_infor_select.rds")

# Extract unique metagenes
metagenes <- unique(all_gene_pathway$gene)
all_genes <- rownames(all_dataset_SCT_filter@assays$SCT@counts)

# Create a raw data matrix for the relevant genes
raw_data <- as.matrix(all_dataset_SCT_filter@assays$SCT@counts[all_genes[all_genes %in% metagenes], ])
merge_metaSCT_raw <- all_dataset_SCT_filter@meta.data[colnames(all_dataset_SCT_filter@assays$SCT@counts), ]

# Calculate mean expression of genes by cell type and dataset
mean_genes_raw <- aggregate(t(raw_data), list(s = as.character(merge_metaSCT_raw$Cell_type), 
                                              p = merge_metaSCT_raw$dataset, 
                                              type = merge_metaSCT_raw$class_new), 
                            function(x) mean(x))
mean_genes_df_raw <- mean_genes_raw[, -c(1:3)]
rownames(mean_genes_df_raw) <- paste(mean_genes_raw[, 1], mean_genes_raw[, 2], mean_genes_raw[, 3], sep = "_")

# Plot correlation heatmap
cor_raw <- plotCorHeatmap2(mean_genes_df_raw, "Correlation")

################################# Plot metabolic activity 
# Aggregate metabolic activity data
meta_data <- all_dataset_SCT_filter@meta.data
meta_aggre <- aggregate(meta_data[, "metabolic_activity"], 
                        list(patient = meta_data$patient, 
                             type = meta_data$class_new, 
                             cell = meta_data$Cell_type, 
                             dataset = meta_data$dataset), 
                        mean)

# Create box plot for total metabolic activity
p_mat <- ggplot(meta_aggre[meta_aggre$cell != "Mast", ], aes(x = type, y = x, col = dataset)) +
  geom_boxplot() +
  theme_hd_plain() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        plot.margin = unit(c(1, 2, 1, 2), "cm"), 
        strip.text = element_text(angle = 60), 
        axis.title.x = element_blank()) + 
  ylab("Total metabolic activity") + 
  facet_wrap(~cell, nrow = 1) + 
  scale_color_manual(name = "Dataset", values = c("#7FDBB3", "#F5AD75", "#85A4E6", "#F593CE", "#CAF77B")) + 
  geom_signif(comparisons = list(c("Normal", "Tumor")), map_signif_level = TRUE, color = "black")

# Combine plots into a single figure layout
pall0 <- ggarrange(p1_flowchart, p2_umap, nrow = 1, ncol = 2, labels = letters[1:2], 
                   font.label = list(size = 18, color = "black", face = "bold"))

pall1 <- ggarrange(p2_dot, p3_cell, nrow = 2, ncol = 1, heights = c(1.4, 1), 
                   labels = letters[c(3, 5)], font.label = list(size = 18, color = "black", face = "bold"))

pall2 <- ggarrange(pall1, cor_raw$pp, labels = c("", letters[4]), 
                   font.label = list(size = 18, color = "black", face = "bold"))

pall3 <- ggarrange(p2_cell, p_mat, labels = letters[6:7], 
                   font.label = list(size = 18, color = "black", face = "bold"), 
                   common.legend = TRUE, legend = "bottom")

# Final combined figure layout
figure1_p <- ggarrange(pall0, pall2, pall3, nrow = 3, ncol = 1, heights = c(1.2, 1.4, 1))

# Save the combined figure as a PDF
pdf(paste("/data/Project/Project1/figures/", "figure1.pdf", sep = ""), height = 20, width = 18)
print(figure1_p)
dev.off()

# Save the combined figure as a TIFF
tiff(paste("/data/Project/Project1/figures/", "figure1.tiff", sep = ""), height = 24, width = 19, res = 300, units = "in", compression = "lzw")
print(figure1_p)
dev.off()

# Create violin plots for RNA features
p_vln1 <- VlnPlot(all_dataset_SCT_filter, features = c("nFeature_RNA"), group.by = "patient", split.by = "class_new", pt.size = 0) +
  scale_fill_manual(values = c("#00BFFF", "#EE7621"))
p_vln2 <- VlnPlot(all_dataset_SCT_filter, features = c("nCount_RNA"), group.by = "patient", split.by = "class_new", pt.size = 0) +
  scale_fill_manual(values = c("#00BFFF", "#EE7621"))
p_vln3 <- VlnPlot(all_dataset_SCT_filter, features = c("percent.mt"), group.by = "patient", split.by = "class_new", pt.size = 0) +
  scale_fill_manual(values = c("#00BFFF", "#EE7621"))

# Define a function to create violin plots for different datasets
p_vln1 <- plotVlnplot(all_dataset_SCT_filter@meta.data, "nFeature_RNA", 'Number of detected genes')
p_vln2 <- plotVlnplot(all_dataset_SCT_filter@meta.data, "nCount_RNA", 'UMI counts')
p_vln3 <- plotVlnplot(all_dataset_SCT_filter@meta.data, "percent.mt", 'Percentage of mitochondrial genes')

# Extract UMAP embeddings and combine with cell type, dataset, and patient metadata
umap <- combined.sct_filter@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = combined.sct_filter@meta.data$Cell_type, 
        dataset = combined.sct_filter@meta.data$dataset, 
        patient = combined.sct_filter@meta.data$patient) 

# Define a color palette for patient groups
colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
           "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffbb78", "#ff9896", 
           "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#f5b0d3", 
           "#c7c7c7", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffbb78", 
           "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", 
           "#f5b0d3", "#c7c7c7", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
           "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", 
           "#ffbb78", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", 
           "#9edae5", "#f5b0d3"
)

# Create a UMAP plot colored by patient
p_patient_umap <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = patient)) +
  geom_point(size = 0.001) +
  theme_hd_plain() +
  facet_wrap(~dataset, nrow = 1) +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(nrow = 4, override.aes = list(size = 4))) +
  xlab('UMAP 1') +
  ylab('UMAP 2')

# Arrange violin plots and UMAP plot into a single figure
ps1_1 <- ggarrange(p_vln1, p_vln2, p_vln3, nrow = 3, ncol = 1, 
                   labels = letters[1:3], 
                   font.label = list(size = 18, color = "black", face = "bold"))

ps1_2 <- ggarrange(ps1_1, p_patient_umap, nrow = 2, ncol = 1, 
                   labels = c("", letters[4]), 
                   font.label = list(size = 18, color = "black", face = "bold"), 
                   heights = c(3, 1))

# Save the combined figure as a PDF
pdf(paste("/data/Project/Project1/figures/", "supfigure1.pdf", sep = ""), height = 20, width = 17)
print(ps1_2)
dev.off()

# Save the combined figure as a TIFF
tiff(paste("/data/Project/Project1/figures/", "supfigure1.tiff", sep = ""), height = 20, width = 17, res = 300, units = "in", compression = "lzw")
print(ps1_2)
dev.off()


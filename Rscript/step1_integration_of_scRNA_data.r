# Load custom plotting functions from an external R script
source("/data/Project/functions/functions_plot.r")

# Specify the path to the HDF5 file containing the data
my_h5_files = "./data/GSE178341_crc10x_full_c295v4_submit.h5"

# Read the count data from the HDF5 file
counts = readH5(my_h5_files)

# Clean up column names by removing the "data|" prefix
colnames(counts) <- gsub("data\\|", "", colnames(counts))

# Read gene features from the HDF5 file
genes <- h5read(my_h5_files, "matrix/features")
# Convert gene features to a data frame with ensembl_id and symbol
genes <- as.data.frame(tibble(ensembl_id = genes$id, symbol = genes$name))

# Set row names for genes and counts
rownames(genes) <- genes$ensembl_id
rownames(counts) <- genes[rownames(counts), "symbol"]

# Create a Seurat object for single-cell analysis
GSE178341 <- CreateSeuratObject(counts = counts, project = "GSE178341", min.cells = 3, min.features = 200)

# Read cluster and metadata information from CSV files
crc10x_full_c295v4_submit_cluster <- read.csv("./data/GSE178341_crc10x_full_c295v4_submit_cluster.csv/crc10x_full_c295v4_submit_cluster.csv")
GSE178341_crc10x_full_c295v4_submit_metatables <- read.csv("./data/GSE178341_crc10x_full_c295v4_submit_metatables.csv/GSE178341_crc10x_full_c295v4_submit_metatables.csv")

# Merge cluster and metadata tables by sampleID and cellID
gse178341_annotation = merge(crc10x_full_c295v4_submit_cluster, GSE178341_crc10x_full_c295v4_submit_metatables, by.x = "sampleID", by.y = "cellID")

# Classify samples as 'Normal' or 'Tumor' based on HistologicTypeSimple
gse178341_annotation$type = ifelse(gse178341_annotation$HistologicTypeSimple == "Normal colon", "Normal", "Tumor")
gse178341_annotation$patient = gse178341_annotation$PID
gse178341_annotation$Cell_type = gse178341_annotation$clMidwayPr

# Set row names of the annotation data frame to sample IDs
rownames(gse178341_annotation) <- gse178341_annotation[, 1]

# Add the metadata to the Seurat object
GSE178341 <- AddMetaData(object = GSE178341, metadata = gse178341_annotation)

# Save the Seurat object to an RDS file
saveRDS(GSE178341, file = "./data/GSE178341_seurat.rds")

# Load the saved Seurat object from the RDS file
GSE178341 <- readRDS("./data/GSE178341_seurat.rds")

# Assign identifiers for original patient types and classes
GSE178341$orig.ident = GSE178341$PatientTypeID
GSE178341$Patient = GSE178341$PID
GSE178341$Class = GSE178341$type

# Calculate the percentage of mitochondrial genes
GSE178341[["percent.mt"]] <- PercentageFeatureSet(GSE178341, pattern = "^MT")

# Create a new column to indicate if the patient has matched normal samples
GSE178341$matched_patient <- ifelse(GSE178341$PID %in% unique(GSE178341$PID[GSE178341$type == "Normal"]), "yes", "no")

# Count the number of cells per patient type
cell_count <- table(GSE178341$PatientTypeID)

# Mark patients with fewer than 1000 cells as 'low_Cell'
GSE178341$low_Cell <- ifelse(GSE178341$PID %in% str_split_fixed(names(cell_count)[cell_count < 1000], "_", 2)[, 1], "yes", "no")

# Subset the Seurat object based on quality control metrics and annotations
GSE178341 <- subset(GSE178341, percent.mt < 25 & nFeature_RNA > 300 & nCount_RNA >= 1000 & 
                      nFeature_RNA < 7500 & matched_patient == "yes" & low_Cell == "no" & Cell_type != "other")

# Split the Seurat object by PatientTypeID
ifnb.list <- SplitObject(GSE178341, split.by = "PatientTypeID")

# Remove duplicates from each split object
ifnb.list = lapply(X = ifnb.list, FUN = DeleteDoubles2)

# Merge the split objects back together
GSE178341_merged <- merge(ifnb.list[[1]], ifnb.list[-1])

# Save the filtered merged Seurat object
saveRDS(GSE178341_merged, file = "./data/GSE178341_seurat_filtered.rds")

# Load the filtered merged Seurat object
GSE178341_merged <- readRDS(file = "./data/GSE178341_seurat_filtered.rds")

# Specify another HDF5 file for additional count data
h5_files = "./data/Count_matrix_allCells_373058_noGeneFiltered.h5"

# Read the count data from the new HDF5 file
counts = Read10X_h5(h5_files)

# Read metadata for epithelial cells and set row names
Epithelial_metadata <- read.csv("./data/singlecell/Epithelial_metadata.csv")
rownames(Epithelial_metadata) <- Epithelial_metadata$cell.ID

# Read metadata for non-epithelial cells and set row names
NonEpithelial_metadata <- read.csv("./data/singlecell/NonEpithelial_metadata.csv")
rownames(NonEpithelial_metadata) <- NonEpithelial_metadata$cell.ID

# Combine epithelial and non-epithelial metadata into one data frame
all_meta <- rbind(Epithelial_metadata, cbind(NonEpithelial_metadata, iCMS = "", msi = ""))

# Create a Seurat object for all cells with specified parameters
all_cells <- CreateSeuratObject(counts = counts, dataset = "all", min.cells = 3, min.features = 200)

# Add combined metadata to the Seurat object
all_cells <- AddMetaData(object = all_cells, metadata = all_meta[colnames(counts), -c(2:3)])

# Save the Seurat object with all metadata
saveRDS(all_cells, file = "./data/all_meta_5cohorts_NGpaper_no_filters.rds")

# Load the saved Seurat object
all_cells <- readRDS(file = "./data/all_meta_5cohorts_NGpaper_no_filters.rds")

# Assign new class and patient identifiers
all_cells$class_new = all_cells$sample.origin
all_cells$patient = all_cells$patient.ID

# Standardize tumor class names
all_cells$class_new[all_cells$class_new == "Tumor-2"] = "Tumor"

# Create a column to indicate if the patient has matched normal samples
all_cells$matched_patient <- ifelse(all_cells$patient.ID %in% unique(all_cells$patient.ID[all_cells$class_new == "Normal"]), "yes", "no")

# Subset the Seurat object to exclude lymph node samples and unmatched patients
all_cells <- subset(all_cells, matched_patient == "yes" & sample.origin != "LymphNode")

# Count the number of cells per patient and class
cell_count <- table(paste(all_cells$patient.ID, all_cells$class_new, sep = "_"))

# Mark patients with fewer than 1000 cells as 'low_Cell'
all_cells$low_Cell <- ifelse(all_cells$patient.ID %in% str_split_fixed(names(cell_count)[cell_count < 1000], "_", 2)[, 1], "yes", "no")

# Create a unique identifier for patient types
all_cells$PatientTypeID <- paste(all_cells$patient.ID, all_cells$class_new, sep = "_")

# Generate a violin plot for mitochondrial percentage by cell type
VlnPlot(all_cells, features = c("percent.mt"), group.by = "cell.type", pt.size = 0, split.by = "class_new") + 
  scale_fill_manual(values = c("#4F94CD", "#EE6AA7")) + 
  geom_hline(yintercept = 10)

# Further subset the Seurat object based on quality control metrics and annotations
all_cells <- subset(all_cells, percent.mt < 25 & nFeature_RNA > 300 & nCount_RNA >= 1000 & 
                      nFeature_RNA < 7500 & matched_patient == "yes" & low_Cell == "no")

# Rename cell type column for consistency
all_cells$Cell_type <- all_cells$cell.type

# Split the Seurat object by PatientTypeID
ifnb.list <- SplitObject(all_cells, split.by = "PatientTypeID")

# Remove duplicates from each split object
ifnb.list = lapply(X = ifnb.list, FUN = DeleteDoubles2)

# Merge the split objects back together
all_cells_merged <- merge(ifnb.list[[1]], ifnb.list[-1])

# Save the merged Seurat object
saveRDS(all_cells_merged, file = "./data/all_meta_5cohorts_NGpaper_filters.rds")

# Load the saved merged Seurat object
all_cells_merged <- readRDS(file = "./data/all_meta_5cohorts_NGpaper_filters.rds")

# Combine the merged datasets from different sources
all_dataset <- merge(all_cells_merged, GSE178341_merged)

# Assign default values for missing dataset entries
all_dataset$dataset[is.na(all_dataset$dataset)] = "Pelka"
all_dataset$patient = all_dataset$patient.ID
all_dataset$patient[is.na(all_dataset$patient)] = all_dataset$PID[is.na(all_dataset$patient)]

# Update sample origin for specific dataset entries
all_dataset$sample.origin[all_dataset$dataset == "Pelka"] = all_dataset$Class[all_dataset$dataset == "Pelka"]

# Standardize class names
all_dataset$class_new = all_dataset$sample.origin
all_dataset$class_new[all_dataset$class_new == "Tumor-2"] = "Tumor"

# Create a column to indicate if the patient has matched normal samples
all_dataset$matched_patient <- ifelse(all_dataset$patient %in% unique(all_dataset$patient[all_dataset$class_new == "Normal"]), "yes", "no")

# Subset the combined dataset to exclude unmatched patients and lymph node samples
all_dataset <- subset(all_dataset, matched_patient == "yes" & sample.origin != "LymphNode")

# Standardize cell type names based on dataset origin
all_dataset$cell.type[all_dataset$dataset == "Pelka"] = all_dataset$Cell_type[all_dataset$dataset == "Pelka"]
all_dataset$cell.type[all_dataset$cell.type == "Epi"] = "Epithelial"
all_dataset$cell.type[all_dataset$cell.type == "Fibro"] = "Fibroblast"
all_dataset$cell.type[all_dataset$cell.type == "PlasmaB"] = "Plasma"
all_dataset$cell.type[all_dataset$cell.type == "Endo"] = "Endothelial"
all_dataset$Cell_type[all_dataset$Cell_type %in% c("Peri", "SmoothMuscle")] = "Fibroblast"

# Preserve old cell type information for reference
all_dataset$Cell_type_old = all_dataset$Cell_type

# Update the Cell_type based on standardized cell type names
all_dataset$Cell_type = all_dataset$cell.type
all_dataset$Cell_type[all_dataset$Cell_type %in% c("T_NK", "TCD4", "ILC", "TCD8", "Tgd", "TZBTB16", "NK")] = "Tcell"
all_dataset$Cell_type[all_dataset$Cell_type %in% c("pDC", "McDC", "DC", "Mono", "Neutrophils", "Macro", "NK", "Granulo")] = "Myloid"

# Save the final dataset to an RDS file
saveRDS(all_dataset, file = "./data/all_meta_6cohorts_NGpaper_no_filter.rds")

# Load the saved dataset
all_dataset <- readRDS(file = "./data/all_meta_6cohorts_NGpaper_no_filter.rds")

# Split the dataset for SCTransform
ifnb.list <- SplitObject(all_dataset, split.by = "dataset")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2")

# Merge the transformed datasets
all_dataset_SCT = merge(ifnb.list[[1]], ifnb.list[-1])

# Prepare for marker finding
all_dataset_SCT <- PrepSCTFindMarkers(all_dataset_SCT)

# Save the prepared dataset
saveRDS(all_dataset_SCT, file = "./data/all_dataset_prepSCT.Rds")

# Load gene pathway information
all_gene_pathway <- readRDS(file = "./data/all_genes_pathway_infor_select.rds")

# Get unique metagenes
metagenes = unique(all_gene_pathway$gene)

# Extract raw expression matrix for genes of interest
all_genes = rownames(all_dataset_SCT@assays$RNA@counts)
subset.matrix <- all_dataset_SCT@assays$RNA@counts[all_genes[all_genes %in% metagenes], ]

# Create a new Seurat object with just the genes of interest
object2 <- CreateSeuratObject(subset.matrix, min.cells = 5, min.features = 100)

# Add metadata to the new object
orig.ident <- all_dataset_SCT@meta.data
object2 <- AddMetaData(object = object2, metadata = orig.ident)

# Run all necessary Seurat analyses
object2 <- runSeuratAll2(object2)

# Save the new Seurat object
saveRDS(object2, file = "./data/all_merged_6_cohort_meta_genes_no_filter.rds")

# Load the saved object
object2 <- readRDS(file = "./data/all_merged_6_cohort_meta_genes_no_filter.rds")

# Split the object for SCTransform
ifnb.list <- SplitObject(object2, split.by = "dataset")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2")

# Select integration features
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)

# Prepare for integration
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

# Run PCA on the datasets
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)

# Find integration anchors
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, reference = c(1), normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca")

# Save the anchors
saveRDS(immune.anchors, file = "./data/all_meta_6cohorts_NGpaper_integrated_result_temp_meta_genes.rds")

# Clean up
rm(ifnb.list)
gc()

# Integrate the data using the anchors
combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)

# Run PCA on the integrated data
combined.sct <- RunPCA(combined.sct, verbose = FALSE)

# Run UMAP and t-SNE for visualization
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- RunTSNE(combined.sct, reduction = "pca", dims = 1:30)

# Find neighbors and clusters
combined.sct <- FindNeighbors(combined.sct, dims = 1:30)
combined.sct <- FindClusters(combined.sct)

# Save the integrated result
saveRDS(combined.sct, file = "./data/all_meta_6cohorts_NGpaper_integrated_result_meta_genes.rds")

# Load the integrated result
combined.sct <- readRDS(file = "./data/all_meta_6cohorts_NGpaper_integrated_result_meta_genes.rds")

# Load the prepared dataset
all_dataset_SCT <- readRDS(file = "./data/all_dataset_prepSCT.Rds")

# Create a cell cluster identifier
all_dataset_SCT$cell_clustered = paste("", as.character(combined.sct@meta.data[rownames(all_dataset_SCT@meta.data), "seurat_clusters"]), sep = "")

# Define markers for different cell types
markers <- as.data.frame(rbind(
  cbind(cell = "Bcell", gene = c("IGHD", "MS4A1", "CD19")),
  cbind(cell = "Plasma", gene = c("CD79A", "MZB1", "CD38")),
  cbind(cell = "T cells", gene = c("CD3E", "CD3D", "CD3G")),
  cbind(cell = "Epithelial cells", gene = c("EPCAM", "MUC1", "PROM1", "CDH1")),
  cbind(cell = "Mast cells", gene = c("HDC", "CPA3", "TPSAB1")),
  cbind(cell = "Myeloids", gene = c("FCGR3A", "CD14", "CD163", "S100A9", "IDO1", "CCR7")),
  cbind(cell = "Endothelial", gene = c("PECAM1", "VWF", "CLDN5")),
  cbind(cell = "Fibroblast", gene = c("FBLN1", "PDGFRA", "FGF7"))
))

# Split markers by cell type
x <- split(markers$gene, markers$cell)
x = x[c(1, 7, 8, 3, 5, 6, 2, 4)]

# Create dot plot and dimensionality plot
p2 = DotPlot(object = all_dataset_SCT, features = x, group.by = "cell_clustered") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) + 
  ggtitle("")

p1 = DimPlot(combined.sct, label = TRUE) + ggtitle("Metabolic Gene Profile")

# Update

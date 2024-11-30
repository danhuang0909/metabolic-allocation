
# Load necessary libraries
library(stringr)
library(writexl)

# Load filtered datasets
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
meta_data <- all_dataset_SCT_filter@meta.data
all_gene_pathway <- readRDS(file = "./data/all_genes_pathway_infor_select.rds")
all_gene_pathway <- all_gene_pathway[all_gene_pathway$gene %in% rownames(all_dataset_SCT_filter@assays$SCT@counts),]

# Split genes by pathway
x <- split(all_gene_pathway$gene, all_gene_pathway$final_name)
allcells_pathwayMeanSCT <- pathwayMean(all_dataset_SCT_filter@assays$SCT@counts, x)
saveRDS(allcells_pathwayMeanSCT, file = "./data/allcells_pathwayMeanSCT.Rds")

# Aggregate by KEGG class
temp <- unique(all_gene_pathway[, c("gene", "kegg_class")])
x <- split(temp$gene, temp$kegg_class)
allcells_pathwayMeanSCTMainClass <- pathwayMean(all_dataset_SCT_filter@assays$SCT@counts, x)
saveRDS(allcells_pathwayMeanSCTMainClass, file = "./data/allcells_pathwayMeanSCTMainClass.Rds")

# Extract metabolic data and normalize
sct_metabolic <- all_dataset_SCT_filter@assays$SCT@counts[rownames(all_dataset_SCT_filter@assays$SCT@counts) %in% all_gene_pathway$gene, ]
sct_metabolic_scale <- apply(sct_metabolic, 2, function(x) x / sum(x))

# Calculate pathway means normalized by metabolic counts
all_gene_pathway <- all_gene_pathway[all_gene_pathway$gene %in% rownames(all_dataset_SCT_filter@assays$SCT@counts), ]
x <- split(all_gene_pathway$gene, all_gene_pathway$final_name)
sct_metabolic_scale_pathwaymean <- pathwayMean(sct_metabolic_scale, x)

temp <- unique(all_gene_pathway[, c("gene", "kegg_class")])
x <- split(temp$gene, temp$kegg_class)
sct_metabolic_scale_pathwaymeanMainClass <- pathwayMean(sct_metabolic_scale, x)

# Save normalized pathway means
saveRDS(sct_metabolic_scale_pathwaymean, file = "./data/allcells_pathwayMeanSCTNormalizedByMetacount.Rds")
saveRDS(sct_metabolic_scale_pathwaymeanMainClass, file = "./data/allcells_pathwayMeanSCTMainClassNormalizedByMetacount.Rds")

# Filter out Mast cells and run Seurat reaction
filter_mast <- meta_data[meta_data$Cell_type != "Mast", ]
pathway_mean_allsct <- RunSeuratReaction(allcells_pathwayMeanSCT[, rownames(filter_mast)], filter_mast)
saveRDS(pathway_mean_allsct, file = "./data/pathway_mean_allsct_filter_mast.Rds")

pathway_mean_allsctClass <- RunSeuratReaction(allcells_pathwayMeanSCTMainClass[, rownames(filter_mast)], filter_mast)
saveRDS(pathway_mean_allsctClass, file = "./data/pathway_mean_allsct_filter_mastClass.Rds")

pathway_mean_allsctClassNormalized <- RunSeuratReaction(sct_metabolic_scale_pathwaymeanMainClass[, rownames(filter_mast)], filter_mast)
saveRDS(pathway_mean_allsctClassNormalized, file = "./data/pathway_mean_allsct_filter_mastClassNormalizedByMetacount.Rds")

# Load pathway means for further analysis
allcells_pathwayMeanSCT <- readRDS(file = "./data/allcells_pathwayMeanSCT.Rds")
allcells_pathwayMeanSCTMainClass <- readRDS(file = "./data/allcells_pathwayMeanSCTMainClass.Rds")
allcells_pathwayMeanSCTMainClassNormalized <- readRDS(file = "./data/allcells_pathwayMeanSCTMainClassNormalizedByMetacount.Rds")
allcells_pathwayMeanSCTNormalized <- readRDS(file = "./data/allcells_pathwayMeanSCTNormalizedByMetacount.Rds")

pathway_mean_allsct_seu <- readRDS(file = "./data/pathway_mean_allsct_filter_mast.Rds")
pathway_mean_allsct_seu$Cell_type <- droplevels(pathway_mean_allsct_seu$Cell_type)

pathway_mean_allsct_seuClass <- readRDS(file = "./data/pathway_mean_allsct_filter_mastClass.Rds")
pathway_mean_allsct_seuClass$Cell_type <- droplevels(pathway_mean_allsct_seuClass$Cell_type)

pathway_mean_allsctClassNormalized <- readRDS(file = "./data/pathway_mean_allsct_filter_mastClassNormalizedByMetacount.Rds")
pathway_mean_allsctClassNormalized$Cell_type <- droplevels(pathway_mean_allsctClassNormalized$Cell_type)

# Function to aggregate data by metadata
getAggregateDate2 <- function(data_x, meta_data) {
  data_y <- aggregate(t(data_x), 
                      list(dataset = meta_data[colnames(data_x), ]$dataset, 
                           type = meta_data[colnames(data_x), ]$class_new, 
                           cell = meta_data[colnames(data_x), ]$Cell_type, 
                           patient = meta_data[colnames(data_x), ]$patient), 
                      mean)
  data_y_df <- data_y[, -c(1:4)]
  rownames(data_y_df) <- paste(data_y[, 1], data_y[, 2], data_y[, 3], data_y[, 4], sep = "_")
  return(data_y_df)
}

# Aggregate normalized pathway means
pathway_classNormalized <- getAggregateDate2(allcells_pathwayMeanSCTMainClassNormalized, pathway_mean_allsctClassNormalized@meta.data)
pathway_class <- getAggregateDate2(allcells_pathwayMeanSCTMainClass, pathway_mean_allsct_seuClass@meta.data)

# Calculate pathway ratios normalized by patient
count_x <- as.matrix(pathway_mean_allsct_seu@assays$RNA@counts)
pathwayRatioNormalized <- c()

for (path_name in unique(all_gene_pathway$kegg_class)) {
  pathway_patient_temp <- count_x[rownames(count_x) %in% all_gene_pathway[all_gene_pathway$kegg_class == path_name, "final_name"], , drop = FALSE]
  temp_x <- apply(pathway_patient_temp, 2, function(x) x / sum(x))
  pathwayRatioNormalized <- rbind(pathwayRatioNormalized, temp_x)
}

pathwayRatioNormalized[is.na(pathwayRatioNormalized)] <- 0
saveRDS(as.data.frame(pathwayRatioNormalized), file = "./pathwayRatioNormalized.rds")

# Aggregate normalized pathway ratios by metadata
samples_x <- colnames(pathwayRatioNormalized)
pathwayRatioNormalizedpatient <- aggregate(t(as.matrix(pathwayRatioNormalized)), 
                                           list(dataset = meta_data[samples_x, "dataset"], 
                                                type = meta_data[samples_x, "class_new"], 
                                                cell = meta_data[samples_x, "Cell_type"], 
                                                patient = meta_data[samples_x, "patient"]), 
                                           mean)

rownames(pathwayRatioNormalizedpatient) <- paste(pathwayRatioNormalizedpatient[, 1], 
                                                 pathwayRatioNormalizedpatient[, 2], 
                                                 pathwayRatioNormalizedpatient[, 3], 
                                                 pathwayRatioNormalizedpatient[, 4], 
                                                 sep = "_")
pathwayRatioNormalizedpatient <- pathwayRatioNormalizedpatient[, -c(1:4)]
saveRDS(as.data.frame(pathwayRatioNormalizedpatient), file = "./pathwayRatioNormalizedpatient.rds")

# Separate Normal and Tumor samples
pathwayRatioNormalized_N <- pathwayRatioNormalizedpatient[grepl("Normal", rownames(pathwayRatioNormalizedpatient)), ]
rownames(pathwayRatioNormalized_N) <- gsub("Normal_", "", rownames(pathwayRatioNormalized_N))
pathwayRatioNormalized_T <- pathwayRatioNormalizedpatient[grepl("Tumor", rownames(pathwayRatioNormalizedpatient)), ]
rownames(pathwayRatioNormalized_T) <- gsub("Tumor_", "", rownames(pathwayRatioNormalized_T))

# Find overlapping samples
overlap_s <- intersect(rownames(pathwayRatioNormalized_N), rownames(pathwayRatioNormalized_T))

# Initialize data frames for distance calculations
diss_nt <- data.frame(matrix(ncol = 0, nrow = length(overlap_s)))
fc_nt <- data.frame(matrix(ncol = 0, nrow = length(overlap_s)))

# Function to calculate column distances
column_distances <- function(matrix1, matrix2, distance_function) {
  n_cols <- ncol(matrix1)
  distances <- c()
  for (i in 1:n_cols) {
    distances <- c(distances, distance_function(matrix1[, i], matrix2[, i]))
  }
  return(distances)
}

# Euclidean distance function
euclidean_distance <- function(x, y) {
  sqrt(sum((x - y)^2))
}

# Calculate distances and fold changes for each KEGG class
for (i in unique(all_gene_pathway$kegg_class)) {
  pathway_normal_temp <- t(pathwayRatioNormalized_N[overlap_s, colnames(pathwayRatioNormalized_N) %in% all_gene_pathway[all_gene_pathway$kegg_class == i, "final_name"]])
  pathway_tumor_temp <- t(pathwayRatioNormalized_T[overlap_s, colnames(pathwayRatioNormalized_T) %in% all_gene_pathway[all_gene_pathway$kegg_class == i, "final_name"]])
  
  temp <- str_split_fixed(colnames(pathway_normal_temp), "_", 4)
  temp_nt <- c()
  temp_fc <- c()
  
  for (cc in unique(temp[, 3])) {
    n_t <- pathway_normal_temp[, temp[, 3] == cc]
    t_t <- pathway_tumor_temp[, temp[, 3] == cc]
    dd <- as.data.frame(column_distances(n_t, t_t, euclidean_distance))
    
    n_t[n_t == 0] <- min(n_t[n_t != 0]) / 2
    t_t[t_t == 0] <- min(t_t[t_t != 0]) / 2
    fc_x <- as.data.frame(abs(log2(apply(t_t / n_t, 2, mean))))
    
    rownames(dd) <- colnames(n_t)
    rownames(fc_x) <- colnames(n_t)
    temp_nt <- rbind(temp_nt, dd)
    temp_fc <- rbind(temp_fc, fc_x)
  }
  
  colnames(temp_nt) <- i
  colnames(temp_fc) <- i
  diss_nt <- cbind(diss_nt, temp_nt)
  fc_nt <- cbind(fc_nt, temp_fc)
}

# Save distance and fold change results
saveRDS(diss_nt, file = "./ratio_distance.rds")

# Save the current workspace
save.image("./data/load_all_files2.rdata")

# Aggregate gene counts for all samples
all_genes <- rownames(all_dataset_SCT_filter@assays$RNA@counts)
subset.matrix <- all_dataset_SCT_filter@assays$RNA@counts
metadata <- all_dataset_SCT_filter@meta.data[colnames(subset.matrix), ]
sum_gene_set <- aggregate(t(subset.matrix), 
                          list(patient = metadata$patient, 
                               class = metadata$class_new, 
                               cell = metadata$Cell_type, 
                               dataset = metadata$dataset), 
                          sum)
sum_gene_set_df <- sum_gene_set[, -c(1:4)]
rownames(sum_gene_set_df) <- paste(sum_gene_set$dataset, sum_gene_set$class, sum_gene_set$cell, sum_gene_set$patient, sep = "_")
saveRDS(sum_gene_set_df, file = "./sum_all_genes_count_rawRNA.rds")

# Write aggregated gene counts to Excel
write_xlsx(cbind(genename = rownames(sum_gene_set_df), as.data.frame(sum_gene_set_df)), 'sum_all_genes_count_rawRNA.xlsx')

# Transpose the sum gene set data frame
sum_gene_set_df <- t(readRDS("./sum_all_genes_count_rawRNA.rds"))
temp_names <- str_split_fixed(colnames(sum_gene_set_df), "_", 4)
cell_types <- unique(temp_names[, 3])

# Initialize result data frame
result_df <- c()
cpm_all <- c()

# Process each cell type
for (cell_type in cell_types) {
  print(cell_type)
  data_df_cell <- sum_gene_set_df[, temp_names[, 3] == cell_type]
  groups <- temp_names[temp_names[, 3] == cell_type, 2]
  temp_cell <- getEdgeRdata(data_df_cell, groups)
  
  cpm_all <- rbind(cpm_all, temp_cell$count)
  result_df <- rbind(result_df, cbind(temp_cell$test$table, genes = rownames(temp_cell$test$table), cell_type = cell_type))
}

# Save results
saveRDS(result_df, file = "./sum_all_genes_count_rawRNA_test.rds")
saveRDS(cpm_all, file = "./sum_all_genes_count_rawRNA_cpm.rds")
write_xlsx(cbind(genename = colnames(cpm_all), as.data.frame(t(cpm_all))), 'sum_all_genes_count_rawRNA_logNormalizedByEdgeR.xlsx')

# Aggregate pseudo data counts
pseudo_all <- aggregate(t(subset.matrix), 
                        list(patient = metadata$patient, 
                             type = metadata$class_new, 
                             dataset = metadata$dataset), 
                        sum)
pseudo_all_df_count <- pseudo_all[, -c(1:3)]
rownames(pseudo_all_df_count) <- paste(pseudo_all$patient, pseudo_all$type, pseudo_all$dataset, sep = "_")
saveRDS(pseudo_all_df_count, file = "./data/pseudo_all_df_count.Rds")

# Load COAD data and process
coad_data <- readRDS(file = "./data/TCGA_COAD_DATA.rds")
count_data <- aggregate(coad_data$count[, -c(1:3)], 
                        list(gene = coad_data$count$gene_name), 
                        mean)
rownames(count_data) <- count_data[, 1]
count_data <- count_data[, -1]

# Process sample names
temp <- as.data.frame(str_split_fixed(colnames(count_data), "-", 5))
temp$sample <- paste(temp[, 1], temp[, 2], temp[, 3], substr(temp[, 4], 1, 2), sep = ".")
temp$patient <- paste(temp[, 1], temp[, 2], temp[, 3], sep = "-")
temp$type <- ifelse(temp[, 4] == "11A", "Normal", ifelse(temp[, 4] %in% c("01A", "01B", "01C"), "Tumor", "others"))
rownames(temp) <- colnames(count_data)

# Get EdgeR data for TCGA
tcga_result <- getEdgeRdata(count_data, temp$type)
tcga_infor <- temp
saveRDS(tcga_result, file = "./data/tcga_result.Rds")

# Extract patient information
patient_infor <- as.data.frame(unique(all_dataset_SCT_filter@meta.data[, c("MMRStatus", "msi", "patient", "TissueSiteSimple", "MMRMLH1Tumor", "TumorStage", "dataset")]))
patient_infor <- patient_infor[!(is.na(patient_infor$MMRStatus) & (is.na(patient_infor$msi) | patient_infor$msi == "")), ]
patient_infor$MSI <- patient_infor$msi
patient_infor$MSI[is.na(patient_infor$MSI)] <- patient_infor$MMRStatus[is.na(patient_infor$MSI)]
patient_infor$MSI[patient_infor$MSI %in% c("MMRd", "MSI-H")] <- "MSI"
patient_infor$MSI[patient_infor$MSI == "MMRp"] <- "MSS"
rownames(patient_infor) <- patient_infor$patient

# Load additional patient information
single_colon_cancer_patient_infor <- as.data.frame(read_excel("/data/Project/Project1/data/single_colon_cancer_patient_infor.xlsx"))
rownames(single_colon_cancer_patient_infor) <- single_colon_cancer_patient_infor$patient.ID


# Merge and process patient information
patient_infor$sideness <- single_colon_cancer_patient_infor[rownames(patient_infor), "Sidedness"]

# Fill in missing sidedness values
patient_infor$sideness[is.na(patient_infor$sideness)] <- patient_infor$TissueSiteSimple[is.na(patient_infor$sideness)]
patient_infor$sideness[patient_infor$sideness %in% c("L", "left")] <- "Left"
patient_infor$sideness[patient_infor$sideness %in% c("R", "right")] <- "Right"

# Process tumor stage information
patient_infor$stage <- single_colon_cancer_patient_infor[rownames(patient_infor), "Group Stage"]
patient_infor$stage[is.na(patient_infor$stage)] <- patient_infor$TumorStage[is.na(patient_infor$stage)]
patient_infor$stage[patient_infor$stage %in% c(2, "I", "II", "IIA", "IIB", "pT1", "pT2")] <- "Early"
patient_infor$stage[patient_infor$stage %in% c(3, "III", "IIIC", "IIIA", "IIIB", "pT3", "pt4a", "pT4a", "pT4b")] <- "Advanced"

# Add transcriptomic information
patient_infor$iCMS.transcriptomic <- single_colon_cancer_patient_infor[rownames(patient_infor), "iCMS.transcriptomic"]

# Save the processed patient information
saveRDS(patient_infor, file = "./data/single_cell_patient_infor.rds")

# Load ICMS gene list
ICMS_genelist <- as.data.frame(read_excel("./data/ICMS_genelist.xlsx"))
rownames(ICMS_genelist) <- ICMS_genelist$GENE

# Classify ICMS states
ICMS_genelist$ICMS_state <- ifelse(ICMS_genelist$ICMS %in% c("iCMS2_Up", "iCMS3_Down"), "iCMS2", "iCMS3")
ICMS_genelist <- as.matrix(ICMS_genelist)

# Extract CPM data for tumor samples
cpm_pseudo_tumor_icms <- t(cpm_pseudo[grepl("Tumor", rownames(cpm_pseudo)), colnames(cpm_pseudo) %in% ICMS_genelist[, "GENE"]])

# Scale and process ICMS data
cpm_pseudo_tumor_icms_data <- getCMSscaledData(cpm_pseudo_tumor_icms, ICMS_genelist, "Pseudo")

# Extract and scale CPM data for epithelial tumor samples
subset.matrix_sct <- all_dataset_SCT_filter@assays$SCT@counts
cpm_epi_tumor_icms <- as.data.frame(subset.matrix_sct[rownames(subset.matrix_sct) %in% ICMS_genelist[, "GENE"], 
                                                      rownames(all_dataset_SCT_filter@meta.data)[all_dataset_SCT_filter@meta.data$class_new == "Tumor" & 
                                                                                                   all_dataset_SCT_filter@meta.data$Cell_type == "Epithelial"]])

# Scale data by dataset
data_x <- cpm_epi_tumor_icms
datasets <- all_dataset_SCT_filter@meta.data[colnames(data_x), "dataset"]

for (dataset in unique(datasets)) {
  dataset_samples <- colnames(data_x)[datasets == dataset]
  data_x[, dataset_samples] <- t(scale(t(data_x[, dataset_samples]), center = TRUE, scale = TRUE))
}

# Replace NA values with 0
data_x[is.na(data_x)] <- 0

# Aggregate mean data by ICMS state
data_x_mean <- aggregate(data_x, list(types = ICMS_genelist[rownames(data_x), "ICMS_state"]), mean)
data_x_mean2 <- t(data_x_mean[, -1])
colnames(data_x_mean2) <- data_x_mean[, 1]

# Further aggregate by patient
data_x_mean_all <- aggregate(data_x_mean2, list(types = all_dataset_SCT_filter@meta.data[rownames(data_x_mean2), "patient"]), mean)
rownames(data_x_mean_all) <- data_x_mean_all[, 1]
data_x_mean_all <- data_x_mean_all[, -1]

# Save aggregated results
saveRDS(data_x_mean_all, file = "./data/icms_SCT_count.rds")
saveRDS(cpm_pseudo_tumor_icms_data, file = "./data/cpm_pseudo_tumor_icms_data.rds")

pre_class<-getCellPreForAllclass(pathway_classNormalized)
 
saveRDS(pre_class,file="./data/prediction_data.rds")

ratioNormalize_pre<-getCellPreForAll(pathwayRatioNormalizedpatient)
saveRDS(ratioNormalize_pre,file="ratioNormalize_prediction_all.rds")


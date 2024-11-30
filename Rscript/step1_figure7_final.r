# Load necessary libraries
library(forestploter)
library(forestplot)
library(readxl)
library(dplyr)
library(survival)
library(forestmodel)

library(survival)

# Load data files
load("./data/load_all_files2.rdata")
setwd("/data/Project/Project1/result_new")
source("./Rscript/functions_plot.r")

# Load patient information and KEGG pathway genes
patient_infor <- readRDS(file="../data/single_cell_patient_infor.rds")
kegg_path_genes <- readRDS(file="./data/all_kegg_pathway_genes.rds")
kegg_signal <- kegg_path_genes[grepl("signal", kegg_path_genes$description),]
all_gene_pathway <- readRDS(file="./data/all_genes_pathway_infor_select.rds")

# Load signal pathway data
signal_path <- readRDS(file="./data/signal_path.rds")

# Create a heatmap for the Spearman correlation coefficients
bk <- c(seq(-0.8, 0.8, by=0.01))
a <- pheatmap(signal_path$all_select$Epithelial,
              col=c(colorRampPalette(colors=c("#00BFFF", "white"))(length(bk)/2),
                    colorRampPalette(colors=c("white", "#EE7621"))(length(bk)/2)),
              legend_breaks=seq(-0.8, 0.8, 0.1),
              breaks=bk,
              main="Spearman's correlation coefficient\nLog2FC(pathway) vs Log2FC(signal genes)",
              fontsize_row=6,
              fontsize_col=4)

# Cluster assignments from the heatmap
labels_x <- c("CGS", "TED", "POS", "ACP", "MIC")
cluster_assignments <- cutree(a$tree_col, k=5)

# Load TCGA results
tcga_result <- readRDS(file="./data/tcga_result.Rds")
cpm_tcga_all <- tcga_result$count

# Process sample identifiers
temp_all <- as.data.frame(str_split_fixed(rownames(cpm_tcga_all), "-", 5))
temp_all$sample <- paste(temp_all[,1], temp_all[,2], temp_all[,3], substr(temp_all[,4], 1, 2), sep=".")
temp_all$patient <- paste(temp_all[,1], temp_all[,2], temp_all[,3], sep="-")
temp_all$type <- ifelse(temp_all[,4] == "11A", "Normal", 
                        ifelse(temp_all[,4] %in% c("01A", "01B", "01C"), "Tumor", "others"))
rownames(cpm_tcga_all) <- paste(temp_all$type, temp_all$patient, sep="_")

# Filter signal genes from TCGA data
cpm_tcga_signal_all <- cpm_tcga_all[, colnames(cpm_tcga_all) %in% kegg_signal$symbol & 
                                      !colnames(cpm_tcga_all) %in% unique(all_gene_pathway$gene)]
temp <- aggregate(cpm_tcga_signal_all, list(samples=rownames(cpm_tcga_signal_all)), mean)
rownames(temp) <- temp[,1]
cpm_tcga_signal_all <- temp[,-1]

# Calculate fold changes for signal data
cpm_tcga_signal_fcAN_all <- getFCdataSinglelogAllnormal(cpm_tcga_signal_all)
cpm_tcga_signal_fcAN_each <- getFCdataSinglelog(cpm_tcga_signal_all)

# Pairing TCGA data
patient_id <- str_split_fixed(rownames(cpm_tcga_signal_all), "_", 2)[,2]
cpm_tcga_signal_all_paired <- cpm_tcga_signal_all[patient_id %in% rownames(cpm_tcga_signal_fcAN_each),]

# PCA analysis
x <- prcomp(cpm_tcga_signal_all_paired)
fviz_pca_ind(x, label="none", col.ind=str_split_fixed(rownames(cpm_tcga_signal_all_paired), "_", 2)[,1])

# PCA for all TCGA data
x <- prcomp(cpm_tcga_signal_all)
p_pca <- fviz_pca_ind(x, label="none", col.ind=str_split_fixed(rownames(cpm_tcga_signal_all), "_", 2)[,1]) +
  scale_color_manual(values=c("#00BFFF", "#EE7621")) + 
  theme_hd_plain() + 
  ggtitle("TCGA cohort")

# Get fold change results for TCGA
cpm_TCGA_fc_allnorm <- GetSignalResultFCBulk(cpm_tcga_signal_fcAN_all, cluster_assignments, paste("TCGA", labels_x))
cpm_TCGA_fc_alleach <- GetSignalResultFCBulk(cpm_tcga_signal_fcAN_each, cluster_assignments, paste("TCGA", labels_x))

# Prepare data for correlation analysis
x1 <- cpm_TCGA_fc_allnorm$mean_x
rownames(x1) <- paste("average", rownames(x1))
all_x <- as.data.frame(t(rbind(cpm_TCGA_fc_alleach$mean_x, x1[, colnames(cpm_TCGA_fc_alleach$mean_x)])))
cor(all_x)

# Scatter plot for pairwise log2FC of genes in modules
p1_tcga_pairwise_average <- ggplot(all_x, aes(x=`average TCGA TED`, y=`TCGA TED`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Pairwise log2FC of genes in modules") +
  xlab("Average log2FC of genes in modules") +
  ggtitle("Log2FC of TED signaling module in TCGA-COAD cohort")

# Calculate Spearman correlation for average log2FC
result_cor_paire_aver <- c()
for(i in labels_x) {
  temp_r <- cor.test(all_x[, paste("average TCGA", i)], all_x[, paste("TCGA", i)], method="spearman")
  result_cor_paire_aver <- rbind(result_cor_paire_aver, c(type=i, rho=round(temp_r$estimate, 2), pvalue=temp_r$p.value))
}

# Convert results to a data frame and format
result_cor_paire_aver <- as.data.frame(result_cor_paire_aver)
result_cor_paire_aver$rho.rho <- as.numeric(result_cor_paire_aver$rho.rho)
result_cor_paire_aver$pvalue <- as.numeric(result_cor_paire_aver$pvalue)
result_cor_paire_aver$type <- factor(result_cor_paire_aver$type, levels=result_cor_paire_aver$type[order(result_cor_paire_aver$rho.rho)])

# Create bar plot for pairwise log2FC vs average log2FC
p_cor_paire_aver <- ggplot(result_cor_paire_aver, aes(x=type, y=rho.rho)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_hd_plain() +
  ggtitle("Pairwise log2FC VS Average log2FC in genes") +
  ylab("Spearman's rho") +
  geom_text(aes(label=rho.rho), vjust=0.1, size=3.5) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Signaling modules")

# Load supplementary clinical data
pan_cancer_TCGA_clinical_data_supplementary <- read_excel("./data/pan_cancer_TCGA_clinical_data_supplementary.xlsx")
colon_clinical <- as.data.frame(pan_cancer_TCGA_clinical_data_supplementary[pan_cancer_TCGA_clinical_data_supplementary$bcr_patient_barcode %in% colnames(cpm_TCGA_fc_allnorm$mean_x),])
colon_clinical_paired <- as.data.frame(pan_cancer_TCGA_clinical_data_supplementary[pan_cancer_TCGA_clinical_data_supplementary$bcr_patient_barcode %in% colnames(cpm_TCGA_fc_alleach$mean_x),])

# Combine clinical data with TCGA fold change data
clinical_all <- cbind(colon_clinical_paired, t(cpm_TCGA_fc_alleach$mean_x)[colon_clinical_paired$bcr_patient_barcode,])
all_re <- c()

# Perform survival analysis
temp_data <- cpm_TCGA_fc_alleach$mean_x[, colnames(cpm_TCGA_fc_alleach$mean_x) %in% colon_clinical$bcr_patient_barcode]
for(path_x in colnames(t(temp_data))) {
  for(i in seq(quantile(clinical_all[, path_x], 0.02), quantile(clinical_all[, path_x], 0.8), 0.05)) {
    for(j in seq(i, quantile(clinical_all[, path_x], 0.8), 0.05)) {
      for(cc in c("OS", "DSS", "DFI", "PFI")) {
        df <- get_surveRaw2(clinical_all, path_x, i, j, cc, paste(cc, ".time", sep=""))
        fit <- survfit(Surv(time, vital_status) ~ FC, data=df$df)
        pp <- as.numeric(str_split_fixed(surv_pvalue(fit)$pval.txt, "= ", 2)[2])
        all_re <- rbind(all_re, c(i, j, cc, path_x, pp))
        if(is.na(pp)) { next }
        if(pp < 0.05) {
          print(c(i, j, cc, surv_pvalue(fit)$pval.txt))
        }
      }
    }
  }
}

# Convert results to a data frame and filter significant results
all_re <- as.data.frame(all_re)
all_re[, 5] <- as.numeric(all_re[, 5])
all_re_sig <- all_re[all_re[, 5] < 0.05,]

# Load ICMS gene list
ICMS_genelist <- as.data.frame(read_excel("./data/ICMS_genelist.xlsx"))
rownames(ICMS_genelist) <- ICMS_genelist$GENE
ICMS_genelist$ICMS_state <- ifelse(ICMS_genelist$ICMS %in% c("iCMS2_Up", "iCMS3_Down"), "iCMS2", "iCMS3")
ICMS_genelist <- as.matrix(ICMS_genelist)

# Prepare TCGA tumor data for ICMS analysis
cpm_tcga_tumor_icms <- t(cpm_tcga_signal_all[grepl("Tumor", rownames(cpm_tcga_signal_all)), 
                                             colnames(cpm_tcga_signal_all) %in% ICMS_genelist[, "GENE"]])
cpm_tcga_tumor_icms_data <- getCMSscaledData(cpm_tcga_tumor_icms, ICMS_genelist, "TCGA")

# Combine clinical data with TCGA fold change and ICMS data
clinical_all <- cbind(colon_clinical, t(cpm_TCGA_fc_allnorm$mean_x)[colon_clinical$bcr_patient_barcode,], 
                      cpm_tcga_tumor_icms_data$mean_x[colon_clinical$bcr_patient_barcode,])

# Perform survival analysis for DFI
path_x <- "TCGA TED"
median_x <- median(clinical_all[, path_x])
df <- get_surveRaw2(clinical_all, path_x, median_x, median_x, "DFI", "DFI.time")
df$df$Log2FC <- df$df$FC
fit <- survfit(Surv(time, vital_status) ~ Log2FC, data=df$df)
p_survival <- ggsurvplot(fit, pval=TRUE)
p_surv <- p_survival$plot + theme_hd_plain() + 
  ggtitle("DFI (Disease free interval)") + 
  scale_color_manual(name="TCGA TED", values=c("#EE7621", "#00BFFF"))


# Load clinical data from a text file
clinical_molecular_public_all <- read.delim("./data/clinical_molecular_public_all.txt")

# Filter clinical data for samples present in TCGA normalized fold change data
clinical_tcga <- clinical_molecular_public_all[clinical_molecular_public_all$sample %in% as.character(cpm_TCGA_fc_allnorm$df_x$Var1), ]
rownames(clinical_tcga) <- clinical_tcga$sample

# Add clinical variables to clinical_all dataframe
clinical_all$cimp <- as.factor(clinical_tcga[rownames(clinical_all), "cimp"])
clinical_all$cms_label <- as.factor(clinical_tcga[rownames(clinical_all), "cms_label"])
clinical_all$MSI <- as.factor(clinical_tcga[rownames(clinical_all), "msi"])

# Process tumor stage information
clinical_all$stage <- clinical_all$ajcc_pathologic_tumor_stage
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA")] <- "Stage I"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage II"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB")] <- "Stage IV"
clinical_all$stage[clinical_all$ajcc_pathologic_tumor_stage %in% c("[Not Available]", "[Discrepancy]")] <- NA

# Define a function for Cox proportional hazards model plotting
PlotCox <- function(clinical_all, time_x, status_x, title_x) {
  # Prepare data for Cox model
  select_Data <- clinical_all %>%
    transmute(time = get(time_x),
              status = get(status_x),
              TED = `TCGA TED`,
              CGS = `TCGA CGS`,
              POS = `TCGA POS`,
              ACP = `TCGA ACP`,
              MIC = `TCGA MIC`,
              iCMS2 = `TCGA iCMS2`,
              iCMS3 = `TCGA iCMS3`,
              Age = age_at_initial_pathologic_diagnosis,
              Sex = gender,
              Stage = stage,
              MSI = MSI,
              CMS_label = cms_label,
              CIMP = cimp)
  
  # Fit the Cox model
  model <- coxph(Surv(time, status) ~ ., select_Data)
  
  # Create a forest plot for the Cox model
  p_multi <- forest_model(model, format_options = forest_model_format_options(text_size = 3.5), 
                          recalculate_width = FALSE, recalculate_height = FALSE)
  
  # Function to prepare data for the Cox plot
  getDataCoxPlot <- function(df_data, title_x) {
    re_c <- c()
    cols <- c("level", "n", "estimate", "std.error", "p.value", "conf.low", "conf.high")
    
    for (i in 1:nrow(df_data)) {
      if (is.na(df_data[i, "variable"]) || df_data[i, "class"] == "numeric") {
        temp <- df_data[i, cols]
        if (df_data[i, "class"] == "numeric") {
          temp[1, 1] = df_data[i, 1]
        } else {
          temp[1, 1] = paste(" ", temp[1, 1])
        }
      } else {
        temp <- rbind(df_data[i, cols], df_data[i, cols])
        temp[1, 1] = df_data[i, 1]
        temp[1, 2:7] = NA
        temp[2, 1] = paste(" ", temp[2, 1])
      }
      re_c <- rbind(re_c, temp)
    }
    
    # Calculate hazard ratios and confidence intervals
    re_c$est = exp(re_c$estimate)
    re_c$low = exp(re_c$conf.low)
    re_c$hi = exp(re_c$conf.high)
    re_c$HR = sprintf("%.2f (%.2f to %.2f)", re_c$est, re_c$low, re_c$hi)
    re_c$HR[is.na(re_c$hi) & is.na(re_c$n)] = NA
    re_c$HR[is.na(re_c$hi) & !is.na(re_c$n)] = "Reference"
    re_c$se <- (log(re_c$hi) - log(re_c$est)) / 1.96
    re_c <- as.data.frame(re_c)
    re_c$` ` <- paste(rep(" ", 38), collapse = " ")
    
    # Format p-values
    re_c$p.value <- ifelse(re_c$p.value < 0.01, "<0.01", round(re_c$p.value, 2))
    
    # Prepare data for forest plot
    x = as.data.frame(re_c[, c(1, 2, 13, 11, 5)])
    if (max(abs(re_c$hi / re_c$low), na.rm = TRUE) > 100) {
      x_t <- "log10"
      t_a <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
    } else {
      x_t <- "log2"
      t_a <- NULL
    }
    
    x[is.na(x)] = ""
    colnames(x) <- c("Variable", "N", "", "Hazard ratio", "P-value")
    
    # Create forest plot
    tm <- forest_theme(base_size = 14, refline_col = "red", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    p = forest(as.data.frame(x),
               est = re_c$est,
               lower = re_c$low,
               upper = re_c$hi,
               ci_column = 3,
               x_trans = x_t, ticks_at = t_a)
    
    # Highlight significant rows in red
    n_sig = which(as.numeric(re_c[, "p.value"]) < 0.05)
    pp <- edit_plot(p, row = n_sig, gp = gpar(col = "red4", fontface = "italic"))
    
    # Add title to the plot
    pp <- insert_text(pp, text = title_x, col = 3, part = "header")
    
    return(list(re_c = re_c, x = x, pp = pp))
  }
  
  re_x <- getDataCoxPlot(p_multi$data, paste(status_x, "(multivariate cox proportional hazards model)"))
  
  # Perform univariate Cox models for each variable
  vars_for_table <- colnames(select_Data)[-c(1:2)]
  univ_formulas <- sapply(vars_for_table, function(x) as.formula(paste('Surv(time, status) ~', x)))
  univ_models <- lapply(univ_formulas, function(x) { coxph(x, data = select_Data) })
  
  # Create forest plot for univariate models
  p_uni = forest_model(model_list = univ_models, covariates = vars_for_table, merge_models = TRUE, 
                       format_options = forest_model_format_options(text_size = 3.5), 
                       recalculate_width = FALSE, recalculate_height = FALSE)
  
  # Additional processing for univariate results can be added here
  # ...
  
  return(list(multi_plot = re_x$pp, uni_plot = p_uni))
}







# Perform Cox proportional hazards model plots for different time metrics
p_dfi <- PlotCox(clinical_all, "DFI.time", "DFI")
p_pfi <- PlotCox(clinical_all, "PFI.time", "PFI")
p_os <- PlotCox(clinical_all, "OS.time", "OS")

# Arrange and display the plots for overall survival (OS)
ggarrange(p_os$p_multi, p_os$p_uni)

# Arrange and display the plots for disease-free interval (DFI)
ggarrange(p_dfi$p_multi, p_dfi$p_uni)

# Arrange and display the plots for progression-free interval (PFI)
ggarrange(p_pfi$p_multi, p_pfi$p_uni)

# Load combined single-cell data
combined_data <- readRDS(file = "./combined_data_single_cell.rds")

# Load gene pathway information
all_gene_pathway <- readRDS(file = "./data/all_genes_pathway_infor_select.rds")

# Load raw count data and filter for epithelial cells
cpm_all <- readRDS(file = "./sum_all_genes_count_rawRNA_cpm.rds")
cpm_epi_all <- cpm_all[grepl("Epithelial", rownames(cpm_all)), ]
cpm_epi_signal_fcAN_each <- getFCdataSinglelog(cpm_epi_all)
rownames(cpm_epi_signal_fcAN_each) <- str_split_fixed(rownames(cpm_epi_signal_fcAN_each), "_", 3)[, 3]

# Get signaling results for epithelial cells
cpm_epi_fc_signal <- GetSignalResultFCBulk(cpm_epi_signal_fcAN_each, cluster_assignments, paste("Epithelial", labels_x))
epi_pathway_fc <- getBulkdataRatioPathway(cpm_epi_all, all_gene_pathway)

# Load tumor epithelial data and process
cpm_epi <- cpm_all[grepl("Tumor_Epithelial", rownames(cpm_all)), ]
rownames(cpm_epi) <- str_split_fixed(rownames(cpm_epi), "_", 4)[, 4]
cpm_epi_tumor_signal <- GetSignalResultFCBulk(cpm_epi, cluster_assignments, paste("Epithelial", labels_x))
epi_pathway <- getBulkdataRatioPathway(cpm_epi, all_gene_pathway)

# Prepare data for correlation analysis
x = cpm_epi_fc_signal$mean_x
rownames(x) <- paste("Log2FC", rownames(x))
all_x_epi <- as.data.frame(t(rbind(x, cpm_epi_tumor_signal$mean_x[, colnames(x)])))

# Calculate correlation
cor(all_x_epi)

# Create scatter plot for TED signaling module
p1_epi_pairwise_average <- ggplot(all_x_epi, aes(y = `Log2FC Epithelial TED`, x = `Epithelial TED`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  xlab("Average TED signaling module\nscore in tumor epithelial cells") +
  ylab("Average Log2FC of TED signaling\nmodule in Epithelial cells") +
  ggtitle("scRNA-seq cohort")

# Create scatter plot for MIC signaling module
p1_epi_pairwise_average_mic <- ggplot(all_x_epi, aes(y = `Log2FC Epithelial MIC`, x = `Epithelial MIC`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  xlab("Average MIC signaling module\nscore in tumor epithelial cells") +
  ylab("Average Log2FC of MIC signaling\nmodule in Epithelial cells") +
  ggtitle("scRNA-seq cohort")

# Perform Spearman correlation tests for each signaling module
result_cor_t_fc <- c()
for (i in labels_x) {
  temp_r <- cor.test(all_x_epi[, paste("Epithelial", i)], all_x_epi[, paste("Log2FC Epithelial", i)], method = "spearman")
  result_cor_t_fc <- rbind(result_cor_t_fc, c(type = i, rho = round(temp_r$estimate, 2), pvalue = temp_r$p.value))
}
result_cor_t_fc <- as.data.frame(result_cor_t_fc)
result_cor_t_fc$rho.rho <- as.numeric(result_cor_t_fc$rho.rho)
result_cor_t_fc$pvalue <- as.numeric(result_cor_t_fc$pvalue)
result_cor_t_fc$type <- factor(result_cor_t_fc$type, levels = result_cor_t_fc$type[order(result_cor_t_fc$rho.rho)])

# Create bar plot for correlation results
p_cor_t_fc <- ggplot(result_cor_t_fc, aes(x = type, y = rho.rho)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_hd_plain() +
  ggtitle("Pairwise log2FC VS Tumor ") +
  ylab("Spearman's rho") +
  geom_text(aes(label = rho.rho), vjust = 0.1, size = 3.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Signaling modules in Epithelial\ncells of scRNA-seq cohort")

# Load drug sensitivity data
Model <- read.csv("./data/drugData/Model.csv")
select_colon_sample <- Model[Model$OncotreePrimaryDisease == "Colorectal Adenocarcinoma" & Model$PrimaryOrMetastasis == "Primary", ]
exp <- read.csv("./data/drugData/Expression_Public_24Q2_subsetted.csv")
rownames(exp) <- exp[, 1]
exp_matrix <- apply(exp[, -1], 1, as.numeric)
rownames(exp_matrix) <- colnames(exp[, -1])
exp_matrix <- exp_matrix[, colnames(exp_matrix) %in% select_colon_sample$ModelID]

# Load drug sensitivity data from different sources
drugall_prism <- read.csv("./data/drugData/Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted_NAsdropped.csv", header = TRUE)
drugall_ctd <- read.csv("./data/drugData/Drug_sensitivity_AUC_(CTD^2)_subsetted.csv", header = TRUE)

drug1 <- read.csv("./data/drugData/Drug_sensitivity_AUC_(Sanger_GDSC1)_subsetted.csv", header = TRUE)
drug2 <- read.csv("./data/drugData/Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted_NAsdropped.csv", header = TRUE)

# Merge drug sensitivity data
drugall_gdsc <- merge(drug1, drug2, by = "X")

# Load gene pathway information and calculate pathway ratios
all_gene_pathway <- readRDS(file = "./data/all_genes_pathway_infor_select.rds")
ccl_pathway <- getBulkdataRatioPathway(t(exp_matrix), all_gene_pathway)

# Get signaling results for CCL
cpm_ccl_signal <- GetSignalResultFCBulk(t(exp_matrix), cluster_assignments[names(cluster_assignments) %in% rownames(exp_matrix)], paste("CCL", labels_x))

# Combine CCL signaling data with pathway ratios
combined_ccl_signal_path <- cbind(t(ccl_pathway$pathwayRatioNormalized), t(cpm_ccl_signal$mean_x[, colnames(ccl_pathway$pathwayRatioNormalized)]))

# Prepare data for boxplots
temp_ccl <- as.data.frame(t(cpm_ccl_signal$mean_x))
Model_select <- Model[Model$ModelID %in% rownames(temp_ccl), ]
rownames(Model_select) <- Model_select$ModelID
temp_ccl$MSI <- Model_select[rownames(temp_ccl), "LegacyMolecularSubtype"]
temp_ccl$MSI[temp_ccl$MSI == ""] <- "MSS"

# Create boxplots for TED and MIC signaling modules
p1_ted <- ggplot(temp_ccl, aes(x = MSI, y = `CCL TED`, fill = MSI)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_hd_plain() +
  scale_fill_manual(name = "Type", values = c("#00BFFF", "#EE7621")) +
  ylab("Average TED signaling module\nscore in cancer cell line")

p1_mic <- ggplot(temp_ccl, aes(x = MSI, y = `CCL MIC`, fill = MSI)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_hd_plain() +
  scale_fill_manual(name = "Type", values = c("#00BFFF", "#EE7621")) +
  ylab("Average MIC signaling module\nscore in cancer cell line")


# Get drug result data for different databases
gdsc_drug <- getDrugResult(drugall_gdsc, exp_matrix, temp_ccl)
ctd_drug <- getDrugResult(drugall_ctd, exp_matrix, temp_ccl)
prism_drug <- getDrugResult(drugall_prism, exp_matrix, temp_ccl)

# Create scatter plots for correlations between TED and MIC signaling module scores
p1_ctd <- ggplot(ctd_drug$all_cor_matrix_all, aes(x = `CCL TED`, y = `CCL MIC`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5, label.y = 1) +
  theme_hd_plain() +
  ggtitle(expression(CTD^2)) +
  xlab("Spearman's rho between the average TED\nsignaling module score and drug AUC values") +
  ylab("Spearman's rho between the average MIC\nsignaling module score and drug AUC values")

p1_prism <- ggplot(prism_drug$all_cor_matrix_all, aes(x = `CCL TED`, y = `CCL MIC`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5, label.y = 1) +
  theme_hd_plain() +
  ggtitle("PRISM") +
  xlab("Spearman's rho between the average TED\nsignaling module score and drug AUC values") +
  ylab("Spearman's rho between the average MIC\nsignaling module score and drug AUC values")

p1_gdsc <- ggplot(gdsc_drug$all_cor_matrix_all, aes(x = `CCL TED`, y = `CCL MIC`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5, label.y = 1) +
  theme_hd_plain() +
  ggtitle("GDSC") +
  xlab("Spearman's rho between the average TED\nsignaling module score and drug AUC values") +
  ylab("Spearman's rho between the average MIC\nsignaling module score and drug AUC values")

# Combine correlation results from different databases
x1 = gdsc_drug$all_cor
x = ctd_drug$all_cor
x2 = prism_drug$all_cor
all_drug_cor <- rbind(cbind(x1, database = "GDSC"), cbind(x, database = "CTD^2"), cbind(x2, database = "PRISM"))
all_drug_cor$drug_names <- str_split_fixed(all_drug_cor$drug, "\\.\\.", 3)[, 1]

# Create Venn diagram for drug overlap
pdrug_all <- plotVennResult(ctd_drug$all_cor, gdsc_drug$all_cor, prism_drug$all_cor, "CCL TED", c("CTD^2", "GDSC", "PRISM"), size_x = 0)

# Customize the Venn diagram
pdrug_overlap <- pdrug_all$p1 + 
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  annotate("text", x = 0, y = 0, label = expression(CTD^2), size = 4, angle = 60, hjust = -5, vjust = -3.5) +
  annotate("text", x = 0, y = 0, label = "GDSC", size = 4, hjust = -4, vjust = -30) +
  annotate("text", x = 0, y = 0, label = "PRISM", size = 4, angle = -60, hjust = 0.9, vjust = -34)

# Analyze significant drug correlations
pdrug_sig <- plotVennResult(ctd_drug$all_cor_sig, gdsc_drug$all_cor_sig, prism_drug$all_cor_sig, "CCL MIC", c("CTD^2", "GDSC", "PRISM"))
x3 = pdrug_sig$x3
xx = table(c(unique(x3$`CTD^2`), unique(x3$GDSC), unique(x3$PRISM)))
sort(xx[xx > 1])

# Analyze significant drug correlations for TED
pdrug_sig <- plotVennResult(ctd_drug$all_cor_sig, gdsc_drug$all_cor_sig, prism_drug$all_cor_sig, "CCL TED", c("CTD^2", "GDSC", "PRISM"))
x3 = pdrug_sig$x3
xx = table(c(unique(x3$`CTD^2`), unique(x3$GDSC), unique(x3$PRISM)))
sort(xx[xx > 1])

# Select drugs that are consistent across databases
select_drug = names(xx[xx > 1])
select_drug = select_drug[select_drug != "BEXAROTENE"]  # Exclude inconsistent results

# Filter and reshape correlation results for selected drugs
all_drug_cor_select <- all_drug_cor[all_drug_cor$drug_names %in% c(select_drug, "FLUOROURACIL", "ELOXATIN", "SN38") & all_drug_cor$type == "CCL TED", ]
all_drug_cor_select_matrix <- reshape2::dcast(all_drug_cor_select[, c("drug_names", "database", "rho")], drug_names ~ database, mean)
rownames(all_drug_cor_select_matrix) <- all_drug_cor_select_matrix[, 1]
all_drug_cor_select_matrix <- all_drug_cor_select_matrix[, -1]

# Aggregate mean correlation results
all_drug_cor_select_mean <- aggregate(all_drug_cor_select[, c("rho", "pvalue")], list(drug_names = all_drug_cor_select$drug_names, database = all_drug_cor_select$database), mean)

# Prepare data for visualization
all_drug_cor_select_mean <- all_drug_cor_select_mean %>%
  mutate(fill_value = ifelse(round(pvalue, 2) > 0.05, NA, rho))
all_drug_cor_select_mean$drug_names <- factor(all_drug_cor_select_mean$drug_names, levels = rev(c("FLUOROURACIL", "SN38", "ELOXATIN", "PEVONEDISTAT", "BORTEZOMIB", "AZ.960", "NVP.TAE.684", "BREFELDIN.A", "DINACICLIB")))

# Create heatmap for drug correlations
pp_cor_drugs <- ggplot(all_drug_cor_select_mean, aes(x = database, y = drug_names, fill = fill_value, label = round(rho, 2))) +
  geom_tile() +
  geom_text(aes(alpha = I(ifelse(is.na(fill_value), 0.25, 1))), show.legend = FALSE) +
  coord_fixed(ratio = 1 / sqrt(2)) +
  theme_hd_minimal_plain() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = ifelse(levels(all_drug_cor_select_mean$drug_names) %in% c("FLUOROURACIL", "ELOXATIN", "SN38"), "orange", "black"))) +
  ggtitle("Spearman's rho of TED signaling\nmodule score and drug AUC value") +
  scale_fill_gradient2(name = "Rho", low = "steelblue", mid = "white", high = "#D15D04", midpoint = 0, na.value = "gray") +
  ylab("Drug names") +
  xlab("Databases") +
  scale_x_discrete(labels = c("CTD^2" = expression(CTD^2)))

# Create scatter plots for specific drugs
p1_drug1 <- ggplot(gdsc_drug$ccl_signal_drug, aes(x = `CCL TED`, y = `PEVONEDISTAT..GDSC2.1529.`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  ylab("AUC value of PEVONEDISTAT") +
  ggtitle("GDSC") +
  xlab("Average TED signaling module\nscore in cancer cell line")

p2_drug1 <- ggplot(ctd_drug$ccl_signal_drug, aes(x = `CCL TED`, y = `PEVONEDISTAT..CTRP.411809.`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  ylab("AUC value of PEVONEDISTAT") +
  ggtitle(expression(CTD^2)) +
  xlab("Average TED signaling module\nscore in cancer cell line")

p1_drug2 <- ggplot(gdsc_drug$ccl_signal_drug, aes(x = `CCL TED`, y = `BORTEZOMIB..GDSC1.104.`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  ylab("AUC value of BORTEZOMIB") +
  ggtitle("GDSC") +
  xlab("Average TED signaling module\nscore in cancer cell line")

p2_drug2 <- ggplot(prism_drug$ccl_signal_drug, aes(x = `CCL TED`, y = `BORTEZOMIB..BRD.BRD.K88510285.001.17.8.`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  ylab("AUC value of BORTEZOMIB") +
  ggtitle("PRISM") +
  xlab("Average TED signaling module\nscore in cancer cell line")

p3_drug1 <- ggplot(gdsc_drug$ccl_signal_drug, aes(x = `CCL TED`, y = `SN38..GDSC1.1494.`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  ylab("AUC value of SN38") +
  ggtitle("GDSC") +
  xlab("Average TED signaling module\nscore in cancer cell line")

# Create scatter plots for NICLOSAMIDE
p5_drug2 <- ggplot(ctd_drug$ccl_signal_drug, aes(x = `CCL MIC`, y = `NICLOSAMIDE..CTRP.32653.`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  ylab("AUC value of NICLOSAMIDE") +
  ggtitle(expression(CTD^2)) +
  xlab("Average MIC signaling module\nscore in cancer cell line")

p6_drug2 <- ggplot(prism_drug$ccl_signal_drug, aes(x = `CCL MIC`, y = `NICLOSAMIDE..BRD.BRD.K35960502.001.20.0.`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col = "black") +
  stat_cor(method = "spearman", size = 5) +
  theme_hd_plain() +
  ylab("AUC value of NICLOSAMIDE") +
  ggtitle("PRISM") +
  xlab("Average MIC signaling module\nscore in cancer cell line")

# Arrange and display all drug plots
ggarrange(p1_drug1, p2_drug1, p1_drug2, p2_drug2, nrow = 1, ncol = 4)

# Combine various plots into a single figure
p7all1 <- ggarrange(p_dfi$p_uni, p_dfi$p_multi, nrow = 1, ncol = 2, labels = c(letters[c(1, 3)]), font.label = list(size = 18, color = "black", face = "bold"))
p7all2 <- ggarrange(p_surv, pdrug_overlap, p1_ctd, pp_cor_drugs, p1_drug1, p2_drug1, nrow = 2, ncol = 3, labels = c(letters[c(2, 4:8)]), font.label = list(size = 18, color = "black", face = "bold"))

# Final arrangement of all plots
p7all <- ggarrange(p7all1, p7all2, nrow = 2, ncol = 1, heights = c(2, 2))

# Save the combined figure to a TIFF file
tiff(paste("./figures/", "figure7new.tiff", sep = ""), height = 20, width = 17, res = 300, units = "in", compression = "lzw")
print(p7all)
dev.off()

# Create and save additional figures for supplementary materials
p_mic_os <- ggarrange(p_os$p_uni, p_os$p_multi, nrow = 2, ncol = 1, labels = c(letters[c(1:2)]), font.label = list(size = 18, color = "black", face = "bold"))
tiff(paste("./figures/", "supfigure10.tiff", sep = ""), height = 17, width = 9, res = 300, units = "in", compression = "lzw")
print(p_mic_os)
dev.off()

# Arrange additional plots for supplementary materials
p_ted_mic <- ggarrange(p1_ted, p1_mic, nrow = 1, ncol = 2, labels = c(letters[c(6:7)]), font.label = list(size = 18, color = "black", face = "bold"), legend = "bottom", common.legend = TRUE)
p7ann_sub <- ggarrange(p_pca, p_cor_paire_aver, p1_tcga_pairwise_average, p_cor_t_fc, p1_epi_pairwise_average, p_ted_mic, p1_prism, p1_gdsc, p1_drug2, p2_drug2, p5_drug2, p6_drug2, nrow = 4, ncol = 3, labels = c(letters[c(1:5)], "", letters[8:13]), font.label = list(size = 18, color = "black", face = "bold"))

# Save the final supplementary figure
tiff(paste("./figures/", "supfigure11.tiff", sep = ""), height = 21, width = 20, res = 300, units = "in", compression = "lzw")
print(p7ann_sub)
dev.off()

# Load required data and libraries
load("./data/load_all_files2.rdata")
setwd("/data/Project/Project1/result_new")
source("./Rscript/functions_plot.r")

# Load patient information and KEGG pathway genes
patient_infor <- readRDS(file="./data/single_cell_patient_infor.rds")
kegg_path_genes <- readRDS(file="./data/all_kegg_pathway_genes.rds")
kegg_signal <- kegg_path_genes[grepl("signal", kegg_path_genes$description), ]

# Load normalized pathway data
pathwayRatioNormalized <- readRDS(file="./pathwayRatioNormalized.rds")
fc_pathwayRatioNormalized <- getFCresult(pathwayRatioNormalized, meta_data)

# Load pathway gene sets
pathway_file_all <- "./data/c5.go.bp.v2023.2.Hs.symbols.gmt"
c2 <- read.gmt(pathway_file_all)

# Load CPM and test result data
cpm_all <- readRDS(file="./data/sum_all_genes_count_rawRNA_cpm.rds")
result_df <- readRDS("./data/sum_all_genes_count_rawRNA_test.rds")
result_df <- result_df[result_df$cell_type != "Mast", ]
result_df_sig <- result_df[result_df$PValue < 0.05, ]
cpm_all <- cpm_all[!grepl("Mast", rownames(cpm_all)), ]

# Select signaling pathway related genes
cpm_all_select_signal <- cpm_all[, colnames(cpm_all) %in% kegg_signal$symbol & !colnames(cpm_all) %in% unique(all_gene_pathway$gene)]
cpm_all_select_signal_fc <- getFCdataSinglelog(cpm_all_select_signal)

# Load signal path data
signal_path <- readRDS(file="./data/signal_path.rds")

# Create a heatmap of Spearman correlation coefficients
bk <- c(seq(-0.8, 0.8, by=0.01))
a <- pheatmap(signal_path$all_select$Epithelial, 
              col=c(colorRampPalette(colors=c("#00BFFF", "white"))(length(bk)/2), 
                    colorRampPalette(colors=c("white", "#EE7621"))(length(bk)/2)),
              legend_breaks=seq(-0.8, 0.8, 0.1), 
              breaks=bk, 
              main="Spearman's correlation coefficient\nLog2FC(pathway) vs Log2FC(signal genes)", 
              fontsize_row=6, 
              fontsize_col=4)

# Cluster assignments for heatmap
labels_x <- c("CGS", "TED", "POS", "ACP", "MIC")
cluster_assignments <- cutree(a$tree_col, k=5)

# Prepare KEGG signaling data for network visualization
kegg_signal_select <- kegg_signal[kegg_signal$symbol %in% names(cluster_assignments), c("description", "symbol")]
kegg_signal_select$description <- gsub(" - Homo sapiens \\(human\\)", "", kegg_signal_select$description)
kegg_signal_select$cluster <- paste("Module", cluster_assignments[kegg_signal_select$symbol], sep=" ")

# Create a graph from the selected KEGG signaling data
g <- graph_from_data_frame(kegg_signal_select, directed=FALSE)

# Set node size and color for the network graph
node_size <- ifelse(V(g)$name %in% kegg_signal_select$symbol, 3, 6)
node_color <- ifelse(V(g)$name %in% kegg_signal_select$symbol, 
                     kegg_signal_select$cluster[match(V(g)$name, kegg_signal_select$symbol)], 
                     "Pathway")
V(g)$name[V(g)$name %in% kegg_signal_select$symbol] <- "" 

# Plot the network graph
p_network <- ggraph(g, layout='graphopt') +
  geom_node_point(aes(size=node_size, color=node_color)) +
  geom_node_text(aes(label=name), size=5, vjust=1.5, check_overlap=TRUE) +
  geom_edge_link(color="gray50", alpha=0.5) +
  scale_size_identity() +
  theme_void() + 
  scale_color_manual(name="", values=c("pink", "#AEEEEE", "#FFBBFF", "#9AFF9A", "#EECFA1", "gray")) + 
  guides(color=guide_legend(label.theme=element_text(size=14), override.aes=list(size=3)))

# Define categories of biological processes
cc <- c("CELL_CYCLE|GROWTH_FACTOR|CELLULAR_RESPONSE_TO_STRESS",
        "REGULATION_OF_TRANSPORT|EPITHELIUM_DEVELOPMENT",
        "PROTEIN_PHOSPHORYLATION|CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND|POSITIVE_REGULATION_OF_INTRACELLULAR_SIGNAL_TRANSDUCTION",
        "APOPTOTIC|POSITIVE_REGULATION_OF_CATALYTIC_ACTIVITY|POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION",
        "MOTILITY|INFLAMMATORY_RESPONSE|RESPONSE_TO_CYTOKINE")

# Initialize data frames for module terms and gene information
gene_infor2 <- c()
module_terms <- c()

# Loop through each cluster to gather gene information
for(i in 1:5){
  x=names(cluster_assignments[cluster_assignments==i])
  c2_dpm=c2[c2$gene %in% x,]
  module_terms<-rbind(module_terms,cbind(module=i,term=names(sort(table(c2_dpm$term),decreasing = T)[1:20]),term_n=sort(table(c2_dpm$term),decreasing = T)[1:20]))
  iterms<-table(as.character(unlist(str_split(c2_dpm$term,"_"))))#
  for(nn in cc){
    c2_dpm_dev<-c2_dpm[grepl(nn,c2_dpm$term),]
    gene_infor2<-rbind(gene_infor2,cbind(key=nn,modules=i,gene=unique(c2_dpm_dev$gene)))
  }
}

# Convert module terms to a data frame and calculate shared terms
module_terms <- as.data.frame(module_terms)
module_terms[, 3] <- as.numeric(module_terms[, 3])
x <- sort(table(module_terms[, 2]))
module_terms$shared <- x[module_terms$term]

# Create a matrix of gene fractions in biological processes
x1 <- as.matrix(table(gene_infor2[, 2], gene_infor2[, 1])) / as.numeric(table(cluster_assignments))
colnames(x1) <- paste(c("ACP", "CGS", "MIC", "POS", "TED"), gsub("_", " ", gsub("\\|", "\n        ", colnames(x1))), sep=":")
rownames(x1) <- paste("Module", rownames(x1), sep=" ")

# Plot heatmap of gene fractions
pheatmap(t(x1)[c(2, 5, 4, 1, 3), ], cluster_cols=FALSE, cluster_rows=FALSE)

# Reshape data for plotting
ratio_data <- reshape2::melt(t(x1))
ratio_data$Var1 <- factor(ratio_data$Var1, levels=colnames(x1)[c(2, 5, 4, 1, 3)])

# Create a heatmap for the fraction of genes in biological processes
p_ratio_heatmap <- ggplot(ratio_data, aes(Var2, Var1, fill=value)) +
  geom_tile(color="transparent") +
  scale_fill_gradient2(low="white", mid="white", high="#EE7621", midpoint=0.4) +
  scale_alpha_identity() +
  geom_text(aes(label=sprintf("%.2f", value)), color="black", size=3) +
  theme_hd_minimal_plain() +
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title=element_blank(), legend.position="none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  ggtitle("Fraction of genes in the biological process") +
  scale_y_discrete(name=NULL, position="right")

# Calculate fold changes for signaling results
cpm_fc_each <- GetSignalResultFC(cpm_all_select_signal_fc, cluster_assignments, labels_x)

# Select epithelial data for further analysis
epi_dpm <- as.data.frame(cpm_all_select_signal_fc[grepl("Epi", rownames(cpm_all_select_signal_fc)), colnames(cpm_all_select_signal_fc) %in% names(cluster_assignments[cluster_assignments == 2])])
rownames(epi_dpm) <- str_split_fixed(rownames(epi_dpm), "_", 3)[, 3]

# Create a heatmap for log2 fold changes of signaling modules
bk <- c(seq(-3, 3, by=0.2))
a <- pheatmap(t(cpm_fc_each$df2_x), cluster_cols=FALSE, 
              col=c(colorRampPalette(colors=c("#0295FD", "white"))(length(bk)/2), 
                    colorRampPalette(colors=c("white", "#FF4A52"))(length(bk)/2)),
              legend_breaks=seq(-3, 3, 0.5), 
              breaks=bk, 
              main="Log2FC(signaling modules)", 
              fontsize_row=8, 
              fontsize_col=8, 
              gaps_col=seq(5, 30, 5))

# Load pseudo-bulk data
pseudo_all_df_count <- readRDS(file="./result_meta/rdata_new/pseudo_all_df_count.Rds")
groups <- str_split_fixed(rownames(pseudo_all_df_count), "_", 4)[, 2]
temp_pseudo <- getEdgeRdata(t(pseudo_all_df_count), groups)
cpm_pseudo <- temp_pseudo$count
result_df_pseudo <- temp_pseudo$test

# Select signaling genes in pseudo-bulk data
cpm_pseudo_signal <- cpm_pseudo[, colnames(cpm_pseudo) %in% kegg_signal$symbol & !colnames(cpm_pseudo) %in% unique(all_gene_pathway$gene)]
cpm_pseudo_signal_fc <- getFCdataSinglelog(cpm_pseudo_signal)
rownames(cpm_pseudo_signal_fc) <- str_split_fixed(rownames(cpm_pseudo_signal_fc), "_", 2)[, 1]

# Calculate fold changes for pseudo-bulk signaling results
cpm_pseudo_fc_each <- GetSignalResultFCBulk(cpm_pseudo_signal_fc, cluster_assignments, paste("Pseudo", labels_x))
pseudo_dpm <- as.data.frame(cpm_pseudo_signal_fc[, colnames(cpm_pseudo_signal_fc) %in% names(cluster_assignments[cluster_assignments == 2])])
colnames(pseudo_dpm) <- paste("Pseudo", colnames(pseudo_dpm), sep="_")

# Combine epithelial and pseudo-bulk data for correlation analysis
combined_data_dpm <- as.data.frame(cbind(epi_dpm, pseudo_dpm[rownames(epi_dpm), ]))
result_cor_dpm <- c()

# Calculate Spearman correlations between epithelial and pseudo-bulk data
for(i in colnames(epi_dpm)) {
  temp_r <- cor.test(combined_data_dpm[, i], combined_data_dpm[, paste("Pseudo", i, sep="_")], method="spearman")
  result_cor_dpm <- rbind(result_cor_dpm, c(type=i, rho=round(temp_r$estimate, 2), pvalue=temp_r$p.value))
}

# Convert results to a data frame
result_cor_dpm <- as.data.frame(result_cor_dpm)
result_cor_dpm$rho <- as.numeric(result_cor_dpm$rho)
result_cor_dpm$pvalue <- as.numeric(result_cor_dpm$pvalue)
result_cor_dpm$type <- factor(result_cor_dpm$type, levels=result_cor_dpm$type[order(result_cor_dpm$rho)])

# Create a bar plot for correlation results
p_cor_ted <- ggplot(result_cor_dpm, aes(x=type, y=rho)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_hd_plain() +
  ylim(0, 1.05) +
  ggtitle("Log2FC in epithelial cells VS pseudo-bulk") +
  ylab("Spearman's rho") +
  geom_text(aes(label=rho), hjust=0.15, size=3.5) +
  xlab("Genes in TED modules") +
  coord_flip()

# Calculate correlations for all genes
all_cor_genes <- GenesCorBulk(cpm_all_select_signal_fc, cpm_pseudo_signal_fc, cluster_assignments, labels_x)

# Create a boxplot for gene correlations
p_cor_genes <- ggplot(all_cor_genes, aes(x=modules, y=rho, fill=modules)) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  scale_fill_manual(name="Class", values=c("pink", "#AEEEEE", "#FFBBFF", "#9AFF9A", "#EECFA1")) +
  theme_hd_plain() +
  theme(legend.position="none") +
  ylab("Spearman's rho of log2FC\n(epithelial cells VS pseudo-bulk)") +
  annotate("text", x=1:length(table(all_cor_genes$modules)), y=0, label=table(all_cor_genes$modules)) +
  xlab("Signaling modules")

# Combine data for further analyses
combined_data <- as.data.frame(t(rbind(cpm_fc_each$df2_x, cpm_pseudo_fc_each$mean_x[, colnames(cpm_fc_each$df2_x)])))

# Save combined data to a file
saveRDS(combined_data, file="./combined_data_single_cell.rds")

# Create correlation plots for TED, MIC, and CGS signaling modules
p1_dpm_cor <- ggplot(combined_data, aes(x=`Pseudo TED`, y=`Epithelial TED`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ggtitle("Log2FC of signaling modules")

p1_mic_cor <- ggplot(combined_data, aes(x=`Pseudo MIC`, y=`Epithelial MIC`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ggtitle("Log2FC of signaling modules")

p1_cgs_cor <- ggplot(combined_data, aes(x=`Pseudo CGS`, y=`Epithelial CGS`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ggtitle("Log2FC of signaling modules")

p1_mic_ted_cor <- ggplot(combined_data, aes(x=`Epithelial TED`, y=`Epithelial MIC`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ggtitle("Log2FC of signaling modules")

p1_acp_pos_cor <- ggplot(combined_data, aes(x=`Epithelial POS`, y=`Epithelial ACP`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ggtitle("Log2FC of signaling modules")





# Initialize correlation results
result_cor <- c()

# Calculate Spearman correlation for each signaling module
for(i in labels_x) {
  temp_r <- cor.test(combined_data[, paste("Epithelial", i)], combined_data[, paste("Pseudo", i)], method="spearman")
  result_cor <- rbind(result_cor, c(type=i, rho=round(temp_r$estimate, 2), pvalue=temp_r$p.value))
}

# Convert results to a data frame
result_cor <- as.data.frame(result_cor)
result_cor$rho.rho <- as.numeric(result_cor$rho.rho)
result_cor$pvalue <- as.numeric(result_cor$pvalue)
result_cor$type <- factor(result_cor$type, levels = result_cor$type[order(result_cor$rho.rho)])

# Create a bar plot for Spearman correlation results
p_cor <- ggplot(result_cor, aes(x=type, y=rho.rho)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_hd_plain() +
  ylim(0, 1.05) +
  ggtitle("Log2FC in epithelial cells\nVS pseudo-bulk") +
  ylab("Spearman's rho") +
  geom_text(aes(label=rho.rho), vjust=0.1, size=3.5) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Signaling modules")

# Calculate correlation matrix for epithelial cells
corr <- cor(combined_data[, grepl("Epi", colnames(combined_data))], method="spearman")
rownames(corr) <- gsub("Epithelial ", "", rownames(corr))
colnames(corr) <- gsub("Epithelial ", "", colnames(corr))

# Find maximum correlation values for rows and columns
max_corr_col <- apply(corr, 2, function(x) {
  max(abs(x[x != 1]))
})
max_corr_row <- apply(corr, 1, function(x) {
  max(abs(x[x != 1]))
})

# Prepare correlation data for plotting
corr_data <- as.data.frame(corr) %>%
  rownames_to_column("variable1") %>%
  pivot_longer(-variable1, names_to="variable2", values_to="correlation") %>%
  mutate(is_max_row = abs(correlation) %in% c(1, abs(max_corr_row[variable1])),
         is_max_col = abs(correlation) %in% c(1, abs(max_corr_col[variable2])),
         is_max = is_max_row | is_max_col,
         alpha = if_else(is_max, 1, 0.5))

# Factor levels for plotting
corr_data$variable1 <- factor(corr_data$variable1, levels=labels_x[c(2, 5, 1, 3, 4)])
corr_data$variable2 <- factor(corr_data$variable2, levels=labels_x[c(2, 5, 1, 3, 4)])

# Create a heatmap for correlation data
p_cor_sig_path <- ggplot(corr_data, aes(variable1, variable2, fill=correlation, alpha=alpha)) +
  geom_tile(color="transparent") +
  scale_fill_gradient2(low="#00BFFF", mid="white", high="#EE7621", midpoint=0) +
  scale_alpha_identity() +
  geom_text(aes(label=sprintf("%.2f", correlation)), color="black", size=3) +
  geom_segment(aes(x=0.5, xend=5.5, y=2.5, yend=2.5), color="black", size=0.5) +
  geom_segment(aes(x=2.5, xend=2.5, y=0.5, yend=5.5), color="black", size=0.5) +
  theme_hd_minimal_plain() +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title = element_blank(), legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colors=c("#00BFFF", "white", "#EE7621"), 
                       breaks=c(-1, 0, 1),
                       guide=guide_colorbar(title="Correlation")) +
  ggtitle("Spearman's rho of log2FC(signaling\nmodules) in Epithelial cells")

# Create a heatmap for log2 fold changes of signaling modules
bk <- c(seq(-2, 2, by=0.2))
pp <- pheatmap(combined_data[, c(1:10, 16:35, 11, 14, 13, 12, 15, 39, 38, 40, 37, 36)],
               col=c(colorRampPalette(colors=c("#00BFFF", "white"))(length(bk)/2), 
                     colorRampPalette(colors=c("white", "#EE7621"))(length(bk)/2)),
               legend_breaks=seq(-2, 2, 0.5), breaks=bk, clustering_method="ward.D2", 
               cluster_cols=FALSE, gaps_col=seq(5, 35, 5), main="Log2FC of signaling modules")

# Select epithelial data for further analysis
data_x = signal_path$all_select$Epithelial
row_km <- cutree(hclust(dist(data_x)), k=5)
epi_select = fc_pathwayRatioNormalized$FC[grepl("Epi", rownames(fc_pathwayRatioNormalized$FC)), rownames(signal_path$all_select$Epithelial)]
rownames(epi_select) <- str_split_fixed(rownames(epi_select), "_", 3)[, 3]

# Initialize result for signal correlations
result_cor_signal <- c()

# Calculate Spearman correlation for pathways against each signaling module
for(i in labels_x) {
  x = t(apply(epi_select, 2, getcorSpearman, combined_data[rownames(epi_select), paste("Epithelial", i)]))
  result_cor_signal <- rbind(result_cor_signal, cbind(type=i, pathway=rownames(x), x))
}

# Convert results to a data frame
result_cor_signal <- as.data.frame(result_cor_signal)
result_cor_signal$rho.rho <- as.numeric(result_cor_signal$rho.rho)
result_cor_signal$pvalue <- as.numeric(result_cor_signal$pvalue)

# Reshape the result for plotting
result_cor_signal_df <- reshape2::dcast(result_cor_signal[, 1:3], pathway ~ type)
rownames(result_cor_signal_df) <- result_cor_signal_df[, 1]
result_cor_signal_df <- result_cor_signal_df[, -1]

# Set factor levels for plotting
result_cor_signal$type <- factor(result_cor_signal$type, levels=labels_x[c(5, 2, 1, 4, 3)])

# Prepare annotation columns for heatmap
ann_col <- combined_data[, c("Epithelial MIC", "Epithelial TED", "Epithelial POS", "Epithelial ACP", "Epithelial CGS")]
colnames(ann_col) <- gsub("Epithelial ", "", colnames(ann_col))
dpm_iir_path <- names(row_km[row_km %in% c(1, 3)])[-c(9, 12, 13, 19, 5, 14, 1, 2, 4, 7)]
x1 = c("Other glycan degradation", "Steroid metabolism", "Purine catabolism")

# Define colors for annotations
ann_colors <- list(CGS=colorRampPalette(c("white", "#458B00"))(100),
                   TED=colorRampPalette(c("white", "#53868B"))(100),
                   POS=colorRampPalette(c("white", "#CD8162"))(100),
                   ACP=colorRampPalette(c("white", "#EE30A7"))(100),
                   MIC=colorRampPalette(c("white", "#FF4A52"))(100))

# Create a heatmap for log2 fold changes
ph1 <- pheatmap(t(epi_select[, c(dpm_iir_path, x1)]),
                col=c(colorRampPalette(colors=c("#00BFFF", "white"))(length(bk)/2 + 20), 
                      colorRampPalette(colors=c("white", "#EE7621"))(length(bk)/2)),
                legend_breaks=seq(-5, 5, 0.5), annotation_col=ann_col, 
                clustering_method="ward.D2", annotation_colors=ann_colors, 
                main="Log2FC in Epithelial(T vs N)", fontsize_col=8)

# Create boxplots for pathways related to MIC and TED modules
p_cor1 <- ggplot(result_cor_signal[result_cor_signal$pathway %in% c(dpm_iir_path, x1), ], 
                 aes(x=type, y=abs(rho.rho), fill=type, col=I("gray"))) +
  geom_boxplot() +
  theme_hd_plain() +
  scale_fill_manual(name="Class", values=c("pink", "#AEEEEE", "#FFBBFF", "#9AFF9A", "#EECFA1")) +
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank()) +
  ylab("Abs(Spearman's rho of log2FC)") +
  ggtitle("Pathways related with\nMIC and TED signaling modules") +
  guides(fill="none")

# Create a heatmap for log2 fold changes for specific pathways
ph2 <- pheatmap(t(epi_select[, names(row_km[row_km %in% c(2, 5)])[-c(7, 6, 8, 18, 11, 14, 23)]]),
                col=c(colorRampPalette(colors=c("#00BFFF", "white"))(length(bk)/2 + 35), 
                      colorRampPalette(colors=c("white", "#EE7621"))(length(bk)/2)),
                legend_breaks=seq(-5, 5, 0.5), annotation_col=ann_col, 
                clustering_method="ward.D2", annotation_colors=ann_colors, 
                main="Log2FC in Epithelial(T vs N)", fontsize_col=8)

# Create boxplots for pathways related to CGS, ACP, and POS modules
p_cor2 <- ggplot(result_cor_signal[result_cor_signal$pathway %in% names(row_km[row_km %in% c(2, 5)])[-c(7, 6, 8, 18, 11, 14, 23)], ], 
                 aes(x=type, y=abs(rho.rho), fill=type, col=I("gray"))) +
  geom_boxplot() +
  theme_hd_plain() +
  scale_fill_manual(name="Class", values=c("pink", "#AEEEEE", "#FFBBFF", "#9AFF9A", "#EECFA1")) +
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.title.x = element_blank()) +
  ylab("Abs(Spearman's rho of log2FC)") +
  ggtitle("Pathways related with CGS,\nACP and POS signaling modules") +
  guides(fill="none")

# Combine data for plotting
combined_data_path <- cbind(combined_data, epi_select[rownames(combined_data), ])

# Create scatter plots for signaling modules against metabolic pathways
p1_signal_path1 <- ggplot(combined_data_path, aes(x=`Epithelial TED`, y=`Leukotriene metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\nleukotriene metabolism in epithelial") +
  xlab("TED signaling module in epithelial") +
  ggtitle("Log2FC (T VS N)")

p1_signal_path2 <- ggplot(combined_data_path, aes(x=`Pseudo TED`, y=`Leukotriene metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\nleukotriene metabolism in epithelial") +
  xlab("TED signaling module in pseudo-bulk data") +
  ggtitle("Log2FC (T VS N)")

# Select epithelial data again for further analysis
epi_select <- fc_pathwayRatioNormalized$FC[grepl("Epi", rownames(fc_pathwayRatioNormalized$FC)), rownames(signal_path$all_select$Epithelial)]
rownames(epi_select) <- str_split_fixed(rownames(epi_select), "_", 3)[, 3]

# Combine signaling data for correlation analysis
epi_signal <- combined_data[, c("Epithelial MIC", "Epithelial TED", "Epithelial POS", "Epithelial ACP", "Epithelial CGS")]
epi_signal_path <- cbind(epi_select, epi_signal[rownames(epi_select), ])

# Calculate Spearman correlation for epithelial signaling pathways
x_epi <- cor(epi_signal_path, method="spearman")
x_epi <- x_epi[!grepl("Epithelial", rownames(x_epi)), grepl("Epithelial", rownames(x_epi))]

# Combine pseudo-bulk signaling data for correlation analysis
pseudo_signal <- combined_data[, c("Pseudo MIC", "Pseudo TED", "Pseudo POS", "Pseudo ACP", "Pseudo CGS")]
pseudo_signal_path <- cbind(epi_select, pseudo_signal[rownames(epi_select), ])

# Calculate Spearman correlation for pseudo-bulk signaling pathways
x_epi_pseudo <- cor(pseudo_signal_path, method="spearman")
x_epi_pseudo <- x_epi_pseudo[!grepl("Pseudo", rownames(x_epi_pseudo)), grepl("Pseudo", rownames(x_epi_pseudo))]

# Combine correlation results for plotting
x_cor_epi_pseudo <- as.data.frame(cbind(x_epi, x_epi_pseudo[rownames(x_epi), ]))

# Initialize result for correlation between epithelial and pseudo-bulk signaling modules
result_cor_epi_pseudo <- c()

# Calculate Spearman correlation for each signaling module
for(i in labels_x) {
  temp_r <- cor.test(x_cor_epi_pseudo[, paste("Pseudo", i)], x_cor_epi_pseudo[, paste("Epithelial", i)], method="spearman")
  result_cor_epi_pseudo <- rbind(result_cor_epi_pseudo, c(type=i, rho=round(temp_r$estimate, 2), pvalue=temp_r$p.value))
}

# Convert results to a data frame
result_cor_epi_pseudo <- as.data.frame(result_cor_epi_pseudo)
result_cor_epi_pseudo$rho.rho <- as.numeric(result_cor_epi_pseudo$rho.rho)
result_cor_epi_pseudo$pvalue <- as.numeric(result_cor_epi_pseudo$pvalue)
result_cor_epi_pseudo$type <- factor(result_cor_epi_pseudo$type, levels=result_cor_epi_pseudo$type[order(result_cor_epi_pseudo$rho.rho)])

# Create a bar plot for epithelial vs pseudo-bulk signaling module correlations
p_cor_epi_pseudo <- ggplot(result_cor_epi_pseudo, aes(x=type, y=rho.rho)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_hd_plain() +
  ggtitle("Epithelial metabolic changes under signaling\nmodule regulation (Epithelial vs pseudo-bulk data)") +
  ylab("Spearman's rho") +
  geom_text(aes(label=rho.rho), vjust=0.1, size=3.5) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Signaling modules")

# Create scatter plots for TED, CGS, and MIC signaling modules
p1_pseudo_epi_ted <- ggplot(x_cor_epi_pseudo, aes(x=`Epithelial TED`, y=`Pseudo TED`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation(epithelial)\nVS TED module(pseudo-bulk data)") +
  xlab("Local metabolic allocation(epithelial)\nVS TED module(epithelial)") +
  ggtitle("Spearman's correlation of log2FC(T vs N)")

p1_pseudo_epi_cgs <- ggplot(x_cor_epi_pseudo, aes(x=`Epithelial CGS`, y=`Pseudo CGS`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation(epithelial)\nVS CGS module(pseudo-bulk data)") +
  xlab("Local metabolic allocation(epithelial)\nVS CGS module(epithelial)") +
  ggtitle("Spearman's correlation of log2FC(T vs N)")

p1_pseudo_epi_mic <- ggplot(x_cor_epi_pseudo, aes(x=`Epithelial MIC`, y=`Pseudo MIC`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation(epithelial)\nVS MIC module(pseudo-bulk data)") +
  xlab("Local metabolic allocation(epithelial)\nVS MIC module(epithelial)") +
  ggtitle("Spearman's correlation of log2FC(T vs N)")

# Create scatter plots for ACP signaling module
p1_signal_path3 <- ggplot(combined_data_path, aes(x=`Epithelial ACP`, y=`Pyruvate metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\npyruvate metabolism in epithelial") +
  xlab("ACP signaling module in epithelial") +
  ggtitle("Log2FC (T VS N)")
# Create a scatter plot for the relationship between ACP signaling module and pyruvate metabolism
p1_signal_path4 <- ggplot(combined_data_path, aes(x=`Pseudo ACP`, y=`Pyruvate metabolism`)) +
  geom_point() +
  geom_smooth(method="lm", formula=y ~ x, se=FALSE, col="black") +
  stat_cor(method="spearman", size=5) +
  theme_hd_plain() +
  ylab("Local metabolic allocation of\npyruvate metabolism in epithelial") +
  xlab("ACP signaling module in pseudo-bulk data") +
  ggtitle("Log2FC (T VS N)")

# Arrange multiple plots into a grid for figure 5
p5all0 <- ggarrange(p_network, p_ratio_heatmap, nrow=1, ncol=2, 
                    labels=letters[c(1:2)], font.label=list(size=18, color="black", face="bold"))

p5all1 <- ggarrange(p_cor_genes, p_cor, p1_dpm_cor, p_cor_sig_path, 
                    nrow=2, ncol=2, labels=letters[c(3:5, 7)], 
                    font.label=list(size=18, color="black", face="bold"))

p5all2 <- ggarrange(p5all1, as_ggplot(pp$gtable), nrow=1, ncol=2, 
                    labels=c("", letters[6]), font.label=list(size=18, color="black", face="bold"), 
                    widths=c(2, 2.5))

p5all4 <- ggarrange(as_ggplot(ph1$gtable), as_ggplot(ph2$gtable), nrow=1, ncol=2, 
                    labels=letters[c(8, 9)], font.label=list(size=18, color="black", face="bold"))

# Combine all parts into one final figure
p5all7 <- ggarrange(p5all0, p5all2, p5all4, nrow=3, ncol=1, heights=c(1.2, 2, 1.6))

# Save the combined figure as a TIFF file
tiff(paste("./figures/", "figure5new.tiff", sep=""), 
     height=24, width=24, res=300, units="in", compression="lzw")
print(p5all7)
dev.off()

# Save the combined figure as a PDF file
pdf(paste("./figures/", "figure5_new.pdf", sep=""), 
    height=24, width=23)
print(p5all7)
dev.off()

# Arrange supplementary plots into a grid for supplementary figure 8
p5all_sub <- ggarrange(p_cor_ted, p1_mic_cor, p1_acp_pos_cor, 
                       p_cor1, p_cor2, p_cor_epi_pseudo, 
                       p1_pseudo_epi_ted, p1_pseudo_epi_mic, 
                       p1_signal_path1, p1_signal_path2, 
                       p1_signal_path3, p1_signal_path4, 
                       nrow=4, ncol=3, labels=letters[c(1:12)], 
                       font.label=list(size=18, color="black", face="bold"), 
                       heights=c(1.7, 1, 1, 1))

# Save the supplementary figure as a TIFF file
tiff(paste("./figures/", "supfigure8.tiff", sep=""), 
     height=22, width=16, res=300, units="in", compression="lzw")
print(p5all_sub)
dev.off()









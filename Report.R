# Load the libraries
library(dplyr)
library(readr)
library(tidyr)
library(R.utils)
library(ggplot2)
library(GEOquery) 
library(data.table)
library(tidyverse)
library(limma)
library(igraph)
library(metafor)
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)



###########################################################################################
###########################################################################################
# GSE59034
# url : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59034
#Extract and annotate using the platform file GPL11532 and separate the samples of the GSE59034 dataset 
# then preform DEG and GO pathway analysis

# Load the dataset and its associated platform file
GSE59034 <- getGEO(filename = 'GSE59034_series_matrix .txt.gz')
GPL11532 <- fread('GPL11532-32230.txt', sep = '\t', header = TRUE)
print(head(GSE59034))
print(head(GPL11532))
# to check how to seperate the samples
GSE59034@phenoData@data[["characteristics_ch1.2"]]
# specifically checking the 'characteristics_ch1.2' column from the phenoData of GSE59034
characteristics_ch1_2 <- GSE59034@phenoData@data[["characteristics_ch1.2"]]
print(characteristics_ch1_2)


#To filter pre and post obesity files
# first remove samples labeled 'obesity status: never-obese'
filtered_samples <- GSE59034[, GSE59034@phenoData@data[["characteristics_ch1.2"]] != "obesity status: never-obese"]
pre_surgery <- filtered_samples[, filtered_samples@phenoData@data[["characteristics_ch1.2"]] == "obesity status: before bariatric surgery"]
post_surgery <- filtered_samples[, filtered_samples@phenoData@data[["characteristics_ch1.2"]] == "obesity status: after bariatric surgery"]

# Check the samples in each group to see if they match 
pre_surgery_samples <- colnames(pre_surgery)
post_surgery_samples <- colnames(post_surgery)
cat("Number of pre-surgery samples: ", ncol(pre_surgery), "\ ")
cat("Number of post-surgery samples: ", ncol(post_surgery), "\ ")

# Extract expression data for pre-surgery and post-surgery samples
exprs_data <- exprs(GSE59034)
pre_surgery_exprs <- exprs_data[, pre_surgery_samples]
post_surgery_exprs <- exprs_data[, post_surgery_samples]
combined_exprs <- cbind(pre_surgery_exprs, post_surgery_exprs)


############ DEG with limma########################
### Create a design matrix to fit for limma
#Create a Design Matrix: Create a design matrix that indicates which samples belong to the pre-surgery and post-surgery groups.
#Fit the Model: Use the limma package to fit a linear model to the combined expression data.
#Perform Differential Expression Analysis: Use the fitted model to perform differential expression analysis and identify differentially expressed genes between the two groups.

conditions <- factor(c(rep("pre_surgery", length(pre_surgery_samples)), rep("post_surgery", length(post_surgery_samples))))
design <- model.matrix(~0 + conditions)
colnames(design) <- levels(conditions)
fit <- lmFit(combined_exprs, design)
contrast_matrix <- makeContrasts(post_surgery - pre_surgery, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust = "fdr", number = Inf)
print(head(results))


# then get the top DEG and then seperate the gene names
GPL11532_subset <- GPL11532[, .(ID, gene_assignment)]
results <- topTable(fit2, adjust = "fdr", number = Inf)
results$ID <- rownames(results)
annotated_results <- merge(results, GPL11532_subset, by.x = 'ID', by.y = 'ID', all.x = TRUE)
print(head(annotated_results))

# some tidying up of gene_assignment column
annotated_results$gene_assignment <- sapply(strsplit(annotated_results$gene_assignment, "//"), `[`, 2)
print(head(annotated_results))
write.csv(annotated_results, 'GSE59034_results_ALL.csv', row.names = FALSE)


###########################################################################################
############# top differential expressed genes#################

# sorthing the genes by adjusted p-value then Filter by logFC and adjusted p-value thresholds
sorted_genes <- annotated_results[order(annotated_results$adj.P.Val), ]
filtered_genes <- subset(sorted_genes, abs(logFC) > 1 & adj.P.Val < 0.05)

top_diff_genes <- head(filtered_genes, 10) # just for comparison
print(top_diff_genes)


# cheking the genes of the absolute value of their log fold change
sorted_genes_by_fc <- annotated_results[order(abs(annotated_results$logFC), decreasing = TRUE), ]
# went with 10 only
top_genes_by_fc <- head(sorted_genes_by_fc, 10)
print(top_genes_by_fc)


#visualizations

#  the volcano plot
volcano_data <- annotated_results
volcano_data$neg_log_pval <- -log10(volcano_data$P.Value)

p <- ggplot(volcano_data, aes(x = logFC, y = neg_log_pval)) +
  geom_point(color = "black") +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  labs(title = "Volcano Plot of GSE59034", 
       x = "Log Fold Change (logFC)", 
       y = "-Log10(P-Value)") +
  theme_minimal() +
  annotate("text", x = max(volcano_data$logFC), y = -log10(0.05), 
           label = "p = 0.05", hjust = 1, vjust = -0.5, color = "red") +
  annotate("text", x = -1, y = max(volcano_data$neg_log_pval), 
           label = "logFC = -1", hjust = -0.1, vjust = 1, angle = 90, color = "red") +
  annotate("text", x = 1, y = max(volcano_data$neg_log_pval), 
           label = "logFC = 1", hjust = -0.1, vjust = 1, angle = 90, color = "red")
print(p)


# Sort the genes by log fold change
upregulated_genes <- annotated_results[order(-annotated_results$logFC), ][1:10, ]
downregulated_genes <- annotated_results[order(annotated_results$logFC), ][1:10, ]
top_genes <- rbind(upregulated_genes, downregulated_genes)
print(top_genes)
top_genes$Regulation <- ifelse(top_genes$logFC > 0, "Upregulated", "Downregulated") # add a column for them


# Barplot
p <- ggplot(top_genes, aes(x = reorder(gene_assignment, logFC), y = logFC, fill = Regulation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 Upregulated and Downregulated Genes for GSE59034", x = "Gene", y = "Log Fold Change (logFC)") +
  theme_minimal() +
  scale_fill_manual(values = c("Upregulated" = "brown2", "Downregulated" = "dodgerblue2"))
print(p)
write.csv(top_genes, 'GSE59034_top_up_down_genes.csv', row.names = FALSE)


##########################################
###########pathway analysis##############
########################################
#  GO pathway enrichment 

top_gene_symbols <- top_genes$gene_assignment
gene_list <- top_genes$gene_assignment
gene_list <- trimws(gene_list)
valid_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL") # Check for valid gene symbols
# Filter gene list to include only valid symbols
gene_list <- gene_list[gene_list %in% valid_genes]
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # gene symbols to Entrez IDs and check
print(entrez_ids)


go_enrich_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)#  for Biological Process (BP)
print(go_enrich_bp)
barplot(go_enrich_bp, showCategory = 10) + ggtitle("GO Biological Process for GSE59034")
go_enrich_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05)#  for Cellular Component (CC)
print(go_enrich_cc)
barplot(go_enrich_cc, showCategory = 10) + ggtitle("GO Cellular Component for GSE59034")
go_enrich_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05)# for Molecular Function (MF)
print(go_enrich_mf)
barplot(go_enrich_mf, showCategory = 10) + ggtitle("GO Pathway results for GSE59034") # Will go for MF because its fits our report aims the most

############################
### EXTRA not included in the report
############################
# Perform KEGG pathway enrichment analysis

top_gene_symbols <- top_genes$gene_assignment
results <- enrichr(top_gene_symbols, databases = c('KEGG_2021_Human', 'GO_Biological_Process_2021', 'Reactome_2021'))
# Re-run the pathway enrichment analysis for KEGG pathways
kegg_results <- enrichr(top_gene_symbols, databases = c('KEGG_2021_Human'))
print(head(kegg_results[['KEGG_2021_Human']], 5)) # Print only the top 10 results
#  a bar plot for KEGG pathways 
p_kegg <- ggplot(head(kegg_results[['KEGG_2021_Human']], 5), aes(x = reorder(Term, -Combined.Score), y = -log10(Adjusted.P.value), fill = Genes)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(title = 'Top Enriched KEGG Pathways', x = 'Pathway', y = '-log10(Adjusted P-value)') +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  ylim(0, max(-log10(head(kegg_results[['KEGG_2021_Human']], 10)$Adjusted.P.value)) + 1) +
  scale_fill_manual(values = c("deeppink", "brown2", "blue", "purple","skyblue" ))
print(p_kegg)

# GO Biological Process pathways
go_results <- enrichr(top_gene_symbols, databases = c('GO_Biological_Process_2021')) # to check all GO pathways
print(head(go_results[['GO_Biological_Process_2021']]))
go_results <- go_results[['GO_Biological_Process_2021']]# Filter 
top_go_results <- go_results[order(go_results$Adjusted.P.value), ][1:10, ]# Select the top pathways 
#  a bar plot for GO Biological Process pathways
p_go <- ggplot(top_go_results, aes(x = reorder(Term, -Combined.Score), y = -log10(Adjusted.P.value), fill = Genes)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(title = 'Top Enriched GO Biological Process Pathways', x = 'Pathway', y = '-log10(Adjusted P-value)') +
  theme_minimal() +
  theme(legend.position = 'bottom')
print(p_go)



############################################################################################
###########################################################################################
###########################################################################################
# GSE66921
# url : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66921
#Extract and annotate using the platform file GPL13607 and seperate the samples of the GSE66921 dataset 
# then preform DEG and GO pathway analysis


# Load the dataset and platform file
GSE66921 <- getGEO(filename = 'GSE66921_series_matrix.txt.gz')
GPL13607 <- read.table('GPL13607-20416.txt', header = TRUE, sep = '\t', quote = '', comment.char = '')
print(head(GSE66921))
print(head(GPL13607))
expr_data <- exprs(GSE66921)
print(head(expr_data))
# will use the metadata to help better separate the samples according to surgery status
pheno_data <- pData(GSE66921)
print(colnames(pheno_data))
selected_columns <- pheno_data[, c("geo_accession", "source_name_ch1")] # to check if before or after
print(selected_columns)

pre_surgery_samples <- c('GSM1634360', 'GSM1634361', 'GSM1634364', 'GSM1634365')
post_surgery_samples <- c('GSM1634358', 'GSM1634359', 'GSM1634362', 'GSM1634363')
combined_expr_data <- expr_data[, c(pre_surgery_samples, post_surgery_samples)]
print(head(combined_expr_data))


### Create a design matrix to fit for limma

#Create a Design Matrix: Create a design matrix that indicates which samples belong to the pre-surgery and post-surgery groups.
#Fit the Model: Use the limma package to fit a linear model to the combined expression data.
#Perform Differential Expression Analysis: Use the fitted model to perform differential expression analysis and identify differentially expressed genes between the two groups.

group <- factor(c(rep('pre_surgery', length(pre_surgery_samples)), rep('post_surgery', length(post_surgery_samples))))
design <- model.matrix(~ group)
print(design)
fit <- lmFit(combined_expr_data, design) # Fit the linear model
print(fit)
fit <- eBayes(fit) # perform DEG
results <- topTable(fit, coef = 2, adjust = "fdr", number = Inf)
print(head(results))


GPL13607 <- read.table('GPL13607-20416.txt', header = TRUE, sep = '\t', quote = '', comment.char = '', skip = 16) #to map the gene IDS 

platform_data <- GPL13607[, c('X6', 'ENST00000322831')]#the relevant columns from the platform file (Probe ID and Gene Name)
colnames(platform_data) <- c('ID', 'gene_symbol')
results_with_genes <- merge(results, platform_data, by.x = 'row.names', by.y = 'ID', all.x = TRUE)
results_with_genes <- results_with_genes[, -1]
print(head(results_with_genes))
write.csv(results_with_genes, 'GSE66921_DEG_ALL.csv', row.names = FALSE)

# Extract the top DEGs based on log fold change (>1) and adjusted p-values (<0.05)
sorted_genes <- results_with_genes[order(results_with_genes$adj.P.Val), ]
filtered_genes <- subset(sorted_genes, abs(logFC) > 1 & adj.P.Val < 0.05)
top_diff_genes <- head(filtered_genes, 10)  # just the top 10
print(top_diff_genes)


#  the top 10 upregulated and top 10 downregulated genes
top_upregulated_genes <- head(filtered_genes[order(filtered_genes$logFC, decreasing = TRUE), ], 10)
top_downregulated_genes <- head(filtered_genes[order(filtered_genes$logFC), ], 10)
top_genes <- rbind(top_upregulated_genes, top_downregulated_genes)
write.csv(top_genes, 'GSE66921_top_10_genes.csv', row.names = FALSE) # Save to a CSV file

top_genes <- top_genes[top_genes$gene_symbol != 'lincRNA:chr7:25978475-25989975_R', ]
print(top_genes)

# Create a bar plot for the top 10 upregulated and top 10 downregulated genes
bar_plot <- ggplot(top_genes, aes(x = reorder(gene_symbol, logFC), y = logFC, fill = logFC > 0)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(title = "Top Upregulated and Downregulated Genes in GSE66921", x = "Gene Symbol", y = "Log Fold Change") +
  coord_flip() +
  scale_fill_manual(values = c("dodgerblue2", "brown2"), labels = c("Downregulated", "Upregulated")) +
  theme(legend.title = element_blank())
print(bar_plot)

#  a volcano plot
volcano_plot <- ggplot(results_with_genes, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_minimal() +
  labs(title = "Volcano Plot of GSE66921", x = "Log Fold Change", y = "-Log10 P-Value") +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed")
print(volcano_plot)


#  an MA plot
ma_plot <- ggplot(results, aes(x = AveExpr, y = logFC)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_minimal() +
  labs(title = "MA Plot", x = "Average Expression", y = "Log Fold Change") +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed")
print(ma_plot)

##################################
#PATHWAY ANALYSIS
##################################


#  GO pathway enrichment 
#visualize the results with dot plots and bar plots.
gene_list <- top_diff_genes$gene_symbol
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) #Convert gene symbols to Entrez IDs using `bitr`.
entrez_df <- data.frame(gene_symbol = names(entrez_ids), entrez_id = entrez_ids, stringsAsFactors = FALSE)


#  GO pathway enrichment 
go_enrich_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)#  for Biological Process (BP)
print(go_enrich_bp)
dotplot(go_enrich_bp, showCategory = 10) + ggtitle("GO Biological Process")
go_enrich_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05)#  for Cellular Component (CC)
print(go_enrich_cc)
dotplot(go_enrich_cc, showCategory = 10) + ggtitle("GO Cellular Component")
go_enrich_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05)# for Molecular Function (MF)
print(go_enrich_mf)
barplot(go_enrich_mf, showCategory = 10) + ggtitle("GO MF results for GSE66921") # MF for the report


#  KEGG pathway enrichment analysis
# had to adjust the pvalues
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.5) # adjusted to 0.5 to get more results 
print(kegg_enrich)
barplot(kegg_enrich, showCategory = 10)
write.csv(kegg_enrich, 'GSE66921_kegg_enrich.csv', row.names = FALSE)



####### EXTRA not included in the report #############

####### network diagram #######
# to define the interactions between the genes
interactions <- data.frame(
  from = c('FOSB', 'FOSB', 'IL1RN', 'IL1RN', 'CCL3', 'CCL3', 'ITGAX', 'ITGAX', 'C5AR1', 'C5AR1', 'EGR2', 'EGR2'),
  to = c('IL1RN', 'CCL3', 'CCL3', 'C5AR1', 'ITGAX', 'C5AR1', 'C5AR1', 'EGR2', 'EGR2', 'FOSB', 'FOSB', 'IL1RN'),
  weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)
g <- graph_from_data_frame(interactions, directed = TRUE)#  a graph object
plot(g, vertex.label.color = 'black', vertex.size = 30, edge.arrow.size = 0.5, main = 'Gene Interaction Network')



###########################################################################################
###########################################################################################
###########################################################################################
#GSE199063
# url:https://www.ncbi.nlm.nih.gov/gds/?term=GSE199063[Accession]


GSE199063 <- getGEO(filename = 'GSE199063_series_matrix.txt.gz')
GPL23126 <- fread('GPL23126-131.txt', sep = '\t', header = TRUE)
str(GPL23126)
str(GSE199063)
expr_data <- exprs(GSE199063) #  the expression data 
str(expr_data)
expr_df <- as.data.frame(expr_data)
expr_df$ID <- rownames(expr_df)
merged_data <- merge(expr_df, GPL23126, by = 'ID')
str(merged_data)

#### Separate the pre and post surgery ###
pheno_data <- pData(GSE199063)
pre_surgery_samples <- pheno_data[pheno_data$`characteristics_ch1.2` == 'time_point: Baseline', ]
post_surgery_samples <- pheno_data[pheno_data$`characteristics_ch1.2` == 'time_point: 2_year', ]
pre_surgery_samples$geo_accession
post_surgery_samples$geo_accession

common_subjects <- intersect(pre_surgery_samples$`characteristics_ch1.3`, post_surgery_samples$`characteristics_ch1.3`)## matching the samples 
pre_surgery_samples_matched <- pre_surgery_samples[pre_surgery_samples$`characteristics_ch1.3` %in% common_subjects, ]
post_surgery_samples_matched <- post_surgery_samples[post_surgery_samples$`characteristics_ch1.3` %in% common_subjects, ]

# whats the samples in each group
head(pre_surgery_samples_matched)
head(post_surgery_samples_matched)
cat('Matched pre-surgery samples:', nrow(pre_surgery_samples_matched), '\ ')
cat('Matched post-surgery samples:', nrow(post_surgery_samples_matched), '\ ')
# sanity check
pre_surgery_samples_matched$geo_accession
post_surgery_samples_matched$geo_accession



###############################################
### Differential gene expression Analysis ####


# macking sure we only use the pre and post-op samples 
expr_data_matched <- expr_data[, c(pre_surgery_samples_matched$geo_accession, post_surgery_samples_matched$geo_accession)] #  subseting the expression data to include only the matched samples
group <- factor(c(rep("Pre", nrow(pre_surgery_samples_matched)), rep("Post", nrow(post_surgery_samples_matched)))) #  indicating pre-surgery and post-surgery samples

design <- model.matrix(~0 + group)
colnames(design) <- c("Pre", "Post")
fit <- lmFit(expr_data_matched, design)
contrast_matrix <- makeContrasts(Post - Pre, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
top_genes <- topTable(fit2, adjust = "fdr", number = Inf)
head(top_genes)

# Extract the gene IDs to map them
gene_ids <- rownames(top_genes)
subset_gpl_data <- GPL23126[, c("ID", "gene_assignment")]
top_genes_df <- data.frame(ID = gene_ids, top_genes)
mapped_genes <- merge(top_genes_df, subset_gpl_data, by = "ID")
print(head(mapped_genes))

#  filter by logFC and adjusted p-value thresholds
sorted_genes <- mapped_genes[order(mapped_genes$adj.P.Val), ]
filtered_genes <- subset(sorted_genes, abs(logFC) > 1 & adj.P.Val < 0.05)
top_diff_genes <- head(filtered_genes, 10)  
print(top_diff_genes)


### Some tidying up to do #####

top_diff_genes$gene_symbol <- sapply(strsplit(top_diff_genes$gene_assignment, " // "), function(x) x[2]) # extract gene symbols from the gene_assignment column
top_diff_genes <- top_diff_genes[!is.na(top_diff_genes$gene_symbol), ]
mapped_genes$gene_symbol <- sapply(strsplit(mapped_genes$gene_assignment, " // "), function(x) x[2])
mapped_genes <- mapped_genes[!is.na(mapped_genes$gene_symbol), ]
top_diff_genes <- top_diff_genes[, -ncol(top_diff_genes) + 1]
mapped_genes <- mapped_genes[, -ncol(mapped_genes) + 1]
print(head(top_diff_genes))
print(head(mapped_genes))
write.csv(mapped_genes, file = "GSE199063_DEG_ALL.csv", row.names = FALSE)

#  a volcano plot
volcano_data <- data.frame(logFC = top_genes$logFC, 
                           negLogP = -log10(top_genes$P.Value))

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = negLogP)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Volcano Plot for GSE199063 ", x = "Log Fold Change", y = "-Log10(P-Value)") +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed")
print(volcano_plot)

#  an MA plot
ma_data <- data.frame(AveExpr = top_genes$AveExpr, 
                      logFC = top_genes$logFC)

ggplot(ma_data, aes(x = AveExpr, y = logFC)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "MA Plot", x = "Average Expression", y = "Log Fold Change") +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed")



# tidying up the gene column
filtered_genes <- filtered_genes[, !names(filtered_genes) %in% c("ID")]
filtered_genes$gene_symbol <- sapply(strsplit(filtered_genes$gene_assignment, " // "), function(x) x[2])
filtered_genes <- filtered_genes[!is.na(filtered_genes$gene_symbol), ]



#for upregulated and downregulated genes
upregulated_genes <- filtered_genes[filtered_genes$logFC > 0, ]
downregulated_genes <- filtered_genes[filtered_genes$logFC < 0, ]
top_10_upregulated <- head(upregulated_genes[order(-upregulated_genes$logFC), ], min(10, nrow(upregulated_genes)))
top_10_downregulated <- head(downregulated_genes[order(downregulated_genes$logFC), ], min(10, nrow(downregulated_genes)))
top_10_combined <- rbind(top_10_upregulated, top_10_downregulated)
top_genes <- top_10_combined

# Create a bar plot for the top 10 upregulated and top 10 downregulated genes
bar_plot <- ggplot(top_genes, aes(x = reorder(gene_symbol, logFC), y = logFC, fill = logFC > 0)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(title = "Top Upregulated and Downregulated Genes in GSE66921", x = "Gene Symbol", y = "Log Fold Change") +
  coord_flip() +
  scale_fill_manual(values = c("dodgerblue2", "brown2"), labels = c("Downregulated", "Upregulated")) +
  theme(legend.title = element_blank())
print(bar_plot)

#################################################
####### Pathway analysis#########################
#################################################

########Pathway Analysis########

# Convert gene symbols to Entrez IDs
top_diff_genes_entrez <- bitr(top_diff_genes$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
mapped_genes_entrez <- bitr(mapped_genes$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#  KEGG pathway enrichment analysis
# adjust the top 10 to pvalue!!!
kegg_enrich_top_diff <- enrichKEGG(gene = top_diff_genes_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
kegg_enrich_mapped <- enrichKEGG(gene = mapped_genes_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
print(head(kegg_enrich_top_diff))
print(head(kegg_enrich_mapped))

# Visualize  KEGG
dotplot(kegg_enrich_mapped, showCategory=20) + ggtitle("KEGG Pathway Enrichment - Mapped Genes")
barplot(kegg_enrich_mapped, showCategory=10, title="KEGG Pathway Enrichment - Mapped Genes")

write.csv(kegg_enrich_top_diff, file = "GSE199063_pathwayTOP10.csv", row.names = FALSE)
write.csv(kegg_enrich_mapped, file = "GSE199063_pathwayALL.csv", row.names = FALSE)

###### GO Analysis########

go_enrich_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)#  for Biological Process (BP)
print(go_enrich_bp)
barplot(go_enrich_bp, showCategory = 10) + ggtitle("GO Biological Process for GSE199063")
go_enrich_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05)#  for Cellular Component (CC)
print(go_enrich_cc)
barplot(go_enrich_cc, showCategory = 10) + ggtitle("GO Cellular Component for GSE199063")
go_enrich_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05)# for Molecular Function (MF)
print(go_enrich_mf)
barplot(go_enrich_mf, showCategory = 10) + ggtitle("GO pathway Molecular Function for GSE199063")


# Another GO enrichment analysis to check if its been consistant
gene_symbols <- top_10_combined$gene_symbol # only for the top genes!!!
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Another GO enrichment analysis to check if its been consistant
go_enrichment <- enrichGO(gene = gene_entrez$ENTREZID, 
                          OrgDb = org.Hs.eg.db, 
                          keyType = "ENTREZID", 
                          ont = "ALL", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05)
barplot(go_enrichment, showCategory = 10) + ggtitle("GO pathway enrichment results for GSE199063")

## Basically all the the GO are involved in inflammation, immunity and lipid metabolism

###########################################################################################
###########################################################################################
############################### Meta-analysis Step ########################################
###########################################################################################
###########################################################################################
# for the meta analysis I will simply load the previous csv files of all the DEG of each data sets that I have saved.
# Since the results of the DEG of each dataset were continuously being filtered through the multiple analysis conducted
# aalso excel can be used for the z-transformation method that will be used for the analysis to help check the results of the code




# Load the saved datasets for the meta-analysis
GSE66921 <- read.csv('GSE66921_DEG_ALL.csv')
GSE59034 <- read.csv('GSE59034_results_ALL.csv')
GSE199063 <- read.csv('GSE199063_DEG_ALL.csv')

# Remove the ID and info columns from GSE59034 and the ID column from GSE199063
GSE199063 <- GSE199063[, -1]
GSE59034 <- GSE59034[, -c(1, 9)]
combined_dataset <- bind_rows(GSE66921, GSE59034, GSE199063)

# Z-transformation

# calculate the mean fold-change for each gene
mean_fold_change <- combined_dataset %>% 
  group_by(gene_symbol) %>% 
  summarise(mean_logFC = mean(logFC, na.rm = TRUE))
# calculate the combined, adjusted P-value through the Z-transform method
combined_p_values <- combined_dataset %>%  
  group_by(gene_symbol) %>% 
  summarise(combined_p_value = pnorm(sum(qnorm(adj.P.Val))/sqrt(n())))
meta_analysis_results <- mean_fold_change %>% 
  inner_join(combined_p_values, by = "gene_symbol")
meta_analysis_results <- meta_analysis_results %>%  # Adjust the combined P-value
  mutate(adj_combined_p_value = p.adjust(combined_p_value, method = "BH"))

# sort by adjusted combined P-value and mean fold-change
meta_analysis_results <- meta_analysis_results %>% 
  arrange(adj_combined_p_value, desc(abs(mean_logFC)))


# the top DEG
top_genes <- meta_analysis_results %>% 
  head(10)


gene_symbols <- meta_analysis_results$gene_symbol
annotations <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = c("SYMBOL", "GENENAME"), keytype = "SYMBOL")
annotated_results <- meta_analysis_results %>% 
  dplyr::left_join(annotations, by = c("gene_symbol" = "SYMBOL")) %>% 
  dplyr::select(gene_symbol, GENENAME, mean_logFC, adj_combined_p_value) %>% 
  arrange(adj_combined_p_value, desc(abs(mean_logFC))) %>% 
  head(10)
print(annotated_results)



# Create a bar plot for the top differentially expressed genes
plot <- ggplot(annotated_results, aes(x = reorder(gene_symbol, -mean_logFC), y = mean_logFC, fill = adj_combined_p_value)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "dodgerblue2", high = "red") +
  labs(title = "Top Differentially Expressed Genes",
       x = "Gene Symbol",
       y = "Mean Log Fold-Change",
       fill = "Adjusted Combined P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot)



# Highlight significant genes
significant_genes <- meta_analysis_results %>%
  filter(adj_combined_p_value < 0.05 & abs(mean_logFC) > 1)

# Create an improved volcano plot
volcano_plot <- ggplot(meta_analysis_results, aes(x = mean_logFC, y = -log10(adj_combined_p_value))) +
  geom_point(aes(color = adj_combined_p_value), alpha = 0.6) +
  geom_point(data = significant_genes, aes(x = mean_logFC, y = -log10(adj_combined_p_value)), color = "red", size = 2) +
  scale_color_gradient(low = "dodgerblue2", high = "brown2") +
  labs(title = "Volcano Plot of Meta-Analysis Results",
       x = "Mean Log Fold-Change",
       y = "-log10(Adjusted Combined P-Value)",
       color = "Adjusted Combined P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(volcano_plot)


# Select the top differentially expressed genes (top 10 upregulated and top 10 downregulated)
top_upregulated_genes <- meta_analysis_results %>%
  filter(mean_logFC > 0) %>%
  head(10)

top_downregulated_genes <- meta_analysis_results %>%
  filter(mean_logFC < 0) %>%
  head(10)

top_genes <- bind_rows(top_upregulated_genes, top_downregulated_genes)

# Create barplot for the top upregulated and downregulated genes
ggplot(top_genes, aes(x = reorder(gene_symbol, mean_logFC), y = mean_logFC, fill = mean_logFC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Upregulated and Downregulated Genes overall" ,
       x = "Gene Symbol",
       y = "Mean Log Fold Change") +
  scale_fill_manual(values = c("dodgerblue2", "brown2"), labels = c("Downregulated", "Upregulated")) +
  theme_minimal()



##############################
### Pathway analysis ####
#############################


########################################
####### GO pathway enrichment  #############

# Convert gene symbols to Entrez IDs
gene_symbols <- c("HYOU1", "MMP9", "FGR", "HN1", "ITGB2", "TFRC", "CCL2", "HMOX1", "TNFRSF11B", "SYK")
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
print(entrez_ids)
go_enrich_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)#  for Biological Process (BP)
print(go_enrich_bp)
barplot(go_enrich_bp, showCategory = 10) + ggtitle("GO Pathway results for Meta-Analysis ")
go_enrich_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05)#  for Cellular Component (CC)
print(go_enrich_cc)
barplot(go_enrich_cc, showCategory = 10) + ggtitle("GO Cellular Component for Meta-Analysis ")
go_enrich_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05)# for Molecular Function (MF)
print(go_enrich_mf)
barplot(go_enrich_mf, showCategory = 10) + ggtitle("GO Molecular Function for Meta-Analysis ")


#  KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(gene = gene_entrez$ENTREZID, 
                           organism = "hsa", 
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05)
head(kegg_results)# View the top KEGG pathways by p-value of 0.05!


# Get the pathway ID for the most significant pathway ro print pathview
most_significant_pathway <- kegg_results@result[1, ]
pathview_id <- most_significant_pathway$ID
print(pathview_id)
pathview(gene.data = entrez_ids$ENTREZID, pathway.id = pathview_id, species = "hsa", kegg.native = TRUE) 
# nitrogen metabolism is the top pathway involved


# KEGG pathway enrichment analysis with a p-value cutoff of 0.5!!!
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa', pvalueCutoff = 0.5)
print(kegg_enrich)
#compare the 0.05 and 0.5 pvalue cutoff
barplot(kegg_enrich, showCategory=10, title="KEGG Pathway Enrichment - Meta-analysis - pvalueCutoff = 0.5")
barplot(kegg_results, showCategory=10, title="KEGG Pathway Enrichment - Meta-analysis - pvalueCutoff = 0.05")

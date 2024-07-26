# in this R script , The 3 datasets where loaded directly from GEO
# The pre- and post op samples where then separated and filtered
# then the ComBat method was used to handle batch effect 
# PCA was used to see how well batch effect was handled between filtered and normalized dataset
# PCA did not show much improvement in handling batch effect so it was not included in this report
# the normalized data can be used in the LIMMAnorm.R Rscript to show the Results of the DEG 
# but not much overlap can be seen between the individual DEG results of and the batch corrected dataset from this R script


# Load the libraries
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(sva)
library(limma)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(visNetwork)
library(RankProd)
library(pathview)

########## start by loading and preprocessing each dataset ############

#######################
###### GSE66921###### 
#################


# Load series and platform data from GEO directly
gset <- getGEO("GSE66921", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL13607", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex_GSE66921 <- exprs(gset) # Extract expression data
print(head(ex_GSE66921))
platform_GSE66921 <- getGEO("GPL13607", AnnotGPL = TRUE)

platform_table <- Table(platform_GSE66921)
print(colnames(platform_table))
#use the ID column for probe IDs and the GeneName column for gene names. 

probe_to_gene <- platform_table[, c("ID", "GeneName")]
mapped_ex <- merge(data.frame(ID = rownames(ex_GSE66921), ex_GSE66921), probe_to_gene, by = "ID")
print(head(mapped_ex))


#################################################
########## Seperating the samples ##############

metadata <- pData(gset)
print(colnames(metadata))
# separate the metadata into pre- and post-surgery samples based on source_name_ch1 and source_name_ch2
pre_surgery_metadata <- metadata[metadata$source_name_ch1 == 'Subcutaneous adipose tissue before bariatric surgery' | metadata$source_name_ch2 == 'Subcutaneous adipose tissue before bariatric surgery', ]
post_surgery_metadata <- metadata[metadata$source_name_ch1 == 'Subcutaneous adipose tissue after bariatric surgery' | metadata$source_name_ch2 == 'Subcutaneous adipose tissue after bariatric surgery', ]

pre_surgery_samples <- rownames(pre_surgery_metadata)
post_surgery_samples <- rownames(post_surgery_metadata)

preop_GSE66921 <- mapped_ex[, c("ID", "GeneName", pre_surgery_samples)]
postop_GSE66921 <- mapped_ex[, c("ID", "GeneName", post_surgery_samples)]


write.csv(preop_GSE66921, file = 'preop_GSE66921.csv', row.names = FALSE)
write.csv(postop_GSE66921, file = 'postop_GSE66921.csv', row.names = FALSE)

######################
###### GSE59034###### 
#################

# Load series and platform data from GEO
gsetGSE59034 <- getGEO("GSE59034", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gsetGSE59034) > 1) idx <- grep("GPL11532", attr(gsetGSE59034, "names")) else idx <- 1
gsetGSE59034 <- gsetGSE59034[[idx]]
exGSE59034 <- exprs(gsetGSE59034)
platform_GPL11532 <- getGEO("GPL11532", AnnotGPL = TRUE)
platform_table <- Table(platform_GPL11532)

probes_GSE59034 <- platform_table[, c("ID", "Gene symbol")]
mapped_GSE59034 <- merge(data.frame(ID = rownames(exGSE59034)), probes_GSE59034, by = "ID")

mapped_GSE59034 <- merge(data.frame(ID = rownames(exGSE59034), exGSE59034), probes_GSE59034, by = "ID")
print(head(mapped_GSE59034))

# rename the 'Gene symbol' column to 'Gene name' and move it to the 2nd column
colnames(mapped_GSE59034)[colnames(mapped_GSE59034) == 'Gene symbol'] <- 'GeneName'
mapped_GSE59034 <- mapped_GSE59034[, c('ID', 'GeneName', setdiff(colnames(mapped_GSE59034), c('ID', 'Gene name')))]
missing_gene_names <- sum(is.na(mapped_GSE59034$GeneName) | mapped_GSE59034$GeneName == "")
print(missing_gene_names) #### 11127 missing gene IDs!!
gse_metadata <- pData(gsetGSE59034)

gse_metadata <- pData(gsetGSE59034)
gse_metadata_obesity_columns <- colnames(gse_metadata)[grepl("obesity", colnames(gse_metadata), ignore.case = TRUE)]
print(gse_metadata_obesity_columns)
gse_metadata_obesity <- gse_metadata[, gse_metadata_obesity_columns]
print(head(gse_metadata_obesity))

###########################################################
# separate the metadata into pre- and post-surgery samples
gse_metadata <- pData(gsetGSE59034)
pre_surgery_metadata <- gse_metadata[gse_metadata$`obesity status:ch1` == "before bariatric surgery", ]
post_surgery_metadata <- gse_metadata[gse_metadata$`obesity status:ch1` == "after bariatric surgery", ]

pre_surgery_sample_ids <- rownames(pre_surgery_metadata)
post_surgery_sample_ids <- rownames(post_surgery_metadata)

pre_surgery_expression_data <- exGSE59034[, pre_surgery_sample_ids]
post_surgery_expression_data <- exGSE59034[, post_surgery_sample_ids]

# by matching the ID in mapped_GSE59034 to the row names in each pre_surgery_expression_data and post_surgery_expression_data

preop_GSE59034 <- merge(mapped_GSE59034[, c('ID', 'GeneName')], pre_surgery_expression_data, by.x = 'ID', by.y = 'row.names')
postop_GSE59034 <- merge(mapped_GSE59034[, c('ID', 'GeneName')], post_surgery_expression_data, by.x = 'ID', by.y = 'row.names')

write.csv(preop_GSE59034, file = "preop_GSE59034.csv", row.names = FALSE)
write.csv(postop_GSE59034, file = "postop_GSE59034.csv", row.names = FALSE)



#########################
###### GSE199063###### 
####################

# Load series and platform data from GEO
gsetGSE199063 <- getGEO("GSE199063", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gsetGSE199063) > 1) idx <- grep("GPL23126", attr(gsetGSE199063, "names")) else idx <- 1
gsetGSE199063 <- gsetGSE199063[[idx]]
exGSE199063 <- exprs(gsetGSE199063)# Extract expression data
print(head(exGSE199063))
platform_GPL23126 <- getGEO("GPL23126", AnnotGPL = TRUE)
platform_table <- Table(platform_GPL23126)# Extract the data table 
print(colnames(platform_table))
print(head(platform_table))

# separate the 'gene_assignment' column based on '//' delimiter
platform_table$gene_assignment_split <- strsplit(as.character(platform_table$gene_assignment), ' // ')

#  a function to extract gene symbol and gene name
extract_gene_info <- function(x) {
  if (length(x) >= 3) {
    gene_symbol <- x[2]
    gene_name <- x[3]
    return(c(gene_symbol, gene_name))
  } else {
    return(c(NA, NA))
  }
}

gene_info <- t(sapply(platform_table$gene_assignment_split, extract_gene_info))
platform_table$GeneName <- gene_info[,1] # add back to platform file 
platform_table$gene_info <- gene_info[,2]
# Print the first few rows of the updated platform_table
print(head(platform_table[, c('ID', 'GeneName', 'gene_info')]))

#########################################
# seperate the pre and post surgery sampless

pheno_data <- pData(gsetGSE199063)

pre_surgery_samples <- pheno_data[pheno_data$`characteristics_ch1.2` == 'time_point: Baseline', ]
post_surgery_samples <- pheno_data[pheno_data$`characteristics_ch1.2` == 'time_point: 2_year', ]

pre_surgery_samples$geo_accession
post_surgery_samples$geo_accession

### matching the samples !!!!

common_subjects <- intersect(pre_surgery_samples$`characteristics_ch1.3`, post_surgery_samples$`characteristics_ch1.3`)
pre_surgery_samples_GSE199063 <- pre_surgery_samples[pre_surgery_samples$`characteristics_ch1.3` %in% common_subjects, ]
post_surgery_samples_GSE199063 <- post_surgery_samples[post_surgery_samples$`characteristics_ch1.3` %in% common_subjects, ]

pre_surgery_samples_GSE199063$geo_accession
post_surgery_samples_GSE199063$geo_accession
pre_surgery_expression_data_matched <- exGSE199063[, pre_surgery_samples_GSE199063$geo_accession]
post_surgery_expression_data_matched <- exGSE199063[, post_surgery_samples_GSE199063$geo_accession]

print(head(pre_surgery_expression_data_matched))
print(head(post_surgery_expression_data_matched))

id_gene_mapping <- platform_table[, c('ID', 'GeneName')]
print(head(id_gene_mapping))

pre_surgery_gene_names <- id_gene_mapping$GeneName[match(rownames(pre_surgery_expression_data_matched), id_gene_mapping$ID)]
print(sum(is.na(pre_surgery_gene_names)))
pre_surgery_with_genes <- cbind(
  ID = rownames(pre_surgery_expression_data_matched),
  GeneName = pre_surgery_gene_names,
  pre_surgery_expression_data_matched
)

pre_surgery_with_genes <- as.data.frame(pre_surgery_with_genes) #  the matrix to a data frame
print(head(pre_surgery_with_genes[, 1:5]))

# repeat the process for post_surgery
post_surgery_gene_names <- id_gene_mapping$GeneName[match(rownames(post_surgery_expression_data_matched), id_gene_mapping$ID)]
print(sum(is.na(post_surgery_gene_names)))
post_surgery_with_genes <- cbind(
  ID = rownames(post_surgery_expression_data_matched),
  GeneName = post_surgery_gene_names,
  post_surgery_expression_data_matched
)
post_surgery_with_genes <- as.data.frame(post_surgery_with_genes)
print(head(post_surgery_with_genes[, 1:5]))

preop_GSE199063 <- pre_surgery_with_genes
postop_GSE199063 <- post_surgery_with_genes


write.csv(preop_GSE199063, 'preop_GSE199063.csv', row.names = FALSE)
write.csv(postop_GSE199063, 'postop_GSE199063.csv', row.names = FALSE)

####################################################################################################################################
####################################################################################################################################


##########
#loading and normalising the datasets
########################


##########
#postop
########################

#  to load and inspect a dataset
load_and_inspect <- function(file_name) {
  data <- read.csv(file_name)
  cat("Dataset:", file_name, "\
")
  cat("Dimensions:", dim(data), "\
")
  cat("Column names:", colnames(data)[1:min(5, ncol(data))], "...\
")
  cat("First few rows:\
")
  print(head(data[, 1:min(5, ncol(data))]))
  cat("\
")
  return(data)
}

postop_GSE59034 <- load_and_inspect("postop_GSE59034.csv")
postop_GSE66921 <- load_and_inspect("postop_GSE66921.csv")
postop_GSE199063 <- load_and_inspect("postop_GSE199063.csv")


#  to load and process each dataset
load_dataset <- function(file_name) {
  data <- read.csv(file_name)
  # Convert ID to character, handling NA values
  data$ID <- as.character(data$ID)
  # Rename columns to include dataset identifier
  names(data)[-(1:2)] <- paste0(file_name, "_", names(data)[-(1:2)])
  return(data)
}

datasets <- list(
  postop_GSE66921 = load_dataset("postop_GSE66921.csv"),
  postop_GSE199063 = load_dataset("postop_GSE199063.csv"),
  postop_GSE59034 = load_dataset("postop_GSE59034.csv")
)

# Display structure of each dataset
for (name in names(datasets)) {
  cat("Structure of", name, ":\
")
  print(str(datasets[[name]]))
  cat("\
")
}

merged_data <- datasets[[1]]
for (i in 2:length(datasets)) {
  merged_data <- full_join(merged_data, datasets[[i]], by = c("ID", "GeneName"))
}

for (name in names(datasets)) {
  col_count <- sum(grepl(name, names(merged_data)))
  cat(name, ":", col_count, "\
")
}


# Save the merged dataset
write.csv(merged_data, "merged_postop_data.csv", row.names = FALSE)

remove_duplicates <- function(df) {
  df <- df[!duplicated(df$GeneName), ]
  return(df)
}

postop_GSE66921_filtered <- remove_duplicates(postop_GSE66921) # Remove duplicates
postop_GSE199063_filtered <- remove_duplicates(postop_GSE199063)
postop_GSE59034_filtered <- remove_duplicates(postop_GSE59034)

common_genes <- Reduce(intersect, list(postop_GSE66921_filtered$GeneName, postop_GSE199063_filtered$GeneName, postop_GSE59034_filtered$GeneName))
print(length(common_genes)) # checking the number of genes after duplicates are removed

# function to check data
check_data <- function(df, dataset_name) {
  print(paste("Checking", dataset_name))
  print(paste("Dimensions:", nrow(df), "rows,", ncol(df), "columns"))
  print("Column names:")
  print(head(colnames(df)))
  print("First few rows:")
  print(head(df[, 1:min(5, ncol(df))]))
  print("Data types:")
  print(sapply(df, class))
  print("Number of NA values:")
  print(colSums(is.na(df)))
}

check_data(postop_GSE66921_filtered, "postop_GSE66921_filtered")
check_data(postop_GSE199063_filtered, "postop_GSE199063_filtered")
check_data(postop_GSE59034_filtered, "postop_GSE59034_filtered")

#prepare data for batch correction
prepare_data <- function(df, common_genes) {
  df <- df[df$GeneName %in% common_genes, ]
  rownames(df) <- df$GeneName
  df <- df[, !(names(df) %in% c("ID", "GeneName"))]
  return(as.matrix(df))
}

expr_data_GSE66921 <- prepare_data(postop_GSE66921_filtered, common_genes)
expr_data_GSE199063 <- prepare_data(postop_GSE199063_filtered, common_genes)
expr_data_GSE59034 <- prepare_data(postop_GSE59034_filtered, common_genes)


write.csv(postop_GSE199063_filtered, "postop_GSE199063_filtered.csv", row.names = FALSE)
write.csv(postop_GSE66921_filtered, "postop_GSE66921_filtered.csv", row.names = FALSE)
write.csv(postop_GSE59034_filtered, "postop_GSE59034_filtered.csv", row.names = FALSE)


#############################################################
###  proceed with batch correction using ComBat !!!
#####################

combined_expr_data <- cbind(expr_data_GSE66921, expr_data_GSE199063, expr_data_GSE59034)
batch <- rep(c("GSE66921", "GSE199063", "GSE59034"), 
             c(ncol(expr_data_GSE66921), ncol(expr_data_GSE199063), ncol(expr_data_GSE59034)))
combat_data <- ComBat(dat = combined_expr_data, batch = batch, par.prior = TRUE, prior.plots = FALSE)

postop_GSE66921_normalized <- combat_data[, 1:ncol(expr_data_GSE66921)]
postop_GSE199063_normalized <- combat_data[, (ncol(expr_data_GSE66921)+1):(ncol(expr_data_GSE66921)+ncol(expr_data_GSE199063))]
postop_GSE59034_normalized <- combat_data[, (ncol(expr_data_GSE66921)+ncol(expr_data_GSE199063)+1):ncol(combat_data)]

# add GeneName column 
postop_GSE66921_normalized <- data.frame(GeneName = rownames(postop_GSE66921_normalized), postop_GSE66921_normalized)
postop_GSE199063_normalized <- data.frame(GeneName = rownames(postop_GSE199063_normalized), postop_GSE199063_normalized)
postop_GSE59034_normalized <- data.frame(GeneName = rownames(postop_GSE59034_normalized), postop_GSE59034_normalized)

write.csv(postop_GSE66921_normalized, "postop_GSE66921_normalized.csv", row.names = FALSE)
write.csv(postop_GSE199063_normalized, "postop_GSE199063_normalized.csv", row.names = FALSE)
write.csv(postop_GSE59034_normalized, "postop_GSE59034_normalized.csv", row.names = FALSE)


####################################################################################
###########
##### Repeating the same exact code for Pre op 


preop_GSE59034 <- load_and_inspect("preop_GSE59034.csv")
preop_GSE66921 <- load_and_inspect("preop_GSE66921.csv")
preop_GSE199063 <- load_and_inspect("preop_GSE199063.csv")

datasets <- list(
  preop_GSE66921 = load_dataset("preop_GSE66921.csv"),
  preop_GSE199063 = load_dataset("preop_GSE199063.csv"),
  preop_GSE59034 = load_dataset("preop_GSE59034.csv")
)

# Display structure of each dataset
for (name in names(datasets)) {
  cat("Structure of", name, ":\
")
  print(str(datasets[[name]]))
  cat("\
")
}

merged_data <- datasets[[1]]
for (i in 2:length(datasets)) {
  merged_data <- full_join(merged_data, datasets[[i]], by = c("ID", "GeneName"))
}

for (name in names(datasets)) {
  col_count <- sum(grepl(name, names(merged_data)))
  cat(name, ":", col_count, "\
")
}

write.csv(merged_data, "merged_preop_data.csv", row.names = FALSE)


# Identify and remove duplicate gene names in each dataset
remove_duplicates <- function(df) {
  df <- df[!duplicated(df$GeneName), ]
  return(df)
}

# remove duplicates
preop_GSE66921_filtered <- remove_duplicates(preop_GSE66921)
preop_GSE199063_filtered <- remove_duplicates(preop_GSE199063)
preop_GSE59034_filtered <- remove_duplicates(preop_GSE59034)


# check the number of common genes across datasets
common_genes <- Reduce(intersect, list(preop_GSE66921_filtered$GeneName, preop_GSE199063_filtered$GeneName, preop_GSE59034_filtered$GeneName))
# check each dataset
check_data(preop_GSE66921_filtered, "preop_GSE66921_filtered")
check_data(preop_GSE199063_filtered, "preop_GSE199063_filtered")
check_data(preop_GSE59034_filtered, "preop_GSE59034_filtered")

# prepare for batch correction
expr_data_GSE66921 <- prepare_data(preop_GSE66921_filtered, common_genes)
expr_data_GSE199063 <- prepare_data(preop_GSE199063_filtered, common_genes)
expr_data_GSE59034 <- prepare_data(preop_GSE59034_filtered, common_genes)


write.csv(preop_GSE199063_filtered, "preop_GSE199063_filtered.csv", row.names = FALSE)
write.csv(preop_GSE66921_filtered, "preop_GSE66921_filtered.csv", row.names = FALSE)
write.csv(preop_GSE59034_filtered, "preop_GSE59034_filtered.csv", row.names = FALSE)


#############################################################
###  proceed with batch correction using ComBat for pre-op

# Combine expression data
combined_expr_data <- cbind(expr_data_GSE66921, expr_data_GSE199063, expr_data_GSE59034)
batch <- rep(c("GSE66921", "GSE199063", "GSE59034"), 
             c(ncol(expr_data_GSE66921), ncol(expr_data_GSE199063), ncol(expr_data_GSE59034)))

combat_data <- ComBat(dat = combined_expr_data, batch = batch, par.prior = TRUE, prior.plots = FALSE)

preop_GSE66921_normalized <- combat_data[, 1:ncol(expr_data_GSE66921)]
preop_GSE199063_normalized <- combat_data[, (ncol(expr_data_GSE66921)+1):(ncol(expr_data_GSE66921)+ncol(expr_data_GSE199063))]
preop_GSE59034_normalized <- combat_data[, (ncol(expr_data_GSE66921)+ncol(expr_data_GSE199063)+1):ncol(combat_data)]

preop_GSE66921_normalized <- data.frame(GeneName = rownames(preop_GSE66921_normalized), preop_GSE66921_normalized)
preop_GSE199063_normalized <- data.frame(GeneName = rownames(preop_GSE199063_normalized), preop_GSE199063_normalized)
preop_GSE59034_normalized <- data.frame(GeneName = rownames(preop_GSE59034_normalized), preop_GSE59034_normalized)

write.csv(preop_GSE66921_normalized, "preop_GSE66921_normalized.csv", row.names = FALSE)
write.csv(preop_GSE199063_normalized, "preop_GSE199063_normalized.csv", row.names = FALSE)
write.csv(preop_GSE59034_normalized, "preop_GSE59034_normalized.csv", row.names = FALSE)


##############
#visualize the pre op
# so we have the filtred which is before batch normalisation 
#and then we have the normalised after combat 

#renaming for simplicity 
preop_GSE66921 <- preop_GSE66921_normalized
preop_GSE199063 <- preop_GSE199063_normalized
preop_GSE59034 <- preop_GSE59034_normalized

#  to plot distributions to ceck for normalisation
plot_distributions <- function(df, title) {
  df_melt <- melt(df, id.vars = "GeneName")
  p <- ggplot(df_melt, aes(x = value, color = variable)) +
    geom_density() +
    labs(title = title, x = "Expression Value", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend for clarity
  return(p)
}

# Plot distributions
p1 <- plot_distributions(preop_GSE66921, "Pre OP GSE66921 Normalized")
p2 <- plot_distributions(preop_GSE199063, "Pre OP GSE199063 Normalized")
p3 <- plot_distributions(preop_GSE59034, "Pre OP GSE59034 Normalized")

# Arrange plots in a grid
grid.arrange(p1, p2, p3, ncol = 3)




#######################################################################################
############### PCA to see if it worked ######################
#######################################################################################



# Remember PCA the filtered is pre ComBat and normalised is Post - Combat !!

# Check for NA/NaN values in the datasets
preop_GSE66921_filtered[is.na(preop_GSE66921_filtered)] <- 0
preop_GSE59034_filtered[is.na(preop_GSE59034_filtered)] <- 0
preop_GSE199063_filtered[is.na(preop_GSE199063_filtered)] <- 0

# Ensure all data columns are numeric
preop_GSE66921_filtered[, -c(1,2)] <- lapply(preop_GSE66921_filtered[, -c(1,2)], as.numeric)
preop_GSE59034_filtered[, -c(1,2)] <- lapply(preop_GSE59034_filtered[, -c(1,2)], as.numeric)
preop_GSE199063_filtered[, -c(1,2)] <- lapply(preop_GSE199063_filtered[, -c(1,2)], as.numeric)

#  to prepare data for PCA
prepare_for_pca <- function(df, name) {
  df_t <- t(df[, -c(1,2)])  # Transpose and remove the first two columns (ID and gene names)
  colnames(df_t) <- df$GeneName  # Set gene names as column names
  df_t <- as.data.frame(df_t)
  df_t$dataset <- name
  return(df_t)
}

# Prepare datasets for PCA
gse66921_pca <- prepare_for_pca(preop_GSE66921_filtered, "GSE66921")
gse59034_pca <- prepare_for_pca(preop_GSE59034_filtered, "GSE59034")
gse199063_pca <- prepare_for_pca(preop_GSE199063_filtered, "GSE199063")

# Find common genes across all datasets
common_genes <- Reduce(intersect, list(colnames(gse66921_pca), colnames(gse59034_pca), colnames(gse199063_pca)))
common_genes <- common_genes[common_genes != "dataset"]

# subset to common genes
gse66921_pca <- gse66921_pca[, c(common_genes, "dataset")]
gse59034_pca <- gse59034_pca[, c(common_genes, "dataset")]
gse199063_pca <- gse199063_pca[, c(common_genes, "dataset")]
combined_data <- rbind(gse66921_pca, gse59034_pca, gse199063_pca)

# Perform PCA
pca_result <- prcomp(combined_data[, !colnames(combined_data) %in% c("dataset")], scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$dataset <- combined_data$dataset

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("GSE66921" = "red", "GSE59034" = "blue", "GSE199063" = "darkviolet")) +
  theme_minimal() +
  labs(title = "PCA of Combined Filtered Datasets",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "% variance)"))

print(pca_plot)
print(summary(pca_result)$importance[2, 1:5]) #  variance explained

print(paste("Number of common genes:", length(common_genes)))
total_samples <- nrow(combined_data)
print(paste("Total number of samples:", total_samples))
samples_per_dataset <- combined_data %>% 
  group_by(dataset) %>% 
  summarise(count = n())
print("Number of samples per dataset:")
print(samples_per_dataset)
#########################################################
#########################################################
########### PCA post batch effect #####################

# Double check to see if the correct data is loaded
gse66921_normalized <- preop_GSE66921_normalized
gse59034_normalized <- preop_GSE59034_normalized
gse199063_normalized <- preop_GSE199063_normalized

# prepare function for PCA
gse66921_pca <- prepare_for_pca(gse66921_normalized, "GSE66921")
gse59034_pca <- prepare_for_pca(gse59034_normalized, "GSE59034")
gse199063_pca <- prepare_for_pca(gse199063_normalized, "GSE199063")

# find common genes and subset
common_genes <- Reduce(intersect, list(colnames(gse66921_pca), colnames(gse59034_pca), colnames(gse199063_pca)))
common_genes <- common_genes[common_genes != "dataset"]
gse66921_pca <- gse66921_pca[, c(common_genes, "dataset")]
gse59034_pca <- gse59034_pca[, c(common_genes, "dataset")]
gse199063_pca <- gse199063_pca[, c(common_genes, "dataset")]
combined_data <- rbind(gse66921_pca, gse59034_pca, gse199063_pca)

#  PCA
pca_result <- prcomp(combined_data[, !colnames(combined_data) %in% c("dataset")], scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$dataset <- combined_data$dataset
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = dataset)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("GSE66921" = "red", "GSE59034" = "blue", "GSE199063" = "darkviolet")) +
  theme_minimal() +
  labs(title = "PCA of Combined Normalized Datasets (After ComBat)",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "% variance)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "% variance)"))

print(pca_plot)


print(summary(pca_result)$importance[2, 1:5]) #  variance explained
print(paste("Number of common genes:", length(common_genes)))
total_samples <- nrow(combined_data)
print(paste("Total number of samples:", total_samples))
samples_per_dataset <- combined_data %>% 
  group_by(dataset) %>% 
  summarise(count = n())
print(samples_per_dataset)


####################################################################
##### next step is to use the normalised data with limma ###########
#### and compare the results with other analysis method ###########
####################################################################
#################################################
##### LIMMA WITH THE NORMALISED DATA ###########
#################################################


### better to put them in a list 


postop_GSE66921 <- postop_GSE66921_normalized
preop_GSE66921 <- preop_GSE66921_normalized
postop_GSE59034 <- postop_GSE59034_normalized
preop_GSE59034 <- preop_GSE59034_normalized
preop_GSE199063 <- preop_GSE199063_normalized
postop_GSE199063 <- postop_GSE199063_normalized

files <- c('postop_GSE66921_normalized.csv', 'preop_GSE66921_normalized.csv', 
           'postop_GSE59034_normalized.csv', 'preop_GSE59034_normalized.csv', 
           'preop_GSE199063_normalized.csv', 'postop_GSE199063_normalized.csv')

# setting gene names as row names and remove the GeneName column
process_data <- function(df) {
  rownames(df) <- df$GeneName
  df$GeneName <- NULL
  return(df) }

postop_GSE66921 <- process_data(postop_GSE66921)
preop_GSE66921 <- process_data(preop_GSE66921)
postop_GSE59034 <- process_data(postop_GSE59034)
preop_GSE59034 <- process_data(preop_GSE59034)
preop_GSE199063 <- process_data(preop_GSE199063)
postop_GSE199063 <- process_data(postop_GSE199063)


########################
### DEG on Pre VS Post
#####################

# combine all pre- and post-op  datasets
preop_combined <- cbind(preop_GSE66921, preop_GSE59034, preop_GSE199063)
postop_combined <- cbind(postop_GSE66921, postop_GSE59034, postop_GSE199063)
combined_data <- cbind(preop_combined, postop_combined) # a combined dataset with a new column
condition <- factor(c(rep('pre', ncol(preop_combined)), rep('post', ncol(postop_combined))))
design <- model.matrix(~condition)

# DEG with limma
fit <- lmFit(combined_data, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, number=Inf)
results$gene_name <- rownames(results) 
results <- results[, c('gene_name', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] # reorder columns
write.csv(results, 'DEG_with_LIMMA', row.names = FALSE)


##############
###VOLCANO PLOT
###################

volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, 'significant', 'not significant')), alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'gray') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'gray') +
  scale_color_manual(values = c('significant' = 'red', 'not significant' = 'black')) +
  labs(title = 'Volcano Plot of Differential Gene Expression',
       x = 'Log2 Fold Change',
       y = '-Log10 Adjusted P-value',
       color = 'Significance') +
  theme_minimal() +
  theme(legend.position = 'bottom')
print(volcano_plot)

#######################
##  TOP AND BOTTOM ##
#######################

top_up <- results[results$adj.P.Val < 0.05 & results$logFC > 0, ]
top_up <- top_up[order(-top_up$logFC), ]
top_down <- results[results$adj.P.Val < 0.05 & results$logFC < 0, ]
top_down <- top_down[order(top_down$logFC), ]
print(head(top_up[, c('gene_name', 'logFC', 'adj.P.Val')], 10))
print(head(top_down[, c('gene_name', 'logFC', 'adj.P.Val')], 10))

##############################
###################################
#  gene symbols to Entrez IDs
gene_list <- results$logFC
names(gene_list) <- results$gene_name
gene_list <- sort(gene_list, decreasing = TRUE)
entrez_ids <- mapIds(org.Hs.eg.db, keys = names(gene_list), keytype = "SYMBOL", column = "ENTREZID")

# remove NA values and sort
gene_list <- gene_list[!is.na(entrez_ids)]
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
sorted_genes <- sort(gene_list, decreasing = TRUE)

# select top 20 upregulated and top 20 downregulated genes
top_up <- head(sorted_genes, 20)
top_down <- tail(sorted_genes, 20)
top_genes <- c(top_up, top_down)

# for plotting
plot_data <- data.frame(
  gene = names(top_genes),
  logFC = top_genes,
  regulation = ifelse(top_genes > 0, "Upregulated", "Downregulated") )
p <- ggplot(plot_data, aes(x = reorder(gene, logFC), y = logFC, fill = regulation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  coord_flip() +
  labs(title = "Top 20 Upregulated and Downregulated Genes",
       x = "Gene",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
print(p)


#######################################################################
############Pathway enrichment  ###########################
###########################################################################

sorted_genes <- results[order(results$logFC, decreasing = TRUE), ]# genes by logFC then select top 10 
top_up <- head(sorted_genes, 10)
top_down <- tail(sorted_genes, 10)
top_genes <- rbind(top_up, top_down)


# gene symbols to ENTREZ IDs
gene_ids <- bitr(top_genes$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO enrichment analysis
go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

#  KEGG pathway analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

print(head(go_enrich@result, 10)) #  top GO terms
print(head(kegg_enrich@result, 10)) #  top KEGG pathways

# view the pathways
barplot(go_enrich, showCategory = 10)
print(barplot(go_enrich))
barplot(kegg_enrich, showCategory = 10)
print(barplot(kegg_enrich))


############################################# not included in the report ####################################

##### Using Pathview 
#### pathview gives a png of that specific pathways 

gene_data <- setNames(gene_list, entrez_ids) # makeing a named vector with Entrez IDs
pathway_id <- "hsa04080"  # the pathway ID to visualize
pathway_name <- sub("^hsa", "", pathway_id)
# Visualize the pathway
pathview(gene.data = gene_data,
         pathway.id = pathway_id,
         species = "hsa",  
         out.suffix = pathway_name,
         kegg.native = TRUE)

########################################################################################################
######################################################################################
# visNetwork for an interactive network 

# prepare the nodes dataframe with the top 5 genes and include downregulated genes
# color column based on logFC
# prepare the edges dataframe 
# -> fully connected network for these 9 genes

nodes <- data.frame(
  id = 1:9,
  label = c("RALB", "CKS2", "OR2D3", "MITF", "OR8K1", "EPHB4", "DMC1", "NFKBIZ", "AIM2"),
  logFC = c(3.091187, 2.772502, 2.665304, 2.471437, 2.477737, -1.911062e-02, -9.181207e-05, -2.399448e-04, -2.202732e-05)
)

nodes$color <- ifelse(nodes$logFC > 0, "pink", "green")

edges <- expand.grid(from = 1:9, to = 1:9)
edges <- edges[edges$from != edges$to, ]

# create the network visualization without addNodeDragMode
visNetwork(nodes, edges) %>%
  visNodes(
    shape = "dot",
    size = 30,
    font = list(size = 14)
  ) %>%
  visEdges(
    arrows = "to",
    smooth = TRUE
  ) %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1),
    nodesIdSelection = TRUE
  ) %>%
  visLayout(randomSeed = 123) %>%
  visPhysics(stabilization = FALSE) %>%
  visInteraction(navigationButtons = TRUE) %>%
  visLegend(
    addNodes = data.frame(
      label = c("Upregulated", "Downregulated"),
      shape = "dot",
      color = c("pink", "green")
    ),
    useGroups = FALSE,
    width = 0.1
  )

####################################################################################
####################################################################################
#################### RANKPROD method to compare results ############################

## -> dosent work -> can't load it all at once : vector memory exhausted (limit reached?)


# Combine all pre- and post- operative datasets
preop_combined <- cbind(preop_GSE66921, preop_GSE59034, preop_GSE199063)
postop_combined <- cbind(postop_GSE66921, postop_GSE59034, postop_GSE199063)

combined_data <- cbind(preop_combined, postop_combined)
condition <- factor(c(rep('pre', ncol(preop_combined)), rep('post', ncol(postop_combined))))

condition_numeric <- as.numeric(condition == 'post')
print(dim(combined_data))

#  RankProd analysis
RP.out <- RankProducts(combined_data, condition_numeric, logged = TRUE, na.rm = TRUE, gene.names = rownames(combined_data))
#  the top upregulated and downregulated genes
RP.top <- topGene(RP.out, cutoff = 0.05, method = "pfp", logged = TRUE)

top_up <- RP.top$Table1
top_down <- RP.top$Table2

print(head(top_up, 10))
print(head(top_down, 10))
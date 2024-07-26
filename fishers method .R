library(ggplot2)
library(dplyr)
library(reshape2)
library(stats)

##############################
# using the fisher's method for meta-analysis
##############################################################
###########################################################################
# first testing the top 10 results of each dataset


# Read the datasets
GSE59034 <- read.csv("GSE59034_top_10_genes.csv")
GSE199063 <- read.csv("GSE199063_DEG.csv")
GSE66921 <- read.csv("GSE66921_top_10_genes.csv")

#  to prepare data
prepare_data <- function(df, gene_col) {
  df <- df[, c(gene_col, "logFC", "adj.P.Val")]
  colnames(df)[1] <- "gene"
  df$gene <- trimws(df$gene)  # Remove leading/trailing whitespace
  return(df) }

data1 <- prepare_data(GSE59034, "gene_assignment")
data2 <- prepare_data(GSE199063, "gene_symbol")
data3 <- prepare_data(GSE66921, "gene_symbol")

combined_data <- rbind(
  cbind(data1, dataset = "GSE59034"),
  cbind(data2, dataset = "GSE199063"),
  cbind(data3, dataset = "GSE66921") )

print(head(combined_data))
# check for NA values in the combined dataset
na_summary_combined <- sapply(combined_data, function(x) sum(is.na(x)))
print(na_summary_combined)
str(combined_data)

#  the Fisher's method function
fisher_method <- function(p_values) {
  #  the log-transformed p-values 
  log_p_values <- log(p_values)
  print(log_p_values)
  
  chi_square <- -2 * sum(log_p_values)
  combined_p <- pchisq(chi_square, df = 2 * length(p_values), lower.tail = FALSE)
  return(combined_p) }

# calculate mean fold-change and combine p-values
results <- aggregate(cbind(logFC, adj.P.Val) ~ gene, data = combined_data, 
                     FUN = function(x) c(mean = mean(x), combined_p = fisher_method(x)))

str(results)
print(head(results))

# Define the Fisher's method function with validation
fisher_method <- function(p_values) {
  # Remove any invalid p-values (e.g., non-positive values)
  valid_p_values <- p_values[p_values > 0 & !is.na(p_values)]
  
  # are there any valid p-values left?
  if (length(valid_p_values) == 0) {
    return(NA) }
  
  chi_square <- -2 * sum(log(valid_p_values))
  combined_p <- pchisq(chi_square, df = 2 * length(valid_p_values), lower.tail = FALSE)
  return(combined_p) }

results <- aggregate(cbind(logFC, adj.P.Val) ~ gene, data = combined_data, 
                     FUN = function(x) c(mean = mean(x), combined_p = fisher_method(x)))
results <- data.frame(
  gene = results$gene,
  mean_logFC = results$logFC[,"mean"],
  combined_p = results$adj.P.Val[,"combined_p"] )

results <- results[order(results$combined_p),] # sort results by combined p-value

#  top 10 up- and down-regulated genes
print("Top 10 up-regulated genes:")
top_up <- head(results[results$mean_logFC > 0,], 10)
print(top_up)
top_down <- head(results[results$mean_logFC < 0,], 10)
print(top_down)
write.csv(results, "combined_gene_expression_resultsTOP10.csv", row.names = FALSE)

top_genes <- rbind(top_up, top_down)
top_genes$gene <- factor(top_genes$gene, levels = top_genes$gene[order(top_genes$mean_logFC)]) #  a factor for gene names to control the order in the plot

# barplot
p <- ggplot(top_genes, aes(x = gene, y = mean_logFC, fill = mean_logFC > 0)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = c('blue', 'red'), guide = 'none') +
  labs(title = 'Top Upregulated and Downregulated Genes overall',
       x = 'Gene',
       y = 'Mean Log2 Fold Change') +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
print(p)


############################################################################
############################################################################
############################################################################
#  Test fisher's method on All DEG results from the 3 datasets

# load the datasets
GSE66921 <- read.csv('GSE66921_DEG_ALL.csv')
GSE59034 <- read.csv('GSE59034_results_ALL.csv')
GSE199063 <- read.csv('GSE199063_DEG_ALL.csv')

#prepare datasets
GSE59034 <- GSE59034[, !names(GSE59034) %in% c('ID')]
names(GSE59034)[names(GSE59034) == 'gene_assignment'] <- 'gene_symbol'
GSE199063 <- GSE199063[, !names(GSE199063) %in% c('ID')]
data1 <- prepare_data(GSE59034, "gene_symbol")
data2 <- prepare_data(GSE199063, "gene_symbol")
data3 <- prepare_data(GSE66921, "gene_symbol")

print(head(data1))
print(head(data2))
print(head(data3))

combined_data <- rbind(
  cbind(data1, dataset = "GSE59034"),
  cbind(data2, dataset = "GSE199063"),
  cbind(data3, dataset = "GSE66921") )

###########################
# calculate mean fold-change and combine p-values
results <- aggregate(cbind(logFC, adj.P.Val) ~ gene, data = combined_data, 
                     FUN = function(x) c(mean = mean(x), combined_p = fisher_method(x)))

# reshape results
data.frame(
  gene = results$gene,
  mean_logFC = results$logFC[,"mean"],
  combined_p = results$adj.P.Val[,"combined_p"] ) -> results

# sort results by combined p-value
results <- results[order(results$combined_p),]

# print the structure of the results
str(results)
print(head(results))


# separate upregulated and downregulated genes and sort them
upregulated <- results[results$mean_logFC > 0, ]
downregulated <- results[results$mean_logFC < 0, ]
upregulated <- upregulated[order(-upregulated$mean_logFC), ]
downregulated <- downregulated[order(downregulated$mean_logFC), ]
#  top 10 genes 
top_upregulated <- head(upregulated, 10)
top_downregulated <- head(downregulated, 10)
top_genes <- rbind(top_upregulated, top_downregulated)


ggplot(top_genes, aes(x = reorder(gene, mean_logFC), y = mean_logFC, fill = mean_logFC > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("blue", "red"), name = "Regulation", labels = c("Down", "Up")) +
  labs(title = "Top 10 Up-regulated and Down-regulated Genes",
       x = "Gene",
       y = "Mean Log Fold Change") +
  theme_minimal()

print(top_upregulated[, c("gene", "mean_logFC", "combined_p")])
print(top_downregulated[, c("gene", "mean_logFC", "combined_p")])

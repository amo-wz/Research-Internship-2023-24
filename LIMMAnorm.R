#################################################
##### LIMMA WITH THE NORMALISED DATA ###########
#################################################

# Load the required libraries
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



### better to put them in a list 


postop_GSE66921 <- read.csv('postop_GSE66921_normalized.csv')
preop_GSE66921 <- read.csv('preop_GSE66921_normalized.csv')
postop_GSE59034 <- read.csv('postop_GSE59034_normalized.csv')
preop_GSE59034 <- read.csv('preop_GSE59034_normalized.csv')
preop_GSE199063 <- read.csv('preop_GSE199063_normalized.csv')
postop_GSE199063 <- read.csv('postop_GSE199063_normalized.csv')

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
###DEG on Pre VS Post
#####################

# combine all pre- and post-op  datasets
preop_combined <- cbind(preop_GSE66921, preop_GSE59034, preop_GSE199063)
postop_combined <- cbind(postop_GSE66921, postop_GSE59034, postop_GSE199063)

# create a combined dataset with a new column
combined_data <- cbind(preop_combined, postop_combined)
condition <- factor(c(rep('pre', ncol(preop_combined)), rep('post', ncol(postop_combined))))
design <- model.matrix(~condition)

# DEG with limma
fit <- lmFit(combined_data, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, number=Inf)
results$gene_name <- rownames(results) 

# reorder columns
results <- results[, c('gene_name', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')]
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
nodes <- data.frame(
  id = 1:9,
  label = c("RALB", "CKS2", "OR2D3", "MITF", "OR8K1", "EPHB4", "DMC1", "NFKBIZ", "AIM2"),
  logFC = c(3.091187, 2.772502, 2.665304, 2.471437, 2.477737, -1.911062e-02, -9.181207e-05, -2.399448e-04, -2.202732e-05)
)

# color column based on logFC
nodes$color <- ifelse(nodes$logFC > 0, "pink", "green")

# prepare the edges dataframe (fully connected network for these 9 genes)
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

# Create a combined dataset with a condition column
combined_data <- cbind(preop_combined, postop_combined)
condition <- factor(c(rep('pre', ncol(preop_combined)), rep('post', ncol(postop_combined))))

# Convert condition to numeric
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
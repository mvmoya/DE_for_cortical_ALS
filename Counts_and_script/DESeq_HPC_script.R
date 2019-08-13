##########################################################################################################
###                                                                                                    ###
### This is the script I used to load and analyze the files from the May-June 2019 Gprin3 and Colgalt2 ###
### IP experiments. It has been optimized to work for DESeq2 input and output.                         ###
###                                                                                                    ###
##########################################################################################################

# Here, we load the libraries that we'll need to run all of the analysis. Some of these can be deprecated

library(DESeq2)
library(edgeR)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(Heatplus)
library(RColorBrewer)
library(pheatmap)
library(ggfortify)
library(dplyr)
library(stringr)
library(ggrepel)

# Here, we set the working directory to wherever the counts files are all housed:

#setwd("/Volumes/HeintzBambi3/WORK/Vicky/DESeq_HPC/")
setwd('/Users/mvmoya/Desktop/Bioinformatics/Seq_Counts_files/New Gprin3 vs Colgalt2/HPC')

# This line gets a list of the files in the directory if they have the word "ID" in them:

sampleFiles <- grep("(Gprin3|Colgalt2).*(Counts)", list.files(getwd()), value = TRUE)

# Here, we create a list of what we want the sample names to be. As written, this has to match the
# order of the files on the server. Sorry for the brute force approach:

sampleCondition <- c('Gprin3_nIP_1', 'Gprin3_nIP_2', 'Gprin3_nIP_3',
                     'Gprin3_Input_1', 'Gprin3_Input_2', 'Gprin3_Input_3',
                     'Colgalt2_nIP_1', 'Colgalt2_nIP_2',
                     'Colgalt2_Input_1', 'Colgalt2_Input_2',
                     'Colgalt2_nIP_3', 'Colgalt2_Input_3')
                     #'Colgalt2_oIP_1', 'Colgalt2_oIP_2', 'Colgalt2_oIP_3')

#sampleCondition <- c('Gprin3_IP_1', 'Gprin3_IP_2', 'Gprin3_IP_3',
#                     'Gprin3_Input_1', 'Gprin3_Input_2', 'Gprin3_Input_3',
#                     'Colgalt2_IP_1', 'Colgalt2_IP_2',
#                     'Colgalt2_Input_1', 'Colgalt2_Input_2')

# Now we make a reference table for the sample names from above and the corresponding filename:

sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)

# Here we read the first file in the list and seed a data frame from it. I've changed the header
# names to be easier to reference, and have removed rows where gene names are NA values:

data.num <- read.delim(toString(sampleTable[1, c('fileName')]), header = TRUE)
data.num <- data.num[,2:3] # For HPC output
#data.num <- data.num[1:32658,] # For Heintz pipeline output
data.num <- data.num[complete.cases(data.num),]
colnames(data.num) <- c('genes', sampleCondition[1])

# Now, we add the remaining sample files to this table by merging the tables by gene name.
# This part loops through every filename specified in our reference table above:

for (i in 2:nrow(sampleTable)) {
  new_cols <- read.delim(toString(sampleTable[i,"fileName"]), header = TRUE)
  new_cols <- new_cols[,2:3] # For HPC output
  #new_cols <- new_cols[1:32658,] # For Heintz pipeline output
  colnames(new_cols) <- c('genes', sampleCondition[i])
  data.num <- merge(data.num, new_cols, by = "genes")
}

# Now we make the rownames for the data frame into the gene names:

rownames(data.num) <- data.num[,1]
data.num[,1] <- NULL

# Here, we create the experimental information table. Using regular expressions, I make a list
# that will tell DESeq if our sample was a Gprin3 or Colgalt2 IP or Input sample:

data.num <- subset(data.num, select=-c(Colgalt2_Input_3)) #Colgalt2_nIP_2, Colgalt2_nIP_3, 

sample_type <- as.vector(sub("(.*_.*)_.*", "\\1", colnames(data.num)))
coldata <- data.frame(sample_type)
colnames(coldata) <- factor(c('sample_type'))

# Now, we create the DESeq object, specifying that we want sample_type to be what groups our samples:

ddsPre <- DESeqDataSetFromMatrix(countData = data.num, colData = coldata, design = ~ sample_type)
dds <- DESeq(ddsPre)

# These lines create a table of normalized counts from the DESeq object above:

rld <- rlog(dds, blind = FALSE)
norm_counts <- 2**assay(rld)

# (Optional) And we can print out counts values for specific genes of interest with this command:

genes_of_interest <- c('Crym', 'Nefh', 'Gprin3', 'Colgalt2', 'Fezf2', 'Bcl11b', 'Pcp4', 
                       'Slco2a1', 'Npnt')
glial_genes <- c('Aldh1l1', 'Aif1', 'Gfap', 'Pvalb')
norm_counts[c(genes_of_interest, glial_genes),grep('IP', colnames(norm_counts))]

# Here, we generate a PCA of our different sample types. The color levels should match the groups
# This is a rather convoluted way of doing it, but gives you more control over which PCs are plotted

vsd_of_interest <- rld

rv <- rowVars(assay(vsd_of_interest))
select <- order(rv, decreasing=TRUE)#[seq_len(500)]
pca <- prcomp(t(assay(vsd_of_interest)[select,]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = sample_type, name = colnames(data.num))
ggplot(data = d, aes_string(x = 'PC1', y = 'PC2', color = 'group', label = 'name')) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 2, box.padding = 0.15, 
                  point.padding = 0.5, segment.color = 'grey50') +
  #geom_label_repel(aes(label = name), size = 2, box.padding = 0.35, 
  #                 point.padding = 0.5, segment.color = 'grey50') +
  xlab(paste('PC1: ', round(percentVar[1] * 100), '% variance')) +
  ylab(paste('PC2: ', round(percentVar[2] * 100), '% variance')) +
  coord_fixed() +
  ggtitle('PCA of IP and Input samples')

View(pca$rotation)

#print comp, rotation, pc loading print

# Here, we perform heirarchical clustering of our samples based on gene expression:

sampleDists <- dist(t(assay(vsd_of_interest)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd_of_interest)
colnames(sampleDistMatrix) <- colnames(vsd_of_interest)
pheatmap(sampleDistMatrix, clustering_distance_cols = sampleDists,
         clustering_distance_rows = sampleDists)

# Here, we perform replicate correlation analysis using Pearson Correlation:

arld <- assay(rld)
logcpm <- cpm(arld[,order(colnames(arld))]) #Creates a matrix of normalized log2 counts per million values from the vsd object
logcpm_df <- as.data.frame(logcpm)
pcor <- rcorr(logcpm)
melted_pcor <- melt(pcor$r)
ggplot(data = melted_pcor, aes(x = Var1, y = Var2, fill = value)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle('Pearson Correlation')+geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.995, limit = c(0.99,1), space = "Lab", 
                       name = "Pearson\nCorrelation\nValue")

# Here, we make a heatmap that clusters our samples by expression of genes of interest from above.
# If you want to specify a new set of these genes, do that below:

#genes_of_interest <- c('whatever you want')
pheatmap(subset(norm_counts[,order(colnames(norm_counts))], rownames(norm_counts) %in% c(genes_of_interest, glial_genes, mito_genes)),
         scale = 'row',
         color = colorRampPalette(c('black', 'mediumvioletred', 'darkorange', 'cornsilk'), space='rgb')(44), 
         cluster_cols = FALSE,
         main = 'Marker genes')

# Now, we finally run DE analysis
# We tell DESeq to compare an experimental sample against a control from our "sample_type"s
# Note that the control sample is listed last in the comparison:

res_gprin_colgalt <- results(dds, contrast = c('sample_type', 'Gprin3_nIP', 'Colgalt2_nIP'))
#res_gprin_old_colgalt <- results(dds, contrast = c('sample_type', 'Gprin3_nIP', 'Colgalt2_oIP'))
res_gprin_ip_input <- results(dds, contrast = c('sample_type', 'Gprin3_nIP', 'Gprin3_Input'))
res_colgalt_ip_input <- results(dds, contrast = c('sample_type', 'Colgalt2_nIP', 'Colgalt2_Input'))

res_gprin_colgalt_drop <- as.data.frame(res_gprin_colgalt)
res_gprin_colgalt_drop[is.na(res_gprin_colgalt_drop)] <- 1

#res_gprin_old_colgalt_drop <- as.data.frame(res_gprin_old_colgalt)
#res_gprin_old_colgalt_drop[is.na(res_gprin_old_colgalt_drop)] <- 1

res_gprin_ip_input_drop <- as.data.frame(res_gprin_ip_input)
res_gprin_ip_input_drop[is.na(res_gprin_ip_input_drop)] <- 1

res_colgalt_ip_input_drop <- as.data.frame(res_colgalt_ip_input)
res_colgalt_ip_input_drop[is.na(res_colgalt_ip_input_drop)] <- 1

# Here, let's specify some genes of interest for highlighting in further analysis/plotting:

gprin_genes <- c('Crym', 'Gprin3', 'Slco2a1')
colgalt_genes <- c('Colgalt2', 'Fezf2', 'Pcp4', 'Npnt')
L5b_genes <- c('Nefh', 'Bcl11b')
striat_genes <- c('Drd1', 'Drd2')
mito_genes <- c('Cox6c', 'Ndufa4', 'Cox6a1', 'Ndufb3', 'Uqcrh', 'Atp5g3')

# First, we specify which results we want to make MA plots for:

of_interest <- res_gprin_colgalt_drop

# Now we plot the MA plot from these results,
# Coloring the genes with an adjusted pvalue less than 0.05 pink:

DESeq2::plotMA(res_gprin_colgalt, cex = 0.6, alpha = 0.0, ylim = c(-10, 10), main = toString())
points(of_interest[(of_interest$padj < 0.05) & (of_interest$log2FoldChange > 1 | of_interest$log2FoldChange < -1),], col = 'lightcoral', pch = 16, cex = 0.6)
abline(0, 0)

# These lines color the points for the genes of interest listed above:

points(of_interest[gprin_genes,], col = 'green', pch = 16, cex = 0.6)
text(of_interest[gprin_genes,], labels = gprin_genes, col = 'green', pos = 4, cex = 0.6, font = 2)

points(of_interest[colgalt_genes,], col='purple', pch=16, cex=0.6)
text(of_interest[colgalt_genes,], labels = colgalt_genes, col = 'purple', pos = 4, cex = 0.6, font = 2)

points(of_interest[L5b_genes,], col = 'black', pch = 16, cex = 0.6)
text(of_interest[L5b_genes,], labels = L5b_genes, col = 'black', pos = 4, cex = 0.6, font = 2)

points(of_interest[mito_genes,], col = 'red', pch = 16, cex = 0.6)
text(of_interest[mito_genes,], labels = mito_genes, col = 'red', pos = 4, cex = 0.6, font = 2)

# If you want to look at the values for the genes of interest, you can run this line:

res_gprin_colgalt_drop[c(gprin_genes, colgalt_genes, L5b_genes),]

# This is an optional bar plot:

rowsoi <- of_interest[c(gprin_genes, colgalt_genes, L5b_genes),]
barplot(rowsoi$log2FoldChange,
        names.arg = c(gprin_genes, colgalt_genes, L5b_genes),
        col = c('purple', 'green')[as.factor(rowsoi$log2FoldChange > 0)],
        main = 'Log2FoldChange of marker genes',
        cex.names = 0.6,
        ylim = c(-7, 7))

# Now we write the txt file with our results. We'll add the individual log counts for each replicate:

ind_counts <- as.data.frame(assay(rld))
relevant_cols <- ind_counts[,grep('IP', colnames(ind_counts))]
all_data <- merge(of_interest, relevant_cols, by ='row.names')

write.table(all_data, file = 'Gprin3_all3reps_vs_Colgalt2_all3reps_080619.txt', sep ='\t', row.names = FALSE, col.names = TRUE)

# The End!


newCol <- strsplit(as.vector(coldata$sample_type),"_") %>% 
  unlist %>% 
  matrix(ncol=2,byrow=TRUE) %>% 
  as.data.frame %>% 
  mutate(MouseLine=V1,IP=V2) %>% 
  select(-V1,-V2)

colData(dds) <- newCol

ddss <- DESeqDataSetFromMatrix(counts(dds,normalized=FALSE),
             colData=newCol,
             design = ~0+MouseLine+IP+MouseLine*IP)
ddss <- DESeq(ddss)

resultsNames(ddss)
res <- results(ddss,contrast = list(c("MouseLineColgalt2","MouseLineGprin3"))

               

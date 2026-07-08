# necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggforce)
library(ggVennDiagram)
library(NbClust)
library(ComplexHeatmap)
library(edgeR)
library(DEGreport)
library(fgsea)
library(purrr)
library(tibble)
library(reshape2)
library(patchwork)
library(stringr)
library(ggpattern)
#(inDAGO)
library(readr)
#library(apeglm)
library(ashr)
# setting working directory 
setwd("/home/gospozha/haifa/hiba/op_align_new/")

#### Preparing necessary files ####

# list files with gene counts for each sample
dir = "/home/gospozha/haifa/hiba/op_align_new/"
files = list.files(paste0(dir, "count.tables"), "*ReadsPerGene.out.tab", full.names = T)
countData = data.frame(fread(files[1]))[c(1,3)]

# looping and reading the 3rd column from the remaining files
for(i in 2:length(files)) {
  countData = cbind(countData, data.frame(fread(files[i]))[3])
}

# skipping the first 4 lines, since count data starts on the 5th line
countData = countData[c(5:nrow(countData)),]

# renaming columns as sample names
#colnames(countData) = c("GeneID", gsub(paste0(dir,"count.tables/"), "", files))
# Get just the file names without the path
fn <- basename(files)

# Strip the suffix to get sample names (adjust pattern if needed)
sample_names <- sub("_ReadsPerGene\\.out\\.tab$", "", fn)

# Now set colnames
colnames(countData) <- c("GeneID", sample_names)
colnames(countData) = gsub("ReadsPerGene.out.tab", "", colnames(countData))
rownames(countData) = countData$GeneID
countData = countData[,c(2:ncol(countData))]
names <- colnames(countData)

# writing count matrix to a file
write.csv(countData, file="CountMatrix.csv")


#### 30 vs 1 ####
# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('Metadata1.csv', header=TRUE, sep=",")

# sample names in both objects
samples_meta  <- MetaData$id
samples_count <- colnames(countData)
# find common samples
common_samples <- intersect(samples_meta, samples_count)
# subset and reorder count matrix
countData<- countData[, common_samples]
# reorder metadata to match countData
MetaData <- MetaData[match(common_samples, MetaData$id), ]
# must be TRUE
all(colnames(countData) == MetaData$id)
MetaData$condition <- as.factor(MetaData$condition)
#MetaData$batch <- as.factor(MetaData$batch)
#MetaData$origin <- as.factor(MetaData$sample)

#### Initial quality check ####

# Convert counts to DGEList
dge <- DGEList(counts = countData)

# > dim(dge)
# [1] 39482    26
# remove low counts
smallestGroupSize <- 4
keep <- rowSums(dge$counts >= 10) >= smallestGroupSize
dge <- dge[keep,]
dim(dge)
#[1] 20367    26

# Calculate FPM (Fragments Per Million)
fpm_values <- cpm(dge, normalized.lib.sizes = TRUE)  # edgeR's CPM is equivalent to FPM

# Convert to long format for plotting
fpm_df <- as.data.frame(fpm_values) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "FPM") %>%
  left_join(MetaData, by = c("Sample" = "id"))  # Merge with metadata

ggplot(fpm_df, aes(x = FPM, color = condition)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  theme_minimal() +
  labs(title="Density Plot of FPM Values per Condition",
       x="FPM (log10 scaled)")

# statistical comparison
anova_res <- aov(FPM ~ condition, data = fpm_df)
summary(anova_res)
TukeyHSD(anova_res)
# they are the same

# reads per samples
library_sizes <- colSums(countData)
library_sizes

# Plot
barplot(library_sizes,
        las=2,
        main="Library sizes (tag-seq)",
        ylab="Total reads")
abline(h = 5e5, col="red", lty=2)  # warning threshold

# low-count genes percentage
low_count_fraction <- apply(countData, 2, function(x) mean(x < 10))
low_count_fraction

barplot(low_count_fraction,
        las=2,
        main="Fraction of low-count genes (<10)",
        ylab="Fraction")

abline(h = 0.7, col="red", lty=2)
abline(h = 0.8, col="darkred", lty=2)

# mean-variance check

gene_means <- rowMeans(countData)
gene_vars  <- apply(countData, 1, var)

qplot(log10(gene_means + 1), log10(gene_vars + 1),
      alpha = 0.3,
      main="Mean–variance distribution (Tag-seq)")

# saturation curve

saturation_plot <- inDAGO:::Saturation(
  matrix = countData, 
  method = "sampling", 
  max_reads = 30000000,  
  palette = "Polychrome::palette36"          
)

# Display the plot
print(saturation_plot)

genes <- function(reads) {
  1.2 * reads^0.42
}

genes(2136)
genes(584)
genes(1068)
genes(5000)


# low count samples and library sizes
lib <- colSums(countData)

# fraction of genes <10
frac_low10 <- apply(countData, 2, function(x) mean(x < 10))

data.frame(
  sample = colnames(countData),
  library_size = lib,
  frac_low10 = round(frac_low10, 3)
)

# Quick plot
par(mfrow = c(1,2))
barplot(lib, las=2, main = "Library sizes")
barplot(frac_low10, las=2, main="Fraction of genes <10", ylim=c(0,1))
abline(h=0.7, col="red", lty=2)


#### DESeq2 model ####
# creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ condition)

smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) 
# 16617     9 30 vs 1

# running a model
dds <- DESeq(dds)
res <- results(dds)

# Plotting histograms of p-values
hist(res$pvalue, breaks=50, col="skyblue", main="~ condition",
     xlab="p-value", xlim=c(0,1), ylim=c(0, max(table(cut(res$pvalue, breaks=50)))))

# saving a DESeq2 model to an R object
saveRDS(dds, file = "dds_30vs1.rds")
dds <- readRDS(file = "dds_30vs1.rds")

#### PCA and sample distances using rlog ####

# estimating size factors to determine if it's better to use rlog
SF <- estimateSizeFactors(dds) 
print(sizeFactors(SF))

# vst transformation
rlog <- vst(dds)
mat <- assay(rlog)
# PCA plot
pcaData <- plotPCA(rlog, intgroup=c("condition"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA.full.pdf",width=7)
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  ggtitle("PCA of gene counts") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
pca
ggsave("pca_30vs1.jpg", pca, width = 6.5, height = 6)
norm.counts <- assay(rlog)
write.csv(norm.counts, file="./30vs1.vst.counts.csv")

#dev.off()

# sample distances
sampleDists <- dist(t(assay(rlog)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlog$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#pdf("Dist.all.pdf",width=7)
dist <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)
dist
ggsave("dist2.jpg", dist, width = 6, height = 6)
#dev.off()


# DE genes in contrasts
alpha <- 0.05
lfc_thr <- 1

# all condition levels
conds <- levels(dds$condition)

# all pairwise combinations
pair_list <- combn(conds, 2, simplify = FALSE)

# containers
summary_list <- list()
de_gene_list <- list()

for (p in pair_list) {
  c1 <- p[1]  # denominator
  c2 <- p[2]  # numerator
  
  # log2FC = log2(c2 / c1)
  res <- results(dds,
                 contrast = c("condition", c2, c1),
                 alpha = alpha)
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # filter DE genes
  keep <- !is.na(res_df$padj) &
    res_df$padj < alpha &
    abs(res_df$log2FoldChange) >= lfc_thr
  
  de_df <- res_df[keep, ]
  
  # add metadata
  de_df$contrast <- paste0(c2, "_vs_", c1)
  de_df$direction <- ifelse(de_df$log2FoldChange > 0, "up", "down")
  
  # summary counts
  n_de   <- nrow(de_df)
  n_up   <- sum(de_df$log2FoldChange > 0)
  n_down <- sum(de_df$log2FoldChange < 0)
  
  summary_list[[paste0(c2, "_vs_", c1)]] <- data.frame(
    contrast = paste0(c2, "_vs_", c1),
    n_DE     = n_de,
    n_up     = n_up,
    n_down   = n_down
  )
  
  # full DE gene table
  if (nrow(de_df) > 0) {
    de_gene_list[[paste0(c2, "_vs_", c1)]] <- de_df[, c(
      "contrast", "gene_id", "log2FoldChange", "padj",
      "pvalue", "baseMean", "lfcSE", "stat", "direction"
    )]
  }
}

# combine all results
de_summary <- bind_rows(summary_list)
de_genes_table <- bind_rows(de_gene_list)

# view
de_summary
head(de_genes_table)

# save
# write.csv(de_summary, "de_summary_30v2.csv", row.names = FALSE)
# write.csv(de_genes_table, "de_genes_with_ids_30v1.csv", row.names = FALSE)


#### files by contrast (not separated) ####
res <- results(dds, contrast = c("condition", "30", "1"), alpha = 0.05)
summary(res)

res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1) ~ "Up",
                                log2FoldChange < -log(1) ~ "Down"))
# write.csv(res.ordered, "de_genes_30v1.csv")

#### unfiltered files by contrast (not separated) ####
res <- results(dds, contrast = c("condition", "30", "1"))
summary(res)

res.ordered <- data.frame(res) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

# annotation table
annot <- read_tsv("gene_annotations_description.tsv", show_col_types = FALSE) %>%
  select(gene_id, ipr, description) %>%
  distinct(gene_id, .keep_all = TRUE) 

res.ordered.annot <- res.ordered %>%
    left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_30v1.unfiltered.csv")

#### lfcShrink ####
resLFC <- lfcShrink(dds, coef="condition_30_vs_1", type="ashr")
summary(resLFC)

res.ordered <- data.frame(resLFC) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

# annotation table
annot <- read_tsv("gene_annotations_description.tsv", show_col_types = FALSE) %>%
  select(gene_id, ipr, description) %>%
  distinct(gene_id, .keep_all = TRUE) 

res.ordered.annot <- res.ordered %>%
  left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_30v1.unfiltered.lfcshrink.csv")
#### 25 vs 5 vs 10 ####
# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('Metadata2.csv', header=TRUE, sep=",")

# sample names in both objects
samples_meta  <- MetaData$id
samples_count <- colnames(countData)
# find common samples
common_samples <- intersect(samples_meta, samples_count)
# subset and reorder count matrix
countData<- countData[, common_samples]
# reorder metadata to match countData
MetaData <- MetaData[match(common_samples, MetaData$id), ]
# must be TRUE
all(colnames(countData) == MetaData$id)
MetaData$condition <- as.factor(MetaData$condition)
#MetaData$batch <- as.factor(MetaData$batch)
#MetaData$origin <- as.factor(MetaData$sample)
#### DESeq2 model ####
# creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ condition)

smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) 
# 16617     9 30 vs 1

# running a model
dds <- DESeq(dds)
res <- results(dds)

# Plotting histograms of p-values
hist(res$pvalue, breaks=50, col="skyblue", main="~ condition",
     xlab="p-value", xlim=c(0,1), ylim=c(0, max(table(cut(res$pvalue, breaks=50)))))

# saving a DESeq2 model to an R object
saveRDS(dds, file = "dds_25vs5vs10.rds")
dds <- readRDS(file = "dds_25vs5vs10.rds")

# estimating size factors to determine if it's better to use rlog
SF <- estimateSizeFactors(dds) 
print(sizeFactors(SF))

# vst transformation
rlog <- vst(dds)

mat <- assay(rlog)
# PCA plot
pcaData <- plotPCA(rlog, intgroup=c("condition"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA.full.pdf",width=7)
pca<-ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  ggtitle("PCA of gene counts") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
pca
ggsave("pca_25vs5vs10.jpg", pca, width = 6.5, height = 6)
norm.counts <- assay(rlog)
write.csv(norm.counts, file="./25vs5vs10.vst.counts.csv")

#### DE genes in contrasts ####
alpha <- 0.05
lfc_thr <- 1

# all condition levels
conds <- levels(dds$condition)

# all pairwise combinations
pair_list <- combn(conds, 2, simplify = FALSE)

# containers
summary_list <- list()
de_gene_list <- list()

for (p in pair_list) {
  c1 <- p[1]  # denominator
  c2 <- p[2]  # numerator
  
  # log2FC = log2(c2 / c1)
  res <- results(dds,
                 contrast = c("condition", c2, c1),
                 alpha = alpha)
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # filter DE genes
  keep <- !is.na(res_df$padj) &
    res_df$padj < alpha &
    abs(res_df$log2FoldChange) >= lfc_thr
  
  de_df <- res_df[keep, ]
  
  # add metadata
  de_df$contrast <- paste0(c2, "_vs_", c1)
  de_df$direction <- ifelse(de_df$log2FoldChange > 0, "up", "down")
  
  # summary counts
  n_de   <- nrow(de_df)
  n_up   <- sum(de_df$log2FoldChange > 0)
  n_down <- sum(de_df$log2FoldChange < 0)
  
  summary_list[[paste0(c2, "_vs_", c1)]] <- data.frame(
    contrast = paste0(c2, "_vs_", c1),
    n_DE     = n_de,
    n_up     = n_up,
    n_down   = n_down
  )
  
  # full DE gene table
  if (nrow(de_df) > 0) {
    de_gene_list[[paste0(c2, "_vs_", c1)]] <- de_df[, c(
      "contrast", "gene_id", "log2FoldChange", "padj",
      "pvalue", "baseMean", "lfcSE", "stat", "direction"
    )]
  }
}

# combine all results
de_summary <- bind_rows(summary_list)
de_genes_table <- bind_rows(de_gene_list)

# view
de_summary
head(de_genes_table)

# save
# 
# write.csv(de_summary, "de_summary_10v25v5.csv", row.names = FALSE)
# write.csv(de_genes_table, "de_genes_with_ids_10v25v5.csv", row.names = FALSE)

#### files by contrast (not separated) ####
res <- results(dds, contrast = c("condition", "25", "5"), alpha = 0.05)
summary(res)

res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1) ~ "Up",
                                log2FoldChange < -log(1) ~ "Down"))
# write.csv(res.ordered, "de_genes_25v5.csv")

res <- results(dds, contrast = c("condition", "25", "10"), alpha = 0.05)
summary(res)

res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1) ~ "Up",
                                log2FoldChange < -log(1) ~ "Down"))
# write.csv(res.ordered, "de_genes_25v10.csv")

res <- results(dds, contrast = c("condition", "10", "5"), alpha = 0.05)
summary(res)

res.ordered <- data.frame(res) %>%
  filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(padj) %>%
  mutate(Expression = case_when(log2FoldChange > log(1) ~ "Up",
                                log2FoldChange < -log(1) ~ "Down"))
# write.csv(res.ordered, "de_genes_10v5.csv")

#### unfiltered files by contrast (not separated) ####

# annotation table
annot <- read_tsv("gene_annotations_description.tsv", show_col_types = FALSE) %>%
  select(gene_id, ipr, description) %>%
  distinct(gene_id, .keep_all = TRUE) 

res <- results(dds, contrast = c("condition", "25", "10"))
summary(res)

res.ordered <- data.frame(res) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

res.ordered.annot <- res.ordered %>%
  left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_25v10.unfiltered.csv")

res <- results(dds, contrast = c("condition", "25", "5"))
summary(res)

res.ordered <- data.frame(res) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

res.ordered.annot <- res.ordered %>%
  left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_25v5.unfiltered.csv")

res <- results(dds, contrast = c("condition", "10", "5"))
summary(res)

res.ordered <- data.frame(res) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

res.ordered.annot <- res.ordered %>%
  left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_10v5.unfiltered.csv")

#### lfcShrink ####
resLFC <- lfcShrink(dds, contrast = c("condition", "25", "10"), type="ashr")

res.ordered <- data.frame(resLFC) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

res.ordered.annot <- res.ordered %>%
  left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_25v10.unfiltered.lfcshrink.csv")

resLFC <- lfcShrink(dds, coef="condition_25_vs_5", type="ashr")

res.ordered <- data.frame(resLFC) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

res.ordered.annot <- res.ordered %>%
  left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_25v5.unfiltered.lfcshrink.csv")


resLFC <- lfcShrink(dds, coef="condition_10_vs_5", type="ashr")

res.ordered <- data.frame(resLFC) %>%
  #filter(padj<.05 & abs(log2FoldChange)>=1)  %>%
  arrange(log2FoldChange) %>%
  rownames_to_column("gene_id")

res.ordered.annot <- res.ordered %>%
  left_join(annot, by="gene_id")

# write.csv(res.ordered.annot, "de_genes_10v5.unfiltered.lfcshrink.csv")

 #### biomineralization 30v1 ####
dds <- readRDS(file = "dds_30vs1.rds")
# the same using rlog transformation
rlog <- vst(dds)
# heatmap of Z-scores
biomineralization_gene_list <- read.csv("biomin/genes.biomin.accessions.txt", sep = " ", header = F) 
bio_genes <- c(biomineralization_gene_list$V1) 
bio_genes <- bio_genes[bio_genes %in% rownames(rlog)]
mat <- assay(rlog)[bio_genes, ]
mat_z <- t(scale(t(mat)))  # Z-score by gene
anno_col <- data.frame(condition = colData(rlog)$condition)
rownames(anno_col) <- colnames(rlog)  

# heatmap of selected genes from rlog
pheatmap(mat_z,
         annotation_col = anno_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE)  # optional for clean plots


# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_30v1.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
write.csv(lrt.biomin, file="./biomin/DE.biomin.genes.30v1.csv")

#### somp 30v1 ####
biomineralization_gene_list <- read.csv("somp/genes.biomin.accessions.txt", sep = " ", header = F) 
bio_genes <- c(biomineralization_gene_list$V1) 
bio_genes <- bio_genes[bio_genes %in% rownames(rlog)]
mat <- assay(rlog)[bio_genes, ]
mat_z <- t(scale(t(mat)))  # Z-score by gene
anno_col <- data.frame(condition = colData(rlog)$condition)
rownames(anno_col) <- colnames(rlog)  

# heatmap of selected genes from rlog
pheatmap(mat_z,
         annotation_col = anno_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE)  # optional for clean plots


# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_30v1.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
write.csv(lrt.biomin, file="./somp/DE.somp.genes.30v1.csv")

#### biomineralization 25v10v5 ####
dds <- readRDS(file = "dds_25vs5vs10.rds")
# the same using rlog transformation
rlog <- vst(dds)
# heatmap of Z-scores
biomineralization_gene_list <- read.csv("biomin/genes.biomin.accessions.txt", sep = " ", header = F) 
bio_genes <- c(biomineralization_gene_list$V1) 
bio_genes <- bio_genes[bio_genes %in% rownames(rlog)]
mat <- assay(rlog)[bio_genes, ]
mat_z <- t(scale(t(mat)))  # Z-score by gene
anno_col <- data.frame(condition = colData(rlog)$condition)
rownames(anno_col) <- colnames(rlog)  

# heatmap of selected genes from rlog
pheatmap(mat_z,
         annotation_col = anno_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE)  # optional for clean plots


# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_25v10.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
write.csv(lrt.biomin, file="./biomin/DE.biomin.genes.25v10.csv")

# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_25v5.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
#empty
#write.csv(lrt.biomin, file="./biomin/DE.biomin.genes.25v5.csv")

# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_10v5.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
write.csv(lrt.biomin, file="./biomin/DE.biomin.genes.10v5.csv")


# summary
files <- c(
  "biomin/DE.biomin.genes.25v10.csv",
  "biomin/DE.biomin.genes.30v1.csv",
  "biomin/DE.biomin.genes.10v5.csv"
)

combined_df <- map_dfr(files, function(f) {
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      contrast = str_extract(f, "\\d+v\\d+")
    )
})

# view
combined_df
write.csv(combined_df, file="./biomin/DE.biomin.genes.summary2.csv")

## plotting
combined_df <- read.csv("./biomin/DE.biomin.genes.summary.csv")
# make a clean gene label: gene_id + annotation
res_long <- combined_df %>%
  mutate(
    gene_id = `...1`,
    gene_name = ifelse(is.na(name) | name == "", "unannotated", name),
    facet_label = paste0(gene_id, " \n (", gene_name, ")")
  )

# order contrasts
contrast_levels <- c("10v5","25v10", "30v1")

# colors and shapes
contrast_colors <- c(
  "30v1"  = "#00A6ED",
  "25v5"  = "#009E73",
  "25v10" = "#ee65aa",
  "10v5"  = "#FFB400"
)

contrast_shapes <- c(
  "10v5"  = 24,
  "25v10" = 23,
  "30v1"  = 21
)

res_long$contrast <- factor(res_long$contrast, levels = contrast_levels)

# reorder genes by mean logFC across contrasts
res_long$facet_label <- factor(
  res_long$facet_label,
  levels = res_long %>%
    group_by(facet_label) %>%
    summarise(mean_logFC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_logFC)) %>%
    pull(facet_label)
)

biomin_bar <- ggplot(
  res_long,
  aes(x = facet_label, y = log2FoldChange, fill = contrast, shape = contrast)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.6) +
  geom_point(size = 2.5, position = position_dodge(width = 0.8), color = "black") +
  coord_flip() +
  scale_y_continuous(name = "log2FC") +
  scale_x_discrete(
    name = "Gene",
    labels = function(x) str_wrap(x, width = 45)
  ) +
  scale_fill_manual(values = contrast_colors) +
  scale_shape_manual(values = contrast_shapes) +
  theme_minimal() +
  #labs(title = "DE biomineralization genes by contrast") +
  theme(
    plot.title = element_text(size = 13, margin = margin(t = 9, b = 6)),
    plot.margin = margin(10, 10, 8, 8),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

biomin_bar
ggsave("biomin_genes.pdf", biomin_bar, width = 6.5, height = 5)


#### somp 25v10v5 ####
biomineralization_gene_list <- read.csv("somp/genes.biomin.accessions.txt", sep = " ", header = F) 
bio_genes <- c(biomineralization_gene_list$V1) 
bio_genes <- bio_genes[bio_genes %in% rownames(rlog)]
mat <- assay(rlog)[bio_genes, ]
mat_z <- t(scale(t(mat)))  # Z-score by gene
anno_col <- data.frame(condition = colData(rlog)$condition)
rownames(anno_col) <- colnames(rlog)  

# heatmap of selected genes from rlog
pheatmap(mat_z,
         annotation_col = anno_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE)  # optional for clean plots


# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_25v10.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
#empty
#write.csv(lrt.biomin, file="./somp/DE.somp.genes.25v10.csv")

# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_25v5.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
write.csv(lrt.biomin, file="./somp/DE.somp.genes.25v5.csv")

# are these genes significantly participate in depth change?
res.ordered <- read.csv("de_genes_10v5.csv", row.names = 1)
lrt.biomin <- res.ordered[bio_genes, ]
lrt.biomin <- na.omit(lrt.biomin)
write.csv(lrt.biomin, file="./somp/DE.somp.genes.10v5.csv")


# summary
files <- c(
  "somp/DE.somp.genes.25v5.csv",
  "somp/DE.somp.genes.30v1.csv",
  "somp/DE.somp.genes.10v5.csv"
)

combined_df <- map_dfr(files, function(f) {
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      contrast = str_extract(f, "\\d+v\\d+")
    )
})

# view
combined_df
write.csv(combined_df, file="./somp/DE.somp.genes.summary.csv")

## plotting somp and biomin
combined_df <- read.csv("./somp/DE.biomin.somp.genes.summary.csv")
# make a clean gene label: gene_id + annotation
res_long <- combined_df %>%
  mutate(
    gene_id = `...1`,
    gene_name = ifelse(is.na(name) | name == "", "unannotated", name),
    facet_label = paste0(gene_name, " \n (", gene_id, ")"),
    contrast = dplyr::recode(
      contrast,
      "10v5"  = "T10vC3",
      "25v5"  = "T25vC3",
      "25v10" = "T25vT10",
      "30v1"  = "N30vN3"
  ))

# order contrasts
contrast_levels <- c("N30vN3", "T10vC3", "T25vC3", "T25vT10" )

# colors and shapes
# contrast_colors <- c(
#   "DvS"  = "#00A6ED",
#   "25v3"  = "#009E73",
#   "25v10" = "#ee65aa",
#   "10v3"  = "#FFB400"
# )
contrast_colors <- c(
  "N30vN3"  = "#00A6ED",
  "T25vC3"  = "#009E73",
  "T25vT10" = "#F393C3",
  "T10vC3"  = "#FFB400"
)
# shapes for category
# change names here if your real category labels are slightly different
category_shapes <- c(
  "toolkit" = 3,        
  "proteome" = 4,  
  "both" = 8             # star
)

res_long$contrast <- factor(res_long$contrast, levels = contrast_levels)
res_long$category <- factor(res_long$category, levels = c("toolkit", "proteome", "both"))

# reorder genes by mean logFC across contrasts
res_long$facet_label <- factor(
  res_long$facet_label,
  levels = res_long %>%
    group_by(facet_label) %>%
    summarise(mean_logFC = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_logFC)) %>%
    pull(facet_label)
)

biomin_bar <- ggplot(
  res_long,
  aes(x = facet_label, y = log2FoldChange, fill = contrast, shape = category)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.6) +
  geom_point(size = 2.8, position = position_dodge(width = 0.8), color = "black") +
  coord_flip() +
  scale_y_continuous(name = "log2FC") +
  scale_x_discrete(name = "") +
  scale_fill_manual(values = contrast_colors) +
  scale_shape_manual(values = category_shapes) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, margin = margin(t = 9, b = 6)),
    plot.margin = margin(10, 10, 8, 8),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 9),
    legend.position = "right"
  )

biomin_bar
ggsave("biomin_somp.genes.pdf", biomin_bar, width = 6.5, height = 5)
saveRDS(biomin_bar, "biomin_somp_barplot.RDS")

#### annotation with gene name ####
deg_files <- list(
  "25v5"  = "de_genes_25v5.csv",
  "25v10" = "de_genes_25v10.csv",
  "10v5"  = "de_genes_10v5.csv",
  "30v1" = "de_genes_30v1.csv"
)
# annotation table
annot <- read_tsv("gene_annotations_description.tsv", show_col_types = FALSE) %>%
  select(gene_id, ipr, description) %>%
  distinct(gene_id, .keep_all = TRUE)

# loop
for (nm in names(deg_files)) {
  
  file_path <- deg_files[[nm]]
  
  deg_df <- read_csv(file_path, show_col_types = FALSE)
  
  deg_annot <- deg_df %>%
    left_join(annot, by = "gene_id")
  
  write_csv(deg_annot, paste0("de_genes_", nm, "_annotated.csv"))
}



#### LRT and patterns ####
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ condition)

dds_LRT <- DESeq(dds, test = "LRT", reduced= ~1)
res_lrt <- results(object = dds_LRT, )
summary(res_lrt)

res.ordered <- res_lrt[order(res_lrt$padj),]
# adding Expression column to show the direction of change in expression, if present.
# here, the cutoff values are 0.1 for padj and 1.5 for log2FC.
res.ordered <- data.frame(res.ordered) %>%
  mutate(Expression = case_when(log2FoldChange >= log(1) & padj <= 0.05 ~ "Upregulated",
                                log2FoldChange <= -log(1) & padj <= 0.05 ~ "Downregulated",
                                TRUE ~ "Unchanged"))
head(res.ordered)
write.csv(res.ordered, file="LRT.DE.genes.csv")

vst <-  read.csv("25vs5vs10.vst.counts.csv", row.names = 1)
topgenes <- head(rownames(res_lrt[order(res_lrt$padj), ]), 50)
mat <- vst[topgenes,]
mat <- na.omit(mat)
mat <- t(scale(t(mat)))
df <- as.data.frame(colData(dds_LRT)[,c("condition")])

pheatmap(mat,    
         annotation_col = df)


# 1. Filter the VST matrix for significant LRT genes only
# This reduces noise and speeds up the clustering
sig_genes <- rownames(res_lrt[which(res_lrt$padj < 0.05), ])
vst_sig <- as.matrix(vst[sig_genes, ])

# 2. Run degPatterns
# 'time' is the column name for your depth (e.g., "condition")
# 'col' is usually NULL or for sample replicates
rownames(MetaData) <- MetaData$id
vst_sig <- vst_sig[, rownames(MetaData)]
vst_sig <- na.omit(vst_sig)

# 3. Double-check the match (should return TRUE)
all(colnames(vst_sig) == rownames(MetaData))
clusters <- degPatterns(vst_sig, metadata = MetaData, time = "condition", summarize = "padj", col = NULL)


# posthoc
# 1. Extract the genes belonging to Cluster 1
# 1. Use the matrix that contains ALL 520 significant genes
# We need to scale it so the 'means' are comparable (Z-scores)
mat_full <- t(scale(t(vst_sig)))

# 2. Now extract Cluster 1 genes
cluster1_genes <- clusters$df %>% filter(cluster == 1) %>% pull(genes)

# 3. Subset the full matrix
c1_mat <- mat_full[cluster1_genes, ]

# 4. Calculate the average Z-score per sample
c1_means <- colMeans(c1_mat)

# 5. Create the testing dataframe
# Make sure MetaData is synced with the columns of c1_mat
test_df <- data.frame(
  Zscore = c1_means, 
  Condition = MetaData$condition
)

# 6. Run the ANOVA and Tukey
model <- aov(Zscore ~ Condition, data = test_df)
summary(model)
TukeyHSD(model)

# cluster 3
# 2. Now extract Cluster 1 genes
cluster3_genes <- clusters$df %>% filter(cluster == 3) %>% pull(genes)

# 3. Subset the full matrix
c3_mat <- mat_full[cluster3_genes, ]

# 4. Calculate the average Z-score per sample
c3_means <- colMeans(c3_mat)

# 5. Create the testing dataframe
# Make sure MetaData is synced with the columns of c1_mat
test_df <- data.frame(
  Zscore = c3_means, 
  Condition = MetaData$condition
)

# 6. Run the ANOVA and Tukey
model <- aov(Zscore ~ Condition, data = test_df)
summary(model)
TukeyHSD(model)

# add annotation
colnames(clusters$df) <- c("gene_id", "cluster")
clusters$df %>% left_join(annot, by = "gene_id") -> clusters_annot
write.csv2(clusters_annot, "degPatterns.clusters.csv")



#### Venn diagrams ####
native <- read.csv("de_genes_30v1.csv")
trans_25v5  <- read.csv("de_genes_25v5.csv")
trans_10v5  <- read.csv("de_genes_10v5.csv")
trans_25v10 <- read.csv("de_genes_25v10.csv")

# control
x <- list(
  A = native$gene_id,
  B = trans_10v5$gene_id, 
  G = trans_25v5$gene_id,
  H = trans_25v10$gene_id)

cont <- ggVennDiagram(x,  
                      category.names = c("30v1", "10v5", '25v5', '25v10'),
                      label_alpha = 0, label = "count", set_size = 3.2, label_size = 4)+
  #scale_fill_gradient(low = "#00A9FF", high = "#f75f55", name = "DEGs count")+
  scale_fill_gradient(low = "deepskyblue3", high = "coral1", name = "DEGs count")+
  guides(fill = guide_colorbar(title.position = "top")) +
  #labs(title = "Venn diagrams of shared DEGs")+
  theme(legend.title = element_text(face = "bold"), 
        plot.margin = margin(10, 5, 5, 5),
        plot.title = element_text(size = 13, margin = margin(t = 9, b = 5) ),
        legend.key.height = unit(0.4, "cm"))

cont

ggsave("ggvenn2.jpg", cont, width = 5, height = 5)
library("tidyverse")
library("DESeq2")
library("RColorBrewer")
library("WGCNA")
library("flashClust")
library("gridExtra")
library("ComplexHeatmap")
#library("goseq")
library("dplyr")
library("clusterProfiler")
#library("simplifyEnrichment")
library(dendsort)
library(edgeR)
library(multcompView)
library(ggforce)
library(gridExtra)
library(tibble)
library(ggpattern)

setwd("/home/gospozha/haifa/hiba/op_align_new/wgcna")
load(file="010426.13.70.RData")

#### 25 v 10 v 5 ####
#### reading files ####
# reading count matrix from a file
countData  <- read.csv2('../CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('../Metadata2.csv', header=TRUE, sep=",")

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
# filtering genes 
# original - 4, but here we will filter by 80% of samples 
smallestGroupSize <- 13
keep <- rowSums(countData >= 10) >= smallestGroupSize
countData <- countData[keep,]

#### vst transform (already done during DESeq2 analysis) ####
dds <- DESeqDataSetFromMatrix(countData = countData,
                                              colData = MetaData,
                                              design = ~ condition)
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)   # genes x samples
# 18795
# 80% samples - 12932
gene_var <- apply(expr, 1, var)
keep <- gene_var >= quantile(gene_var, 0.7)
expr_filt <- expr[keep, ]
#expr_filt <- expr
# 5639
# 3880 with 80% samples
#vst transformation using deseq2
#vsd <- read.csv(file="../25vs5vs10.vst.counts.csv", check.names = F, header = T, row.names = 1)
#### check for outliers ####
sampleTree <- hclust(dist(t(expr_filt)), method = "average")
# Plot the sample tree
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree,
     main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2
)
#Plot a line showing the cut-off
abline(h = 220, col = "red")

# remove two outliers since WGCNA is more sensitive to them
# I wont remove op_14 at first
datExpr <- as.data.frame(t(expr_filt))
dim(datExpr) #only 3880 genes have left
#keepSamples <- !(rownames(datExpr) %in% c("op_14"))  
#datExpr <- datExpr[keepSamples, ]
#MetaData <- MetaData[keepSamples, ]

#Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #Should return TRUE

#### pick soft power threshold ####
powers = c(c(1:10), seq(from = 9, to=25, by=1))
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

# mean connectivity plot
sft_df <- data.frame(sft$fitIndices) %>%
  mutate(model_fit = -sign(slope) * SFT.R.sq)

names(sft_df)
ggplot(sft_df, aes(x = Power, y = mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  xlab("Soft Threshold (power)") +
  ylab("Mean connectivity") +
  ggtitle("Mean connectivity") +
  theme_classic()
# manual calculations
#softPower=13 # 70% variation, 4 samples
softPower=9 # 70% variation, 14 samples
#### get modules ####
# automatic modules
bwnet <- blockwiseModules(datExpr,
                          maxBlockSize = 16000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = softPower, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
                          minModuleSize = 30,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          mergeCutHeight = 0.25,
                          deepSplit = 2
)
readr::write_rds(bwnet,
file = file.path("wgcna_results.25v5v10.RDS"))
readr::write_rds(bwnet,
                 file = file.path("wgcna_results.25v5v10.70.14.RDS"))
readr::write_rds(bwnet,
                 file = file.path("wgcna_results.25v5v10.14.RDS"))
bwnet <- read_rds("wgcna_results.25v5v10.70.14.RDS")

# n of genes in each module
table(bwnet$colors) #4119 are not assigned
# 1208 with 14. 18 modules
# plotting dendrogram
mergedColors = labels2colors(bwnet$colors)
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)
# We can also see if our eigengenes relate to our metadata labels. 
# First we double check that our samples are still in order.

all.equal(MetaData$id, rownames(module_eigengenes))

#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# making binary traits vector
MetaData_bin <- MetaData %>%
  dplyr::select(id, condition)  %>%
  mutate(value = 1) %>%  # Create a helper column for filling binary values
  pivot_wider(names_from = condition, values_from = value, values_fill = list(value = 0)) %>%
  column_to_rownames("id")

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, MetaData_bin, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(15,8)

#Displaying correlations and its p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3))

#Displaying the correlation values in a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MetaData_bin),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Plot as clustered Heatmap
#add bold sigignificant p-values, dendrogram with WGCNA MEtree cut-off, module clusters

#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)

#Make list of heatmap row colors
htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))


#### complex heatmap ####
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

ht=Heatmap(moduleTraitCor, name = "Eigengene", column_title = "Module-trait eigengene correlation", 
           col = blueWhiteRed(50), 
           row_names_side = "left", row_dend_side = "left",
           width = unit(4.2, "cm"), height = unit(14, "cm"), 
           column_order = 1:3, column_dend_reorder = FALSE, cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), column_split = 6, column_dend_height = unit(0.5, "in"),
           #cluster_rows = METree, row_split = 6, row_gap = unit(2.5, "mm"), border = TRUE,
           cell_fun = function(j, i, x, y, w, h, col) {
             if(heatmappval[i, j] <= 0.05) {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 6, fontface = "bold"))
             }
             else {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 6, fontface = "plain"))
             }},
           column_names_gp =  gpar(fontsize = 6),
           row_names_gp = gpar(fontsize = 6, alpha = 0.75, border = TRUE, fill = htmap.colors))
ht <- draw(ht)

library(ggplotify)
ht_gg <- as.ggplot(function() {
  draw(ht, heatmap_legend_side = "right")
})


# horizontal heatmap
moduleTraitCor_t <- t(moduleTraitCor)
heatmappval_t <- t(heatmappval)
ht <- Heatmap(
  moduleTraitCor_t,
  name = "Eigengene",
  #column_title = "Module-trait eigengene correlation",
  col = blueWhiteRed(50),
  
  row_names_side = "left",
  column_names_rot = 45,
  
  width = unit(14, "cm"),
  height = unit(1.8, "cm"),
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(
      sprintf("%s", heatmappval_t[i, j]),
      x, y,
      gp = gpar(
        fontsize = 6,
        fontface = ifelse(heatmappval_t[i, j] <= 0.05, "bold", "plain")
      )
    )
  },
  
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 7, fontface = "bold"),
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  )
)

draw(ht)

## try colored design
module_cols <- gsub("^ME", "", colnames(moduleTraitCor_t))

top_anno <- HeatmapAnnotation(
  Module = anno_simple(
    module_cols,
    col = setNames(module_cols, module_cols),
    height = unit(1.5, "mm")
  ),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  moduleTraitCor_t,
  name = "Eigengene",
  col = blueWhiteRed(50),
  
  top_annotation = top_anno,
  
  row_names_side = "left",
  column_names_rot = 45,
  
  width = unit(14, "cm"),
  height = unit(1.8, "cm"),
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.text(
      sprintf("%s", heatmappval_t[i, j]),
      x, y,
      gp = gpar(
        fontsize = 6,
        fontface = ifelse(heatmappval_t[i, j] <= 0.05, "bold", "plain")
      )
    )
  },
  
  column_names_gp = gpar(
    fontsize = 7,
    col = "black"
  ),
  
  row_names_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 7, fontface = "bold"),
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  )
)

draw(ht)

###  Gene relationship to trait and important modules: Gene Significance and Module Membership

#Colors of the modules

modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

traitData <- as.data.frame(MetaData_bin)

geneTraitSignificance <- as.data.frame(cor(datExpr, traitData, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste0("GS.", names(traitData))
names(GSPvalue) <- paste0("p.GS.", names(traitData))

## Summary output of network analysis results
# Make a dataframe that connects traits, genes, and gene annotation

probes <- colnames(datExpr)

geneInfo0 <- data.frame(
  gene_id = probes,
  moduleColor = paste0("ME", mergedColors),
  geneTraitSignificance,
  GSPvalue
)

for (mod in 1:ncol(geneModuleMembership)) {
  geneInfo0[[paste0("MM.", modNames[mod])]] <- geneModuleMembership[, mod]
  geneInfo0[[paste0("p.MM.", modNames[mod])]] <- MMPvalue[, mod]
}

geneInfo <- geneInfo0[order(geneInfo0$moduleColor), ]

geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")

write.csv(geneInfo, file = "geneInfo.14.70.csv", row.names = FALSE)


#### hub genes ####

# This gives correlation between each gene and each module eigengene
kME <- as.data.frame(cor(datExpr, MEs, use = "p"))
kME_pval <- as.data.frame(corPvalueStudent(as.matrix(kME), nSamples = nrow(datExpr)))

topHubGenes <- list()

allColors <- unique(mergedColors)
allColors <- allColors[allColors != "grey"]  # grey = unassigned

for (color in allColors) {
  # Gene indices for the current module
  inModule <- mergedColors == color
  
  # Column of the module eigengene in kME table
  column <- paste0("ME", color)
  
  # Subset kME and p-values for genes in the module
  kME_col <- kME[inModule, column]
  pval_col <- kME_pval[inModule, column]
  
  genes <- colnames(datExpr)[inModule]
  
  # Filter for significant genes only (e.g. p < 0.05)
  sig_idx <- which(pval_col < 0.05 & !is.na(pval_col) & !is.na(kME_col))
  
  if(length(sig_idx) == 0) {
    warning(paste0("No significant genes found in module ", color))
    next
  }
  
  kME_sig <- kME_col[sig_idx]
  genes_sig <- genes[sig_idx]
  
  # Order significant genes by absolute kME (can also use signed kME if preferred)
  topGenes <- genes_sig[order(abs(kME_sig), decreasing = TRUE)]
  
  # Save top 10 significant genes, or fewer if less than 10 significant genes
  topHubGenes[[color]] <- topGenes[1:min(200, length(topGenes))]
}

topHubGenes


hub_df <- do.call(rbind, lapply(names(topHubGenes), function(mod) {
  data.frame(module = mod,
             gene = topHubGenes[[mod]],
             rank = 1:length(topHubGenes[[mod]]))
}))
write.csv(hub_df, "top_hub_genes_long_format.csv", row.names = FALSE)


hubs = chooseTopHubInEachModule(datExpr, mergedColors)


#### analysis of specific modules ####

# single module

# pairs and trios

# Long format
cor_long <- as.data.frame(moduleTraitCor) %>%
  rownames_to_column("module") %>%
  pivot_longer(-module, names_to = "trait", values_to = "correlation")

pval_long <- as.data.frame(moduleTraitPvalue) %>%
  rownames_to_column("module") %>%
  pivot_longer(-module, names_to = "trait", values_to = "pvalue")

# Join
trait_assoc <- left_join(cor_long, pval_long, by = c("module", "trait"))


# modules and traits correlation
trait_assoc %>%
  group_by(module) %>%
  summarise(
    n_strong = sum(abs(correlation) > 0.2 & pvalue < 0.05),
    traits = paste(trait[abs(correlation) > 0.2 & pvalue < 0.05], collapse = ", ")
  ) %>%
  filter(n_strong > 0) -> traitcor


write.csv(traitcor, "module.trait.corr.csv", row.names = FALSE)


save.image(file="010426.13.70.RData") 
#load("rin/070725.vst.RData")

#### bar plot for modules ####
# modules - black, greenyellow, blue, yellow
target_modules <- c("MEblack", "MEgreenyellow", "MEblue", "MEyellow")

# subset MEs directly instead of averaging them into clusters
# transpose immediately because limma expects features as rows and samples as columns
expressionProfile_data <- as.data.frame(t(MEs[, target_modules]))

# set up the design matrix using your numeric conditions
# ensure condition is a factor so limma handles it correctly
MetaData$condition <- factor(MetaData$condition, levels = c("5", "10", "25"))

des_mat <- model.matrix(~ 0 + MetaData$condition)
colnames(des_mat) <- c("Group5", "Group10", "Group25") # avoids R syntax bugs with pure numbers

# run the limma pipeline with explicit contrasts
fit <- limma::lmFit(expressionProfile_data, design = des_mat)

# define the exact pairwise comparisons you want to test
contrast.matrix <- makeContrasts(
  Group10_vs_Group5  = Group10 - Group5,
  Group25_vs_Group10 = Group25 - Group10,
  Group25_vs_Group5  = Group25 - Group5,
  levels = des_mat
)

fit2 <- limma::contrasts.fit(fit, contrast.matrix)
fit2 <- limma::eBayes(fit2)

# loop over each WGCNA module and get adjusted p-values for all contrasts
conditions <- c("5", "10", "25")
contrasts <- combn(conditions, 2, simplify = FALSE)

module_pval_lists <- vector("list", length(target_modules))
names(module_pval_lists) <- target_modules

for (i in seq_along(target_modules)) {
  mod <- target_modules[i]
  pvals <- c()
  
  for (con in contrasts) {
    # generate the contrast name to match our matrix above (e.g., "Group10_vs_Group5")
    sorted_con <- sort(as.numeric(con), decreasing = TRUE)
    contrast_name <- paste0("Group", sorted_con[1], "_vs_Group", sorted_con[2])
    
    # extract the p-value for the specific module safely
    tt <- topTable(fit2, coef = contrast_name, number = Inf, adjust.method = "BH")
    pval <- tt[mod, "adj.P.Val"]
    
    label <- paste(con, collapse = "-")
    pvals[label] <- pval
  }
  
  module_pval_lists[[mod]] <- pvals
}

# generate significance letters for each module
module_letters <- lapply(module_pval_lists, function(pv) {
  multcompLetters(pv, threshold = 0.05)$Letters
})

# pivot data for plotting/visualization
eigengenes_long <- expressionProfile_data %>%
  rownames_to_column("Module") %>%
  pivot_longer(-Module, names_to = "Sample", values_to = "Eigengene") %>%
  left_join(MetaData %>% select(Sample = id, Condition = condition), by = "Sample")


#### joint plots for all clusters ####
#  modules to include 
modules_to_plot <- c("MEblack", "MEgreenyellow", "MEblue", "MEyellow")

# filter long data to keep only the selected 4 modules
eigengenes_subset <- eigengenes_long %>% 
  filter(Module %in% modules_to_plot)

# convert module_letters list into a long dataframe for ggplot annotations
annotation_letters <- bind_rows(
  lapply(names(module_letters), function(mod) {
    data.frame(
      Module = mod,
      Condition = names(module_letters[[mod]]),
      Letter = module_letters[[mod]]
    )
  }),
  .id = "id"
) %>% 
  filter(Module %in% modules_to_plot)

# calculate the highest y point to place the text safely above the boxplots
y_pos <- eigengenes_subset %>%
  group_by(Module, Condition) %>%
  summarise(y = max(Eigengene) + 0.05, .groups = "drop")

# merge the calculated y coordinates with your significance letters
annotation_letters <- annotation_letters %>%
  left_join(y_pos, by = c("Module", "Condition"))

# set order for conditions so the numeric values display chronologically
annotation_letters$Condition <- factor(annotation_letters$Condition, levels = c("5", "10", "25"))
eigengenes_subset$Condition <- factor(eigengenes_subset$Condition, levels = c("5", "10", "25"))

# generate the multi-panel pattern boxplot
bars <- ggplot(eigengenes_subset, aes(x = Condition, y = Eigengene)) +
  geom_boxplot(
    aes(fill = Condition),
    width = 0.4, # slightly widened from 0.1 to look better with 3 groups instead of 6
    alpha = 1,
    lwd = 0.4
  ) +
  geom_text(
    data = annotation_letters,
    aes(x = Condition, y = y, label = Letter),
    nudge_y = 0.05 * (max(eigengenes_subset$Eigengene) - min(eigengenes_subset$Eigengene)),
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Module, ncol = 2, scales = "free_y") + # set to 2 columns for a balanced 2x2 grid
  labs(
    #title = "Eigengene expression across modules",
    y = "Eigengene expression",
    x = "Condition"
  ) +
  guides(fill = "none", pattern = "none", color = "none") + # updated guide syntax from FALSE to 'none'
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"), 
    title = element_text(size = 11)
  )

print(bars)
ggsave("barplots.jpg", bars)
ggsave("barplots.pdf", bars, width = 5, height = 5)
save.image(file="010426.13.70.RData") 
#load("rin/070725.vst.RData")


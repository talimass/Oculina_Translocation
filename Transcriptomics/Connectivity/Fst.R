library(adegenet)
library(StAMPP)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(SNPRelate)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggforce)
library(dplyr)
library(cowplot)

setwd("/home/gospozha/haifa/hiba/op_align_new/snp/pruned")


# PLINK to GDS
snpgdsBED2GDS("Ocupat_vcf_pruned_filt_0.3.bed", 
              "Ocupat_vcf_pruned_filt_0.3.fam", 
              "Ocupat_vcf_pruned_filt_0.3.bim", 
              "temp.gds", cvt.chr="char")

genofile <- snpgdsOpen("temp.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
geno <- snpgdsGetGeno(genofile) # Rows = Samples, Cols = SNPs
snpgdsClose(genofile)

# Load Metadata 
meta <- read.csv("pop.csv")

meta <- meta[match(sample.id, meta$id), ]

# Genotype conversion 
# convert 0,1,2 directly to AA, AB, BB
genotype_matrix <- matrix(NA, nrow=nrow(geno), ncol=ncol(geno))
genotype_matrix[geno == 0] <- "AA"
genotype_matrix[geno == 1] <- "AB"
genotype_matrix[geno == 2] <- "BB"

# StAMPP 'r' format dataframe
# StAMPP format 'r' requires: Sample, Pop, Ploidy, Format, SNP1, SNP2...
st_df <- data.frame(
  Sample = sample.id,
  Pop = meta$condition,  
  Ploidy = 2,
  Format = "BiA",
  stringsAsFactors = FALSE
)

# combine metadata with genotype Matrix
final_df <- cbind(st_df, genotype_matrix)

# filter out NAs in Pop if any
final_df <- final_df[!is.na(final_df$Pop), ]

# convert to StAMPP format
genotype.st <- stamppConvert(final_df, "r")

# run Fst
genotype.fst <- stamppFst(genotype.st, nboots = 1000, percent = 95, nclusters = 3)

print(genotype.fst$Fsts)
print(genotype.fst$Pvalues)


# Calculate Nei's Distance 
neis_dist <- stamppNeisD(genotype.st, pop = TRUE)


#### Fst plot ####

# prepare the data from stampp object
fst_mat <- genotype.fst[["Fsts"]]
p_mat <- genotype.fst[["Pvalues"]]

# ensure P-values are symmetric 
p_mat[upper.tri(p_mat)] <- t(p_mat)[upper.tri(p_mat)]

# combined matrix for display
combined_mat <- fst_mat
#  upper triangle with the P-values
combined_mat[upper.tri(combined_mat)] <- p_mat[upper.tri(p_mat)]

# label matrix
n <- nrow(combined_mat)
label_matrix <- matrix("", nrow = n, ncol = n)

for(i in 1:n) {
  for(j in 1:n) {
    if(i > j) {
      label_matrix[i, j] <- sprintf("%.3f", fst_mat[i, j]) 
    } else if (i < j) {
      label_matrix[i, j] <- ifelse(
        p_mat[i, j] < 0.001,
        "p<0.001",
        paste0("p=", sprintf("%.3f", p_mat[i, j]))
      )
    }
  }
}

# color 
breaks <- pretty(range(fst_mat, na.rm = TRUE), n = 5)

col_fun <- colorRamp2(
  breaks = breaks,
  colors = colorRampPalette(
    c("#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#5a5eaa")
  )(length(breaks))
)


ht <- Heatmap(
  fst_mat, 
  name = "Fst",
  col = col_fun,
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7)
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    if(i < j) {
      grid.rect(x, y, width, height, gp = gpar(fill = "white", col = NA))
    }
    
    label <- label_matrix[i, j]
    if (label != "") {
      text_col <- ifelse(i < j, "grey40", "black")
      grid.text(label, x, y, gp = gpar(fontsize = 6, col = text_col))
    }
    
    if (i == j) {
      grid.lines(
        x = unit.c(x - 0.5 * width, x + 0.5 * width),
        y = unit.c(y + 0.5 * height, y - 0.5 * height),
        gp = gpar(col = "grey80")
      )
    }
  }
)

draw(ht)

#### PCA plot ####
genotype.st  -> genotype.st2
geno.nas <- genotype.st2[, !(names(genotype.st2) %in% c("Sample", "Pop", "pop.num", "ploidy", "format"))]
pop.vector <- as.factor(genotype.st2[[2]])
convert_to_genind_format <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 1, "11",
                ifelse(x == 0.5, "12",
                       ifelse(x == 0, "22", NA))))
}

geno.char <- as.data.frame(lapply(geno.nas, convert_to_genind_format), stringsAsFactors = FALSE)
geno.obj <- df2genind(geno.char, pop = pop.vector, ploidy = 2, NA.char = NA, sep = "")
dist.mtx <- dist(geno.obj)
pca <- dudi.pca(tab(geno.obj, NA.method = "mean"), scannf = F, nf = 2)
s.class(pca$li, fac = pop.vector, col = rainbow(length(unique(pop.vector))))


# Extract individual scores (coordinates)
scores <- as.data.frame(pca$li)
scores$Pop <- pop.vector
scores$SampleID <- rownames(scores)

# Variance explained
eig <- 100 * pca$eig / sum(pca$eig)

# Calculate centroids per population
centroids <- scores %>%
  group_by(Pop) %>%
  dplyr::summarise(
    Axis1 = mean(Axis1),
    Axis2 = mean(Axis2)
  )

# Add centroid coords to scores
scores <- scores %>%
  left_join(centroids, by = "Pop", suffix = c("", ".centroid"))

# Plot

pt_to_mm <- function(pt) pt / 2.845

condition_colors <- c(
  "N3"  = "#FF8900",
  "C3"        = "#FE6100",
  "T10&T25" = "#AA42B8",
  "N30" = "#3f51b6"
)

pca_plot2 <- ggplot(scores, aes(x = Axis1, y = Axis2, color = Pop)) +
  geom_segment(
    aes(xend = Axis1.centroid, yend = Axis2.centroid),
    alpha = 0.4,
    linewidth = 0.3
  ) +
  geom_text(
    data = centroids,
    aes(x = Axis1, y = Axis2, label = Pop, color = Pop),
    size = pt_to_mm(8),
    show.legend = FALSE,
    alpha = 0.6
  ) +
  geom_point(size = 1) +
  labs(
    x = paste0("PC1 (", round(eig[1], 1), "%)"),
    y = paste0("PC2 (", round(eig[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )+
  scale_color_manual(values = condition_colors)

pca_plot2

#### save pics ####
ht_grob <- grid.grabExpr(draw(ht))

plotA <- ggdraw(ht_grob)
plotB <- pca_plot2 

final_combined_plot <- plot_grid(
  plotA, plotB, 
  ncol = 2, 
  labels = "AUTO",       
  label_size = 9,       
  rel_widths = c(1, 1) 
)

fig_width_cm <- 14
fig_height_cm <- 6

save_plot(
  "Fst_PCA_Combined_Plot.pdf",
  final_combined_plot,
  base_width = fig_width_cm / 2.54,
  base_height = fig_height_cm / 2.54
)

save_plot(
  "Fst_PCA_Combined_Plot.png",
  final_combined_plot,
  base_width = fig_width_cm / 2.54,
  base_height = fig_height_cm / 2.54,
  dpi = 600
)

plot(final_combined_plot)


save.image("050426.Rdata")
save.image("090726.Rdata")
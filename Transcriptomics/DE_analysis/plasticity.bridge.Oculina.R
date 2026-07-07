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
library(readr)
library(emmeans)
library(multcomp)
library(multcompView)

# setting working directory 
setwd("/home/gospozha/haifa/hiba/op_align_new/")

# reading count matrix from a file
countData  <- read.csv2('CountMatrix.csv', header=TRUE, row.names=1, sep=',', check.names = F)
# reading metadata file
MetaData <- read.csv2('Metadata.csv', header=TRUE, sep=",")


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
MetaData$condition <- as.factor(MetaData$Condition)
MetaData$year <- as.factor(MetaData$Batch)

#### DESeq2 model ####
# creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = MetaData,
                              design = ~ year + condition)

smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) 
# 20367


#### PCA ####

# Run a global DESeq2/VST just to get normalized log-transformed data
dds_all <- DESeqDataSetFromMatrix(countData = countData, colData = MetaData, design = ~ year + condition)
vsd_all <- vst(dds_all, blind = FALSE)
mat_all <- assay(vsd_all)

# run PCA on uncorrected data
#pca_clean <- prcomp(t(mat_all), scale. = TRUE)
gene_var <- apply(mat_all, 1, var, na.rm = TRUE)

mat_all_pca <- mat_all[
  is.finite(gene_var) & gene_var > 0,
]

top_n <- min(5000, nrow(mat_all_pca))

mat_all_pca <- mat_all_pca[
  order(apply(mat_all_pca, 1, var), decreasing = TRUE)[1:top_n],
]

pca_clean <- prcomp(t(mat_all_pca), scale. = TRUE)
var_exp <- (pca_clean$sdev)^2 / sum((pca_clean$sdev)^2) * 100

# Plotting dataframe
pca_clean_df <- data.frame(
  PC1 = pca_clean$x[,1],
  PC2 = pca_clean$x[,2],
  Condition = MetaData$condition,
  Year = MetaData$year
)

# PERMANOVA check
uncor <- t(mat_all)
# running PERMANOVA on uncorrected dataset 
permanova_terms <- adonis2(
  uncor ~ year + condition,
  data = MetaData,
  method = "euclidean",
  permutations = 999,
  by = "margin"
)
permanova_terms
p_year <- permanova_terms["year", "Pr(>F)"]
p_condition <- permanova_terms["condition", "Pr(>F)"]
p_text <- paste0("P_Year = ", p_year,"\nP_Condition = ", p_condition) 
# ggplot code
pt_to_mm <- function(pt) pt / 2.845

pca_uncor <- ggplot(pca_clean_df, aes(x = PC1, y = PC2, color = Condition, shape = Year)) +
  geom_point(size = 1.2, alpha = 0.8) +
  theme_minimal() +
  labs(#title = "Unified PCA (limma Batch-Corrected)",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"))+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = p_text,
    hjust = -0.05,
    vjust = 1.2,
    size = pt_to_mm(6)
  )+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )


pca_uncor
# now remove batch effect
# Define the biological design we want to PROTECT while removing batch
biological_design <- model.matrix(~ condition, data = MetaData)

# Clean the batch effect! 
# This subtracts the Year1 vs Year2 difference calculated from the Shallow bridge
mat_clean <- removeBatchEffect(mat_all_pca, batch = MetaData$year, design = biological_design)

# Run standard unconstrained PCA on the cleaned matrix
pca_clean <- prcomp(t(mat_clean), scale. = TRUE)
var_exp <- (pca_clean$sdev)^2 / sum((pca_clean$sdev)^2) * 100

# Plotting dataframe
pca_clean_df <- data.frame(
  PC1 = pca_clean$x[,1],
  PC2 = pca_clean$x[,2],
  Condition = MetaData$condition,
  Year = MetaData$year
)

# permanova check

correct <- t(mat_clean)
permanova_terms<-adonis2(
  correct ~ year + condition,
  data = MetaData,
  method = "euclidean",
  permutations = 999,
  by = "margin"
)
permanova_terms
p_year <- permanova_terms["year", "Pr(>F)"]
p_condition <- permanova_terms["condition", "Pr(>F)"]
p_text <- paste0("P_Year = ", p_year,"\nP_Condition = ", p_condition) 

# ggplot code
pca_cor <- ggplot(pca_clean_df, aes(x = PC1, y = PC2, color = Condition, shape = Year)) +
  geom_point(size = 1.2, alpha = 0.8) +
  theme_minimal() +
  labs(#title = "Unified PCA (limma Batch-Corrected)",
       x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)"))+
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = p_text,
    hjust = -0.05,
    vjust = 1.2,
    size = pt_to_mm(6)
  )+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )


pca_cor

#### PCA plot ####
combined_pca <- pca_uncor + pca_cor +
  plot_annotation(tag_levels = 'A')

ggsave("combined_pca_corrected.jpg", combined_pca, width = 10, height = 4)

combined_pca <- pca_uncor + pca_cor +
  plot_annotation(tag_levels = "A") &
  scale_color_discrete(
    labels = c("N3 & C3", "T10", "T25", "N30")
  ) &
  labs(
    shape = "Experiment"
  )

ggsave("combined_pca_corrected.jpg", combined_pca, width = 10, height = 4)

#### plasticity analysis ####
setwd("/home/gospozha/haifa/hiba/op_align_new/plasticity/")
pca <- prcomp(t(mat_clean), scale. = FALSE)

pc_scores <- as.data.frame(pca$x)
pc_scores$sample <- rownames(pc_scores)
pc_scores$condition <- MetaData$condition
pc_scores$year <- MetaData$year

# variance explained by each PC
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var_explained)


n_pcs <- 2

pcs_use <- paste0("PC", seq_len(n_pcs))
weights <- var_explained[seq_len(n_pcs)]

n_pcs
cum_var[n_pcs]

shallow_group <- "3"
deep_group <- "30"

shallow_centroid <- colMeans(
  pc_scores[pc_scores$condition == shallow_group, pcs_use, drop = FALSE]
)

deep_centroid <- colMeans(
  pc_scores[pc_scores$condition == deep_group, pcs_use, drop = FALSE]
)

table(pc_scores$condition)

weighted_dist <- function(x, centroid, weights) {
  sqrt(sum((x - centroid)^2 * weights))
}

pc_scores$dist_to_shallow <- apply(
  pc_scores[, pcs_use, drop = FALSE],
  1,
  weighted_dist,
  centroid = shallow_centroid,
  weights = weights
)

pc_scores$dist_to_deep <- apply(
  pc_scores[, pcs_use, drop = FALSE],
  1,
  weighted_dist,
  centroid = deep_centroid,
  weights = weights
)

pc_scores$deep_like_index <- pc_scores$dist_to_shallow /
  (pc_scores$dist_to_shallow + pc_scores$dist_to_deep)

S <- shallow_centroid
D <- deep_centroid

deep_axis <- D - S

project_deep_axis <- function(x, S, deep_axis) {
  sum((x - S) * deep_axis) / sum(deep_axis^2)
}

pc_scores$deep_axis_score <- apply(
  pc_scores[, pcs_use, drop = FALSE],
  1,
  project_deep_axis,
  S = S,
  deep_axis = deep_axis
)


ggplot(pc_scores, aes(x = condition, y = dist_to_shallow, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Weighted PCA distance from shallow centroid"
  )

ggplot(pc_scores, aes(x = condition, y = deep_axis_score, fill = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Position along shallow-to-deep transcriptomic axis"
  )


#### position plot ####

# choose groups to plot
# change the levels if your exact condition names differ
deep_axis_df <- pc_scores %>%
  filter(condition %in% c("3", "10", "25", "30")) %>%
  mutate(
    condition = factor(condition, levels = c("3", "10", "25", "30"))
  )

#### Summary statistics ####

deep_axis_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean = mean(deep_axis_score, na.rm = TRUE),
    sd = sd(deep_axis_score, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

#### Statistics ####

fit_lm_axis <- lm(deep_axis_score ~ condition, data = deep_axis_df)
summary(fit_lm_axis)

# Assumption checks
shapiro.test(residuals(fit_lm_axis))
car::leveneTest(deep_axis_score ~ condition, data = deep_axis_df)

#### Multcomp letters ####

emm_axis <- emmeans(fit_lm_axis, ~ condition)

letters_axis <- cld(
  emm_axis,
  adjust = "BH",      # or "tukey"; BH is often less conservative
  Letters = letters
) %>%
  as.data.frame() %>%
  mutate(
    .group = gsub(" ", "", .group)
  ) %>%
  dplyr::select(condition, letters = .group)

#### Plot statistics ####

# Extract p-value for condition effect from ANOVA
anova_axis <- anova(fit_lm_axis)
p_lm_axis <- anova_axis["condition", "Pr(>F)"]

p_label_axis <- paste0("LM P = ", signif(p_lm_axis, 3))

deep_axis_stats <- deep_axis_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_p = mean(deep_axis_score, na.rm = TRUE),
    sd_p = sd(deep_axis_score, na.rm = TRUE),
    se_p = sd_p / sqrt(n),
    ci_val = 1.96 * se_p,
    lower = mean_p - ci_val,
    upper = mean_p + ci_val,
    max_group = max(deep_axis_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(letters_axis, by = "condition") %>%
  mutate(
    letter_y = max_group + 0.08 * diff(range(deep_axis_df$deep_axis_score, na.rm = TRUE))
  )

y_label_axis <- max(deep_axis_df$deep_axis_score, na.rm = TRUE) +
  0.18 * diff(range(deep_axis_df$deep_axis_score, na.rm = TRUE))

#### Plot ####
#### Overlapping distribution plot: shallow-to-deep axis score ####


#### Data ####

deep_axis_df <- pc_scores %>%
  filter(condition %in% c("3", "10", "25", "30")) %>%
  mutate(
    condition = factor(condition, levels = c("3", "10", "25", "30"))
  )

#### Statistics ####

fit_lm_axis <- lm(deep_axis_score ~ condition, data = deep_axis_df)
summary(fit_lm_axis)

shapiro.test(residuals(fit_lm_axis))
car::leveneTest(deep_axis_score ~ condition, data = deep_axis_df)

# non normal

#### Kruskal-Wallis test ####

#### Kruskal-Wallis test ####

kruskal_axis <- kruskal.test(
  deep_axis_score ~ condition,
  data = deep_axis_df
)

kruskal_axis

p_kw_axis <- kruskal_axis$p.value
p_label_axis <- paste0("Kruskal-Wallis P = ", signif(p_kw_axis, 1))

#### Pairwise Wilcoxon tests ####

pairwise_res_axis <- pairwise.wilcox.test(
  deep_axis_df$deep_axis_score,
  deep_axis_df$condition,
  p.adjust.method = "BH"
)

pairwise_res_axis

#### Convert pairwise p-value matrix to multcomp letters ####

library(multcompView)

pw_p_axis <- pairwise_res_axis$p.value

# multcompLetters needs a named vector of pairwise p-values
pw_vec_axis <- pw_p_axis[!is.na(pw_p_axis)]

names(pw_vec_axis) <- apply(
  which(!is.na(pw_p_axis), arr.ind = TRUE),
  1,
  function(i) {
    paste(
      rownames(pw_p_axis)[i[1]],
      colnames(pw_p_axis)[i[2]],
      sep = "-"
    )
  }
)

cld_axis <- multcompView::multcompLetters(
  pw_vec_axis,
  threshold = 0.05
)

letter_df_axis <- data.frame(
  condition = names(cld_axis$Letters),
  Letter = cld_axis$Letters
)

letter_df_axis

deep_axis_stats <- deep_axis_df %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_p = mean(deep_axis_score, na.rm = TRUE),
    median_p = median(deep_axis_score, na.rm = TRUE),
    sd_p = sd(deep_axis_score, na.rm = TRUE),
    se_p = sd_p / sqrt(n),
    .groups = "drop"
  ) %>%
  left_join(letter_df_axis, by = "condition")


#### Plot ####

p_axis_density <- ggplot(
  deep_axis_df,
  aes(x = deep_axis_score, fill = condition, color = condition)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey40"
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey40"
  ) +
  geom_density(
    alpha = 0.35,
    linewidth = 0.4,
    adjust = 1
  ) +
  geom_rug(
    aes(color = condition),
    sides = "b",
    alpha = 0.6,
    linewidth = 0.25
  ) +
  geom_vline(
    data = deep_axis_stats,
    aes(xintercept = median_p, color = condition),
    linetype = "solid",
    linewidth = 0.35,
    alpha = 0.8,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = mean(range(deep_axis_df$deep_axis_score, na.rm = TRUE)),
    y = Inf,
    label = p_label_axis,
    vjust = 1.5,
    size = pt_to_mm(6),
    color = "black"
  ) +
  theme_classic() +
  labs(
    x = "Position along shallow-to-deep transcriptomic axis",
    y = "Density",
    fill = "Condition",
    color = "Condition"
  ) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )


density_ymax <- max(ggplot_build(p_axis_density)$data[[3]]$density, na.rm = TRUE)

deep_axis_stats <- deep_axis_stats %>%
  mutate(
    letter_y = density_ymax * 1.05
  )

p_axis_density_letters <- p_axis_density +
  geom_text(
    data = deep_axis_stats,
    aes(x = median_p, y = letter_y, label = Letter),
    color = "black",
    size = pt_to_mm(7),
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.18))
  )

p_axis_density_letters

condition_colors <- c(
  "3"  = "#ff5e4c",
  "10" = "#DC267F",
  "25" = "#785EF0",
  "30" = "#3f51b6"
)


condition_labels <- c(
  "3"  = "N3 & C3",
  "10" = "T10",
  "25" = "T25",
  "30" = "N30"
)

p_axis_density_letters <- p_axis_density_letters +
  scale_fill_manual(
    labels = condition_labels,
    values = condition_colors
  ) +
  scale_color_manual(
    values = condition_colors,
    labels = condition_labels
  )+
  scale_x_continuous(
    breaks = c(0, 1),
    labels = c("Shallow", "Deep"),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.28))
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 10, r = 5, b = 5, l = 5)
  )


pca_uncor <- pca_uncor +
  scale_color_manual(
    values = condition_colors,
    labels = condition_labels
  )&
  labs(
    shape = "Experiment"
  )

pca_cor <- pca_cor +
  scale_color_manual(
    values = condition_colors,
    labels = condition_labels
  )&
  labs(
    shape = "Experiment"
  )


combined_pca_plast <- (pca_uncor + pca_cor) / p_axis_density_letters +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        plot.margin = margin(t = 10, r = 5, b = 5, l = 5))

combined_pca_plast

ggsave(
  "/home/gospozha/haifa/hiba/op_align_new/plasticity/combined_pca_plast2.pdf",
  plot = combined_pca_plast,
  width = 14,
  height = 11,
  units = "cm",
  device = "pdf",
  useDingbats = FALSE
)
ggsave(
  "/home/gospozha/haifa/hiba/op_align_new/plasticity/combined_pca_plast2.jpg",
  plot = combined_pca_plast,
  width = 14,
  height = 11,
  units = "cm",
  device = "jpg"
)


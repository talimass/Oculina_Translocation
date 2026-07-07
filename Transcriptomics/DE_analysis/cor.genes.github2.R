library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
# setting working directory 
setwd("/home/gospozha/haifa/hiba/op_align_new/correlation/")

# # read DEseq2 results

native <- read.csv("../de_genes_30v1.unfiltered.lfcshrink.csv")
trans_25v5  <- read.csv("../de_genes_25v5.unfiltered.lfcshrink.csv")
trans_10v5  <- read.csv("../de_genes_10v5.unfiltered.lfcshrink.csv")
trans_25v10 <- read.csv("../de_genes_25v10.unfiltered.lfcshrink.csv")

# merge in one table
merged <- native %>%
  inner_join(trans_10v5, by = "gene_id", suffix = c("_native", "_10v5")) %>%
  inner_join(trans_25v5, by = "gene_id") %>%
  dplyr::rename(log2FoldChange_25v5 = log2FoldChange, 
         padj_25v5 = padj) %>%
  inner_join(trans_25v10, by = "gene_id") %>%
  dplyr::rename(log2FoldChange_25v10 = log2FoldChange, 
         padj_25v10 = padj)

# interesting gene sets
annot <- read_tsv("../gene_annotations_description.tsv", show_col_types = FALSE) %>%
  dplyr::select(gene_id, ipr, description) %>%
  distinct(gene_id, .keep_all = TRUE) 

biomineralization_gene_list <- read.csv("../biomin/genes.biomin.accessions.txt", sep = " ", header = F) 
somp_gene_list <- read.csv("../somp/genes.biomin.accessions.txt", sep = " ", header = F) 
biomin <- unique(c(biomineralization_gene_list$V1, somp_gene_list$V1))

# subsets
biomin_subset <- merged[merged$gene_id %in% biomin, ]
somp_subset <- merged[merged$gene_id %in% somp_gene_list$V1, ]
toolkit_subset <- merged[merged$gene_id %in% biomineralization_gene_list$V1, ]

#### spearman correlation ####

# --- 30v1 and 25v5 ---
b.model.25v5 <- lm(log2FoldChange_25v5 ~ log2FoldChange_native, data = biomin_subset)
summary(b.model.25v5)
cor.test(biomin_subset$log2FoldChange_native, biomin_subset$log2FoldChange_25v5, method = "spearman")

s.model.25v5 <- lm(log2FoldChange_25v5 ~ log2FoldChange_native, data = somp_subset)
summary(s.model.25v5)
cor.test(somp_subset$log2FoldChange_native, somp_subset$log2FoldChange_25v5, method = "spearman")

t.model.25v5 <- lm(log2FoldChange_25v5 ~ log2FoldChange_native, data = toolkit_subset)
summary(t.model.25v5)
cor.test(toolkit_subset$log2FoldChange_native, toolkit_subset$log2FoldChange_25v5, method = "spearman")


# --- 30v1 and 10v5 ---
b.model.10v5 <- lm(log2FoldChange_10v5 ~ log2FoldChange_native, data = biomin_subset)
summary(b.model.10v5)
cor.test(biomin_subset$log2FoldChange_native, biomin_subset$log2FoldChange_10v5, method = "spearman")

s.model.10v5 <- lm(log2FoldChange_10v5 ~ log2FoldChange_native, data = somp_subset)
summary(s.model.10v5)
cor.test(somp_subset$log2FoldChange_native, somp_subset$log2FoldChange_10v5, method = "spearman")

t.model.10v5 <- lm(log2FoldChange_10v5 ~ log2FoldChange_native, data = toolkit_subset)
summary(t.model.10v5)
cor.test(toolkit_subset$log2FoldChange_native, toolkit_subset$log2FoldChange_10v5, method = "spearman")


# --- 30v1 and 25v10 ---
b.model.25v10 <- lm(log2FoldChange_25v10 ~ log2FoldChange_native, data = biomin_subset)
summary(b.model.25v10)
cor.test(biomin_subset$log2FoldChange_native, biomin_subset$log2FoldChange_25v10, method = "spearman")

s.model.25v10 <- lm(log2FoldChange_25v10 ~ log2FoldChange_native, data = somp_subset)
summary(s.model.25v10)
cor.test(somp_subset$log2FoldChange_native, somp_subset$log2FoldChange_25v10, method = "spearman")

t.model.25v10 <- lm(log2FoldChange_25v10 ~ log2FoldChange_native, data = toolkit_subset)
summary(t.model.25v10)
cor.test(toolkit_subset$log2FoldChange_native, toolkit_subset$log2FoldChange_25v10, method = "spearman")


#### plot ####

# GeneSet column 
somp_plot    <- somp_subset %>% mutate(GeneSet = "proteome")
toolkit_plot <- toolkit_subset %>% mutate(GeneSet = "toolkit")

# Combine and pivot
plot_data <- bind_rows(somp_plot, toolkit_plot) %>%
  dplyr::select(gene_id, GeneSet, log2FoldChange_native, log2FoldChange_10v5, log2FoldChange_25v5, log2FoldChange_25v10) %>%
  pivot_longer(cols = c(log2FoldChange_10v5, log2FoldChange_25v5, log2FoldChange_25v10), 
               names_to = "Contrast", 
               values_to = "L2FC_Translocation") %>%
  mutate(Contrast = factor(dplyr::recode(Contrast, 
                                  "log2FoldChange_10v5" = "T10vC3", 
                                  "log2FoldChange_25v5" = "T25vC3",
                                  "log2FoldChange_25v10" = "T25vT10"),
                           levels = c("T10vC3", "T25vC3", "T25vT10")))
# plot

full_model <- ggplot(plot_data, aes(x = log2FoldChange_native, y = L2FC_Translocation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  
  # toolkit first = background
  geom_point(
    data = plot_data %>% filter(GeneSet == "toolkit"),
    aes(color = GeneSet),
    alpha = 0.45,
    size = 1.2
  ) +
  
  # proteome second = printed on top
  geom_point(
    data = plot_data %>% filter(GeneSet == "proteome"),
    aes(color = GeneSet),
    alpha = 0.45,
    size = 1.2
  ) +
  
  stat_cor(
    data = plot_data,
    aes(
      color = GeneSet,
      label = paste(..r.label.., ..p.label.., sep = "~`,`~")
    ),
    method = "spearman",
    cor.coef.name = "rho",
    label.x.npc = "left",
    size = 2,
    p.accuracy = 0.001
  ) +
  
  facet_wrap(~Contrast) +
  scale_color_manual(values = c("proteome" = "#BF44D2", "toolkit" = "#95a5a6")) +
  coord_cartesian(ylim = c(-1.1, 1.5)) +
  theme_minimal() +
  labs(
    x = "Native depth response (log2FC)",
    y = "Translocation response (log2FC)",
    color = "Gene set"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
# Render plot
print(full_model)

hidden_dots <- plot_data %>%
  filter(L2FC_Translocation < -1.1| L2FC_Translocation > 1.5) %>%
  select(gene_id, GeneSet, Contrast, log2FoldChange_native, L2FC_Translocation)

# View the results
print(hidden_dots)


ggsave("toolkit_proteome.spearman.jpg", full_model, width = 8, height = 5)

saveRDS(full_model, "cor.biomin.plot.RDS")


#### combined plot ####

setwd("/home/gospozha/haifa/hiba/op_align_new/correlation/")
biomin_bar <- readRDS("../biomin_somp_barplot.RDS")
gsea_lollipop <- readRDS("../cell atlas/gsea_lollipop.RDS")
cor_biomin <- readRDS("cor.biomin.plot.RDS")

compact_theme <- theme(
  text = element_text(size = 7),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.text.x = element_text(size = 6),
  legend.title = element_text(size = 7, face = "bold"),
  legend.text = element_text(size = 6),
  strip.text = element_text(size = 7, face = "bold"),
)

# Panel A: gene-level biomineralization plot
biomin_bar2 <- biomin_bar +
  compact_theme +
  theme(
    #legend.position = "bottom",
    axis.title.x = element_text(size = 7),
    axis.title.y = element_blank(),
    plot.margin = margin(2, 4, 2, 2)
  )

# Panel B: correlation plot

cor_biomin2 <- cor_biomin +
  compact_theme +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = -6, r = 0, b = 0, l = 0),
    legend.margin = margin(t = -4, r = 0, b = 0, l = 0),
    legend.spacing.y = unit(0, "mm"),
    plot.margin = margin(2, 2, 2, 2),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7)
  )

# Panel C: GSEA lollipop
gsea_lollipop2 <- gsea_lollipop +
  compact_theme +
  theme(
    #legend.position = "bottom",
    axis.title.x = element_text(size = 7),
    axis.title.y = element_blank(),
    theme(plot.margin = margin(8, 2, 2, 2))
  )

# Right column: B above C
right_col <- wrap_elements(full = cor_biomin2) / gsea_lollipop2 +
  plot_layout(heights = c(0.9, 1.1))

# Full combined figure
combined_biomin <- biomin_bar2 | right_col

combined_biomin <- combined_biomin +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 10, face = "bold"),
    plot.tag.position = c(0, 1)
  )

combined_biomin

# Save: final width 18 cm
ggsave(
  "../biomin/combined_biomineralization_figure.pdf",
  combined_biomin,
  width = 18,
  height = 13.5,
  units = "cm",
  device = cairo_pdf
)

ggsave(
  "../biomin/combined_biomineralization_figure.png",
  combined_biomin,
  width = 18,
  height = 13.5,
  units = "cm",
  dpi = 600
)


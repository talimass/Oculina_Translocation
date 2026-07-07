library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(fgsea)
library(purrr)
library(tidyverse)
library(forcats)
library(stringr)
library(scales)

# setting working directory 
setwd("/home/gospozha/haifa/hiba/op_align_new/cell atlas/")

# read DEseq2 results
native <- read.csv("../de_genes_30v1.unfiltered.lfcshrink.csv")
trans_25v5  <- read.csv("../de_genes_25v5.unfiltered.lfcshrink.csv")
trans_10v5  <- read.csv("../de_genes_10v5.unfiltered.lfcshrink.csv")
trans_25v10 <- read.csv("../de_genes_25v10.unfiltered.lfcshrink.csv")

annot <- read_tsv("../gene_annotations_description.tsv", show_col_types = FALSE) %>%
  dplyr::select(gene_id, ipr, description) %>%
  distinct(gene_id, .keep_all = TRUE) 

#### GSEA analysis ####
res_list <- list(
  native = native,
  trans_25v5 = trans_25v5,
  trans_25v10 = trans_25v10,
  trans_10v5 = trans_10v5
)


gene_set_dir <- "./gene_sets"  
id_col <- "id"                 

auto_gene_sets <- list.files(
  path = gene_set_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
) %>%
  set_names(~ tools::file_path_sans_ext(basename(.x))) %>%
  map(~ read.delim(.x, stringsAsFactors = FALSE) %>%
        mutate(
          id = sub("^.*_(G[0-9]+)$", "ACROYT_\\1", `gene.ID`)
        )) %>%
  map(~ unique(na.omit(.x$id)))

names(auto_gene_sets)
lengths(auto_gene_sets)

biomineralization_gene_list <- read.csv("../biomin/genes.biomin.accessions.txt", sep = " ", header = F) 
somp_gene_list <- read.csv("../somp/genes.biomin.accessions.txt", sep = " ", header = F) 
biomin <- unique(c(biomineralization_gene_list$V1, somp_gene_list$V1))

symbiotic <- read.csv("symbiotic.csv", sep=",")
aposymbiotic <- read.csv("aposymbiotic.csv", sep=",")


manual_gene_sets <- list(
  "biomin proteome" = somp_gene_list$V1,
  "biomin toolkit" = biomineralization_gene_list$V1,
  "symbiotic" = unique(symbiotic$id),
  "aposymbiotic" = unique(aposymbiotic$id)
)

bio_gene_set <- c(auto_gene_sets, manual_gene_sets)

names(bio_gene_set)
lengths(bio_gene_set)

# if we need to select certain sets
# sets_to_keep <- c(
#   "oocytes",
#   "calicoblasts",
#   "gastro_algae",
#   "biomin toolkit",
#   "biomin proteome"
# )
# 
# bio_gene_set <- bio_gene_set[sets_to_keep]


# Function to run fgsea on a DESeq2 result
run_fgsea <- function(res, gene_set, name) {
  res_df <- as.data.frame(res)
  
  # Remove rows with NA values in LFC or padj
  res_df <- res_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
  
  # Compute ranking score - rank according to pval and logfc
  res_df$rank_score <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue + 1e-300 )

  # Create named vector
  ranked <- res_df$rank_score
  names(ranked) <- res_df$gene_id
  
  # Remove infinite or NaN values (can happen if padj = 0)
  ranked <- ranked[is.finite(ranked)]
  
  # Sort in decreasing order
  ranked <- sort(ranked, decreasing = TRUE)
  
  # Run fgsea
  fgsea_res <- fgsea(pathways = gene_set, stats = ranked)
  fgsea_res$contrast <- name
  fgsea_res
}


# 4. Run fgsea on all DE results and bind results into one table
all_fgsea <- imap_dfr(res_list, ~ run_fgsea(.x, bio_gene_set, .y))
# 5. Filter significant results (FDR < 0.05)

all_results_final <- all_fgsea %>%
  mutate(global_padj = p.adjust(pval, method = "BH"))

signif_fgsea <- all_results_final %>%
  filter(global_padj < 0.05)


all_fgsea_long <- signif_fgsea %>%
 unnest(cols = c(leadingEdge)) %>%
 dplyr::rename(gene_id = leadingEdge) %>%
 left_join(annot, by = "gene_id")# Ends up near 5m levels

# 2. Save as a standard CSV
write.csv(all_fgsea_long, "fgsea_selected.lfcshrink.pval.csv", row.names = FALSE)

## lollipop plot

# Clean and prepare data
fgsea_plot_df <- signif_fgsea %>%
  mutate(
    contrast_clean = case_when(
      contrast == "native"     ~ "N30vN3",
      contrast == "trans_25v5" ~ "T25vC3",
      contrast == "trans_10v5" ~ "T10vC3",
      contrast == "trans_25v10" ~ "T25vT10", 
      TRUE ~ contrast
    ),
    contrast_clean = factor(
      contrast_clean,
      levels = c("N30vN3", "T10vC3", "T25vC3", "T25vT10")
    ),
    pathway_label = str_wrap(pathway, width = 45)
  )

# Order pathways by strongest absolute NES
pathway_order <- fgsea_plot_df %>%
  group_by(pathway_label) %>%
  summarise(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs_NES) %>%
  mutate(pathway_y = row_number())

fgsea_plot_df <- fgsea_plot_df %>%
  left_join(pathway_order, by = "pathway_label")

# Manual vertical offsets for contrasts
offset_df <- tibble(
  contrast_clean = factor(c("N30vN3", "T10vC3", "T25vC3", "T25vT10"),
                          levels = c("N30vN3", "T10vC3", "T25vC3", "T25vT10")),
  offset = c(-0.27, -0.09, 0.09, 0.27)
)

fgsea_plot_df <- fgsea_plot_df %>%
  left_join(offset_df, by = "contrast_clean") %>%
  mutate(y_pos = pathway_y + offset)

contrast_colors <- c(
  "N30vN3"  = "#00A6ED",
  "T25vC3"  = "#009E73",
  "T25vT10" = "#F393C3",
  "T10vC3"  = "#FFB400"
)


setwd("/home/gospozha/haifa/hiba/op_align_new/cell atlas/")
# Plot
gsea_lollipop <- ggplot(fgsea_plot_df) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.4
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      y = y_pos,
      yend = y_pos,
      color = contrast_clean
    ),
    linewidth = 0.6,
    alpha = 0.65
  ) +
  
  geom_point(
    aes(
      x = NES,
      y = y_pos,
      size = size,
      fill = contrast_clean
    ),
    shape = 21,
    color = "white",
    stroke = 0.2,
    alpha = 0.95,
    show.legend = TRUE
  ) +
  
  scale_y_continuous(
    breaks = pathway_order$pathway_y,
    labels = pathway_order$pathway_label,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  
  scale_color_manual(values = contrast_colors, name = "Comparison") +
  scale_fill_manual(values = contrast_colors, name = "Comparison") +
  
  scale_size_continuous(
    name = "Genes in set",
    range = c(1, 2.5),
    breaks = c(50, 150, 250),
    limits = c(0, max(fgsea_plot_df$size, na.rm = TRUE))
  ) +
  
  scale_x_continuous(
    breaks = seq(-3, 3, by = 1),
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL
  ) +
  
  guides(
    color = "none",
    fill = guide_legend(
      title = "contrast",
      override.aes = list(size = 5),
      order = 1
    ),
    size = guide_legend(
      title = "genes in set",
      order = 2,
      override.aes = list(fill = "grey70", color = "grey40")
    )
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.box = "vertical"
  )

gsea_lollipop

ggsave ("fgsea.lollipop2.jpg", gsea_lollipop, width = 8, height = 4)
ggsave ("fgsea.lollipop2.full.jpg", gsea_lollipop, width = 5, height = 6)


saveRDS(gsea_lollipop, "gsea_lollipop.RDS")



# leading edge genes - check if they overlap

leading_edge_table <- signif_fgsea %>%
  dplyr::select(pathway, contrast, NES, pval, padj, global_padj, leadingEdge) %>%
  unnest(cols = c(leadingEdge)) %>%
  dplyr::rename(gene_id = leadingEdge) %>%
  left_join(annot, by = "gene_id")

leading_edge_summary <- leading_edge_table %>%
  dplyr::count(pathway, contrast, name = "n_leading_edge_genes") %>%
  arrange(contrast, desc(n_leading_edge_genes))

leading_edge_summary

shared_leading_genes <- leading_edge_table %>%
  filter(
    contrast == "10v3",
    pathway %in% c("oocyte", "cnido_1", "cnido_2")
  ) %>%
  dplyr::count(gene_id, ipr, sort = TRUE) %>%
  filter(n >= 2)

shared_leading_genes

leading_edge_table %>%
  filter(contrast == "DvS",
         pathway %in% c("gastro_algae", "aposymbiotic")) %>%
  dplyr::count(gene_id, ipr, sort = TRUE) %>%
  filter(n >= 2)


gene_set_overlap <- combn(names(bio_gene_set), 2, simplify = FALSE) %>%
  map_dfr(function(x) {
    tibble(
      set1 = x[1],
      set2 = x[2],
      overlap = length(intersect(bio_gene_set[[x[1]]], bio_gene_set[[x[2]]])),
      set1_size = length(bio_gene_set[[x[1]]]),
      set2_size = length(bio_gene_set[[x[2]]])
    )
  }) %>%
  mutate(
    overlap_fraction_set1 = overlap / set1_size,
    overlap_fraction_set2 = overlap / set2_size
  ) %>%
  arrange(desc(overlap))

gene_set_overlap

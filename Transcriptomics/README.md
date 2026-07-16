# RNA-seq analysis of *O. patagonica* native depth comparisons and downward translocation 

## Samples

To analyze the transcriptional response of *O. patagonica*, we sequenced 26 samples belonging to different conditions and experiments. In the native depth comparison, we had 5 replicates of native deep (N30) and 4 replicates of native shallow (N3) colonies (the samples collected and physiologically described in Martinez et al., 2021). In the translocation experiment, we had 6 replicates of native shallow control samples (C3), 6 replicates of 10 m transplants (T10) and 5 replicates of 25 m transplants (T25). Paired 150-bp from 3'-end libraries were obtained using Illumina NovaSeq X.

## Computational resources, tools and scripts

All computational-intensive tasks were performed on Hive2 computer cluster of the Faculty of Natural Sciences at University of Haifa. Necessary tools were installed locally through miniconda. For each task, a slurm script was created and executed via sbatch.

## Quality control, filtering, mapping

The initial QC analysis of raw data and the following mapping is detailed in [QC_Mapping](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/QC_Mapping). Filtered reads were aligned to *O. patagonica* reference genome [GCF_052425735.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_052425735.1/) (Levy et al., 2025) 

## Connectivity (Fst) analysis

The connectivity of coral populations from different depths was measured through a Fst (fixation index) by analyzing RNA-seq-derived SNPs. The analysis was performed using GATK4 pipeline. The detailed description can be found in the corresponding [Connectivity](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/Connectivity) folder.

## Zooxantellae identification

Reads, that did not align to *O. patagonica genome*, were extracted with --outReadsUnmapped Fastx argument and presudoaligned to available Symbiodinaceae transcripts. Details are provided in [Zoox_identification](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/Zoox_identification) folder.

## Differential expression analysis

The folder [DE_analysis](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/DE_analysis) contains the details of DESeq2 DE analysis, The behaviour of candidate biomineralization-related genes was also analysed via spearman correlation. The combined transcriptomic divergence analyses allowed to access the plasticity of the transplanted colonies. 

## WGCNA analysis

VST-transformed counts were analyzed with WGCNA v.1.73 R package to identify modules of co-expressed genes and their association with experimental treatments. The details are present in the [WGCNA](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/WGCNA) folder.

## GO enrichment analysis

Significant DE genes and significant genes from WGCNA module clusters were used as an input for GO enrichment analysis via ClusterProfiler, as described in the [GO_enrichment](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/GO_enrichment) folder. GO annotation of *O. patagonica* genome was obtained using InterPro 5 and eggNOG 5 databases (detail in the [Genome_annotation](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/Genome_annotation) folder. 

## Gene set enrichment analysis

GSEA analysis was performed on selected *O. patagonica* gene sets: [GSEA](https://github.com/talimass/Oculina_Translocation/tree/main/Transcriptomics/GSEA)

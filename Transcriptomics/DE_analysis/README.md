## Differential expression, biomineralization genes, correlation analysis, plasticity analysis

[DESeq2.Oculina.R](DESeq2.Oculina.R) is used to complete two independent DE analyses (one on native depth comparison, the second one on translocation experiment) and two biomineralization genes analyses with X and X gene sets. [Metadata1.csv](Metadata1.csv) and [Metadata2.csv](Metadata2.csv) are used for separate analyses. 

[cor.genes.github2.R](cor.genes.github2.R) is the script for correlation analysis between conditions and candidate biomineralization gene sets expression.

[plasticity.bridge.Oculina.R](plasticity.bridge.Oculina.R) is the analysis of broad transcriptomic divergence in the combined experiment design where native shallow samples from two experiments are used as bridging samples. This analysis is using combined [Metadata.csv](Metadata.csv)

[proteome.genes.biomin.accessions.txt](proteome.genes.biomin.accessions.txt) represents genome orthologs of previously identifyed genes from SOMP *O. patagonica* proteome (Zaquin et al., 2022). [toolkit.genes.biomin.accessions.txt](toolkit.genes.biomin.accessions.txt) represents orthologs of known coral biomineralization genes (Drake et al., 2013; Mass et al., 2013; Ramos-Silva et al., 2013; Takeuchi et al., 2016; Peled et al., 2020; Mummadisetti et al., 2021). 

Here and in the following scripts, condition labels were updated for final plotting to more accurately reflect the experimental design after most analyses had been completed. Samples previously labeled “1” were collected at 1–3 m and are now designated N3. Samples previously labeled “5” were collected at the same depth and are now designated C3:

1 → N3  
5 → C3  
10 → T10  
25 → T25  
45 → T45  

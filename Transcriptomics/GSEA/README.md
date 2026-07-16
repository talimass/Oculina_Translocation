## Gene set enrichment analysis

The script [GSEA.github.R](GSEA.github.R) uses several gene sets:

1. [proteome.genes.biomin.accessions.txt](proteome.genes.biomin.accessions.txt) represents genome orthologs of previously identifyed genes from SOMP *O. patagonica* proteome (Zaquin et al., 2022).  
2. [toolkit.genes.biomin.accessions.txt](toolkit.genes.biomin.accessions.txt) represents orthologs of known coral biomineralization genes (Drake et al., 2013; Mass et al., 2013; Ramos-Silva et al., 2013; Takeuchi et al., 2016; Peled et al., 2020; Mummadisetti et al., 2021).  
3. symbiotic and aposymbiotic gene sets - genes, significantly enriched in naturally symbiotic or aposymbiotic colonies (Levy et al., 2025)
4. gene sets from *O. patagonica* cell atlas (Levy et al., 2025; https://sebelab.crg.eu/multicoral-sc-atlas/): oocytes (cell type = germline_oocytes), calicoblasts (cell type = calicoblast), digestive filaments (cell type = digestive_filaments), glands (cell type = gland_1_Xbp – gland_8) and immune-related genes (cell type = immune_like_1, immune_like_2, immune_macrophage_like)

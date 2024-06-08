setwd("~/Desktop/PhD_Project_related/GENERATION R DATA RESULTS/Data/Data/ONLY THESE RESULTS TO USE WITH CORRECT COEFFICIENT/RESULTS FROM IRENE/All the models Irene ran")


pm25_adj=read.csv(file="GenR_M3_pm25_glmFit_18787.csv",header=T,sep=" ",stringsAsFactors = F)
pm10_adj=read.table(file="GenR_M3_pm10_glmFit_18787.csv",header=T,sep=" ",stringsAsFactors = F)
no2_adj=read.table(file="GenR_M3_no2_glmFit_18787.csv",header=T,sep=" ",stringsAsFactors = F)

pm25_pos=subset(pm25_adj,pm25_adj$logFC > 0)
pm25_neg= subset(pm25_adj,pm25_adj$logFC < 0)

library(tidyr)

pm25_adj$Gene=rownames(pm25_adj)
pm10_adj$Gene=rownames(pm10_adj)
no2_adj$Gene=rownames(no2_adj)




library(stringr)
pm25_adj[c('Gene', 'Dot')] <- str_split_fixed(pm25_adj$Gene, '\\.', 2)
pm10_adj[c('Gene', 'Dot')] <- str_split_fixed(pm10_adj$Gene, '\\.', 2)
no2_adj[c('Gene', 'Dot')] <- str_split_fixed(no2_adj$Gene, '\\.', 2)

pm25_adj$p.adjusted=p.adjust(pm25_adj$PValue,method="BH")
pm10_adj$p.adjusted=p.adjust(pm10_adj$PValue,method="BH")
no2_adj$p.adjusted=p.adjust(no2_adj$PValue,method="BH")


library(org.Hs.eg.db)

pm25_adj$Entrez <- mapIds(org.Hs.eg.db, pm25_adj$Gene,keytype="ENSEMBL", column="ENTREZID")
pm25_adj$Symbol <- mapIds(org.Hs.eg.db, pm25_adj$Entrez,keytype="ENTREZID", column="SYMBOL")

pm10_adj$Entrez <- mapIds(org.Hs.eg.db, pm10_adj$Gene,keytype="ENSEMBL", column="ENTREZID")
pm10_adj$Symbol <- mapIds(org.Hs.eg.db, pm10_adj$Entrez,keytype="ENTREZID", column="SYMBOL")

no2_adj$Entrez <- mapIds(org.Hs.eg.db, no2_adj$Gene,keytype="ENSEMBL", column="ENTREZID")
no2_adj$Symbol <- mapIds(org.Hs.eg.db, no2_adj$Entrez,keytype="ENTREZID", column="SYMBOL")

pm25_adj$Entrez=as.character(pm25_adj$Entrez)
pm10_adj$Entrez=as.character(pm10_adj$Entrez)
no2_adj$Entrez=as.character(no2_adj$Entrez)

pm25_adj$Symbol=as.character(pm25_adj$Symbol)
pm10_adj$Symbol=as.character(pm10_adj$Symbol)
no2_adj$Symbol=as.character(no2_adj$Symbol)


pm25_sig=subset(pm25_adj,pm25_adj$p.adjusted < 0.25)



library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(magrittr)

msigdbr_species()

hs_msigdb_df <- msigdbr(species = "Homo sapiens")

head(hs_msigdb_df)


#Filter the human data frame to the KEGG pathways that are included in the
# curated gene sets
hs_GO_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C5", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("GO:BP","GO:CC","GO:MF") # This is because we only want KEGG pathways
  )

hs_KEGG_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("CP:KEGG") # This is because we only want KEGG pathways
  )


hs_reactome_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("CP:REACTOME") # This is because we only want KEGG pathways
  )




## FIRST PERFORMING ENRICHMENT ANALYSIS
GO_ora_results_pm25 <- enricher(
  gene = pm25_sig$Gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = pm25_adj$Gene,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_KEGG_df,
    gs_name,
    ensembl_gene
  )
)

View(GO_ora_results_pm25@result)
enrich_plot <- enrichplot::dotplot(GO_ora_results_pm25, showCategory=15,font.size=12,title="Enrichment using MsigDB for genes significantly associated with PM2.5",orderBy= "p.adjust")
enrich_plot ## Dis



KEGG_ora_results_pm25 <- enricher(
  gene = pm25_sig$Gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = pm25_adj$Gene,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_KEGG_df,
    gs_name,
    ensembl_gene
  )
)

enrich_plot <- enrichplot::dotplot(KEGG_ora_results_pm25, showCategory=15,font.size=12,title="Enrichment using MsigDB for genes significantly associated with PM2.5",orderBy= "p.adjust")
enrich_plot ## Dis





GO_ora_results_pm10 <- enricher(
  gene = pm10_sig$Gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = pm10_adj$Gene,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)


enrich_plot <- enrichplot::dotplot(GO_ora_results_pm10, showCategory=15,font.size=12,title="Enrichment using MsigDB for genes significantly associated with PM10",orderBy= "p.adjust")
enrich_plot ## Dis






GO_ora_results_no2 <- enricher(
  gene = no2_sig$Gene, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = no2$Gene,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_GO_df,
    gs_name,
    ensembl_gene
  )
)


enrich_plot <- enrichplot::dotplot(GO_ora_results_pm25, showCategory=15,font.size=12,title="Enrichment using MsigDB for genes significantly associated with NO2",orderBy= "p.adjust")
enrich_plot ## Dis


## Now performing GSEA
# Let's create a named vector ranked based on the t-statistic values
no2_g_vector <- no2_adj$logFC
names(no2_g_vector) <- no2_adj$Gene

pm10_g_vector <- pm10_adj$logFC
names(pm10_g_vector) <- pm10_adj$Gene

pm25_g_vector <- pm25_adj$logFC
names(pm25_g_vector) <- pm25_adj$Gene


# We need to sort the t-statistic values in descending order here
no2_g_vector <- sort(no2_g_vector, decreasing = TRUE)
pm10_g_vector <- sort(pm10_g_vector, decreasing = TRUE)
pm25_g_vector <- sort(pm25_g_vector, decreasing = TRUE)


# Set the seed so our results are reproducible:
set.seed(2020)


# Attach the library
library(clusterProfiler)

# Package that contains MSigDB gene sets in tidy format
library(msigdbr)

# We will need this so we can use the pipe: %>%
library(magrittr)


dr_hallmark_df <- msigdbr(species = "Homo sapiens",  category = "H")
dr_immune_df <- msigdbr(species = "Homo sapiens",  category = "C7",subcategory = "IMMUNESIGDB")
dr_celltype_df <- msigdbr(species = "Homo sapiens",  category = "C8")
dr_canonicalpathway_df <- msigdbr(species = "Homo sapiens",  category = "C2",subcategory = "CP")
dr_keggpathway_df <- msigdbr(species = "Homo sapiens",  category = "C2",subcategory = "CP:KEGG")
dr_regulatory_df <- msigdbr(species = "Homo sapiens",  category = "C3",subcategory = "TFT:GTRD")
dr_reactome_df <- msigdbr(species = "Homo sapiens",  category = "C2",subcategory = "CP:REACTOME")
dr_biocarta_df <- msigdbr(species = "Homo sapiens",  category = "C2",subcategory = "CP:BIOCARTA")
dr_wikipathway_df <- msigdbr(species = "Homo sapiens",  category = "C2",subcategory = "CP:WIKIPATHWAYS")


gsea_no2 <- GSEA(
  geneList = no2_g_vector, # Ordered ranked gene list
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    dr_wikipathway_df,
    gs_name,
    ensembl_gene
  )
)
View(gsea_no2@result)
enriched_plot <- enrichplot::dotplot(gsea_no2,font.size=12,title="GSEA using MsigDB Reactome DB for NO2 association results",orderBy= "p.adjust") + facet_grid(.~.sign)
enriched_plot



gsea_pm25 <- GSEA(
  geneList = pm25_g_vector, # Ordered ranked gene list
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    dr_immune_df,
    gs_name,
    ensembl_gene
  )
)

View(gsea_pm25@result)
library(ggplot2)
enriched_plot <- enrichplot::dotplot(gsea_pm25,font.size=15,title="",orderBy= "p.adjust",decreasing=FALSE,showCategory=10) + facet_grid(.~.sign)
enriched_plot

pm25_tab_GenR= gsea_pm25@result
genR_pos=subset(pm25_tab_GenR,pm25_tab_GenR$NES > 0)
genR_neg=subset(pm25_tab_GenR,pm25_tab_GenR$NES < 0)

pm25_top_10=pm25_tab_GenR[1:10,]



p_pm25=ggplot(pm25_top_10, aes(y=reorder(Description,NES), x=NES,fill=p.adjust)) +
  theme_classic(base_size = 16)+
  geom_bar(stat="identity")+
  scale_fill_continuous(low="blue", high="red")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12))+
  xlab("Normalised enrichment score")+ylab("Immunologic DB terms")
p_pm25 

tiff(filename="PM2.5 Immunologic GSEA barplot.tiff",width = 6000,height = 2000,res = 300)
p_pm25
dev.off()





#write.table(gsea_pm25@result,file="PM2.5_Model3_Immunologic_GSEA_results.txt",col.names = T,row.names = F,sep="\t",quote = F)
write.table(gsea_pm25@result,file="PM2.5_Model3_Hallmark_GSEA_results.txt",col.names = T,row.names = F,sep="\t",quote = F)


View(gsea_pm25@result) ## Interferon alha and beta, TNF alpha suppressed
## with 15 min, no of sig terms = 

gsea_pm10 <- GSEA(
  geneList = pm10_g_vector, # Ordered ranked gene list
  minGSSize = 5, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p-value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    dr_wikipathway_df,
    gs_name,
    ensembl_gene
  )
)

View(gsea_pm10@result) 

enriched_plot <- enrichplot::dotplot(gsea_pm10,font.size=12,title="GSEA using MsigDB Reactome DB for PM10 association results",orderBy= "p.adjust",showCategory=15) + facet_grid(.~.sign)
enriched_plot


View(gsea_pm10@result)
View(gsea_pm25@result)
View(gsea_no2@result)


library(clusterProfiler)
library(enrichplot)
library(ggplot2)

## Nothing for canonical pathways, kegg pathways, TFT:GTRD
## Celltype and immune signature is done




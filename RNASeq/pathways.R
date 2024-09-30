# BiKE app modules

#########################
## Load data
#########################

setwd("datasets/KI/BiKE_Omics_App/BiKE_Lookup_App/")

RNA_Seq_M1_Peripheral = read.delim("app_data/S.v.AS.in.PBMC.Peripheral.tsv")
RNA_Seq_M2_Local = read.delim("app_data/S.v.AS.in.PBMC.Local.tsv")
RNA_Seq_M3_Plaque = read.delim("app_data/S.v.AS.in.Plaque.tsv")
RNA_Seq_M4_Plasma = read.delim("app_data/rnaseq.model4.PBMC.Peripheral.and.PBMC.Local.tsv", row.names = 1)



# fix model 4 formatting
head(RNA_Seq_M4_Plasma)

RNA_Seq_M4_Plasma_Symptomatic = RNA_Seq_M4_Plasma[,1:3]
colnames(RNA_Seq_M4_Plasma_Symptomatic) = c("Gene.ID","log2FoldChange", "pvalue")
RNA_Seq_M4_Plasma_Symptomatic$padj = p.adjust(RNA_Seq_M4_Plasma_Symptomatic$pvalue, method = "BH")

RNA_Seq_M4_Plasma_Location = RNA_Seq_M4_Plasma[,c(1,4,5)]
colnames(RNA_Seq_M4_Plasma_Location) = c("Gene.ID","log2FoldChange", "pvalue")
RNA_Seq_M4_Plasma_Location$padj = p.adjust(RNA_Seq_M4_Plasma_Location$pvalue, method = "BH")


RNASeqModels = lapply(list(Symptomatic_Peripheral = RNA_Seq_M1_Peripheral, Symptomatic_Local = RNA_Seq_M2_Local, Symptomatic_Plaque = RNA_Seq_M3_Plaque, Symptomatic_Plasma = RNA_Seq_M4_Plasma_Symptomatic, Location_Plasma = RNA_Seq_M4_Plasma_Location), function(DS){
  newDF = DS[,c("log2FoldChange","pvalue","padj")]
  rownames(newDF) = DS$Gene.ID
  colnames(newDF) = c("Log2FC", "p-value", "adjusted p-value")
  return(newDF)
})

mergedRNASeq = read.delim("app_data/all_RNA_Seq.txt")
colnames(mergedRNASeq) = gsub("_|\\."," ", colnames(mergedRNASeq))

# exclude the non-adjusted pvalues for display reasons
mergedRNASeq = mergedRNASeq[-grep("^p value.*",colnames(mergedRNASeq))]
colnames(mergedRNASeq) = gsub("p value", "p-value", colnames(mergedRNASeq))

#######################

require(fgsea)
require(msigdbr)

## get pathway information from https://www.gsea-msigdb.org/gsea/msigdb/
genesets = msigdbr()

## we will use H and C2 (CP:KEGG)
splitGS = split(genesets, paste(genesets$gs_cat, genesets$gs_subcat))

KEGG_GS = split(splitGS$`C2 CP:KEGG`$gene_symbol, splitGS$`C2 CP:KEGG`$gs_name)
HALLMARK_GS = split(splitGS$`H `$gene_symbol, splitGS$`H `$gs_name)
combinedGS = c(KEGG_GS, HALLMARK_GS)
GOBP = split(splitGS$`C5 GO:BP`$gene_symbol, splitGS$`C5 GO:BP`$gs_name)
# 20-200size gene sets
GOBP_filt = GOBP[unlist(lapply(GOBP, length)) > 19 & unlist(lapply(GOBP, length)) < 201]
length(GOBP_filt)

## use the model coefficients for ranks - it combined the pvalue and FC
geneFCs = mergedRNASeq[,grep("Log2FC", colnames(mergedRNASeq))]

## run fast GSEA
fgseaRes = lapply(geneFCs, function(GR){
  names(GR) = rownames(geneFCs)
  GR = sort(GR)
  print(head(GR))
  
  hallmarkKEGGfgseaRes <- fgsea(combinedGS, GR, maxSize=100)
  GOBPfgseaRes <- fgsea(GOBP_filt, GR, maxSize=100)
  ## check the results
  GO.ordered.pathways.res = data.frame(hallmarkKEGGfgseaRes[order(hallmarkKEGGfgseaRes$padj)])
  HK.ordered.pathways.res = data.frame(GOBPfgseaRes[order(GOBPfgseaRes$padj)])
  GO.sig.pathways = GO.ordered.pathways.res[GO.ordered.pathways.res$padj < 0.05,]
  HK.sig.pathways = HK.ordered.pathways.res[HK.ordered.pathways.res$padj < 0.05,]
  sig.pathways = list(GOres = GO.sig.pathways,hallmarkKEGGres = HK.sig.pathways)
  return(sig.pathways)
})


## plot the results
require(ggplot2)

combinedGSForPlotting = c(combinedGS, GOBP_filt)

lapply(1:length(fgseaRes), function(I){
  modelResults = fgseaRes[[I]]
  modelGeneRanks = geneRanks[[I]]
  modelName = names(fgseaRes)[I]
  
  merged.res = do.call(rbind, modelResults)
  merged.res = merged.res[order(merged.res$padj),]
  
  print(merged.res[,c(1,3,5)])
  
  if(!nrow(merged.res)>0){return(NA)}
  
  #pdf(file = paste0("datasets/KI/Proteomics/",modelName, ".pdf"))
  
  lapply(head(merged.res$pathway, 5), function(H) {
    
    EP = plotEnrichment(combinedGSForPlotting[[H]], modelGeneRanks) + 
      labs(title= paste0(modelName, " ", H, "\n-  BH Adj. p-value = ", 
                         signif(merged.res$padj[merged.res$pathway == H], 4)) )
    EP
  })
})



# write the results to excel


list.for.export = lapply(fgseaRes, function(RES){
  merged.res = do.call(rbind, RES)
  merged.res$leadingEdge = unlist(lapply(merged.res$leadingEdge, paste, collapse = " "))
  merged.res = merged.res[order(merged.res$padj),]
})

writexl::write_xlsx(list.for.export, "~/datasets/KI/RNA-Seq/FGSEA_results.xlsx")



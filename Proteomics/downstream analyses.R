## input model results
require(tidyverse)

olink.model.results = read.delim("datasets/KI/Proteomics/Olink_models_results.txt")
olink_symbols = readxl::read_excel("datasets/KI/Proteomics/Olink symbol mapping.xlsx")

OlinkSymbols = unlist(lapply(strsplit(rownames(olink.model.results), " \\|"), '[[', 1))
mapped_gene_symbols = olink_symbols$To[match(OlinkSymbols, olink_symbols$Assay)]

olink.model.results$GeneSymbol = mapped_gene_symbols
olink.model.results$m1.pval.adj = p.adjust(olink.model.results$m1.pval, method = "BH")
olink.model.results$m2.pval.adj = p.adjust(olink.model.results$m2.pval, method = "BH")
olink.model.results$m3.pval.adj = p.adjust(olink.model.results$m3.pval, method = "BH")
olink.model.results$m4.sym.pval.adj = p.adjust(olink.model.results$m4.sym.pval, method = "BH")
olink.model.results$m4.loc.pval.adj = p.adjust(olink.model.results$m4.loc.pval, method = "BH")

plot(olink.model.results$m1.coef, -log10(olink.model.results$m1.pval.adj))
m1.sig = olink.model.results$m1.pval.adj < 0.05
m2.sig = olink.model.results$m2.pval.adj < 0.05
m3.sig = olink.model.results$m3.pval.adj < 0.05
m4.sym.sig = olink.model.results$m4.sym.pval.adj < 0.05
m4.loc.sig = olink.model.results$m4.loc.pval.adj < 0.05 & abs(olink.model.results$m4.loc.coef) > 0.5

symptomatic.sig = olink.model.results[m1.sig|m4.sym.sig,]

write.table(unique(symptomatic.sig$GeneSymbol), file = "datasets/KI/Proteomics/SigProteins_Symptomatic.txt", sep = "\t", quote = F, row.names = F, col.names = F)

location.sig = olink.model.results[m4.loc.sig,]

write.table(unique(location.sig$GeneSymbol), file = "datasets/KI/Proteomics/SigProteins_Location.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#combined
write.table(unique(c(symptomatic.sig$GeneSymbol, location.sig$GeneSymbol)), file = "datasets/KI/Proteomics/SigProteins_Combined.txt", sep = "\t", quote = F, row.names = F, col.names = F)


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
protCoefs = olink.model.results[,grep("coef", colnames(olink.model.results))]

OlinkSymbols = unlist(lapply(strsplit(rownames(olink.model.results), " \\|"), '[[', 1))
mapped_gene_symbols = olink_symbols$To[match(OlinkSymbols, olink_symbols$Assay)]

#rownames(geneRanks) = mapped_gene_symbols

geneRanks = lapply(protCoefs, function(GR){
  names(GR) = mapped_gene_symbols
  GR = sort(GR[!duplicated(names(GR))])
  return(GR)
})

## run fast GSEA
fgseaRes = lapply(geneRanks, function(GR){
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

writexl::write_xlsx(list.for.export, "datasets/KI/Proteomics/FGSEA_results.xlsx")



#####################################################################################
#####################################################################################
#####################################################################################

########### clustering

# load("datasets/KI/Proteomics/Clean_OLink_Data.RData")
load("datasets/KI/Proteomics/Clean_OLink_Data.RData")
#head(Clean_OLink_Data)
require(OlinkAnalyze)
# my_NPX_data <- read_NPX(filename = "~/datasets/KI/Proteomics/20200468_Matic_NPX PseudoID_for_analysis.xlsx")
# 
# opw = pivot_wider(data = my_NPX_data, id_cols = c("SampleID"), values_from = "NPX", names_from = c("Assay","Panel"))
# clusterMat = data.matrix(opw[,-1])
# rownames(clusterMat) = opw$SampleID
# clusterMat = clusterMat[-grep("CONTROL.*", rownames(clusterMat)),]


#########################################
## PCA Analyses
#########################################

npx_df = data.frame(readxl::read_xlsx("datasets/KI/Proteomics/20200468_Matic_NPX PseudoID_for_analysis.xlsx", 3))

#plate effects
npx_plateIDs = npx_df[,466:470]
require(gplots)
heatmap.2(na.omit(do.call(cbind, lapply(npx_plateIDs, function(X){ as.numeric(as.factor(X)) }))), scale = "none", trace = "none", main = "Plate IDs", cexCol = 0.5)

# metadata
npx_annots = npx_df[,c(1:5,466)]
rownames(npx_annots) = paste0(npx_annots$PID, "_", npx_annots$BiKE.ID.prefix)
npx_annots$Olink_outliers = rownames(npx_annots)%in%names(OL_OLs)

# actual data
npx_df2 = npx_df[,-c(1:5,466:ncol(npx_df))]
npx_mat = data.matrix(do.call(cbind,lapply(npx_df2, as.numeric)))
rownames(npx_mat) = rownames(npx_annots)

non_samples = grep("control|na", rownames(npx_annots), ignore.case = T) 

pca_mat = t(scale(npx_mat[-non_samples,]))
pca_mat[is.na(pca_mat)] = 0

##############
PCA_res = prcomp(pca_mat)

PCAcols = c("red","blue")[grepl("EP",rownames(PCA_res[[2]]))+1]
plot(PCA_res$rotation[,c(1,2)], col = PCAcols, pch = 16, main = "PCA plot")
legend('topright', legend = c("STP - local","EP - peripheral"), fill = c("red","blue"))

plot(PCA_res)
plot(PCA_res[[1]][1:100], main = "variance captured by each PC")
# use first 13 PCs
numpcs = 13
abline(v = numpcs, col = "grey")
PCmat = PCA_res[[2]][,1:numpcs]

### peripheral only
PCmat = na.omit(PCmat[sample_annots[rownames(PCmat),"BiKE.ID.prefix"] == "EP",])

### UMAP

sample.sex = sample_annots[rownames(PCmat),"Sex"]
sample.type = sample_annots[rownames(PCmat),"BiKE.ID.prefix"]
sample.plate = sample_annots[rownames(PCmat),"Plate.ID...466"]
sample.outliers = sample_annots[rownames(PCmat),"Olink_outliers"]
sample.outliers = sample_annots[rownames(PCmat),"Symptomatic"]
require(M3C)
umap(t(PCmat), colvec=c("red","blue","green"), labels = sample.type, dotsize = 2)
umap(t(PCmat), colvec=c("red","blue","green"), labels = sample.sex, dotsize = 2)
umap(t(PCmat), colvec=c("red","blue","green"), labels = sample.plate, dotsize = 2)
#umap(t(PCmat), colvec=c("red","blue","green"), labels = as.numeric(sample.outliers), dotsize = 2)
umap(t(PCmat), colvec=c("red","blue","green"), labels = sample.outliers, dotsize = 2)

#umap(t(PCmat), colvec=rainbow(length(unique(all.sample.annots$Som.cluster))), labels = as.factor(all.sample.annots$Som.cluster))

m3c_clusters <- M3C(data.frame(t(PCmat)), method=2)
umap(t(PCmat), colvec=rainbow(length(unique(m3c_clusters$assignments))), labels = as.character(m3c_clusters$assignments), dotsize = 2)

heatmap(m3c_clusters$realdataresults[[2]]$consensus_matrix, scale = "none", Rowv = F, Colv = F)

PCmat_sample.annots = sample_annots[rownames(PCmat),]
PCmat_sample.annots$M3Cclusters = m3c_clusters$assignments

###########################################
###########################################


##################################################
## compare clusters with limma

require(limma)

annots = PCmat_sample.annots
annots$Cluster = as.factor(annots$M3Cclusters)

#########

#gene.annots = read.delim("datasets/ENS_Uniprot_ID.txt")

work.eset <- npx_mat[-non_samples,]
work.eset <- na.omit(work.eset[sample_annots[rownames(work.eset),"BiKE.ID.prefix"] == "EP",])

OlinkSymbols = unlist(lapply(strsplit(colnames(work.eset), "\\.\\.\\."), '[[', 1))
mapped_gene_symbols = olink_symbols$To[match(OlinkSymbols, gsub("[ -/]","\\.",olink_symbols$Assay))]
colnames(work.eset) = mapped_gene_symbols

work.annots <- annots[rownames(work.eset),]

design <- model.matrix(~ Cluster + Sex + 0, work.annots)
colnames(design) <- make.names(colnames(design))
#design

fit <- lmFit(t(work.eset), design)

list.of.contrasts <- list(
  c1_vs_All = "Cluster1-(Cluster2+Cluster3+Cluster4+Cluster5+Cluster6+Cluster7)/6",
  c2_vs_All = "Cluster2-(Cluster1+Cluster3+Cluster4+Cluster5+Cluster6+Cluster7)/6",
  c3_vs_All = "Cluster3-(Cluster1+Cluster2+Cluster4+Cluster5+Cluster6+Cluster7)/6",
  c2and3_vs_All = "(Cluster3+Cluster2)/2-(Cluster1+Cluster4+Cluster5+Cluster6+Cluster7)/5",
  c2vs3 = "Cluster2-Cluster3",
  c4_vs_All = "Cluster4-(Cluster1+Cluster3+Cluster2+Cluster5+Cluster6+Cluster7)/6",
  c5_vs_All = "Cluster5-(Cluster1+Cluster3+Cluster4+Cluster2+Cluster6+Cluster7)/6",
  c6_vs_All = "Cluster6-(Cluster1+Cluster3+Cluster4+Cluster5+Cluster2+Cluster7)/6",
  c7_vs_All = "Cluster7-(Cluster1+Cluster3+Cluster4+Cluster5+Cluster6+Cluster2)/6"
)

all_results <- lapply(list.of.contrasts, function(X){
  
  cont.matrix <- makeContrasts(contrasts=X,levels=design) ##
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  
  tT <- topTable(fit2, adjust.method = "BH", sort.by="p", number=Inf)
  #tT$Gene <- trimws(gene.annots[match(rownames(tT), gene.annots$UniProtKB.Swiss.Prot.ID),'HGNC.symbol'])
  
  volcanoplot(fit = fit2, highlight = 50, main = names(which(list.of.contrasts == X)), names = rownames(fit2), xlim = c(min(tT$logFC) - 0.1, max(tT$logFC) + 0.1)) 
  
  return(tT)
  
})



save(all_results, work.eset, work.annots, file="/novo/users/drdj/datasets/KI/Proteomics/clusterComparisonResults.Rdata")


#########


#sanity check boxplots
GOIs=c("CCL21","SPON1","ANGPT1","EGF")
GOIs=c("FGF2","HGF","TNFSF12","PTN")
GOIs=head(all_results[[3]]$ID, 4)
GOIs=head(all_results[[4]]$ID, 4)
GOIs=head(all_results[[5]]$ID, 4)
GOIs=head(all_results[[6]]$ID, 4)
GOIs=head(all_results[[7]]$ID, 4)
GOIs=head(all_results[[8]]$ID, 4)
GOIs=head(all_results[[9]]$ID, 4)

par(mfrow = c(2,2))
par(mar = c(5,5,1,1))
lapply(GOIs, function(GOI){
  print(GOI)
  
  probeID = na.omit(rownames(all_results[[1]])[all_results[[1]]$ID == GOI])
  splittedvals = split(unlist(t(work.eset)[GOI, ]), work.annots$Cluster)
  boxplot(splittedvals, col = "grey", main = GOI, las = 2, xlab = "cluster", ylab = "prot level")
  #points(x = c(rep(1, length(splittedvals$))+rnorm(length(splittedvals$A), sd = 0.1), rep(2, length(splittedvals$B)))+rnorm(length(splittedvals$B), sd = 0.1), y = unlist(splittedvals[c("A","B")]), pch = 16)
})

#######################################
#######################################

## use the model t stat
protCoefs = lapply(all_results, function(X){
  tstats = X$t
  names(tstats) = X$ID
  return(tstats[!duplicated(names(tstats))])
})


## run fast GSEA
fgseaRes = lapply(protCoefs, function(GR){
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

writexl::write_xlsx(list.for.export, "datasets/KI/Proteomics/FGSEA_results.xlsx")



# 
# my_palette <- colorRampPalette(c("darkblue","blue", "white", "red","firebrick"))(n = 41)
# 
# source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
# 
# hcl = hclust(dist(clusterMat))
# hclusters = cutree(hcl, k = 10)
# 
# rowColors = data.frame(sampleType = c("blue","black")[as.numeric(sample_annots[rownames(clusterMat),"BiKE.ID.prefix"])], 
#                        sex = c("pink","dodgerblue")[as.numeric(sample_annots[rownames(clusterMat),"Sex"])],
#                        symptomatic = c("forestgreen","firebrick1")[as.numeric(sample_annots[rownames(clusterMat),"Symptomatic"])],
#                        cluster = rainbow(length(unique(hclusters)))[hclusters])
# 
# heatmap.3(clusterMat, trace = "none", scale = "col", col = my_palette,
#           RowSideColors = t(as.matrix(rowColors)),
#           RowSideColorsSize=4, margins = c(12,5))



# 
# ## run GAGE
# gageRes = lapply(geneRanks, function(GR){
#   #gageres = gage(geneRanks, gsets = combinedGS, ref = NULL, samp = NULL)
#   hallmarkKEGGgageRes <- gage(GR, gsets = combinedGS, ref = NULL, samp = NULL)
#   GOBPgageRes <- gage(GR, gsets = GOBP, ref = NULL, samp = NULL)
#   ## check the results
#   
#   GOMERGED = merge(GOBPgageRes$greater, GOBPgageRes$less, suffixes = c("GREATER","LESS"))
#   HKMERGED = merge(hallmarkKEGGgageRes$greater, hallmarkKEGGgageRes$less, suffixes = c("GREATER","LESS"))
#   
#   GO.ordered.pathways.res = data.frame(hallmarkKEGGgageRes[order(hallmarkKEGGgageRes$padj)])
#   HK.ordered.pathways.res = data.frame(GOBPgageRes[order(GOBPgageRes$padj)])
#   
#   GO.sig.pathways = GO.ordered.pathways.res[GO.ordered.pathways.res$padj < 0.05,]
#   HK.sig.pathways = HK.ordered.pathways.res[HK.ordered.pathways.res$padj < 0.05,]
#   
#   sig.pathways = list(GOres = GO.sig.pathways,hallmarkKEGGres = HK.sig.pathways)
#   return(sig.pathways)
# })## run GAGE
# 




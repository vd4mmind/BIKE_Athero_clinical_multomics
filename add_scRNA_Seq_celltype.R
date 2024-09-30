#################################################################
#################################################################
#################################################################
require(Seurat)
## WIRKA et al *patient) Jias path read.delim("~/../hjiw/data/TTAorta/public/PMID_31359001/")
## patient
load("~/../hjiw/data/TTAorta/public/PMID_31359001/PMID_31359001_patient.RData")
DimPlot(patient, label = T)

Wirka_clusters_list = list(SMC_clusters = c("modSMC", "SMC1", "SMC2"),
                           SMC_plus_fib_clusters = c("modSMC", "SMC1", "SMC2","Fibro","Pericyte"),
                           EC_clusters = c("EC1","EC2","EC3"),
                           MP_clusters = c("Macro","Mono"),
                           #TCell_clusters = c("NK","TCell1","TCell2"),
                           #Plasma_clusters = c("Plasma1","Plasma2"),
                           FB_clusters = c("Fibro"))

## Alsaigh et al
## seurat.obj
load("datasets/AS_scRNAseq_paired_plaques/seurat.obj.RData")
DimPlot(seurat.obj, label = T)


# need to label the SMC
DotPlot(seurat.obj, features = c("ACTA2","MYH11","TAGLN"))
Alsaigh_SMC_clusters = c(5,9,11,13,15,20)
renamed.seurat.obj = RenameIdents(seurat.obj, '5' = "SMC", '9' = "SMC", '11'= "SMC", '13'= "SMC",'15'= "SMC",'20'= "SMC")
DimPlot(renamed.seurat.obj, label = T)

# and fibroblasts
DotPlot(renamed.seurat.obj, features = c("COL1A1","COL1A2","DCN","FBLN1","FBN1","OGN"))
FeaturePlot(renamed.seurat.obj, features = c("COL1A1","COL1A2","DCN","FBLN1","FBN1","OGN"))
renamed.seurat.obj = RenameIdents(renamed.seurat.obj, '14' = "FB")
DimPlot(renamed.seurat.obj, label = T)


# and ECs
DotPlot(renamed.seurat.obj, features = c("VWF","PECAM1"))
FeaturePlot(renamed.seurat.obj, features =  c("VWF","PECAM1"))
renamed.seurat.obj = RenameIdents(renamed.seurat.obj, '12' = "EC", '8' = "EC", '7'= "EC")
DimPlot(renamed.seurat.obj, label = T)

# and MP/MCs
DotPlot(renamed.seurat.obj, features = c("IL1B","CD14","CD86","FCGR1A"))
FeaturePlot(renamed.seurat.obj, features =  c("IL1B","CD14","CD86","FCGR1A"))
renamed.seurat.obj = RenameIdents(renamed.seurat.obj, '2' = "MP", '3' = "MP", '4'= "MP", '16'= "MP")
DimPlot(renamed.seurat.obj, label = T)

FeaturePlot(renamed.seurat.obj, features = c("IL1B","VWF","ACTA2","FBLN1"))

# and T cells
FeaturePlot(renamed.seurat.obj, features = c("NKG7","PTPRC","CCL5","GZMB"), order = T, pt.size = 0.5, label = F)
DotPlot(renamed.seurat.obj, features = c("NKG7","PTPRC","CCL5","GZMB"))
renamed.seurat.obj = RenameIdents(renamed.seurat.obj, '0' = "T-Cells", '1' = "T-Cells", '6'= "T-Cells")
DimPlot(renamed.seurat.obj, label = T)



#######
Alsaigh_clusters_list = list(SMC_clusters = c("SMC"),
                           SMC_plus_fib_clusters = c("SMC","FB"),
                           EC_clusters = c("EC"),
                           MP_clusters = c("MP"),
                           FB_clusters = c("FB"))
#########


## BiKE genes
## 
TARLIST = read.delim("datasets/KI/Single_Omics_Targets_Enriched.txt")
all.genes = sort(unique(c(TARLIST$gene_name)))
length(all.genes)
# 268


#################################################################
#################################################################
#################################################################
#wirka et al ##############################
gc = DotPlot(patient, features = all.genes)
# 14 were not found in the dataset

#alsaigh et al ############################
gc = DotPlot(renamed.seurat.obj, features = all.genes)
# 2 were not found in the dataset

gcsplit = split(gc$data, gc$data$features.plot)

gcDF = reshape(data = gc$data, timevar = "id", idvar = "features.plot", v.names = c("avg.exp", "pct.exp", "avg.exp.scaled"), direction = "wide")

#write.table(file = "datasets/ODIN/Wirka_dotplot_table.txt", gcDF, quote = F, row.names = F, col.names = T, sep = "\t")
#write.table(file = "datasets/ODIN/Alsaigh_dotplot_table.txt", gcDF, quote = F, row.names = F, col.names = T, sep = "\t")


#################################################################
## and non expressed genes

max.pct.exps = apply((gcDF[,grepl("pct.exp", colnames(gcDF))]), 1, max)
hist(max.pct.exps, breaks = 100, main = "max % exp")
abline(v = 5, col = "red")

#write.table(file = "datasets/ODIN/Wirka_max_perc_cells.txt", data.frame(max.pct.exps), quote = F, row.names = T, col.names = F, sep = "\t")
#write.table(file = "datasets/ODIN/Alsaigh__max_perc_cells.txt", data.frame(max.pct.exps), quote = F, row.names = T, col.names = F, sep = "\t")

#################################################################


#################################################################
#################################################################
#################################################################

## deg approach ##

load("~/StemCells/ScParams.RData")

find_DEG <- function(sel_ident1, sel_ident2, seurObj){
  
  cluster_markers = FindMarkers(object = seurObj,
                                ident.1 = sel_ident1,
                                ident.2 = sel_ident2,
                                min.pct = 0.25)
  
  require(data.table)
  cluster_markers = data.table(cluster_markers, keep.rownames = T)
  cluster_markers = merge.data.frame(cluster_markers,
                                     #seurObj$params$gtf[, c('gene_id','gene_name')], ### changed this line for this dataset
                                     ScParams[, c('gene_id','gene_name')], ### changed this line for this dataset
                                     by.x ='rn', by.y = 'gene_id', all.x =T)
  colnames(cluster_markers)[1] = 'ens_id'
  
  #####
  
  # calculate mean exBCession for selected idents
  all_sel_idents = c(sel_ident1, sel_ident2)
  subsets_seurat_idents = subset(seurObj, idents = all_sel_idents)
  subsets_rna_mat = as.matrix(GetAssayData(subsets_seurat_idents, assay = 'RNA', slot = 'data'))
  idents_rna_exp = data.table(ens_id = rownames(subsets_rna_mat),
                              mean_exp = apply(subsets_rna_mat, 1, mean))
  # add mean exp value
  cluster_markers = merge(cluster_markers, idents_rna_exp, by='ens_id')
  # save to reactive value
  #all_data$cluster_markers = cluster_markers
  
  #print(head(cluster_markers))
  
  require(dplyr)
  showed_markers = as.tbl(cluster_markers) %>%
    dplyr::select(gene_name, mean_exp, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ens_id) %>%
    mutate_at(2:7, signif, 3)
  
  colnames(showed_markers) = c('UNKNOWN', 'Mean expr. value', 'average logFC', 'P_value', 'Adjusted P_value', 'Cluster1 cell frac expr.', 'Cluster2 cell frac expr.', 'Gene Symbol')
  showed_markers$percDiff = showed_markers$`Cluster1 cell frac expr.` - showed_markers$`Cluster2 cell frac expr.`
  return(showed_markers)
}


#################################################################
#################################################################
#################################################################

Wirka_cluster_DEG = lapply(Wirka_clusters_list, function(activeClusters){
  
  print(activeClusters)
  all.clusters = as.character(unique(patient@active.ident))
  
  activemarkers = find_DEG(sel_ident1 = activeClusters, sel_ident2 = setdiff(all.clusters, activeClusters), seurObj = patient)
  #print(head(data.frame(activemarkers)))
  
  upPercCutoff = 0.2
  upFcCutoff = 0.5
  pvalcutoff = 0.05
  
  # {
  #   plot(activemarkers$`average logFC`, activemarkers$percDiff, pch = 16, main = paste0("Wirka ", paste(activeClusters),"\nSignature definition"))
  #   
  #   abline(v = 0, h = 0, lty = 1, col = "grey")
  #   abline(v = upFcCutoff, h = upPercCutoff, lty = 2, col = "red")
  # }
  
  activemarkers$singleMeasure = activemarkers$`average logFC` * abs(activemarkers$percDiff)
  activemarkers$include = activemarkers$`Adjusted P_value` < pvalcutoff & activemarkers$percDiff > upPercCutoff & activemarkers$`average logFC`> upFcCutoff 
  
  return(activemarkers)
  
})

#################################################################

Alsaigh_cluster_DEG = lapply(Alsaigh_clusters_list, function(activeClusters){
  
  print(activeClusters)
  all.clusters = as.character(unique(renamed.seurat.obj@active.ident))
  
  activemarkers = find_DEG(sel_ident1 = activeClusters, sel_ident2 = setdiff(all.clusters, activeClusters), seurObj = renamed.seurat.obj)
  
  upPercCutoff = 0.2
  upFcCutoff = 0.5
  pvalcutoff = 0.05
  
  # {
  #   plot(activemarkers$`average logFC`, activemarkers$percDiff, pch = 16, main = paste0("Alsaigh ", paste(activeClusters),"\nSignature definition"))
  #   
  #   abline(v = 0, h = 0, lty = 1, col = "grey")
  #   abline(v = upFcCutoff, h = upPercCutoff, lty = 2, col = "red")
  # }
  
  activemarkers$singleMeasure = activemarkers$`average logFC` * abs(activemarkers$percDiff)
  activemarkers$include = activemarkers$`Adjusted P_value` < pvalcutoff & activemarkers$percDiff > upPercCutoff & activemarkers$`average logFC`> upFcCutoff 
  
  return(activemarkers)
  
})

#################################################################

Alsaigh_clean_res = lapply(names(Alsaigh_cluster_DEG), function(nDEG){

  DEG = Alsaigh_cluster_DEG[[nDEG]]
  subDEG = DEG[, c("Gene Symbol", "include", "singleMeasure")]
  colnames(subDEG) = c("Gene", paste0(gsub("_clusters","",nDEG), c("_marker", "_metric")))
  return(subDEG)
  
  })

names(Alsaigh_clean_res) = names(Alsaigh_cluster_DEG)

Alsaigh_merged = plyr::join_all(Alsaigh_clean_res, type = "full")

#FeaturePlot(renamed.seurat.obj, features =  c("ACTA2","ACTR2","ADIRF","A2M"))

colnames(Alsaigh_merged) = c("Gene", paste0("Alsaigh_", colnames(Alsaigh_merged[-1])))

#################################################################

Wirka_clean_res = lapply(names(Wirka_cluster_DEG), function(nDEG){
  
  DEG = Wirka_cluster_DEG[[nDEG]]
  subDEG = DEG[, c("Gene Symbol", "include", "singleMeasure")]
  colnames(subDEG) = c("Gene", paste0(gsub("_clusters","",nDEG), c("_marker", "_metric")))
  return(subDEG)
  
})

names(Wirka_clean_res) = names(Wirka_cluster_DEG)

Wirka_merged = plyr::join_all(Wirka_clean_res, type = "full")

#FeaturePlot(patient, features =  c("AIF1","C1QC","TYROBP","FCER1G"))

colnames(Wirka_merged) = c("Gene", paste0("Wirka_", colnames(Wirka_merged[-1])))

#################################################################

both.merged = merge(Wirka_merged, Alsaigh_merged, all = T)

######

write.table(both.merged, file = "datasets/KI/Wirka_Alsaigh_cellTypes.txt", quote = F, row.names = F, col.names = T, sep = "\t")



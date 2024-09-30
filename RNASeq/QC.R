# KI RNA-Seq 

#install.packages("aws.s3")
#install.packages("SASxport")

library("SASxport")
library("aws.s3")

# #from NNEDL credentials profile
# [dl_user]
# aws_access_key_id = ASIATIJCHR3KJUKEICQX
# aws_secret_access_key = UWAuqHwDW6gK1gQOP+boGuDTXVlaeyMTmXmT1ijx
# aws_session_token = FwoGZXIvYXdzELr//////////wEaDECfzz2M7Jeq2xYa3CL0A/JYBFmkw0XPL0tE+lSc7NKQxwECksY5wc3FiBXbCxzuJydJkVM+7HxFtqDxGSDY6cK7Ebaqhp9XLTreBh+6oBh80UKPMGQ6M1oJaLveEdzJKXaWZCF80o2mVqGkQsk+vbuVoq8Zwn2xkSwfkW06M99NrhpY9JchZQnSIBEEI5fDTmtPLN8OK5TJWu9qkctWCSNS2Ml2zVXSfqAYJ0yAA3L4gXBEVaE6bprce7xfytZ5xecYGodX/GOuTywhpkbYSkLDkcJ+RQdU7CuA3S6UvIXdad70mxCBUrkzzMsLLWFLubtzoqQep5RFoimaLv3qpg77z5wA5edEc4PsNVnp2pwVWUsdlOku7xudN0EX5dxu+o7BsHpVgJ9OgBGMU8kVHYVdDacP7ECEhD9vkGt3HGKZOe/ELl/YvH+0IARtq6lWqQD0c5aec7FsPXnyztGyzPXVmPA7dGTE9f6YR7oHxRSJ7fl5vDtMPFnF2HGtQoN68KIvi+7u7um35yBl7247dM414OPA/ZNOqZywmRay5V+FZZAYTV+3LhIFC/6XutZ68TL8Bp91unFUwOIQ0FDqT/HSpFPc0XvZcIu8pNScFRnMd+0QA0OauHQnGQ1nRKlWpwVZWcGBiDHUH6ecZHJpSroWxAOd3p8EOfqWh1WfmVlSIbH0KIqnsYYGMpQB7yyqfy3LMZtOxJQbolOWa2kqM3mxBqks1b0t7tOvAHd2wIYRRzKim+Kat0nS26/wRUxlzQM/YpUGk+rga07f5/6pCKS52Zjn7hBlclzucRSzSZ18kFeC6KIjwYWus3wjZaPnMiI4ZolP/LWnrGZg59Z0QhhQc0aSDRJaxksGY0w7J7d6e+a/7oUM6ZtHAuoXGSOJaA==

Sys.setenv(
  "AWS_ACCESS_KEY_ID"  = "ASIATIJCHR3KJUKEICQX",
  "AWS_SECRET_ACCESS_KEY" = "UWAuqHwDW6gK1gQOP+boGuDTXVlaeyMTmXmT1ijx",
  "AWS_SESSION_TOKEN" = "FwoGZXIvYXdzELr//////////wEaDECfzz2M7Jeq2xYa3CL0A/JYBFmkw0XPL0tE+lSc7NKQxwECksY5wc3FiBXbCxzuJydJkVM+7HxFtqDxGSDY6cK7Ebaqhp9XLTreBh+6oBh80UKPMGQ6M1oJaLveEdzJKXaWZCF80o2mVqGkQsk+vbuVoq8Zwn2xkSwfkW06M99NrhpY9JchZQnSIBEEI5fDTmtPLN8OK5TJWu9qkctWCSNS2Ml2zVXSfqAYJ0yAA3L4gXBEVaE6bprce7xfytZ5xecYGodX/GOuTywhpkbYSkLDkcJ+RQdU7CuA3S6UvIXdad70mxCBUrkzzMsLLWFLubtzoqQep5RFoimaLv3qpg77z5wA5edEc4PsNVnp2pwVWUsdlOku7xudN0EX5dxu+o7BsHpVgJ9OgBGMU8kVHYVdDacP7ECEhD9vkGt3HGKZOe/ELl/YvH+0IARtq6lWqQD0c5aec7FsPXnyztGyzPXVmPA7dGTE9f6YR7oHxRSJ7fl5vDtMPFnF2HGtQoN68KIvi+7u7um35yBl7247dM414OPA/ZNOqZywmRay5V+FZZAYTV+3LhIFC/6XutZ68TL8Bp91unFUwOIQ0FDqT/HSpFPc0XvZcIu8pNScFRnMd+0QA0OauHQnGQ1nRKlWpwVZWcGBiDHUH6ecZHJpSroWxAOd3p8EOfqWh1WfmVlSIbH0KIqnsYYGMpQB7yyqfy3LMZtOxJQbolOWa2kqM3mxBqks1b0t7tOvAHd2wIYRRzKim+Kat0nS26/wRUxlzQM/YpUGk+rga07f5/6pCKS52Zjn7hBlclzucRSzSZ18kFeC6KIjwYWus3wjZaPnMiI4ZolP/LWnrGZg59Z0QhhQc0aSDRJaxksGY0w7J7d6e+a/7oUM6ZtHAuoXGSOJaA==",
  "AWS_DEFAULT_REGION" = "eu-central-1"
)


# bike <- get_bucket(bucket = 'nnedl-core-prd-eu-central-1-source', prefix='atherobike', max = Inf)
# 
# length(bike)
# 
# # Print the files included in the bucket
# AllBikeFiles = unlist(lapply(1:length(bike), function(i){
#   bike[i]$Contents$Key
# }))
# 
# 
# str(bike)

main.counts.file = "atherobike/KI-cohort/PseudoID.Export/RNAseq/200819_A00605_0144_AH5YKTDSXY/results/featureCounts/gene_counts"
main.counts.folder = "atherobike/KI-cohort/PseudoID.Export/RNAseq/200819_A00605_0144_AH5YKTDSXY/results/"

# #testFile = s3read_using(FUN = read.csv, bucket='nnedl-core-prd-eu-central-1-source', object="grohfpef/300436_Iziah_Sama_P1/rseqc/2460003/rseqc.txt" )
# table(unlist(lapply(strsplit(AllBikeFiles[grep("RNA",AllBikeFiles)], "_"), "[[", 1)))
# table(AllBikeFiles[grep(main.counts.file, AllBikeFiles)])
# table(unlist(lapply(strsplit(AllBikeFiles[grep(main.counts.file,AllBikeFiles)], "/"), "[[", 9)))
# head(table(unlist(lapply(strsplit(AllBikeFiles[grep(main.counts.folder,AllBikeFiles)], "_"), "[[", 1))))
sum(grepl(main.counts.folder, AllBikeFiles))

counts.files = AllBikeFiles[grep(main.counts.file,AllBikeFiles)]
sampleIDs = gsub(".*TA-2453-(.*)Aligned.*","\\1",counts.files)

counts.list = lapply(counts.files, function(COUNTFILE){
  NNEDLfile = s3read_using(FUN = read.delim, bucket='nnedl-core-prd-eu-central-1-source', object= COUNTFILE, skip = 1)
  return(NNEDLfile[,c(1,6,8)])
})
names(counts.list) = sampleIDs

counts.matrix = do.call(cbind, lapply(counts.list, "[[", 3))
rownames(counts.matrix) = counts.list[[1]]$Geneid

write.table(counts.matrix, file = "datasets/KI/RNA-Seq/counts_matrix.txt", sep = "\t", quote = F)

checksums.list = lapply(AllBikeFiles[grep("checksum",AllBikeFiles)], function(CHECKFILE){
  NNEDLfile = s3read_using(FUN = read.delim, bucket='nnedl-core-prd-eu-central-1-source', object= CHECKFILE, header = F)
  print(head(NNEDLfile))
  return(NNEDLfile)
})
names(checksums.list) = AllBikeFiles[grep("checksum",AllBikeFiles)]

lapply(checksums.list, dim)

head(checksums.list$`atherobike/KI-cohort/PseudoID.Export/RNAseq/checksums.md5`)

pontusChecksum = read.delim("datasets/KI/RNA-Seq/checksums.md5", header = F)

#####################################################
# ---------------------------------
#   -what are labels: ST (is this local PBMCs?), CPT (is this peripheral PBMCs?) and PL (is this plaque?) 
#   ST=Local PBMC(RNA from cells from the local plasma taken in EDTA tubes) CPT=periperal PBMC(RNA from cells from the periperal plasma taken in CPT tubes) PL=RNA from plaque
# 
# -which samples have no prefix?
#   No prefix=RNA from cells from the periperal plasma taken in EDTA tubes
# 
# -as you can see from the attached slide it seems that individual 846 appears with several labels: 846, 846E and ST846. what is what?
#   846 chould be the RNA from cells from periperal EDTA plasma but becouse I found a mistake by mixing up tubes I was not completely shore that it was right I put a new sample there that was right  and named it 846E(846 should be excl).
# 846E= RNA from cells from periperal EDTA plasma ST846= RNA from cells from local EDTA plasma
# 
# -similar for individual 445, we have samples labelled PL 445 and PL445B, what is what?
#   PL445 and 445B is the same dubbel pipetted by mistake (PL=RNA from plaque)
# 
# -for individual 1011 there are many labels, what is that?
#   PL1011 A-Q is the same sample  and we put it there to have it as a controll on each plate becouse we send in the plates in three different batches and they run them on different times.
# ---------------------------------
#   
#####################################################
#####################################################

# 
# phenos = read.delim("datasets/KI/RNA-Seq/BiKE phenotypic qualifiers RNAseq.txt")
# 
# sample.annots = data.frame(row.names = colnames(counts.matrix), sampleID = colnames(counts.matrix), sampleType = unlist(lapply(strsplit(colnames(counts.matrix), "[0-9]"), "[[", 1)), PID = gsub("[A-Z]+","",colnames(counts.matrix)))
# sample.annots$sampleType[sample.annots$sampleType == ""] = "CPT-EDTA"
# sample.annots = merge(sample.annots, phenos, by.x = 'PID', by.y = 'BiKE.ID')
# rownames(sample.annots) = sample.annots$sampleID
# 
# write.table(sample.annots, file = "datasets/KI/RNA-Seq/sample_annots.txt", sep = "\t", quote = F)


######################################################
######################################################

counts.matrix = read.delim("datasets/KI/RNA-Seq/counts_matrix.txt", check.names = F)
sample.annots = read.delim("datasets/KI/RNA-Seq/sample_annots.txt", header = T, row.names = 1)
sample.annots = na.omit(sample.annots)

counts.norm = edgeR::cpm(counts.matrix, log=T)

########################################################
# 
# ██████╗░░█████╗░░█████╗░
# ██╔══██╗██╔══██╗██╔══██╗
# ██████╔╝██║░░╚═╝███████║
# ██╔═══╝░██║░░██╗██╔══██║
# ██║░░░░░╚█████╔╝██║░░██║
# ╚═╝░░░░░░╚════╝░╚═╝░░╚═╝
#
########################################################

genes.to.test = rownames(counts.norm)[rowSums(counts.norm)>100]
rowvars = apply(counts.norm, 1, var)
mostvar = names(head(rev(sort(rowvars)), 2000))
## pca

pca.mat = counts.norm[intersect(mostvar, genes.to.test), ]

##########


PCA_res = prcomp(t(scale(t(pca.mat))))
plot(PCA_res)

PCAcols = rainbow(6)[as.factor(sample.annots$sampleType)]
names(PCAcols) = sample.annots$sampleID
  
plot(PCA_res$rotation[,c(1,2)], col = PCAcols[ rownames(PCA_res$rotation[,c(1,2)]) ], pch = 16, main = "PCA plot - scaled")
legend('topright', legend = levels(as.factor(sample.annots$sampleType)), fill = rainbow(6))


###########

## PCA correlation heatmap
to.exclude = !rownames(PCA_res[[2]]) %in% rownames(sample.annots)

matchedAnnots = data.matrix(sample.annots[rownames(PCA_res[[2]])[!to.exclude],])
#matchedAnnots[unlist(apply(matchedAnnots, 1, function(X){return(sum(is.na(X)))}))>0,]

pcaCors = cor(PCA_res$rotation[,1:20][!to.exclude,], matchedAnnots)
require(gplots)
heatmap.2(pcaCors, trace = "none", scale = "none", col = colorRampPalette(c("blue", "white", "red"))(100), Colv = F, Rowv = F, margins = c(10,10), main = "correlation of PCs with metadata")

barplot(cor(PCA_res$rotation[,1:20][!to.exclude,], matchedAnnots[,"Symtomatic"])[,1], main = "correlation with Symptomatic status", las = 2)

#############

PCAcols = rainbow(6)[as.factor(sample.annots$Sex)]
names(PCAcols) = sample.annots$sampleID

plot(PCA_res$rotation[,c(3,5)], col = PCAcols[ rownames(PCA_res$rotation[,c(3,5)]) ], pch = 16, main = "PCA plot - scaled")
legend('topright', legend = levels(as.factor(sample.annots$Sex)), fill = rainbow(6))


#############
# 
# █▀ █░█ █▄▄ █▀ █▀▀ ▀█▀   █▀█ █▀▀ ▄▀█
# ▄█ █▄█ █▄█ ▄█ ██▄ ░█░   █▀▀ █▄▄ █▀█
# ## subset males only

males.only.edta.plasma = rownames(sample.annots)[sample.annots$Sex == "M" & sample.annots$sampleType %in% c("CPT-EDTA","ST")]
sub.annots = sample.annots[males.only.edta.plasma,]

sub.pca.mat = counts.norm[intersect(mostvar, genes.to.test), males.only.edta.plasma]

PCA_res = prcomp(t(scale(t(sub.pca.mat))))
plot(PCA_res)

PCAcols = rainbow(6)[as.factor(sub.annots$sampleType)]
names(PCAcols) = sub.annots$sampleID

plot(PCA_res$rotation[,c(1,2)], col = PCAcols[ rownames(PCA_res$rotation[,c(1,2)]) ], pch = 16, main = "PCA plot - scaled - men EDTA plasma only")
legend('topleft', legend = levels(as.factor(sub.annots$sampleType)), fill = rainbow(6))

###

to.exclude = !rownames(PCA_res[[2]]) %in% rownames(sub.annots)

matchedAnnots = data.matrix(sub.annots[rownames(PCA_res[[2]])[!to.exclude],])
#matchedAnnots[unlist(apply(matchedAnnots, 1, function(X){return(sum(is.na(X)))}))>0,]

pcaCors = cor(PCA_res$rotation[,1:20][!to.exclude,], matchedAnnots)
require(gplots)
heatmap.2(pcaCors, trace = "none", scale = "none", col = colorRampPalette(c("blue", "white", "red"))(100), Colv = F, Rowv = F, margins = c(10,10), main = "correlation of PCs with metadata")

barplot(cor(PCA_res$rotation[,1:20][!to.exclude,], matchedAnnots[,"Symtomatic"])[,1], main = "correlation with Symptomatic status", las = 2)

#############

PCAcols = rainbow(3)[as.factor(sample.annots$Symtomatic)]
names(PCAcols) = sample.annots$sampleID

plot(PCA_res$rotation[,c(4,3)], col = PCAcols[ rownames(PCA_res$rotation[,c(4,3)]) ], pch = 16, main = "PCA plot - scaled")
legend('topright', legend = levels(as.factor(sample.annots$Symtomatic)), fill = rainbow(3))

require(rgl)

plot3D::scatter3D(x = PCA_res$rotation[,1], y = PCA_res$rotation[,3], z = PCA_res$rotation[,4], col = PCAcols[ rownames(PCA_res$rotation[,c(1,2)]) ], pch = 16, 
                  bty = "u", colkey = FALSE, 
                  main ="PCs: 1 3 4", col.panel ="lightgrey",  
                  col.grid = "white",
                  theta = 20, phi = 20
                  )

################################################
# ################################################
# 
# ░█████╗░████████╗██╗░░██╗███████╗██████╗░
# ██╔══██╗╚══██╔══╝██║░░██║██╔════╝██╔══██╗
# ██║░░██║░░░██║░░░███████║█████╗░░██████╔╝
# ██║░░██║░░░██║░░░██╔══██║██╔══╝░░██╔══██╗
# ╚█████╔╝░░░██║░░░██║░░██║███████╗██║░░██║
# ░╚════╝░░░░╚═╝░░░╚═╝░░╚═╝╚══════╝╚═╝░░╚═╝

GOI = "SVEP1"
genes.to.test = c("SVEP1", "NTN4", "MFGE8", "CXCL6", "HHIP")
names(genes.to.test) = genes.to.test

ttestres = lapply(genes.to.test, function(GOI){
  #print(GOI)
  Gvals = counts.norm[GOI,]
  split.vals = split(Gvals, paste(sample.annots[names(Gvals),]$sampleType, sample.annots[names(Gvals),]$Symtomatic) )
  tres = try(t.test(split.vals$`PL AS`, split.vals$`PL S`), silent = T)
  #plot = boxplot(split.vals[!names(split.vals)=="NA NA"], main = paste0(GOI, " - KI data :)"))
  boxplot(split.vals[!names(split.vals)=="NA NA"], main = paste0(GOI, " - KI data :)"))
  if(inherits(tres, "try-error"))
  {
    return(NA)
  }
  return(list(pval = tres$p.value, tstat = tres$statistic))#, boxplot = plot))
})

all.pvals = lapply(ttestres, "[[", "pval")  
all.tstats = lapply(ttestres, "[[", "tstat")  


plot(unlist(all.tstats), -log10(unlist(all.pvals)))
tail(sort(abs(unlist(all.tstats))))

countsNormForPlotting = log2(counts.matrix / (colSums(counts.matrix) / 1000000) + 1)


topGOIS = gsub("\\.t","",names(tail(sort(abs(unlist(all.tstats))), 20)))
topGOIS = c("MFGE8")
topGOIS = c("SLC16A14")
topGOIS = rownames(counts.norm)[grep("SMPD",rownames(counts.norm))]
topGOIS = rownames(counts.norm)[grep("SGMS",rownames(counts.norm))]
topGOIS = genes.to.test

lapply(topGOIS, function(GOI){
  Gvals = unlist(countsNormForPlotting[GOI,])
  split.vals = split(Gvals, paste(sample.annots[names(Gvals),]$sampleType, sample.annots[names(Gvals),]$Symtomatic) )
  #names(split.vals)
  names(split.vals) = c("Plasma\n(CPT)\nAsymptomatic", "Plasma\n(CPT)\nSymptomatic", "Plasma\nAsymptomatic", "Plasma\nSymptomatic", "NA NA", "Plaque\nAsymptomatic", "Plaque\nSymptomatic", "Local Plasma\nAsymptomatic", "Local Plasma\nSymptomatic" )
  
  tres = try(t.test(split.vals$`Plaque\nAsymptomatic`, split.vals$`Plaque\nSymptomatic`), silent = T)
  print(tres$p.value)
  #plot = boxplot(split.vals[!names(split.vals)=="NA NA"], main = paste0(GOI, " : KI data : p = ", signif(tres$p.value,3)))
  par(mgp=c(2.8,1.5,0))
  #plot = boxplot(split.vals[!names(split.vals)=="NA NA"], main = paste0(GOI, " : BiKE RNA-Seq - Human Carotid Atherosclerotic Plaques" ), ylab = "cpm")
  plot = boxplot(split.vals[!names(split.vals) %in% c("NA NA","Plasma\n(CPT)\nAsymptomatic","Plasma\n(CPT)\nSymptomatic")][c(1,2,5,6,3,4)], 
                 main = paste0(GOI, " : BiKE RNA-Seq - Human Carotid Atherosclerotic Plaques" ), 
                 ylab = "Gene Expression: Log2(CPM+1)", 
                 #col =c("firebrick1","firebrick3","firebrick2","firebrick4","tomato1", "tomato3")
                 #col =c("tomato1","firebrick2","tomato2","firebrick3","tomato3", "firebrick3")
                 col =c("dodgerblue1","firebrick2","dodgerblue2","firebrick3","dodgerblue3", "firebrick3")
                 #col =c("dodgerblue1","dodgerblue4","magenta1","magenta4","firebrick1", "firebrick4")
  )
})


##########################
# check il18 vs il21r

CPTsamps = as.character(sample.annots$sampleID[sample.annots$sampleType == "CPT-EDTA"])
corres = cor.test(counts.norm["IL18",CPTsamps], counts.norm["IL21R",CPTsamps])
plot(counts.norm["IL18",CPTsamps], counts.norm["IL21R",CPTsamps], pch = 16, xlab = "IL18", ylab = "IL21R", main = paste0("IL18 vs. IL21R in BiKE RNA-Seq CPT EDTA\n cor = ",signif(corres$estimate, 3)," pval = ",signif(corres$p.value, 3)))

##########################
# check il18 vs il21r

CPTsamps = as.character(sample.annots$sampleID[sample.annots$sampleType == "CPT-EDTA"])
corres = cor.test(counts.norm["IL18",CPTsamps], counts.norm["IL21R",CPTsamps])
plot(counts.norm["IL18",CPTsamps], counts.norm["IL6",CPTsamps], pch = 16, xlab = "IL18", ylab = "IL6", main = paste0("IL18 vs. IL21R in BiKE RNA-Seq CPT EDTA\n cor = ",signif(corres$estimate, 3)," pval = ",signif(corres$p.value, 3)))


#########################################
## Basic boxplots
#########################################

basic.boxplot <- function(GOI, to.plot = T){
  Gvals = unlist(countsNormForPlotting[GOI,])
  split.vals = split(Gvals, paste(sample.annots[names(Gvals),]$sampleType, sample.annots[names(Gvals),]$Symtomatic) )
  #names(split.vals)
  names(split.vals) = c("Plasma\n(CPT)\nAsymptomatic", "Plasma\n(CPT)\nSymptomatic", "Plasma\nAsymptomatic", "Plasma\nSymptomatic", "NA NA", "Plaque\nAsymptomatic", "Plaque\nSymptomatic", "Local Plasma\nAsymptomatic", "Local Plasma\nSymptomatic" )
  
  #tres = try(t.test(split.vals$`Plaque\nAsymptomatic`, split.vals$`Plaque\nSymptomatic`), silent = T)
  #print(tres$p.value)
  #plot = boxplot(split.vals[!names(split.vals)=="NA NA"], main = paste0(GOI, " : KI data : p = ", signif(tres$p.value,3)))
  par(mgp=c(2.8,1.5,0))
  #plot = boxplot(split.vals[!names(split.vals)=="NA NA"], main = paste0(GOI, " : BiKE RNA-Seq - Human Carotid Atherosclerotic Plaques" ), ylab = "cpm")
  bbpl = boxplot(split.vals[!names(split.vals) %in% c("NA NA","Plasma\n(CPT)\nAsymptomatic","Plasma\n(CPT)\nSymptomatic")][c(1,2,5,6,3,4)], 
                 main = paste0(GOI, " : BiKE RNA-Seq - Human Carotid Atherosclerotic Plaques" ), 
                 ylab = "Gene Expression: Log2(CPM+1)", 
                 col =c("dodgerblue1","firebrick2","dodgerblue2","firebrick3","dodgerblue3", "firebrick3"),
                 plot = to.plot
  )
  return(bbpl)
}

#########################################
## Precompute all boxplots
#########################################

all.gene.boxplots = lapply(rownames(countsNormForPlotting), basic.boxplot, to.plot = F)
names(all.gene.boxplots) = rownames(countsNormForPlotting)
# 
save(all.gene.boxplots, file = "datasets/KI/RNA-Seq/All_RNASeq_Boxplots.RData")
save(all.gene.boxplots, file = "datasets/KI/BiKE_Omics_App/BikE_Omics_App/All_RNASeq_Boxplots.RData")


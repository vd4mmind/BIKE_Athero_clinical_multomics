#########################################
## 
## BiKE Omics - Olink data - QC
## Script by Djordje Djordjevic drdj@novonordisk.com
##
#########################################

#########################################
## Install packages
#########################################

# install.packages("devtools")
#devtools::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze')

#old_path <- Sys.getenv("LD_LIBRARY_PATH") 
#Sys.setenv(LD_LIBRARY_PATH = paste(old_path, "/opt/R/4.0.3/lib64/R/library", sep = ":"))

require(OlinkAnalyze)

#########################################
## Read in the data
#########################################
# I had to modify the file back to the original format to be able to use the read NPX function

# reading Olink NPX data 
my_NPX_data <- read_NPX(filename = "~/datasets/KI/Proteomics/20200468_Matic_NPX PseudoID_for_analysis.xlsx")
my_NPX_header <- read_excel("~/datasets/KI/Proteomics/20200468_Matic_NPX PseudoID_for_analysis.xlsx", 2)

panels = unique(my_NPX_data$Panel)

#########################################
## Plot distributions
#########################################

plot(c(-5,15),c(0,0.3), main = "Distributions of values", xlab = "NPX val", ylab = "Density")

lapply(1:length(panels), function(P){
  PANEL = panels[P]
  a = my_NPX_data %>% filter(Panel == PANEL) %>% select(NPX)
  #  hist(a$NPX, col = "grey", xlim = c(-5,15))
  lines(density(na.omit(a$NPX)), col = rainbow(5)[P], lwd = 2)
  
  text(x = -3, y = 0.05*P, PANEL, col = rainbow(5)[P]) 
})

# visualize the NPX distribution per sample per panel, 
lapply(panels, function(PANEL){
  olink_dist_plot(my_NPX_data %>% filter(Panel == PANEL)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c('turquoise3', 'red'))
})

lapply(panels, function(PANEL){
  # visualize potential outliers by IQR vs. sample median per panel, example for one panel
  olink_qc_plot(my_NPX_data %>% filter(Panel == PANEL)) +
    scale_color_manual(values = c('turquoise3', 'red'))
})

#########################################
## Extract OLink defined outliers
#########################################

OLINK_outliers = lapply(panels, function(PANEL){
  # visualize potential outliers by IQR vs. sample median per panel, example for one panel
  OLs = olink_qc_plot(my_NPX_data %>% filter(Panel == PANEL))$data$Outlier
  names(OLs) = olink_qc_plot(my_NPX_data %>% filter(Panel == PANEL))$data$SampleID
  return(OLs)
})
names(OLINK_outliers) = panels
multi_OLs = rowSums(do.call(cbind,OLINK_outliers))
table(multi_OLs)
OL_OLs = rev(sort(na.omit(multi_OLs[multi_OLs > 0])))

#########################################
## Check duplicate proteins
#########################################

duplicate.proteins = names(which(table(my_NPX_data[,c("Assay")]) > 900))
columns = c("SampleID","Assay","Panel","NPX")
subd = my_NPX_data %>% select(columns) 

par(mfrow=c(3,6))

lapply(duplicate.proteins, function(DUPLICATE){
  subm = subd %>% filter(Assay == DUPLICATE)
  split.data = split(subm$NPX, subm$Panel)
  dup.cor = cor(na.omit(cbind(scale(split.data[[1]]), scale(split.data[[2]]))))[2]
  #plot(scale(split.data[[1]]), scale(split.data[[2]]), xlab = names(split.data)[1], ylab = names(split.data)[2], main = paste0(DUPLICATE, " : ", signif(dup.cor,2)), pch = 16, col = "black")
  plot((split.data[[1]]), (split.data[[2]]), xlab = names(split.data)[1], ylab = names(split.data)[2], main = paste0(DUPLICATE, " : ", signif(dup.cor,2)), pch = 16, col = "black")
  lines(x=c(-5,20), y =c(-5,20), col = "red")  
})

par(mfrow=c(1,1))


#########################################
## Remove failing samples
#########################################


rows.failing.QC = my_NPX_data$QC_Warning == "Warning"
table(rows.failing.QC)
#rows.failing.QC
#FALSE   TRUE 
#373520   3680 
#1% failure

samples.fail.QC = unique(my_NPX_data$SampleID[rows.failing.QC])
assays.fail.QC = unique(my_NPX_data$Assay[rows.failing.QC])

Clean_Data = my_NPX_data[!rows.failing.QC,]


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
plot(PCA_res)

PCAcols = c("red","blue")[grepl("EP",rownames(PCA_res[[2]]))+1]
plot(PCA_res$rotation[,c(1,2)], col = PCAcols, pch = 16, main = "PCA plot")
legend('topright', legend = c("STP - local","EP - peripheral"), fill = c("red","blue"))

PCAcols = c("green", "purple")[npx_annots$Olink_outliers[match(rownames(PCA_res[[2]]), rownames(npx_annots))]+1]
plot(PCA_res$rotation[,c(1,2)], col = PCAcols, pch = 16, main = "PCA plot")
legend('topright', legend = c("O-Link outlier","Within range"), fill = c("purple","green"))

#########################################
## PCA based outliers
#########################################

## +3 SD away from mean
pc1sd = sd(PCA_res$rotation[,1])
pc1mean = mean(PCA_res$rotation[,1])
#pc1mean = median(PCA_res$rotation[,1])
pc1outliers = which(abs(PCA_res$rotation[,1]) > pc1mean+3*pc1sd)

pc2sd = sd(PCA_res$rotation[,2])
pc2mean = mean(PCA_res$rotation[,2])
#pc2mean = median(PCA_res$rotation[,2])
pc2outliers = which(abs(PCA_res$rotation[,2]) > pc2mean+3*pc2sd)

PCAcols = c("green", "purple")[rownames(PCA_res[[2]]) %in% c(names(pc1outliers), names(pc2outliers)) +1]
plot(PCA_res$rotation[,c(1,2)], col = PCAcols, pch = 16, main = "PCA plot")
legend('topright', legend = c("PCA outlier","Within range"), fill = c("purple","green"))


#########################################
## PCA correlation heatmap
#########################################

pcaCors = cor(PCA_res$rotation[,1:20], data.matrix(npx_annots[match(rownames(PCA_res[[2]]), rownames(npx_annots)),]))
require(gplots)
heatmap.2(pcaCors, trace = "none", scale = "none", col = colorRampPalette(c("blue", "white", "red"))(100), Colv = F, Rowv = F, margins = c(10,10), main = "correlation of PCs with metadata")

##
PCAcols = c("green","purple")[grepl("A",npx_annots$Symtomatic[match(rownames(PCA_res[[2]]), rownames(npx_annots))])+1]

barplot(cor(grepl("A",npx_annots$Symtomatic[match(rownames(PCA_res[[2]]), rownames(npx_annots))])+1, PCA_res$rotation[,1:20]), main = "correlation PCs with symptomatic status", las = 2)

plot(PCA_res$rotation[,c(3,5)], col = PCAcols, pch = 16, main = "PCA plot")
legend('topright', legend = c("S - Symptomatic","AS - Asymptomatic"), fill = c("green","purple"))
plot(PCA_res$rotation[,c(5,10)], col = PCAcols, pch = 16, main = "PCA plot")
legend('topright', legend = c("S - Symptomatic","AS - Asymptomatic"), fill = c("green","purple"))

#########################################
## Remove outliers
#########################################

PCA_outliers = unique(c(names(pc1outliers), names(pc2outliers)))
combined_outliers = unique(c(PCA_outliers, names(OL_OLs)))
length(combined_outliers)  
#44 outliers

### remove OLink outliers in a panel specific way
OLINK_outliers_names = lapply(OLINK_outliers, function(X){names(which(na.omit(X)>0))})

naCl = na.omit(Clean_Data)
dim(naCl)

naCl$Olink_outlier = F

lapply(names(OLINK_outliers_names), function(Panel){
  naCl[naCl$Panel == Panel & naCl$SampleID %in% OLINK_outliers_names[[Panel]], "Olink_outlier"] <<- TRUE 
})

dim(naCl[naCl$Olink_outlier == F,])

#Clean_OLink_Data = naCl[naCl$Olink_outlier == F & !naCl$SampleID %in% PCA_outliers,]
Clean_OLink_Data = naCl[naCl$Olink_outlier == F,]

#########################################
## Prepare to save the files
#########################################

Clean_OLink_Data$ProtPan = paste(Clean_OLink_Data$Assay, "|", gsub("Olink Target 96 ","",Clean_OLink_Data$Panel))

##

npxA = npx_annots[!grepl("control|na", rownames(npx_annots), ignore.case = T),]
npxA = npxA[!npxA$Symtomatic == "NA",]

npxA$Age = as.numeric(npxA$Age)
npxA$Sex = as.factor(npxA$Sex)
npxA$BiKE.ID.prefix = as.factor(npxA$BiKE.ID.prefix)
npxA$PID = as.factor(npxA$PID)

npxA$Symptomatic = as.factor(npxA$Symtomatic) # note the fix of spelling

sample_annots = npxA[,!colnames(npxA)=="Symtomatic"]

Clean_OLink_Data = Clean_OLink_Data[,!colnames(Clean_OLink_Data) == "Olink_outlier"]
  
#########################################
## Save the files
#########################################

write.table(Clean_OLink_Data, file = "datasets/KI/Proteomics/Clean_OLink_Data.txt", sep="\t", row.names = F, col.names = T, quote = F)
write.table(sample_annots, file = "datasets/KI/Proteomics/Clean_OLink_Sample_Annots.txt", sep="\t", row.names = F, col.names = T, quote = F)

save(sample_annots, Clean_OLink_Data, file = "datasets/KI/Proteomics/Clean_OLink_Data.RData")

#########################################
## END SCRIPT
#########################################

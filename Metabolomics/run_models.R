#########################################
## 
## BiKE Omics - Metabolon data - Run the models template
## Script by Djordje Djordjevic drdj@novonordisk.com
##
#########################################

#########################################
## Install packages
#########################################

require(lme4)
require(lmerTest)
require(dplyr)
#########################################
## Read in the data
#########################################
## two data files
## Metabolomics QC
metabolite_annots = readxl::read_xlsx("datasets/KI/Metabolomics/KARO-05-20MD DATA TABLES Sample Meta Data PseudoID.xlsx", 2)
sample_annots = readxl::read_xlsx("datasets/KI/Metabolomics/KARO-05-20MD DATA TABLES Sample Meta Data PseudoID.xlsx", 3)
log2vals = readxl::read_xlsx("datasets/KI/Metabolomics/KARO-05-20MD DATA TABLES Sample Meta Data PseudoID.xlsx", 7)


#########################################
## Prepare data
#########################################
## prepare data if neccessary
l2v = signif(data.matrix(log2vals[,-1]),4)
rownames(l2v) = unlist(log2vals[,1])
metabolon.data = t(l2v)

rownames(sample_annots) = sample_annots$PARENT_SAMPLE_NAME
colnames(sample_annots)[colnames(sample_annots) == "Symtomatic"] = "Symptomatic"

bad.samples = sample_annots[sample_annots$Symptomatic == 'NA', ]$PARENT_SAMPLE_NAME
metabolon.data = metabolon.data[,!colnames(metabolon.data) %in% bad.samples]
sample_annots = sample_annots[!rownames(sample_annots) %in% bad.samples,]

bad.metabolites = which(rowSums(metabolon.data) == 0)
metabolon.data = metabolon.data[-bad.metabolites,]

list.of.metabolites = rownames(metabolon.data)

#########################################
## Run the models
#########################################
## for each protein, extract that data  and run the 4 models

lmvals = lapply(list.of.metabolites, function(Pr){
  #print(Pr)
  print(paste0(Pr, ": ", metabolite_annots[metabolite_annots$CHEM_ID == Pr,"CHEMICAL_NAME"]))
  
  # just this metabolite
  SubC = metabolon.data[Pr,]
  
  ############################
  # model 1
  # symptomatic vs asymptomatic - peripheral only
  m1data = SubC[sample_annots[sample_annots$BiKE.ID.prefix == "EP",]$PARENT_SAMPLE_NAME]
  m1annots = na.omit(sample_annots[match(names(m1data), sample_annots$PARENT_SAMPLE_NAME),])
  m1DF = data.frame(val = m1data, m1annots)
  
  m1 <- lm(val ~ Symptomatic + Age + Sex, data = m1DF)
  
  m1.coef = summary(m1)$coefficients["SymptomaticS","Estimate"]
  m1.pval = summary(m1)$coefficients["SymptomaticS","Pr(>|t|)"]
  
  
  ############################
  # model 2
  # symptomatic vs asymptomatic - peripheral  only
  m2data = SubC[sample_annots[sample_annots$BiKE.ID.prefix == "STP",]$PARENT_SAMPLE_NAME]
  m2annots = na.omit(sample_annots[match(names(m2data), sample_annots$PARENT_SAMPLE_NAME),])
  m2DF = data.frame(val = m2data, m2annots)
  
  m2 <- lm(val ~ Symptomatic + Age + Sex, data = m2DF)
  
  m2.coef = summary(m2)$coefficients["SymptomaticS","Estimate"]
  m2.pval = summary(m2)$coefficients["SymptomaticS","Pr(>|t|)"]
  
  
  # ############################
  # # model 3
  # # protdiff: local - peripheral (only those with both)
  # m3data = SubC
  # pids = gsub("(.*)_.*","\\1",m3data$SampleID)
  # multiple_pids = names(which(table(pids)>1))
  # protdiff = lapply(multiple_pids, function(PID){
  #   STP_ID = grep(paste0(PID, "_STP"), m3data$SampleID)
  #   EP_ID = grep(paste0(PID, "_EP"), m3data$SampleID)
  #   return(unlist(m3data[STP_ID, "NPX"] - m3data[EP_ID, "NPX"]))
  # })
  # names(protdiff) = multiple_pids
  # m3annots = sample_annots[sample_annots$PID %in% names(protdiff) & !duplicated(sample_annots$PID),]
  # rownames(m3annots) = m3annots$PID
  # m3DF = data.frame(val = as.numeric(unlist(protdiff)), m3annots[names(protdiff),])
  # 
  # m3 <- lm(val ~ Symptomatic + Age + Sex, data = m3DF)
  # 
  # m3.coef = summary(m3)$coefficients["SymptomaticS","Estimate"]
  # m3.pval = summary(m3)$coefficients["SymptomaticS","Pr(>|t|)"]

  m3.coef = 0
  m3.pval = 1
  
  ############################
  # model 4
  # integrated model with random effect - all data points
  # two covariates 
  m4data = SubC
  m4annots = sample_annots[match(names(m4data), sample_annots$PARENT_SAMPLE_NAME),]
  m4DF = data.frame(val = m4data, m4annots)
  
  m4 <- lmer(val ~ Symptomatic + Age + Sex + BiKE.ID.prefix + (1|PseudoID.metabolomics), data = m4DF)
  
  m4.sym.coef = summary(m4)$coefficients["SymptomaticS","Estimate"]
  m4.sym.pval = summary(m4)$coefficients["SymptomaticS","Pr(>|t|)"]
  
  m4.loc.coef = summary(m4)$coefficients["BiKE.ID.prefixSTP","Estimate"]
  m4.loc.pval = summary(m4)$coefficients["BiKE.ID.prefixSTP","Pr(>|t|)"]
  
  
  res = data.frame(m1.coef = m1.coef, m1.pval = m1.pval,
                   m2.coef = m2.coef, m2.pval = m2.pval,
                   m3.coef = m3.coef, m3.pval = m3.pval,
                   m4.sym.coef = m4.sym.coef, m4.sym.pval = m4.sym.pval, m4.loc.coef = m4.loc.coef, m4.loc.pval = m4.loc.pval)
  return(res)  
})
names(lmvals) = list.of.metabolites

resDF = do.call(rbind, lmvals)

metaboLongNames = metabolite_annots[match(rownames(resDF), metabolite_annots$CHEM_ID), ]$CHEMICAL_NAME

rownames(resDF) = metaboLongNames

#########################################
## Save results
#########################################

write.table(resDF, file = "datasets/KI/Metabolomics/Metabolon_models_results.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(resDF, file = "datasets/KI/BiKE_Omics_App/BikE_Omics_App/Metabolon_models_results.txt", sep = "\t", quote = F, row.names = T, col.names = T)

#########################################
## Plot results
#########################################

sig.metabolites.plot.and.list = lapply(c("m1","m2","m3","m4.sym", "m4.loc"), function(M){
  
  coefs = resDF[,paste0(M,".coef")]  
  pvals = resDF[,paste0(M,".pval")] 
  metaboNames = rownames(resDF)
  #metaboLongNames = metabolite_annots[match(metaboNames, metabolite_annots$CHEM_ID), ]$CHEMICAL_NAME
  metaboLongNames = metaboNames
  minusLog10AdjPvals = -log10(p.adjust(pvals, "BH"))
  mostsig = which(minusLog10AdjPvals > -log10(0.05))
  
  plot(coefs, minusLog10AdjPvals, col = "lightgreen", pch = 16, main = paste0("Metabolon Symptomatic Coefficient: ", M))
  
  abline(h = -log10(0.05), lty = 2)
  abline(v = 0, lty = 2, col = "grey")
  
  if(length(mostsig) > 0){
    text(coefs[mostsig], minusLog10AdjPvals[mostsig]+0.02+0.01*(rnorm(length(mostsig))), labels = metaboLongNames[mostsig], cex = 0.7)
    return(mostsig)
  }
  
  return(NULL)
})
names(sig.metabolites.plot.and.list) = c("m1","m2","m3","m4.sym", "m4.loc")



# how many sig results
lapply(sig.metabolites.plot.and.list, length)

#########################################
## Split the data for boxplots
#########################################

# symptomatic coefficient
sig.metabolites = sort(rownames(resDF)[sig.metabolites.plot.and.list$m4.sym])
sig.metabolites = sort(rownames(resDF)[sig.metabolites.plot.and.list$m1])
sig.metabolites = sort(rownames(resDF))

# location coefficient
#sig.metabolites = rownames(resDF)[sig.metabolites.plot.and.list$m4.loc]
#sig.metabolites = sort(rownames(resDF)[head(order(resDF$m4.loc.pval), 9)])

valsU = lapply(sig.metabolites, function(Mc){
  Pr = as.character(metabolite_annots$CHEM_ID[metabolite_annots$CHEMICAL_NAME == Mc])
  subC = metabolon.data[Pr,]
  subA = sample_annots[match(names(subC), sample_annots$PARENT_SAMPLE_NAME),]
  groups = paste(as.character(subA$Symptomatic), as.character(subA$BiKE.ID.prefix), sep = ".")
  
  vals_split = split(subC, groups)
  
  symp_peri = vals_split$S.EP
  symp_local = vals_split$S.STP
  asymp_peri = vals_split$AS.EP
  asymp_local = vals_split$AS.STP
  
  valsList = list(symp_peri =symp_peri, asymp_peri = asymp_peri, symp_local = symp_local, asymp_local = asymp_local)
})

names(valsU) = sig.metabolites


#########################################
## Basic boxplots
#########################################

basic.boxplot <- function(G, to.plot = T){
  bbpl = boxplot(valsU[[G]], main = G, col = c("firebrick1","dodgerblue","firebrick4","steelblue"), plot = to.plot)
  return(bbpl)
}

#########################################
## Precompute all boxplots
#########################################

 all.metabolite.boxplots = lapply(rownames(resDF), basic.boxplot, to.plot = F)
 names(all.metabolite.boxplots) = rownames(resDF)
# 
 save(all.metabolite.boxplots, file = "datasets/KI/Metabolomics/All_Metabolon_Boxplots.RData")
 save(all.metabolite.boxplots, file = "datasets/KI/BiKE_Omics_App/BikE_Omics_App/All_Metabolon_Boxplots.RData")

#########################################
## Volcano plots
#########################################





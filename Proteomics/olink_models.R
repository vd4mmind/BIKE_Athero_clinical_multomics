#########################################
## 
## BiKE Omics - Olink data - models
## Script by Djordje Djordjevic drdj@novonordisk.com
##
#########################################

#########################################
## Install packages
#########################################

require(lme4)
require(lmerTest)
#########################################
## Read in the data
#########################################

load("datasets/KI/Proteomics/Clean_OLink_Data.RData")
ls()
#"Clean_OLink_Data" "sample_annots" 

#########################################
## Prepare data
#########################################
list.of.proteins.by.panel = na.omit(unique(Clean_OLink_Data$ProtPan))


#########################################
## Run the models
#########################################
## for each protein, extract that data  and run the 4 models

lmvals = lapply(list.of.proteins.by.panel, function(Pr){
  print(Pr)
  
  # just this protein * panel
  subC = Clean_OLink_Data %>% 
    filter(ProtPan == Pr) %>% 
    filter(!grepl("CONTROL|NA", SampleID, ignore.case = T)) 
  
  ############################
  # model 1
  # symptomatic vs asymptomatic - peripheral only
  m1data = subC %>% 
    filter(grepl("EP", SampleID)) 
  m1annots = sample_annots[m1data$SampleID,]
  m1DF = data.frame(val = m1data$NPX, m1annots)
  
  m1 <- lm(val ~ Symptomatic + Age + Sex, data = m1DF)
  
  m1.coef = summary(m1)$coefficients["SymptomaticS","Estimate"]
  m1.pval = summary(m1)$coefficients["SymptomaticS","Pr(>|t|)"]
  
  
  ############################
  # model 2
  # symptomatic vs asymptomatic - peripheral  only
  m2data = subC %>% 
    filter(grepl("STP", SampleID))
  m2annots = sample_annots[m2data$SampleID,]
  m2DF = data.frame(val = m2data$NPX, m2annots)
  
  m2 <- lm(val ~ Symptomatic + Age + Sex, data = m2DF)
  
  m2.coef = summary(m2)$coefficients["SymptomaticS","Estimate"]
  m2.pval = summary(m2)$coefficients["SymptomaticS","Pr(>|t|)"]
  
  
  ############################
  # model 3
  # protdiff: local - peripheral (only those with both)
  m3data = subC
  pids = gsub("(.*)_.*","\\1",m3data$SampleID)
  multiple_pids = names(which(table(pids)>1))
  protdiff = lapply(multiple_pids, function(PID){
    STP_ID = grep(paste0(PID, "_STP"), m3data$SampleID)
    EP_ID = grep(paste0(PID, "_EP"), m3data$SampleID)
    return(unlist(m3data[STP_ID, "NPX"] - m3data[EP_ID, "NPX"]))
  })
  names(protdiff) = multiple_pids
  m3annots = sample_annots[sample_annots$PID %in% names(protdiff) & !duplicated(sample_annots$PID),]
  rownames(m3annots) = m3annots$PID
  m3DF = data.frame(val = as.numeric(unlist(protdiff)), m3annots[names(protdiff),])
  
  m3 <- lm(val ~ Symptomatic + Age + Sex, data = m3DF)
  
  m3.coef = summary(m3)$coefficients["SymptomaticS","Estimate"]
  m3.pval = summary(m3)$coefficients["SymptomaticS","Pr(>|t|)"]
  
  ############################
  # model 4
  # integrated model with random effect - all data points
  # two covariates 
  m4data = subC
  m4annots = sample_annots[m4data$SampleID,]
  m4DF = data.frame(val = m4data$NPX, m4annots)
  
  m4 <- lmer(val ~ Symptomatic + Age + Sex + BiKE.ID.prefix + (1|PID), data = m4DF)
  
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
names(lmvals) = list.of.proteins.by.panel

resDF = do.call(rbind, lmvals)

#########################################
## Save results
#########################################

write.table(resDF, file = "datasets/KI/Proteomics/Olink_models_results.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(resDF, file = "datasets/KI/BiKE_Omics_App/BikE_Omics_App/Olink_models_results.txt", sep = "\t", quote = F, row.names = T, col.names = T)

#########################################
## Plot results
#########################################
## nominal P value
nominal.sig.prots.plot.and.list = lapply(c("m1","m2","m3","m4.sym", "m4.loc"), function(M){
  
  coefs = resDF[,paste0(M,".coef")]  
  pvals = resDF[,paste0(M,".pval")] 
  protnames = rownames(resDF)
  minusLog10Pvals = -log10(pvals)
  minusLog10AdjPvals = -log10(p.adjust(pvals, "BH"))
  mostsig = which(minusLog10Pvals > -log10(0.05))
  
  plot(coefs, minusLog10Pvals, col = "lightgreen", pch = 16, main = paste0("Olink Symptomatic Coefficient: ", M))
  
  abline(h = -log10(0.05), lty = 2)
  abline(v = 0, lty = 2, col = "grey")
  
  if(length(mostsig) > 0){
    text(coefs[mostsig], minusLog10Pvals[mostsig]+0.02+0.01*(rnorm(length(mostsig))), labels = gsub("(.*) \\|.*","\\1",protnames[mostsig]), cex = 0.7)
    return(mostsig)
  }
  
  return(NULL)
})
names(nominal.sig.prots.plot.and.list) = c("m1","m2","m3","m4.sym", "m4.loc")

## adjusted p-value
sig.prots.plot.and.list = lapply(c("m1","m2","m3","m4.sym", "m4.loc"), function(M){

  coefs = resDF[,paste0(M,".coef")]
  pvals = resDF[,paste0(M,".pval")]
  protnames = rownames(resDF)
  minusLog10AdjPvals = -log10(p.adjust(pvals, "BH"))
  mostsig = which(minusLog10AdjPvals > -log10(0.05))

  plot(coefs, minusLog10AdjPvals, col = "lightgreen", pch = 16, main = paste0("Olink Symptomatic Coefficient: ", M))

  abline(h = -log10(0.05), lty = 2)
  abline(v = 0, lty = 2, col = "grey")

  if(length(mostsig) > 0){
    text(coefs[mostsig], minusLog10AdjPvals[mostsig]+0.02+0.01*(rnorm(length(mostsig))), labels = gsub("(.*) \\|.*","\\1",protnames[mostsig]), cex = 0.7)
    return(mostsig)
  }

  return(NULL)
})
names(sig.prots.plot.and.list) = c("m1","m2","m3","m4.sym", "m4.loc")

# how many sig results
lapply(sig.prots.plot.and.list, length)

#########################################
## Split the data for boxplots
#########################################

# symptomatic coefficient
#sig.prots = sort(rownames(resDF)[sig.prots.plot.and.list$m4.sym])
#sig.prots = sort(rownames(resDF)[sig.prots.plot.and.list$m1])

# location coefficient
#sig.prots = rownames(resDF)[sig.prots.plot.and.list$m4.loc]
#sig.prots = sort(rownames(resDF)[head(order(resDF$m4.loc.pval), 9)])

valsU = lapply(sort(rownames(resDF)), function(Pr){
  subC = Clean_OLink_Data %>% filter(ProtPan == Pr) %>% filter(! SampleID %in% rownames(sample_annots)[grep("control|na", rownames(sample_annots), ignore.case = T) ])
  subA = sample_annots[subC$SampleID,]
  groups = paste(as.character(subA$Symptomatic), as.character(subA$BiKE.ID.prefix), sep = ".")
  
  vals_split = split(subC$NPX, groups)
  
  symp_peri = vals_split$S.EP
  symp_local = vals_split$S.STP
  asymp_peri = vals_split$AS.EP
  asymp_local = vals_split$AS.STP
  
  valsList = list(symp_peri =symp_peri, asymp_peri = asymp_peri, symp_local = symp_local, asymp_local = asymp_local)
})

names(valsU) = sort(rownames(resDF))

#########################################
## Fancy boxplots 
#########################################

multiplots = function(G){
  print(G)
  ggdf = data.frame(OLink_Value = unlist(valsU[[G]]))
  ggdf$group = gsub("[0-9]+","",rownames(ggdf))
  
  # box plot
  p <- ggplot(ggdf, aes(x=group, y=OLink_Value, fill=group)) + 
    #geom_violin()
    geom_jitter(shape=1, position=position_jitter(0.2), alpha = 0.1)+
    geom_boxplot(outlier.alpha = 0, alpha = 0.85) +
    theme_classic()+ 
    scale_fill_manual(values=c("dodgerblue","steelblue","firebrick1","firebrick4",0.4)) + 
    labs(title=G) + 
    theme(legend.position = "none")
  
  return(p)
}

prots.to.plot = lapply(sig.prots, multiplots)

batchedPlots = split(prots.to.plot, ceiling(seq_along(prots.to.plot)/9))
lapply(batchedPlots, function(BB){
  gridExtra::grid.arrange( grobs = BB, nrow = 3, ncol = 3 )
  return("Done")
})



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

 all.protein.boxplots = lapply(rownames(resDF), basic.boxplot, to.plot = F)
 names(all.protein.boxplots) = rownames(resDF)
# 
 save(all.protein.boxplots, file = "datasets/KI/Proteomics/All_Olink_Boxplots.RData")
 save(all.protein.boxplots, file = "datasets/KI/BiKE_Omics_App/BikE_Omics_App/All_Olink_Boxplots.RData")

#########################################
## Volcano plots
#########################################




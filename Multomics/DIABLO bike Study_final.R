#########################################
## 
## BiKE Multi Omics -script
## Script by Xiang Zhang and Vivek Das (vvda@novnordisk.com)
##
#########################################

#########################################
## Install packages
#########################################

library(mixOmics) # import the mixOmics library
library(dplyr)
library(DescTools)
library(ggplot2)
library(showtext)
library(psych)
library(corrplot)
library(stringr)
library (VennDiagram)
library(reshape2)
library(gdata)
library(pheatmap)
library(PCAtools)
library(kmodR)
library(ropls)
library(limma)
library(statmod)
library(variancePartition)
require(lmerTest)
library(openxlsx)
library(randomcoloR)
#require(lme4)
library(mlbench)
library(caret)
library(dendextend)
library(ggrepel)
library(ggforce)
library(factoextra)
library(NbClust)
library(cluster)
library(reshape2)
library(Mfuzz)
library(Hmisc)
library(corrplot)
library(DescTools)
library(BiocParallel)

set.seed(123) # for reproducibility, remove for normal use

#########################################
## Read in the data
#########################################

## -------------------------------------------------------------------------------------------------------------------
#data(.bike) # load in the data
####
### input matrix of metabolomics proteomics and transcriptomics

data = list(meta = dmc, # set a list of all the X dataframes
            prot = dpc,
	    mRNA = dtc)

lapply(data, dim) # check their dimensions

Y = factor(clini$S.AS) # set the response variable as the Y dataframe
summary(Y)


####Data-driven option: we could perform regression analyses 
#with PLS to further understand the correlation between data sets. 
pls1 <- pls(data$meta, data$prot, ncomp = 1)
cor(pls1$variates$X, pls1$variates$Y)

pls2 <- pls(data$meta, data$mRNA, ncomp = 1)
cor(pls2$variates$X, pls2$variates$Y)

pls3 <- pls(data$prot, data$mRNA, ncomp = 1)
cor(pls3$variates$X, pls3$variates$Y)

m1 <- cor(pls1$variates$X, pls1$variates$Y)[1,1]
m2 <- cor(pls2$variates$X, pls2$variates$Y)[1,1]
m3 <- cor(pls3$variates$X, pls3$variates$Y)[1,1]


#indicating that a design with weights could be chosen.


## ---- fig.show = "hold", out.width = "33%", 
#fig.cap = "FIGURE 1: Circle Correlation Plots for pairwise PLS models on the  bike data. 
#Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5. "----
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)


spls1 <- spls(data[["meta"]], data[["prot"]], keepX = list.keepX, keepY = list.keepY)
spls2 <- spls(data[["meta"]], data[["mRNA"]], keepX = list.keepX, keepY = list.keepY) # generate three pairwise PLS models
spls3 <- spls(data[["mRNA"]], data[["prot"]], keepX = list.keepX, keepY = list.keepY)

plotVar(spls1, cutoff = 0.5, title = "(a) meta vs mRNA", legend = c("meta", "mRNA"), # plot features of first PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(spls2, cutoff = 0.3, title = "(a) meta vs prot", legend = c("meta", "prot"), # plot features of second PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(spls3, cutoff = 0.3, title = "(a) mRNA vs prot", legend = c("mRNA", "prot"), # plot features of third PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))


## ---- echo = FALSE--------------------------------------------------------------------------------------------------
#pls1 <- spls(data[["meta"]], data[["mRNA"]], ncomp = 1, keepX = 25, keepY = 25)
#pls2 <- spls(data[["meta"]], data[["prot"]], ncomp = 1, keepX = 25, keepY = 25)
#pls3 <- spls(data[["mRNA"]], data[["prot"]], ncomp = 1, keepX = 25, keepY = 25)


## -------------------------------------------------------------------------------------------------------------------
cor(spls1$variates$X, spls1$variates$Y) # calculate correlation of meta and mRNA
cor(spls2$variates$X, spls2$variates$Y) # calculate correlation of meta and proteins
cor(spls3$variates$X, spls3$variates$Y) # calculate correlation of mRNA and proteins


ms1 <- cor(spls1$variates$X, spls1$variates$Y)[1,1]
ms2 <- cor(spls2$variates$X, spls2$variates$Y)[1,1]
ms3 <- cor(spls3$variates$X, spls3$variates$Y)[1,1]


## -------------------------------------------------------------------------------------------------------------------
#values above 0.5 will cause a reduction in predictive ability of the model – and prediction is 
#what is desired in this context. Hence, a value of 0.1 will be used to prioritise the discriminative ability of the model.
#Weighted design matrix can be 10%, 30% or no contribution that is set as 0.1, 0.3 or 0 
#across each omics layers assuming their contribution to predict outcome of Asympto-matic vs Symptomatic
#Weights can be unequal as well. 
#This is not always easy to assume. 
#Needs a few itera-tions. In my experience I have tried 0, 0.1 and 0.3 
#while integrating 3-4 different omics layers and they have worked well. However, this is very much data dependent. 
#design = matrix(0, ncol = length(data), nrow = length(data), # for square matrix filled with 0.1s
#                dimnames = list(names(data), names(data)))
#diag(design) = 0 # set diagonal to 0s

#design

#MyDesign <- matrix(c(0, 0.1, 0.3,
#                    0.1, 0, 0.9,
#                    0.3, 0.9, 0),
#                   byrow=TRUE,
#                   ncol = length(X), nrow = length(X),
#                 dimnames = list(names(X), names(X)))
#MyDesign

#Step 1: with initial 0.1 across each omics and  diagonal you get ncomp = 2. 
#Use that for tune.block.splsda().

m1=0.1
m2=0.1
m3=0.1

design <- matrix(c(0, m1, m2,  ###m or ms or 0 or 0.1 or 0.3
                   m1, 0, m3,
                   m2, m3, 0),
                   byrow=TRUE,
                   ncol = length(data), nrow = length(data),
                 dimnames = list(names(data), names(data)))
design

## -------------------------------------------------------------------------------------------------------------------
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design) # form basic DIABLO model


## ---- fig.cap = "FIGURE 2: Choosing the number of components in 
#`block.plsda` using `perf()` with 10 × 10-fold CV function in the `.bike` study. 
#Classification error rates (overall and balanced, see Section 7.3) are represented on 
#the y-axis with respect to the number of components on the x-axis for each prediction distance presented in PLS-DA"----
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', folds = 10, nrepeat = 10) # run component number tuning with repeated CV

plot(perf.diablo) # plot output of tuning

##test ncomp
## -------------------------------------------------------------------------------------------------------------------
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # set the optimal ncomp value
perf.diablo$choice.ncomp$WeightedVote # show the optimal choice for ncomp for each dist metric

ncomp

#####  longtime ####paremeter test for M(0 , 0.1 , 0.2, 0.3, 0.4)


##################  M0.2
#ncomp = 3
set.seed(123) # for reproducibility, remove for normal use
m1=0.2
m2=0.2
m3=0.2

design <- matrix(c(0, m1, m2,  ###m or ms or 0 or 0.1 or 0.3
                   m1, 0, m3,
                   m2, m3, 0),
                   byrow=TRUE,
                   ncol = length(data), nrow = length(data),
                 dimnames = list(names(data), names(data)))
design


test.keepX = list (meta = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   prot = c(5:9, seq(10, 18, 2), seq(20,30,5)), 
                   mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)))


tune.bike = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
#                              dist = "centroids.dist",BPPARAM = BPPARAM)
                             dist = "centroids.dist")
list.keepX = tune.bike$choice.keepX # set the optimal values of features to retain
list.keepX

save(tune.bike,list.keepX,design, file = 'result-bike-diablo_design_0.2.RData')



list.keepX


## -------------------------------------------------------------------------------------------------------------------
step.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, # set the optimised DIABLO model
                          keepX = list.keepX, design = design)


## -------------------------------------------------------------------------------------------------------------------
step.diablo.model$design # design matrix 


## ---- fig.cap = "FIGURE 3: Diagnostic plot from multiblock sPLS-DA applied on the `bike` study. 
#Samples are represented based on the specified component (here `ncomp = 1`) for 
#each data set (mRNA, meta and protein). Samples are coloured by   subtype 
#and 95% confidence ellipse plots are represented."----
plotDiablo(step.diablo.model, ncomp = 1)
plotDiablo(step.diablo.model, ncomp = 2)
plotDiablo(step.diablo.model, ncomp = 3)


 
## final model ##test which is the best
load('result-bike-diablo_design_0.2.RData') 
ncomp
list.keepX
design

## -------------------------------------------------------------------------------------------------------------------
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, # set the optimised DIABLO model
                          keepX = list.keepX, design = design)


## -------------------------------------------------------------------------------------------------------------------
final.diablo.model$design # design matrix for the final model



#results

selectVar(final.diablo.model, block = 'mRNA', comp = 1)
selectVar(final.diablo.model, block = 'mRNA', comp = 2)
selectVar(final.diablo.model, block = 'mRNA', comp = 3)

selectVar(final.diablo.model, block = 'meta', comp = 1)
selectVar(final.diablo.model, block = 'meta', comp = 2)
selectVar(final.diablo.model, block = 'meta', comp = 3)

selectVar(final.diablo.model, block = 'prot', comp = 1)
selectVar(final.diablo.model, block = 'prot', comp = 2)
selectVar(final.diablo.model, block = 'prot', comp = 3)
## -------------------------------------------------------------------------------------------------------------------

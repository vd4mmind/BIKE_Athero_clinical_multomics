#spare code 

##############
# simple box plots

par(mfrow = c(3,3))
lapply(sig.prots, function(G){
  boxplot(valsU[[G]], main = G, col = c("firebrick1","dodgerblue","firebrick4","steelblue"))
  return(G)
})
par(mfrow = c(1,1))

#############
# fancy violin plots

multiplots = lapply(sig.prots, function(G){
  ggdf = data.frame(OLink_Value = unlist(valsU[[G]]))
  ggdf$group = gsub("[0-9]+","",rownames(ggdf))
  
  # Basic violin plot
  p <- ggplot(ggdf, aes(x=group, y=OLink_Value, fill=group)) + 
    geom_violin()
  
  pp = p + geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.2) + 
    theme_classic()+ 
    scale_fill_manual(values=c("dodgerblue","steelblue","firebrick1","firebrick4",0.4)) + 
    geom_boxplot(width = 0.4, outlier.alpha = 0, fill = "white", alpha = 0.75) +
    labs(title=G) + 
    theme(legend.position = "none")
  
  return(pp)
})


#######################################################
#######################################################
#######################################################


#svm

library(e1071)

subdata = na.omit(Clean_OLink_Data[,c(1,12,14)])
formattedData = pivot_wider(data = subdata, names_from = SampleID, values_from = NPX, values_fill = NA)


valuesDF = data.frame(t(formattedData[,-1]))
colnames(valuesDF) =  formattedData$ProtPan

par(mfrow = c(3,2))
for(i in 1:6){
  hist(as.numeric(valuesDF[,i]), breaks = 100, main = paste0("distribution ", colnames(valuesDF)[i]))
}

sample_annots_formatted = sample_annots
sample_annots_formatted[sample_annots_formatted == "NA"] = NA
sample_annots_formatted = na.omit(sample_annots_formatted)
sample_annots_formatted$Age = as.numeric(sample_annots_formatted$Age)
sample_annots_formatted$Sex = as.factor(sample_annots_formatted$Sex)
sample_annots_formatted$BiKE.ID.prefix = as.factor(sample_annots_formatted$BiKE.ID.prefix)
sample_annots_formatted$PID = as.factor(sample_annots_formatted$PID)


par(mfrow = c(3,2))
hist(sample_annots_formatted$Age, main = "Age", breaks = 100)
barplot(table(sample_annots_formatted$Sex), main = "Sex")
barplot(table(sample_annots_formatted$Symtomatic), main = "Symptomatic")
barplot(table(sample_annots_formatted$BiKE.ID.prefix), main = "Sample type")
barplot(table(table(sample_annots_formatted$PID)), main = "PID")


mergedDF = merge(na.omit(sample_annots_formatted), valuesDF, by = "row.names", all.y = F)

mainDF = mergedDF[,-1]
rownames(mainDF) = mergedDF[,1]

## split data into a train and test set
index <- 1:(nrow(mainDF))
testindex <- sample(index, trunc(length(index)/3))
testset <- na.omit(mainDF[testindex,])
trainset <- na.omit(mainDF[-testindex,])

## svm


run.svm = function(cols.to.keep, trainset, testset, VOI){
  svm.model <- svm(x = data.matrix(trainset[,cols.to.keep]), y = as.factor(trainset[,VOI]), probability = T)
  svm.pred <- predict(svm.model, data.matrix(testset[,cols.to.keep]), probability = T)
  print(table(pred = svm.pred, true = testset[,VOI]))
  return(svm.pred)
}

##

run.svm.with.plot = function(cols.to.keep, trainset, testset, VOI){
  svm.res = run.svm(cols.to.keep, trainset = trainset, testset = testset, VOI = VOI)
  resROC = pROC::roc(response = testset[,"Symtomatic"], predictor = attributes(svm.res)$probabilities[,1])
  resAuc = pROC::auc(response = testset[,"Symtomatic"], predictor = attributes(svm.res)$probabilities[,1])
  plot(resROC)                 
  return(list(roc = resROC, auc = resAuc, svm.res = svm.res))
}


###

predictor.list = list(AgeSex = c("Age", "Sex"), 
                      top10 = head(sig.prots, 10), 
                      sigProts = sig.prots, 
                      allProts = colnames(testset)[grep("\\|",colnames(testset))])

cols.to.keep = sig.prots

##########

allBinsResults = lapply(predictor.list, run.svm.with.plot, trainset = trainset, testset = testset, VOI = "Symtomatic")

plot(allBinsResults[[1]]$roc, col = "white", main = "Prediction of Symptomatic status")
lapply(1:length(allBinsResults), function(RES){
  lines(allBinsResults[[RES]]$roc, col = rainbow(length(allBinsResults))[RES])
})

legend('bottomright', border = rainbow(length(allBinsResults)), fill = rainbow(length(allBinsResults)), legend = paste(names(allBinsResults), "| auc =", round(unlist(lapply(allBinsResults, "[[", "auc")), 3)))

##########
# just peripheral
allBinsResults = lapply(predictor.list, run.svm.with.plot, trainset = trainset[gsub(".*_(.*)","\\1",rownames(trainset)) == "EP",], testset = testset[gsub(".*_(.*)","\\1",rownames(testset)) == "EP",], VOI = "Symtomatic")

plot(allBinsResults[[1]]$roc, col = "white", main = "Prediction of Symptomatic status - peripheral only")
lapply(1:length(allBinsResults), function(RES){
  lines(allBinsResults[[RES]]$roc, col = rainbow(length(allBinsResults))[RES])
})

legend('bottomright', border = rainbow(length(allBinsResults)), fill = rainbow(length(allBinsResults)), legend = paste(names(allBinsResults), "| auc =", round(unlist(lapply(allBinsResults, "[[", "auc")), 3)))


## Metabolomics QC
metabolite_annots = readxl::read_xlsx("datasets/KI/Metabolomics/KARO-05-20MD DATA TABLES Sample Meta Data PseudoID.xlsx", 2)
sample_annots = readxl::read_xlsx("datasets/KI/Metabolomics/KARO-05-20MD DATA TABLES Sample Meta Data PseudoID.xlsx", 3)
log2vals = readxl::read_xlsx("datasets/KI/Metabolomics/KARO-05-20MD DATA TABLES Sample Meta Data PseudoID.xlsx", 7)
l2v = signif(data.matrix(log2vals[,-1]),4)
rownames(l2v) = unlist(log2vals[,1])
l2v = t(l2v)


### boxplot
data.frame(unlist(metabolite_annots[grep("sphing",metabolite_annots$PLOT_NAME),"PLOT_NAME"]))

#query_metabolites = "sphingomyelin (d18:1/22:2, d18:2/22:1, d16:1/24:2)*"
query_metabolites = unlist(metabolite_annots[grep("sphingomyelin \\(",metabolite_annots$PLOT_NAME),"PLOT_NAME"])
query_metabolites = unlist(metabolite_annots[grep("ceramide",metabolite_annots$PLOT_NAME),"PLOT_NAME"])

par(mfrow = c(3,3))

lapply(query_metabolites, function(query_metabolite){
  
  subMA = metabolite_annots[metabolite_annots$PLOT_NAME == query_metabolite,]
  ChID = subMA$CHEM_ID
  query_vals = l2v[rownames(l2v) == ChID,]
  
  subSA = sample_annots[match(names(query_vals), sample_annots$PARENT_SAMPLE_NAME), c(1, 3, 6)]
  subSA$groups = paste0(subSA$BiKE.ID.prefix, " ", subSA$Symtomatic)
  
  split.vals = split(query_vals, subSA$groups)
  split.vals = split.vals[!grepl("NA", names(split.vals))]
  
  EP.t.test.p = t.test(split.vals$`EP AS`, split.vals$`EP S`)$p.value
  STP.t.test.p = t.test(split.vals$`STP AS`, split.vals$`STP S`)$p.value
  min.p = signif(min(EP.t.test.p, STP.t.test.p),3)
  
  #if(min.p < 0.0001)
    
  boxplot(split.vals, col = "grey", main = paste0(query_metabolite, "\n T-test p-value = ", min.p))
})

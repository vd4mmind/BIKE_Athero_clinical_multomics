#KI metadata
metabolon_meta = read.delim("datasets/KI/Metabolomics/Metabolon_meta.txt")
olink_meta = read.delim("datasets/KI/Proteomics/Olink_meta.txt", header = F)
colnames(olink_meta) = as.character(unlist(olink_meta[1,]))
olink_meta = olink_meta[-1,]

head(metabolon_meta)

mms = split(metabolon_meta, metabolon_meta$BiKE.ID.prefix)
ols = split(olink_meta, olink_meta$BiKE.ID.prefix)

mmuni = lapply(list(mmEP = mms$EP, mmSTP = mms$STP), function(sub){ 
  sss = sub[,c(3,4,5,11)]
  uni_combos = apply(sss, 1, function(x){paste(as.character(x), collapse =  " ")})
  
  print(table(table(uni_combos)))
  return(uni_combos)
})

omuni = lapply(list(olsEP = ols$EP, olsSTP = ols$STP), function(sub){ 
  sss = sub[,c(1,2,3,7)]
  uni_combos = apply(sss, 1, function(x){paste(as.character(x), collapse = " ")})
  
  print(table(table(uni_combos)))
  return(uni_combos)
})

sum(omuni$olsEP %in% mmuni$mmEP)
sum(omuni$olsSTP %in% mmuni$mmSTP)

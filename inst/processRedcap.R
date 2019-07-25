#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

.oldf<<-args[2]



system(paste0(" awk '{if(NR==1)sub(/^\\xef\\xbb\\xbf/,\"\");print}' ",.oldf," > ",.oldf,"nobomb.R"))
source(paste0(args[2],"nobomb.R"))

annot=data






f=gsub("R$","RData",gsub("R","DATA",.oldf))
tmp=file.remove(paste0(.oldf,"nobomb.R"))


annot=annot[,which(!paste0(colnames(annot) ,".factor")%in% colnames(annot))]

colnames(annot)=gsub(".factor$","",colnames(annot))
save(annot,file=paste0(f))
# print(scan(args[2],what="character"))




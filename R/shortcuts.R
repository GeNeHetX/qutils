qlimma = function(expressiondata,samplesgrp1,samplesgrp2,thresh=0.01,namesonly=T,addWillcox=F){
  library(limma)
  samplesgrp1=intersect(samplesgrp1,colnames(expressiondata))
  samplesgrp2=intersect(samplesgrp2,colnames(expressiondata))
  expressiondata = expressiondata[,c(samplesgrp2,samplesgrp1)]
  classfactor= factor(c(rep("a",length(samplesgrp2)),rep("b",length(samplesgrp1))))
  design=model.matrix(~classfactor)
  colnames(design)[2] ="avsb"
  contrast.matrix <- makeContrasts("a-b" ,levels=classfactor)
  fit <- lmFit(expressiondata,design)
  fit2 <- eBayes(fit)
  fulltab=topTable(fit2,adjust="BH",number=nrow(expressiondata))
  if(namesonly){
    return(rownames(fulltab)[which(fulltab$adj.P.Val < thresh)])
  }else{

    X=fulltab[which(fulltab$adj.P.Val < thresh),]
    if(addWillcox){

      X$Wilcoxon.Pvalue=unlist(mclapply(rownames(X),function(prob){
        wilcox.test(as.numeric(expressiondata[prob,samplesgrp1]),as.numeric(expressiondata[prob,samplesgrp2]))$p.value}))
    }
    return(X)
  }

}

